
library("MASS")
library("mixAK")
library("lme4")
source("rd_np_funs.R")

## output file
filename <- "rd_np.txt"

data("PBCseq")
mypbc <- PBCseq[,c("id","sex","drug","age","month","lalbumin","lbili", "hepatom")]
# dim(mypbc) 1945 x 17

# Delete cases with NA value
na_obs <- which(is.na(mypbc$hepatom))
na_id <- unique(mypbc$id[sort(unique(na_obs))])
mypbc_na_ex <- mypbc[!(mypbc$id %in% na_id),]

n <- as.numeric(table(mypbc_na_ex$id))  # measure times for each subject
N <- length(n)    # total 252 patients
cn <- c(0, cumsum(n))

# years after enrollment 
years <- mypbc_na_ex$month/12
mypbc_na_ex <- cbind(mypbc_na_ex, years)


#------------------------------------------------------------------
# margina model

fit1 <- lmer(lbili ~ sex + drug + age + years + (1 | id), data = mypbc_na_ex, REML = FALSE)

fit2 <- lmer(lalbumin ~ sex + drug + age + years +  (1 | id), data = mypbc_na_ex, REML = FALSE)

fit3 <- glmer(hepatom ~ sex + drug + age + years + (1 | id), data = mypbc_na_ex, family = binomial,
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))



#------------------------------------------------------------------
# joint model

Y <- mypbc_na_ex[,c("lbili","lalbumin","hepatom")]
Y <- as.matrix(Y)
# cor(Y, method = "kendall")
X <- mypbc_na_ex[,c("sex","drug","age","years")]
X <- as.matrix(X)
cat("Y:", colnames(Y), "\n")
cat("X:", colnames(X), "\n")

# initial value
sgm <- attr(summary(fit1)$varcor$id,"stddev")
lmd <- c(-attr(summary(fit2)$varcor$id,"stddev"), attr(summary(fit3)$varcor$id,"stddev"))/sgm
phi <- c(summary(fit1)$sigma, summary(fit2)$sigma)
bta1 <- fixef(fit1)
bta2 <- fixef(fit2)
bta2[1] <- bta2[1] - bta1[1]
bta3 <- fixef(fit3)
bta3[1] <- bta3[1] - bta1[1]
rho <- c(0,0,0)

K <- 20  # number of mass points
cat("K:", K, "\n")
pii <- rep(1/K, K)
b <- seq(-3*sgm, 3*sgm, length.out = K) + bta1[1]

para <- c(lmd, phi, rho, bta1, bta2, bta3, pii, b)
lpara <- length(para)
lpnob <- length(c(lmd, phi, rho, bta1, bta2, bta3))
cat("start values:", para, "\n")

para0 <- para + 1
mdif <- max(abs(para-para0))
iter <- 0
ml <- 0
while (mdif > 1e-4) {
  para0 <- para
  lmd0 <- para0[1:2]
  phi0 <- para0[3:4]
  rho0 <- para0[5:7]
  bta10 <- para0[8:11]
  bta20 <- para0[12:16]
  bta30 <- para0[17:21]
  pii0 <- para0[22:(21+K)]
  b0 <- para0[(22+K):lpara]
  pnop0 <- para0[-(22:(21+K))]
  
  #  compute Omg
  Omg <- matrix(NA, N, K)
  for (i in 1:N) {
    tt <- rep(1, K)
    for (j in (cn[i]+1):cn[i+1]) {
      tt <- tt*sapply(b0, lik_ijk, y=Y[j,], x=X[j,], bta1=bta10, bta2=bta20, bta3=bta30, lmd=lmd0, phi=phi0, rho=rho0)
    }
    Omg[i,] <- (pii0*tt)/(sum(pii0*tt))
  }
  
  # update pi
  pii <- colSums(Omg)/N
  
  # update other parameters
  pnop0[3:4] <- log(pnop0[3:4])
  pnop0[5:7] <- tran1(pnop0[5:7])
  
  mle <- optim(par=pnop0, loglik_EM, Y=Y, X=X, Omg=Omg, n=n, method = "BFGS", control = list(fnscale=-1))
  pnop <- mle$par
  pnop[3:4] <- exp(pnop[3:4])
  pnop[5:7] <- tran2(pnop[5:7])
  
  iter <- iter + 1
  cat("iteration number:", iter, "\n")
  
  b <- pnop[(lpnob+1):(lpara-K)]
  pii <- pii[order(b)]
  b <- sort(b)

  # combine mass points
  cpiib <- comb(pii,b)
  pii <- cpiib$pii
  b <- cpiib$b

  para <- c(pnop[1:lpnob], pii, b)
  lpara <- length(para)
  if(length(b) == K){
    mdif <- max(abs(para-para0))
    mdif_index <- which.max(abs(para-para0))
    cat("max difference:", mdif, mdif_index, "\n")
  }
  else{
    mdif <- 1
    K <- length(b)
    cat("K:", K, "\n")
  }
  
  ml <- loglik_np(para[-(lpnob+K)], Y, X, n)
  cat("loglik:", ml, "\n")
  # cat("para:", para, "\n")
  cat("\n")
} # EM algorithm finished

write.table(paste("loglik:", ml), file=filename, append=T, row.names=F, col.names=F)
write.table(t(para), file=filename, append=T, row.names=F, col.names=F)

# compute the observed information matrix (Louis, 1982)
M1 <- -optimHess(para[-(lpnob+K)], Ellik, Y=Y, X=X, Omg=Omg, n=n)
delta <- 1e-5
Q <- array(NA, dim=c(N, K, lpara-1))
for(i in 1:N){
  for(k in 1:K){
    for(l in 1:(lpara-1)){
      para1 <- para[-(lpnob+K)]
      para1[l] <- para1[l] + delta
      Q[i,k,l] <- (g(Y[(cn[i]+1):cn[i+1],],X[(cn[i]+1):cn[i+1],],para1,k)-g(Y[(cn[i]+1):cn[i+1],],X[(cn[i]+1):cn[i+1],],para[-(lpnob+K)],k))/delta
    }
  }
}

M2 <- matrix(0, lpara-1, lpara-1)
for(i1 in 1:N){
  for(i2 in 1:N){
    if(i2 == i1) next
    for(k1 in 1:K){
      for(k2 in 1:K){
        M2 <- M2 + Omg[i1,k1]*Omg[i2,k2]*(Q[i1,k1,]%*%t(Q[i2,k2,]))
      }
    }
  }
}
for(i in 1:N){
  for(k in 1:K){
    M2 <- M2 + Omg[i,k]*(Q[i,k,]%*%t(Q[i,k,]))
  }
}
H_L <- M1 - M2
std_L <- sqrt(diag(ginv(H_L)))
write.table(t(std_L), file=filename, append=T, row.names=F, col.names=F)
