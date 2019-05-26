
library("foreach")
library("doParallel")
library("MASS")
library("lme4")
source("funs_np.R")

N <- 100 #  subjects numbers
cat("N:", N, "\n")
maxm_sub <- 10  # maximum of measurements for each subject
cat("maxm_sub:", maxm_sub, "\n")

## true values
rbta1 <- c(-1.5, 2.0, 1.0)
rbta2 <- c(-1.5, 2.0, 1.0)
rbta3 <- c(-0.5, 2.0, 1.0)
lbta <- length(rbta1)
rphi <- c(2.0, 2.0) # standard devition of Gaussian reponses
rrho <- c(0.5, 0.5, 0.3)

# the true random effect distribution is Gaussian distribution
rsgm <- sqrt(1.25) 
# the true random effect distribution is a mixture of two Gaussian distributions
# rmu <- c(-1, 1)
# rsgm1 <- c(0.5, 0.5)
# rpii <- rmu[2]/(rmu[2]-rmu[1])

rpara <- c(rsgm, rphi, rbta1, rbta2, rbta3, rrho)
cat("rpara:",rpara,"\n")

write.table(t(rpara), file="res.txt", append=T, row.names=F, col.names=F)

# parallel
num_ite <- 100 # simulations times
cat("num_ite:",num_ite,"\n")
cores <- 20
cl <- makeCluster(cores, outfile = "pro.txt")
registerDoParallel(cl, cores=cores)
chunk_size <- num_ite/cores

res_all <- foreach(core=1:cores, .combine="rbind", .inorder=FALSE, .packages = c("MASS","lme4")) %dopar%
{  # local data for results
  res <- matrix(0, nrow=chunk_size, ncol=200)
  
  for(x in ((core-1)*chunk_size+1):(core*chunk_size)) {
    set.seed(x)
    
    # generate the number of repeated measurements for each subject
    n <- rep(0, N)
    id <- numeric()
    tt <- numeric()  # time covariates
    
    s <- rep(0, N)  # group covariates
    for(i in 1:N){ 
      n[i] <- rbinom(1, maxm_sub-1, 0.8) + 1
      tt <- c(tt, 1:n[i])
      id <- c(id, rep(i,n[i]))
      if(i <= N/2) s[i] <- 1 
    }
    M <- sum(n)
    cn <- c(0, cumsum(n))
    
    # generate the covariates X
    X <- matrix(NA, M, 2, dimnames=list(NULL, c('x1','x2')))
    X[,1] <-  (tt-5)/10
    X[,2] <- rep(s, n)
    
    Y <- matrix(NA, M, 3, dimnames=list(NULL, c('y1','y2','y3')))
    # sample
    for(i in 1:N){
      bi <- rnorm(1, 0, rsgm)
      # if(runif(1) < rpii) bi <- rnorm(1, rmu[1], rsgm1[1])
      # else bi <- rnorm(1, rmu[2], rsgm1[2])
      for(j in (cn[i]+1):cn[i+1]){
        mu1 <- drop(c(1,X[j,])%*%rbta1) + bi
        mu2 <- drop(c(1,X[j,])%*%rbta2) + bi 
        mu3 <- exp(drop(c(1,X[j,])%*%rbta3)+bi)/(1+exp(drop(c(1,X[j,])%*%rbta3)+bi))
        
        w <- runif(3)
        Y[j,1] <- rphi[1]*qnorm(w[1]) + mu1
        
        Y[j,2] <- rphi[2]*(qnorm(w[2])*sqrt(1-rrho[1]^2)+rrho[1]*(Y[j,1]-mu1)/rphi[1]) + mu2
        
        # P(y3==0|y2,y1)
        temp1 <- (-qnorm(mu3)-rrho[2]*(Y[j,2]-mu2)/rphi[2])/sqrt(1-rrho[2]^2)
        temp2 <- rrho[3]*((Y[j,1]-mu1)/rphi[1]-rrho[1]*(Y[j,2]-mu2)/rphi[2])/sqrt(1-rrho[1]^2)
        temp <- pnorm((temp1-temp2)/sqrt(1-rrho[3]^2))
        
        if(w[3] < temp) Y[j,3] <- 0
        else Y[j,3] <- 1
      } # j end
    } # i end
    
    
    data0 <- data.frame(id, Y, X)
    # margina model, used as initial valuse for joint estimation
    fit1 <- lmer(y1 ~ x1 + x2 + (1 | id), data = data0, REML = FALSE)
    
    fit2 <- lmer(y2 ~ x1 + x2 + (1 | id), data = data0, REML = FALSE)
    
    fit3 <- glmer(y3 ~ x1 + x2 + (1 | id), data = data0, family = binomial,
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
    
    phi <- c(summary(fit1)$sigma, summary(fit2)$sigma)
    bta1 <- as.numeric(fixef(fit1))
    bta2 <- as.numeric(fixef(fit2))
    bta3 <- as.numeric(fixef(fit3))
    
    # initial values for EM algorithm
    bta10 <- bta1[-1]
    bta20 <- bta2
    bta20[1] <- bta2[1] - bta1[1]
    bta30 <- bta3
    bta30[1] <- bta3[1] - bta1[1]
    phi0 <- phi
    rho0 <- c(0,0,0)
    
    K <- 20  # number of mass points
    pii0 <- rep(1/K, K)
    b0 <- seq(-3, 3, length.out = K) + bta1[1]
    
    para <- c(phi0, rho0, bta10, bta20, bta30, pii0, b0)
    lpara <- length(para)
    lpnob <- length(c(phi0, rho0, bta10, bta20, bta30))  
    
    # EM algorithm
    para0 <- para + 1
    mdif <- max(abs(para-para0))
    while (mdif > 1e-4) {
      para0 <- para
      phi0 <- para0[1:2]
      rho0 <- para0[3:5]
      bta10 <- para0[6:7]
      bta20 <- para0[8:10]
      bta30 <- para0[11:13]
      pii0 <- para0[14:(13+K)]
      b0 <- para0[(14+K):lpara]
      pnop0 <- para0[-(14:(13+K))]  # parameters without pi
      
      #  compute Omg
      Omg <- matrix(NA, N, K)
      for (i in 1:N) {
        tt <- rep(1, K)
        for (j in (cn[i]+1):cn[i+1]) {
          tt <- tt*sapply(b0, lik_ijk, y=Y[j,], x=X[j,], bta1=bta10, bta2=bta20, bta3=bta30, phi=phi0, rho=rho0)
        }
        Omg[i,] <- (pii0*tt)/(sum(pii0*tt))
      }
      
      # update pi
      pii <- colSums(Omg)/N
      
      # update other parameters
      pnop0[1:2] <- log(pnop0[1:2])
      pnop0[3:5] <- tran1(pnop0[3:5])
      
      mle <- optim(par=pnop0, loglik_EM, Y=Y, X=X, Omg=Omg, n=n, method = "BFGS", control = list(fnscale=-1))
      pnop <- mle$par
      pnop[1:2] <- exp(pnop[1:2])
      pnop[3:5] <- tran2(pnop[3:5])
      
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
      }
      else{
        mdif <- 1
        K <- length(b)
      }
    } # EM algorithm finished
    
    ml <- loglik_np(para[-(13+K)], Y, X, n)
    
    # compute the observed information matrix (Louis, 1982)
    M1 <- -optimHess(para[-(13+K)], Ellik, Y=Y, X=X, Omg=Omg, n=n)
    
    delta <- 1e-5
    Q <- array(NA, dim=c(N, K, lpara-1))
    for(i in 1:N){
      for(k in 1:K){
        for(l in 1:(lpara-1)){
          para1 <- para[-(13+K)]
          para1[l] <- para1[l] + delta
          Q[i,k,l] <- (g(Y[(cn[i]+1):cn[i+1],],X[(cn[i]+1):cn[i+1],],para1,k)-g(Y[(cn[i]+1):cn[i+1],],X[(cn[i]+1):cn[i+1],],para[-(13+K)],k))/delta
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
    
    print(x)
    print(ml)
    print(para)
    print(std_L)
    res[x-(core-1)*chunk_size, 1:lpara] <- para
    res[x-(core-1)*chunk_size, 101:(99+lpara)] <- std_L
    res[x-(core-1)*chunk_size, 200] <- ml
  } # x end
  res
}
stopCluster(cl)
stopImplicitCluster()

write.table(res_all,file="res.txt",append=T,row.names=F,col.names=F)
