
source("rd_funs.R")
library("fastGHQuad")
library("MASS")
library("mixAK")
library("lme4")

## output file
filename <- "rd.txt"

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
X <- mypbc_na_ex[,c("sex","drug","age","years")]
X <- as.matrix(X)
cat("Y:", colnames(Y), "\n")
cat("X:", colnames(X), "\n")

num_qua <- 80
cat("num_qua:", num_qua, "\n")
rule <- gaussHermiteData(num_qua) # number of quadrature points

# initial value
sgm0 <- attr(summary(fit1)$varcor$id,"stddev")
lmd0 <- c(-attr(summary(fit2)$varcor$id,"stddev"), attr(summary(fit3)$varcor$id,"stddev"))/sgm0
phi0 <- c(summary(fit1)$sigma, summary(fit2)$sigma)
bta0 <- c(fixef(fit1), fixef(fit2), fixef(fit3))
rho0 <- c(0,0,0)

para0 <- c(sgm0, lmd0, phi0, rho0, bta0)
names(para0) <- NULL
# loglik(para0, Y, X, rule, n)
cat("para0:", para0, "\n")
write.table(t(para0),file=filename,append=T,row.names=F,col.names=F)

para0[4:5] <- log(para0[4:5])
para0[6:8] <- tran1(para0[6:8])
mle <- optim(par=para0, loglik_t, Y=Y, X=X, rule=rule, n=n, method = "BFGS", 
             control = list(fnscale=-1, trace=1, REPORT=1), hessian = T)
para <- mle$par
para[4:5] <- exp(para[4:5])
para[6:8] <- tran2(para[6:8])

write.table(paste("loglik:", mle$value), file=filename, append=T, row.names=F, col.names=F)
write.table(t(round(para,4)), file=filename, append=T, row.names=F, col.names=F)

para_hess <- optimHess(para, loglik, Y=Y, X=X, rule=rule, n=n)
para_std <- sqrt(diag(ginv(-para_hess)))

write.table(t(round(para_std,4)), file=filename, append=T, row.names=F, col.names=F)




