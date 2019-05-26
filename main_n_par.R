
library("fastGHQuad")
library("foreach")
library("doParallel")
library("MASS")
library("lme4")
source("funs.R")

N <- 100 #  subjects numbers
cat("N:", N, "\n")
maxm_sub <- 10  # maximum of measurements for each subject
cat("maxm_sub:", maxm_sub, "\n")
num_qua <- 20
cat("num_qua:",num_qua,"\n")
rule <- gaussHermiteData(num_qua) # number of quadrature points

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
lpara <- length(rpara)

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
  res <- matrix(0, nrow=chunk_size, ncol=2*lpara+1)
  
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
    
    sgm0 <- attr(summary(fit1)$varcor$id,"stddev")
    phi0 <- c(summary(fit1)$sigma, summary(fit2)$sigma)
    bta0 <- c(fixef(fit1), fixef(fit2), fixef(fit3))
    rho0 <- c(0,0,0)
    
    lower = as.numeric(c(sgm0-0.5, phi0-1, bta0-1, -0.8, -0.8, -0.8))
    upper = as.numeric(c(sgm0+0.5, phi0+1, bta0+1, 0.8, 0.8, 0.8))
    start = as.numeric(c(sgm0, phi0, bta0, rho0))  
    
    # loglik(para=rpara, X1=X1, X2=X2, X3=X3, Y=Y, rule=rule)
    mle <- optim(par=start, loglik, X=X, Y=Y, rule=rule, n=n, method = "L-BFGS-B",
                 lower = lower, upper = upper, control = list(fnscale = -1), hessian = T)
    para <- mle$par
    para_std <- sqrt(diag(ginv(-mle$hessian)))
    
    print(x)
    print(para)
    print(para_std)
    res[x-(core-1)*chunk_size, 1:lpara] <- para
    res[x-(core-1)*chunk_size, (lpara+1):(2*lpara)] <- para_std
    res[x-(core-1)*chunk_size, 2*lpara+1] <- mle$value
  } # x end
  res
}
stopCluster(cl)
stopImplicitCluster()

write.table(res_all,file="res.txt",append=T,row.names=F,col.names=F)
