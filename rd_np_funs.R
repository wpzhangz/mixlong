# density function for Gaussian copula
d_Gaus <- function(w1, w2, theta){
  i1 <- theta^2*(w1^2+w2^2) - 2*theta*w1*w2
  i1 <- i1/(2*(1-theta^2))
  i2 <- 1/sqrt(1-theta^2)
  return(i2*exp(-i1))
}

# transform rho from limited values to unlimited values
tran1 <- function(x){
  x <- as.vector(x)
  z <- log((1+x)/(1-x))
  z
}

# the inverse function of tran1
tran2 <- function(x){
  x <- as.vector(x)
  z <- (exp(x)-1)/(exp(x)+1)
  z
}


# combine mass points if the distance between them falls below a specified value
comb <- function(pii, b){
  lb <- length(b)
  d <- 0.05
  
  k <- 1
  idx <- rep(1, lb)
  for(i in 2:lb){
    if((b[i]-b[i-1]) > d) k <- k + 1
    idx[i] <- k 
  }
  
  t <- as.numeric(table(idx))
  lt <- length(t)
  cpii <- rep(0, lt)
  cb <- rep(0, lt)
  ct <- c(0, cumsum(t))
  for(i in 1:lt){
    cb[i] <- sum(b[(ct[i]+1):ct[i+1]])/t[i]
    cpii[i] <- sum(pii[(ct[i]+1):ct[i+1]])
  }
  
  aa <- (cpii > 1e-3)
  cpii <- cpii[aa]
  cpii <- cpii/sum(cpii)
  cb <- cb[aa]
  return(list(pii=cpii, b=cb))
}

# compute the value of D-vine for three-dimensional mixed case
# y1, y2 ~ normal distribution; y3 ~ binary distribution
lik_ijk <- function(y, x, bta1, bta2, bta3, lmd, phi, rho, b){
  mu1 <- sum(x*bta1) + b
  mu2 <- sum(c(1,x)*bta2) + lmd[1]*b 
  mu3 <- exp(sum(c(1,x)*bta3)+lmd[2]*b)/(1+exp(sum(c(1,x)*bta3)+lmd[2]*b))
  
  # P(Y3==0|y2,y1)
  if(is.nan(mu3) | mu3 == 1){
    # 1-mu3
    mu3_n <- 1/(1+exp(sum(c(1,x)*bta3)+lmd[2]*b))
    temp1 <- (qnorm(mu3_n)-rho[2]*(y[2]-mu2)/phi[2])/sqrt(1-rho[2]^2)
  }
  else temp1 <- (-qnorm(mu3)-rho[2]*(y[2]-mu2)/phi[2])/sqrt(1-rho[2]^2)
  temp2 <- rho[3]*((y[1]-mu1)/phi[1]-rho[1]*(y[2]-mu2)/phi[2])/sqrt(1-rho[1]^2)
  temp <- pnorm((temp1-temp2)/sqrt(1-rho[3]^2))
  if(y[3] == 1) temp <- pnorm(-(temp1-temp2)/sqrt(1-rho[3]^2))
  
  w1 <- (y[1]-mu1)/phi[1]
  w2 <- (y[2]-mu2)/phi[2]
  z <- dnorm(y[1], mu1, phi[1])*dnorm(y[2], mu2, phi[2])*d_Gaus(w2, w1, rho[1])*temp
  return(z)
}

# used in EM algorithm
loglik_EM <- function(para, Y, X, Omg, n){
  N <- length(n)
  cn <- c(0, cumsum(n))
  lpara <- length(para)
  
  lmd <- para[1:2]
  phi <- para[3:4]
  rho <- para[5:7]
  bta1 <- para[8:11]
  bta2 <- para[12:16]
  bta3 <- para[17:21]
  b <- para[22:lpara]
  
  phi <- exp(phi)
  rho <- tran2(rho)
  
  z <- 0
  for (i in 1:N) {
    for (j in (cn[i]+1):cn[i+1]) {
      z <- z + sum(Omg[i,]*log(sapply(b, lik_ijk, y=Y[j,], x=X[j,], bta1=bta1, bta2=bta2, bta3=bta3, lmd=lmd, phi=phi, rho=rho)))
    }
  }
  z
}

# log-likelihood function for Nonparametric model
loglik_np <- function(para, Y, X, n){
  N <- length(n)
  cn <- c(0, cumsum(n))
  lpara <- length(para)
  
  lmd <- para[1:2]
  phi <- para[3:4]
  rho <- para[5:7]
  bta1 <- para[8:11]
  bta2 <- para[12:16]
  bta3 <- para[17:21]
  K <- (lpara-20)/2
  pii <- para[22:(20+K)]
  pii <- c(pii, 1-sum(pii))
  b <- para[(21+K):lpara]
  
  llik <- 0
  for (i in 1:N) {
    lik_i <- rep(1, K)
    for (j in (cn[i]+1):cn[i+1]) {
      lik_i <- lik_i*sapply(b, lik_ijk, y=Y[j,], x=X[j,], bta1=bta1, bta2=bta2, bta3=bta3, lmd=lmd, phi=phi, rho=rho)
    }
    llik_i <- log(sum(pii*lik_i))
    llik <- llik + llik_i
  }
  return(llik)
}

# the expectation of the complete log-likelihood function
Ellik <- function(para, Y, X, Omg, n){
  N <- length(n)
  cn <- c(0, cumsum(n))
  lpara <- length(para)
  
  lmd <- para[1:2]
  phi <- para[3:4]
  rho <- para[5:7]
  bta1 <- para[8:11]
  bta2 <- para[12:16]
  bta3 <- para[17:21]
  K <- (lpara-20)/2
  pii <- para[22:(20+K)]
  pii <- c(pii, 1-sum(pii))
  b <- para[(21+K):lpara]
  
  z <- 0
  for(i in 1:N){
    for(k in 1:K){
      llik_ik <- 0
      for(j in (cn[i]+1):cn[i+1]){
        llik_ik <- llik_ik + log(lik_ijk(y=Y[j,], x=X[j,], bta1=bta1, bta2=bta2, bta3=bta3, lmd=lmd, phi=phi, rho=rho, b=b[k]))
      }
      z <- z + Omg[i,k]*(llik_ik+log(pii[k]))
    }
  }
  return(z)
}

# compute the log-likelihood value for subject i under the kth mass point
g <- function(Y, X, para, k){
  lmd <- para[1:2]
  phi <- para[3:4]
  rho <- para[5:7]
  bta1 <- para[8:11]
  bta2 <- para[12:16]
  bta3 <- para[17:21]
  K <- (lpara-20)/2
  pii <- para[22:(20+K)]
  pii <- c(pii, 1-sum(pii))
  b <- para[(21+K):lpara]
  
  z <- 0
  if(is.null(dim(Y))) z <- log(lik_ijk(y=Y, x=X, bta1=bta1, bta2=bta2, bta3=bta3, lmd=lmd, phi=phi, rho=rho, b=b[k]))
  else{
    for(j in 1:dim(Y)[1]){
      z <- z + log(lik_ijk(y=Y[j,], x=X[j,], bta1=bta1, bta2=bta2, bta3=bta3, lmd=lmd, phi=phi, rho=rho, b=b[k]))
    }
  }
  z <- z + log(pii[k])
  return(z)
}


