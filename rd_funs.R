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

# compute the value of D-vine for three-dimensional mixed case
# y1, y2 ~ normal distribution; y3 ~ binary distribution
lik_ijk <- function(y, x, bta1, bta2, bta3, lmd, phi, rho, b){
  mu1 <- sum(c(1,x)*bta1) + b
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

# funtion to be maximized in the 'optim' function
loglik_t <- function(para, Y, X, rule, n){
  N <- length(n)
  cn <- c(0, cumsum(n))
  num_qua <- length(rule$x)
  
  sgm <- para[1]
  lmd <- para[2:3]
  phi <- para[4:5]
  rho <- para[6:8]
  bta1 <- para[9:13]
  bta2 <- para[14:18]
  bta3 <- para[19:23]
  
  phi <- exp(phi)
  rho <- tran2(rho)
  
  # compute loglik
  llik <- 0
  for(i in 1:N){
    # f(yi)
    lik_ik <- rep(1, num_qua)
    bi <- sqrt(2)*sgm*rule$x
    for(j in (cn[i]+1):cn[i+1]){
      lik_ik <- lik_ik*sapply(bi, lik_ijk, y=Y[j,], x=X[j,], bta1=bta1, bta2=bta2, bta3=bta3, lmd=lmd, phi=phi, rho=rho)
    }
    lik_i <- sum(rule$w*lik_ik)
    if(is.nan(lik_i) | lik_i == 0) lik_i <- 1e-300
    llik <- llik + log(lik_i) - 1/2*log(pi)
  } # i end
  return(llik)
}

# log-likelihood function
loglik <- function(para, Y, X, rule, n){
  N <- length(n)
  cn <- c(0, cumsum(n))
  num_qua <- length(rule$x)
  
  sgm <- para[1]
  lmd <- para[2:3]
  phi <- para[4:5]
  rho <- para[6:8]
  bta1 <- para[9:13]
  bta2 <- para[14:18]
  bta3 <- para[19:23]
  
  # compute loglik
  llik <- 0
  for(i in 1:N){
    # f(yi)
    lik_ik <- rep(1, num_qua)
    bi <- sqrt(2)*sgm*rule$x
    for(j in (cn[i]+1):cn[i+1]){
      lik_ik <- lik_ik*sapply(bi, lik_ijk, y=Y[j,], x=X[j,], bta1=bta1, bta2=bta2, bta3=bta3, lmd=lmd, phi=phi, rho=rho)
    }
    lik_i <- sum(rule$w*lik_ik)
    if(is.nan(lik_i) | lik_i == 0) lik_i <- 1e-300
    llik <- llik + log(lik_i) - 1/2*log(pi)
  } # i end
  return(llik)
}




