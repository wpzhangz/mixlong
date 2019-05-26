# density function for Gaussian copula
d_Gaus <- function(w1, w2, theta){
  i1 <- theta^2*(w1^2+w2^2) - 2*theta*w1*w2
  i1 <- i1/(2*(1-theta^2))
  i2 <- 1/sqrt(1-theta^2)
  return(i2*exp(-i1))
}

# h function for Gaussian copula
h <- function(u1, u2, theta){
  i1 <- qnorm(u1) - theta*qnorm(u2)
  i2 <- 1/sqrt(1-theta^2)
  return(pnorm(i1*i2))
}

# h inverse function for Gaussian copula
h_inv <- function(u1, u2, theta){
  i1 <- qnorm(u1)*sqrt(1-theta^2)
  i2 <- theta*qnorm(u2)
  return(pnorm(i1+i2))
}

# compute the value of D-vine for three-dimensional mixed case
# y1, y2 ~ normal distribution; y3 ~ binary distribution
lik_ijk <- function(y, x, bta1, bta2, bta3, phi, rho, b){
  mu1 <- sum(c(1,x)*bta1) + b
  mu2 <- sum(c(1,x)*bta2) + b 
  mu3 <- exp(sum(c(1,x)*bta3)+b)/(1+exp(sum(c(1,x)*bta3)+b))
  
  # P(Y3==0|y2,y1)
  if(is.nan(mu3) | mu3 == 1){
    # 1-mu3
    mu3_n <- 1/(1+exp(sum(c(1,x)*bta3)+b))
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

# log-likelihood function for Gaussian model
loglik <- function(para, X, Y, rule, n){
  N <- length(n)
  cn <- c(0, cumsum(n))
  lpara <- length(para)
  num_qua <- length(rule$x)
  
  sgm <- para[1]
  phi <- para[2:3]
  bta1 <- para[4:6]
  bta2 <- para[7:9]
  bta3 <- para[10:12]
  rho <- para[13:15]
  
  # compute loglik
  llik <- 0
  for(i in 1:N){
    # f(yi)
    lik_ik <- rep(1, num_qua)
    bi <- sqrt(2)*sgm*rule$x
    for(j in (cn[i]+1):cn[i+1]){
      lik_ik <- lik_ik*sapply(bi, lik_ijk, y=Y[j,], x=X[j,], bta1=bta1, bta2=bta2, bta3=bta3, phi=phi, rho=rho)
    }
    lik_i <- sum(rule$w*lik_ik)
    
    if(is.nan(lik_i) | lik_i == 0) lik_i <- 1e-300
    llik <- llik + log(lik_i) - 1/2*log(pi)
  } # i end
  return(llik)
}


