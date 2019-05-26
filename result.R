

setwd("D:/study/copula/simu/new-2019-03")
file <- "res.txt"  
rpara <- read.table(file, nrows=1)
a <- read.table(file, skip=1)

##---------------------------------------------------------
## Gaussian model

a1 <- a[,1:(dim(a)[2]-1)]
k <- dim(a1)[2]/2  # number of parameters

# MEAN
MEAN <- apply(a1[,1:k], 2, mean)

# SE
SE <- apply(a1[,(k+1):(2*k)], 2, mean)

# SD
SD <- apply(a1[,1:k], 2, sd)

# CP
u <- qnorm(0.025)
CP <- rep(0, k)
for(j in 1:k){
  for(i in 1:dim(a1)[1]){
    lower <- rpara[j] + u*a1[i,j+k]
    upper <- rpara[j] - u*a1[i,j+k]
    if(a1[i,j]>lower & a1[i,j]<upper) CP[j] <- CP[j] + 1
  }
}
CP <- CP/dim(a1)[1]

# MAB
sum <- 0
for(i in 1:dim(a1)[1])
  sum <- sum + abs(a1[i,1:k]-rpara)
MAB <- sum/dim(a1)[1]

# RMSE
sum <- 0
for(i in 1:dim(a1)[1])
  sum <- sum + (a1[i,1:k]-rpara)^2
RMSE <- sqrt(sum/dim(a1)[1])

# output
summ <- round(rbind(rpara, MEAN, SE, SD, CP, MAB, RMSE), 4)
write.table(summ, file="summ_n.txt")


##--------------------------------------------------------
## nonparametric model

rhi <- rpara[2:3]
rbta1 <- rpara[4:6]
rbta2 <- rpara[7:9]
rbta2[1] <- rbta2[1] - rbta1[1]
rbta3 <- rpara[10:12]
rbta3[1] <- rbta3[1] - rbta1[1]
rrho <- rpara[13:15]

rpara <- c(rhi, rbta1[-1], rbta2, rbta3, rrho)
k <- length(rpara)  # number of parameters

a1 <- a[,c(1:k,101:(100+k))]

# MEAN
MEAN <- apply(a1[,1:k], 2, mean)

# SE
SE <- apply(a1[,(k+1):(2*k)], 2, mean)

# SD
SD <- apply(a1[,1:k], 2, sd)

# CP
u <- qnorm(0.025)
CP <- rep(0, k)
for(j in 1:k){
  for(i in 1:dim(a1)[1]){
    lower <- rpara[j] + u*a1[i,j+k]
    upper <- rpara[j] - u*a1[i,j+k]
    if(a1[i,j]>lower & a1[i,j]<upper) CP[j] <- CP[j] + 1
  }
}
CP <- CP/dim(a1)[1]

# MAB
sum <- 0
for(i in 1:dim(a1)[1])
  sum <- sum + abs(a1[i,1:k]-rpara)
MAB <- sum/dim(a1)[1]

# RMSE
sum <- 0
for(i in 1:dim(a1)[1])
  sum <- sum + (a1[i,1:k]-rpara)^2
RMSE <- sqrt(sum/dim(a1)[1])

# output
summ <- round(rbind(rpara, MEAN, SE, SD, CP, MAB, RMSE), 4)
write.table(summ, file="summ_np.txt")


