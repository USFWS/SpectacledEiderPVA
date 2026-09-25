##Find equilibrium with optim

# 1.26.2022 updated some input values from IPM_MS, added comments
# 1.26.2022 also rescaled N in Atest to Btest to rescale betaN for better conver-
# gence diagnostics in JAGS

Atest = function(N=10000, betaN=0.01){
  ilogit <- function(x){exp(x)/(1+exp(x))}
  BP <- ilogit(2.2) # assumed mean BP was 0.9
  mean.phi0 <- 0.325 # phi.0 + (1-phi.0)*omegaJ, from cjs and EE
  mean.logit.phi0 <- log(mean.phi0/(1-mean.phi0))
  mean.phiA <- 0.9 # phi.A + (1-phi.A)*omegaA, from cjs and EE
  mean.logit.phiA <- log(mean.phiA/(1-mean.phiA))
  mean.log.F <- -0.47 # log(nest success*clutch size*duckling survival, Kig)
  phi0 <- ilogit(mean.logit.phi0)
  logit.phiA <- mean.logit.phiA - betaN*(N/1000) # rescaling N
  phiA <- ilogit(logit.phiA)
  F <- exp(log(phi0) + mean.log.F - betaN*(N/1000) + log(BP) ) # rescaling N
  phi2 <- phiA
  phi1 <- phi2
  alpha <- 1.95/7.1 # from EE
  
  A <- matrix(c(0, 0, F, F, 
                (1-alpha)*phi1, 0, 0, 0, 
                alpha*phi1, 0, 0, 0, 
                0, phi2, phiA, phiA), ncol = 4, byrow = T)
  
  return( (max(Re(eigen(A)$values))-1)^2 )
}

Btest = function(betaN=0.01, N=12500){
  ilogit <- function(x){exp(x)/(1+exp(x))}
  BP <- ilogit(2.2) # assumed mean BP was 0.9
  mean.phi0 <- 0.325 # phi.0 + (1-phi.0)*omegaJ, from cjs and EE
  mean.logit.phi0 <- log(mean.phi0/(1-mean.phi0))
  mean.phiA <- 0.9 # phi.A + (1-phi.A)*omegaA, from cjs and EE
  mean.logit.phiA <- log(mean.phiA/(1-mean.phiA))
  mean.log.F <- -0.47# log(nest success*clutch size*duckling survival, Kig)
  phi0 <- ilogit(mean.logit.phi0)
  logit.phiA <- mean.logit.phiA - betaN*(N/1000) # rescaling N
  phiA <- ilogit(logit.phiA)
  F <- exp(log(phi0) + mean.log.F - betaN*(N/1000) + log(BP) ) # rescaling N
  phi2 <- phiA
  phi1 <- phi2
  alpha <- 1.95/7.1 # from EE
  
  A <- matrix(c(0, 0, F, F, 
                (1-alpha)*phi1, 0, 0, 0, 
                alpha*phi1, 0, 0, 0, 
                0, phi2, phiA, phiA), ncol = 4, byrow = T)
  
  return( (max(Re(eigen(A)$values))-1)^2 )
}

optimize(f=Atest, interval = c(0, 1e7), betaN = 0.0001)

optimize(f=Btest, interval = c(0, 0.1), N=12500)
#-------------------------------------------------------------------------------
##relate this to the EE for the ACP
#high K = 50657 ibb = 50657/2 females = 25328.5
optimize(f=Btest, interval = c(0, 0.1), tol=1e-8, N=25329)$minimum
# [1] 0.008039082
#best K = 25000 ibb = 12500 females
optimize(f=Btest, interval = c(0, 0.1), tol=1e-8, N=12500)$minimum
# [1] 0.01628975
#low K = 10743 ibb = 5371.5 females
optimize(f=Btest, interval = c(0, 0.1), tol=1e-8, N=5372)$minimum
# [1] 0.0379043

##relate this to the EE for the YKD
#high females = 38294
optimize(f=Btest, interval = c(0, 0.1), tol=1e-8, N=38294)$minimum
# [1] 0.005317332
#best  = 19300 females
optimize(f=Btest, interval = c(0, 0.1), tol=1e-8, N=19300)$minimum
# [1] 0.01055036
#low = 8284  females
optimize(f=Btest, interval = c(0, 0.1), tol=1e-8, N=8284)$minimum
# [1] 0.02458014
#-------------------------------------------------------------------------------

#now find a distribution that matches the EE
qmatch.gamma <- function(par=c(3,10000), eeq=c(0.1, 0.7, 2.3), qtarget=c(0.1, 0.5, 0.9)){
  qdist <- qgamma(qtarget, shape=par[1], rate=par[2])
  return( sum((10000*(qdist-eeq))^2) )
}

#ACP
fit <- optim(par=c(3, 10000), fn=qmatch.gamma, method="L-BFGS-B", lower=c(1,1), 
             eeq=c(0.008039082, 0.01628975, 0.0379043))
fit
# $par
# [1]   2.5295 124.3611
# 
# $value
# [1] 430.4592
# 
# $counts
# function gradient 
# 58       58 
# 
# $convergence
# [1] 0
# 
# $message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

qgamma(c(0.1, 0.5, 0.9), shape=fit$par[1], rate=fit$par[2])
# [1] 0.006610512 0.017731524 0.037474102
#YKD
fit <- optim(par=c(3, 10000), fn=qmatch.gamma, method="L-BFGS-B", lower=c(1,1), 
             eeq=c(0.005317332, 0.01055036, 0.02458014))
fit
# $par
# [1]   2.562276 193.844387
# 
# $value
# [1] 203.4094
# 
# $counts
# function gradient 
# 78       78 
# 
# $convergence
# [1] 0
# 
# $message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
qgamma(c(0.1, 0.5, 0.9), shape=fit$par[1], rate=fit$par[2])
# [1] 0.004317903 0.011322344 0.023611622
hist(rgamma(10000, fit$par[1], fit$par[2]))

# Check results
# Plot results from draws of prior:
Nsim <- 1000
DD <- rgamma(Nsim, fit$par[1], fit$par[2])
K <- numeric(Nsim)
for(i in 1:Nsim){
  K[i] <- optimize(f=Atest, interval = c(0, 1e7), betaN = DD[i])$minimum
}
hist(K)
################################################################################
#try fitting distribution with the package rriskDistributions
#-------------------------------------------------------------------------------
#YKD
library(rriskDistributions)
fit <- get.gamma.par(p=c(0.1, 0.5, 0.9), q=c(0.005317332, 0.01055036, 0.02458014))
qgamma(c(0.1, 0.5, 0.9), shape=fit[1], rate=fit[2])
hist(rgamma(10000, fit[1], fit[2]), breaks = 100)
# > fit
# shape       rate 
# 3.109273 251.648023 
# Plot results from draws of prior:
Nsim <- 1000
DD <- rgamma(Nsim, fit[1], fit[2])
K <- numeric(Nsim)
for(i in 1:Nsim){
  K[i] <- optimize(f=Atest, interval = c(0, 1e7), betaN = DD[i])$minimum
}
hist(K)
summary(K)
summary(DD)
plot(DD, K)
plot(DD[K<2e6], K[K<2e6])
max(DD[K<2e6])
hist(DD[K<2e6])
#-------------------------------------------------------------------------------
# ACP
fit <- get.gamma.par(p=c(0.1, 0.5, 0.9), q=c(0.008039082, 0.01628975, 0.0379043))
fit
# shape       rate 
# 3.023973 158.296488 
qgamma(c(0.1, 0.5, 0.9), shape=fit[1], rate=fit[2])
hist(rgamma(10000, fit[1], fit[2]))
################################################################################
#humm, for strong DD, optimization doesn't work. Is there a maximum DD?
lambda = function(N=10000, betaN=0.01){
  ilogit <- function(x){exp(x)/(1+exp(x))}
  BP <- ilogit(2.2) # assumed mean BP was 0.9
  mean.phi0 <- 0.325 # phi.0 + (1-phi.0)*omegaJ, from cjs and EE
  mean.logit.phi0 <- log(mean.phi0/(1-mean.phi0))
  mean.phiA <- 0.9 # phi.A + (1-phi.A)*omegaA, from cjs and EE
  mean.logit.phiA <- log(mean.phiA/(1-mean.phiA))
  mean.log.F <- -0.47 # log(nest success*clutch size*duckling survival, Kig)
  phi0 <- ilogit(mean.logit.phi0)
  logit.phiA <- mean.logit.phiA - betaN*(N/1000) # rescaling N
  phiA <- ilogit(logit.phiA)
  F <- exp(log(phi0) + mean.log.F - betaN*(N/1000) + log(BP) ) # rescaling N
  phi2 <- phiA
  phi1 <- phi2
  alpha <- 1.95/7.1 # from EE
  
  A <- matrix(c(0, 0, F, F, 
                (1-alpha)*phi1, 0, 0, 0, 
                alpha*phi1, 0, 0, 0, 
                0, phi2, phiA, phiA), ncol = 4, byrow = T)
  
  return( max(Re(eigen(A)$values)) )
}

growth <- numeric(Nsim)
DD <- rgamma(Nsim, fit[1], fit[2])
for(i in 1:Nsim){
  growth[i] = lambda(N=1, betaN=DD[i])
}
summary(growth)
plot(DD, growth)
#Ah ha, basic pop bio!, the max DD needs to be bound so that r_max => 0 when N = 0
#In fact DD should be < lambda at zero population size
lambda(N=0, betaN=0.04)
#[1] 1.042269
#Take home, use results from quantile matching in rriskDistributions, above
# shape       rate 
# 3.109273 251.648023 
#-------------------------------------------------------------------------------
# Do population projections
A = function(N=10000, betaN=0.01){
  ilogit <- function(x){exp(x)/(1+exp(x))}
  BP <- ilogit(2.2) # assumed mean BP was 0.9
  mean.phi0 <- 0.325 # phi.0 + (1-phi.0)*omegaJ, from cjs and EE
  mean.logit.phi0 <- log(mean.phi0/(1-mean.phi0))
  mean.phiA <- 0.9 # phi.A + (1-phi.A)*omegaA, from cjs and EE
  mean.logit.phiA <- log(mean.phiA/(1-mean.phiA))
  mean.log.F <- -0.47 # log(nest success*clutch size*duckling survival, Kig)
  phi0 <- ilogit(mean.logit.phi0)
  logit.phiA <- mean.logit.phiA - betaN*(N/1000) # rescaling N
  phiA <- ilogit(logit.phiA)
  F <- exp(log(phi0) + mean.log.F - betaN*(N/1000) + log(BP) ) # rescaling N
  phi2 <- phiA
  phi1 <- phi2
  alpha <- 1.95/7.1 # from EE
  
  A <- matrix(c(0, 0, F, F, 
                (1-alpha)*phi1, 0, 0, 0, 
                alpha*phi1, 0, 0, 0, 
                0, phi2, phiA, phiA), ncol = 4, byrow = T)
  
  return( A )
}
T = 1000
#YKD
#low = ~ 8200; beta = 0.025
n=matrix(0, 4,T)
n[,1] <- c(1, 1, 1, 1)
for(i in 2:T){
  n[,i] = A(N=sum(n[c(3,4),i-1]), betaN=0.005)%*%n[,i-1]
}
plot(1:T, colSums(n[3:4,]))
# Best: N = 19300; beta = 0.01
n=matrix(0, 4,T)
n[,1] <- c(1, 1, 1, 1)
for(i in 2:T){
  n[,i] = A(N=sum(n[c(3,4),i-1]), betaN=0.01)%*%n[,i-1]
}
plot(1:T, colSums(n[3:4,]))
# High: N = 38000; beta = 0.005
n=matrix(0, 4,T)
n[,1] <- c(1, 1, 1, 1)
for(i in 2:T){
  n[,i] = A(N=sum(n[c(3,4),i-1]), betaN=0.005)%*%n[,i-1]
}
plot(1:T, colSums(n[3:4,]))
# Now explore the pop size over the whole prior
nreps = 1000
#beta <- rgamma(nreps, 3.109273, 251.648023) #from above
#beta <- (1/exp(rnorm(nreps, 7.95, 3.31)))/1000 #from original SSA
beta <- (1/rgamma(nreps, 100000, 1))/1000 # from IMP_YKD on GitHub, was this used in 2021?
hist(beta, breaks = 100)
pop <- c()
for( j in 1:nreps){
  n <- matrix(c(1, 1, 1, 1), 4, 1)
  for(i in 2:T){
    n = A(N=sum(n[c(3,4),1]), betaN=beta[j])%*%n
  }
  pop[j] <- sum(n[c(3,4),1])
}

hist(pop)
hist(log(pop))
plot(beta, pop)

