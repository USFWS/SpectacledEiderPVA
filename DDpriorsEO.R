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
#high females = 37.5K
optimize(f=Btest, interval = c(0, 0.1), tol=1e-8, N=37500)$minimum
# [1] 0.005429917
#best  = 20K females
optimize(f=Btest, interval = c(0, 0.1), tol=1e-8, N=20000)$minimum
# [1] 0.0101811
#low = 8.5K  females
optimize(f=Btest, interval = c(0, 0.1), tol=1e-8, N=8500)$minimum
# [1] 0.02395552


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
             eeq=c(0.005429917, 0.0101811, 0.02395552))
fit
# $par
# [1]   2.623867 203.060233
# 
# $value
# [1] 265.7278
# 
# $counts
# function gradient 
# 53       53 
# 
# $convergence
# [1] 0
# 
# $message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
qgamma(c(0.1, 0.5, 0.9), shape=fit$par[1], rate=fit$par[2])
# [1] 0.004317903 0.011322344 0.023611622
hist(rgamma(10000, fit$par[1], fit$par[2]))
