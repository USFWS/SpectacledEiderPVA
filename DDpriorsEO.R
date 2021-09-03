##Find equilibrium with optim

Atest = function(N=10000, betaN=0.0001){
  ilogit <- function(x){exp(x)/(1+exp(x))}
  BP <- ilogit(1.73) # from EE
  mean.phi0 <- 0.36
  mean.logit.phi0 <- log(mean.phi0/(1-mean.phi0))
  mean.phiA <- 0.9
  mean.logit.phiA <- log(mean.phiA/(1-mean.phiA))
  mean.log.F <- -0.2
  phi0 <- ilogit(mean.logit.phi0)
  logit.phiA <- mean.logit.phiA - betaN*N
  phiA <- ilogit(logit.phiA)
  F <- exp(log(phi0) + mean.log.F - betaN*N + log(BP) )
  phi2 <- phiA
  phi1 <- phi2
  alpha <- 1.95/7.1 # from EE
  
  A <- matrix(c(0, 0, F, F, 
                (1-alpha)*phi1, 0, 0, 0, 
                alpha*phi1, 0, 0, 0, 
                0, phi2, phiA, phiA), ncol = 4, byrow = T)
  
  return( (max(Re(eigen(A)$values))-1)^2 )
}
Btest = function(betaN=0.0001, N=12500){
  ilogit <- function(x){exp(x)/(1+exp(x))}
  BP <- ilogit(1.73) # from EE
  mean.phi0 <- 0.36
  mean.logit.phi0 <- log(mean.phi0/(1-mean.phi0))
  mean.phiA <- 0.9
  mean.logit.phiA <- log(mean.phiA/(1-mean.phiA))
  mean.log.F <- -0.2
  phi0 <- ilogit(mean.logit.phi0)
  logit.phiA <- mean.logit.phiA - betaN*N
  phiA <- ilogit(logit.phiA)
  F <- exp(log(phi0) + mean.log.F - betaN*N + log(BP) )
  phi2 <- phiA
  phi1 <- phi2
  alpha <- 1.95/7.1 # from EE
  
  A <- matrix(c(0, 0, F, F, 
                (1-alpha)*phi1, 0, 0, 0, 
                alpha*phi1, 0, 0, 0, 
                0, phi2, phiA, phiA), ncol = 4, byrow = T)
  
  return( (max(Re(eigen(A)$values))-1)^2 )
}

optimize(f=Atest, interval = c(0, 1e6), betaN = 0.0001)

optimize(f=Btest, interval = c(0, 0.001), N=12500)

##relate this to the EE for the ACP
#high K = 50657 ibb = 50657/2 females = 25328.5
optimize(f=Btest, interval = c(0, 0.001), tol=1e-8, N=25329)$minimum
#[1] 1.39755e-05
#best K = 25000 ibb = 12500 females
optimize(f=Btest, interval = c(0, 0.001), tol=1e-8, N=12500)$minimum
#[1] 2.831876e-05
#low K = 10743 ibb = 5371.5 females
optimize(f=Btest, interval = c(0, 0.001), tol=1e-8, N=5372)$minimum
#[1] 6.589424e-05

##relate this to the EE for the YKD
#high females = 37.5K
optimize(f=Btest, interval = c(0, 0.001), tol=1e-8, N=37500)$minimum
#[1] 9.439188e-06
#best  = 20K females
optimize(f=Btest, interval = c(0, 0.001), tol=1e-8, N=20000)$minimum
#[1] 1.769929e-05
#low = 8.5K  females
optimize(f=Btest, interval = c(0, 0.001), tol=1e-8, N=8500)$minimum
#[1] 4.164514e-05


#now find a distribution that matches the EE
qmatch.gamma <- function(par=c(3,10000), eeq=c(0.1, 0.7, 2.3), qtarget=c(0.1, 0.5, 0.9)){
  qdist <- qgamma(qtarget, shape=par[1], rate=par[2])
  return( sum((10000*(qdist-eeq))^2) )
}

#ACP
fit <- optim(par=c(3, 10000), fn=qmatch.gamma, method="L-BFGS-B", lower=c(1,1), 
             eeq=c(0.00001398, 0.00002832, 0.00006589))
fit
# $par
# [1]     2.53054 71564.04780
# 
# $value
# [1] 0.001301901
# 
# $counts
# function gradient 
# 34       34 
# 
# $convergence
# [1] 0
# 
# $message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
qgamma(c(0.1, 0.5, 0.9), shape=fit$par[1], rate=fit$par[2])
# [1] 1.149584e-05 3.082759e-05 6.514170e-05
#YKD
fit <- optim(par=c(3, 10000), fn=qmatch.gamma, method="L-BFGS-B", lower=c(1,1), 
             eeq=c(0.000009439, 0.00001770, 0.00004165))
fit
# $par
# [1] 2.623284e+00 1.167722e+05
# 
# $value
# [1] 0.0008031142
# 
# $counts
# function gradient 
# 39       39 
# 
# $convergence
# [1] 0
# 
# $message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
qgamma(c(0.1, 0.5, 0.9), shape=fit$par[1], rate=fit$par[2])
# [1] 7.505676e-06 1.968395e-05 4.105220e-05