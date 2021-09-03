ilogit <- function(x){exp(x)/(1+exp(x))}

T <- 10000
N <- matrix(0, nrow=6, ncol=T)

N[1,1] <- runif(1, 100, 300) # n1, 1YO 
N[2,1] <- runif(1, 100, 350) # n2nB, 2YO non-breeding
N[3,1] <- runif(1, 100, 250) # n2B, 2YO breeding
N[5,1] <- runif(1, 100, 350) # n3+ from n2B and n3+
N[6,1] <- runif(1, 100, 2500) # n3+ from n2nB
N[4,1] <- N[5,1] + N[6,1] # n3+, 3+YO breeding


#omegaA <- 1/(1+35) # from EE
#omegaJ <- 1/(1+9)
  
BP <- ilogit(1.73) # from EE
mean.phi0 <- 0.36
mean.logit.phi0 <- log(mean.phi0/(1-mean.phi0))
mean.phiA <- 0.9
mean.logit.phiA <- log(mean.phiA/(1-mean.phiA))
mean.log.F <- -.2

#betaN <- c(0.00033, 0.00033) #high K = 50657 ibb = 50657/2 females = 25328.5
#[1] 25348.15
#betaN <- c(0.00090, 0.00090) #best K = 25000 ibb = 12500 females
#[1] 12532
betaN <- c(0.00162, 0.00162) #low K = 10743 ibb = 5371.5 females
#[1] 5362.5

for(t in 1:(T-1)){
  logit.phi0 <- mean.logit.phi0 - betaN[2]*(N[3,t] + N[4,t]) 
  phi0 <- ilogit(logit.phi0)
  logit.phiA <- mean.logit.phiA - betaN[2]*(N[3,t] + N[4,t])
  phiA <- ilogit(logit.phiA)
  F <- exp(log(phi0 + (1-phi0)) + mean.log.F - betaN[1]*(N[3,t] + N[4,t]) 
              + log(BP) )
  phi2 <- phiA
  phi1 <- phi2
  alpha <- 1.95/7.1 # from EE
  
  N[1,t+1] <- F * (N[3,t] + N[4,t]) 
  N[2,t+1] <- (phi1 + (1-phi1))*(1-alpha)*N[1,t]
  N[3,t+1] <- (phi1 + (1-phi1))*alpha*N[1,t]
  N[5,t+1] <- (phi2 + (1-phi2))*N[2,t]
  N[6,t+1] <- (phiA + (1-phiA))*(N[3,t] + N[4,t])
  N[4,t+1] <- N[5,t+1] + N[6,t+1]
  
}

plot(1:T, colSums(N[3:4,]))
sum(N[3:4,T])
#betaN <- c(0.00044, 0.00044) #high K = 50657 ibb = 50657/2 females = 25328.5
#[1] 25348.15
#betaN <- c(0.00089, 0.00089) #best K = 25000 ibb = 12500 females
#[1] 12532
#betaN <- c(0.00208, 0.00208) #low K = 10743 ibb = 5371.5 females
#[1] 5362.5
x <- rgamma(100000, 3, 3300) #use rgamma(10000, 3, 3300)
hist(x)
mean(x)
sd(x)
quantile(x, probs=c(0.1, 0.5, 0.9))

log.invDDfem <- rnorm(10000, 7.2, 0.7)
invDDfem <- exp(log.invDDfem)
DDfem <- 1/invDDfem
mean(DDfem)
round(quantile(DDfem, c(0.1, 0.5, 0.9)), 6)

