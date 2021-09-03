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
  
BP <- 0.9 
mean.phi0 <- 0.25
mean.logit.phi0 <- log(mean.phi0/(1-mean.phi0))
mean.phiA <- 0.9
mean.logit.phiA <- log(mean.phiA/(1-mean.phiA))
mean.log.F <- -.47

betaN <- c(0.00055, 0.00055)

for(t in 1:(T-1)){
  logit.phi0 <- mean.logit.phi0 - betaN[2]*(N[3,t] + N[4,t]) 
  phi0 <- ilogit(logit.phi0)
  logit.phiA <- mean.logit.phiA - betaN[2]*(N[3,t] + N[4,t])
  phiA <- ilogit(logit.phiA)
  F <- exp(log(phi0 + (1-phi0)) + mean.log.F - betaN[1]*(N[3,t] + N[4,t]) 
              + log(BP) )
  phi2 <- phiA
  phi1 <- phi2
  alpha <- 0.3
  
  N[1,t+1] <- F * (N[3,t] + N[4,t]) 
  N[2,t+1] <- (phi1 + (1-phi1))*(1-alpha)*N[1,t]
  N[3,t+1] <- (phi1 + (1-phi1))*alpha*N[1,t]
  N[5,t+1] <- (phi2 + (1-phi2))*N[2,t]
  N[6,t+1] <- (phiA + (1-phiA))*(N[3,t] + N[4,t])
  N[4,t+1] <- N[5,t+1] + N[6,t+1]
  
}

plot(1:T, colSums(N[1:4,]))
sum(N[1:4, T])

# 0.00016, ~ 40000, inv = 6250
# 0.00038, ~ 17000, inv = 2631
# 0.000085, ~ 75000, inv = 11764
log.invDD <- rnorm(10000, 8.7, 0.55)
invDD <- exp(log.invDD)
DD <- 1/invDD
round(quantile(DD, c(0.1, 0.5, 0.9)), 6)

# if we need half, because this is only females (THINK THIS IS TRUE)
# then 0.00032 ~ 20K
# 0.00076 ~ 8.5K
# 0.00017 ~ 37.5K
log.invDDfem <- rnorm(10000, 7.95, 0.55)
invDDfem <- exp(log.invDDfem)
DDfem <- 1/invDDfem
round(quantile(DDfem, c(0.1, 0.5, 0.9)), 6)



