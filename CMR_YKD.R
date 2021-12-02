###This is the CMR component pulled out of IPM_YKD.R in order to add a GAM to the ice functional response

library(jagsUI)

# Data
# Count Data, 1988 - 2019
counts <- read.csv("input_data/YKD_SPEI.csv", header = T)
counts <- counts$Nibb[counts$Year>=1992 & counts$Year<2016]
counts[20] <- (counts[19]+counts[21])/2 #imput missing 2011 count
observer <- c(1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
              5, 5, 5, 5, 5)
# Fecundity Data, 1992 - 2015
fecund.param <- read.csv("input_data/fecundity.csv", header = T)
# convert mean and se of NS estimates to alpha and beta (Beta distribution) for
# input as covariates to the detection model (mark-recap)
fecund.param$ns.obs.alpha <- 
  round((fecund.param$ns.obs.gs-fecund.param$ns.obs.gs^2-fecund.param$sigma.ns.obs^2)*
          fecund.param$ns.obs.gs/fecund.param$sigma.ns.obs^2, 0)
fecund.param$ns.obs.beta <- round(fecund.param$ns.obs.alpha*(1-fecund.param$ns.obs.gs)
                                  /fecund.param$ns.obs.gs, 0)

# Mark-Recap Data, 1992 - 2015
mcr.kig <- read.csv("input_data/mcr.kig.csv", header = T)
# for multistate cap-recap analysis, observations recoded to reflect states
ch <- as.matrix(mcr.kig[, -c(1,26,27)]) # remove column with band number and columns indicating age (juve/adult)
#####
# MS Array Function
marray <- function(ch, unobs = 2){
  ns <- length(table(ch)) - 1 + unobs
  no <- ncol(ch)
  out <- matrix(0, ncol = ns*(no-1)+1, nrow = ns*(no-1))
  # Remove capture histories of individuals that are marked at last occasion
  get.first <- function(x) min(which(x!=0))
  first <- apply(ch, 1, get.first)
  last.only <- which(first==no)
  if (length(last.only) > 0) ch <- ch[-last.only,]
  # Compute m-array
  for (i in 1:nrow(ch)){
    cap.occ <- which(ch[i,]!=0)
    state <- ch[i,cap.occ]
    if (length(state) == 1) {
      out[state[1]+ns*(cap.occ[1]-1), ns*(no-1)+1] <- out[state[1]+ns*(cap.occ[1]-1), ns*(no-1)+1] + 1
    }
    if (length(state) > 1) {
      for (t in 2:length(cap.occ)){
        out[(cap.occ[t-1]-1)*ns+state[t-1], (cap.occ[t]-2)*ns+state[t]] <- out[(cap.occ[t-1]-1)*ns+state[t-1], (cap.occ[t]-2)*ns+state[t]] + 1
      } # t
      if (max(cap.occ) < no){
        out[(cap.occ[t]-1)*ns+state[t], ns*(no-1)+1] <- out[(cap.occ[t]-1)*ns+state[t], ns*(no-1)+1] + 1
      } # if
    } # if
  } # t
  return(out)
}    
ms.arr <- marray(ch)
# number of states for the multistate mark-recapture analysis
ns <- 5

# Survival Covariates, 1992 - 2016

ext.sea.ice <- read.csv("input_data/extreme.sea.ice.csv", header = T)
# year is winter year (Nov - Apr)
ice.data <- ext.sea.ice[ext.sea.ice$year >= 1992 & ext.sea.ice$year <= 2016,]
ice.data <- ice.data$ice.obs
ice.data <- (ice.data-mean(ice.data)) / sd(ice.data)

### Lead Constant
cat(file = "CMR.jags", "
model {

## Priors

# survival and breeding propensity
mean.alpha.inv ~ dgamma(7.1, 1.95) T(1,) # from EE
mean.logit.alpha <- logit(1/mean.alpha.inv)

mean.phi0 ~ dbeta(15,45) # from cjs model
mean.logit.phi0 <- logit(mean.phi0)
tau.phi0 <- pow(sigma.phi0, -2)
sigma.phi0 ~ dunif(0, 1)

mean.phiA ~ dbeta(90,10) # from cjs model
tau.phiA <- pow(sigma.phiA, -2)
sigma.phiA ~ dunif(0, 1)

mean.p ~ dbeta(1,1)
mean.logit.p <- logit(mean.p)
tau.p <- pow(sigma.p, -2)
sigma.p ~ dunif(0, 1)
    
for (i in 1:4){
  beta[i] ~ dnorm(0, 100)
} # i
betaN.inv ~ dgamma(100000, 1) # from EE
betaN <- 1/betaN.inv # DD for phi and F

beta.nest ~ dnorm(0, 100)

## Multistate survival model
# process model
for (t in 1:(n.occasions-1)){
nest[t] ~ dbeta(nest.obs.a[t], nest.obs.b[t])
}

nest.s <- (nest - mean(nest))/sd(nest)

for (t in (1):(n.occasions-1)){
   logit.p[t] <- mean.logit.p + beta.nest*nest.s[t] + eps.p[t]
   eps.p[t] ~ dnorm(0, tau.p)
   p[t] <- ilogit(logit.p[t])
} # t

for (t in 1:(n.occasions - 1)){
    logit.phi0[t] <- mean.logit.phi0 + beta[3]*ice[t+1] 
                     + beta[4]*ice[t+1]*ice[t+1] + eps.phi0[t] 
    eps.phi0[t] ~ dnorm(0, tau.phi0)
    phi0[t] <- ilogit(logit.phi0[t])
    
    logit.phiA[t] <- logit(mean.phiA) + beta[1]*ice[t+1] 
                     + beta[2]*ice[t+1]*ice[t+1] 
                     - betaN*(N[t]) + eps.phiA[t] 
    logit.phi2[t] <- logit(mean.phiA) + beta[1]*ice[t+1] 
                     + beta[2]*ice[t+1]*ice[t+1] 
                     - betaN*(N[t]) + eps.phiA[t]
    eps.phiA[t] ~ dnorm(0, tau.phiA)
    phiA[t] <- ilogit(logit.phiA[t])
    phi2[t] <- ilogit(logit.phi2[t])
    
    phi1[t] <- phi2[t]
    
    logit.alpha[t] <- mean.logit.alpha
    alpha[t] <- ilogit(logit.alpha[t])
} # t

for (t in 1:(n.occasions - 1)){
  # state transition and reencounter probabilities
  psi[1,t,1] <- 0
  psi[1,t,2] <- phi0[t]
  psi[1,t,3] <- 0
  psi[1,t,4] <- 0
  psi[1,t,5] <- 0
  psi[2,t,1] <- 0
  psi[2,t,2] <- 0
  psi[2,t,3] <- phi1[t]*(1-alpha[t])
  psi[2,t,4] <- phi1[t]*alpha[t]
  psi[2,t,5] <- 0
  psi[3,t,1] <- 0
  psi[3,t,2] <- 0
  psi[3,t,3] <- 0
  psi[3,t,4] <- 0
  psi[3,t,5] <- phi2[t]
  psi[4,t,1] <- 0
  psi[4,t,2] <- 0
  psi[4,t,3] <- 0
  psi[4,t,4] <- 0
  psi[4,t,5] <- phiA[t]
  psi[5,t,1] <- 0
  psi[5,t,2] <- 0
  psi[5,t,3] <- 0
  psi[5,t,4] <- 0
  psi[5,t,5] <- phiA[t]
  po[1,t] <- 0
  po[2,t] <- 0
  po[3,t] <- 0
  po[4,t] <- p[t]
  po[5,t] <- p[t]
  
# non-encounter probabilities, dq, and reshape the array for the 
  # encounter probabilities
  for (s in 1:ns){
    dp[s,t,s] <- po[s,t]
    dq[s,t,s] <- 1-po[s,t]
  } # s
  for (s in 1:(ns-1)){
    for (m in (s+1):ns){
      dp[s,t,m] <- 0
      dq[s,t,m] <- 0
    } # m
  } # s
  for (s in 2:ns){
    for (m in 1:(s-1)){
      dp[s,t,m] <- 0
      dq[s,t,m] <- 0
    } # m
  } # s
} # t

# multinomial likelihood
for (t in 1:((n.occasions-1)*ns)){
   marr[t,1:(n.occasions*ns-(ns-1))] ~ dmulti(pr[t,], rel[t])
   } # t

# Define the cell probabilities of the multistate m-array   
# Define matrix U: product of probabilities of state-transition and non-encounter (this is just done because there is no product function for matrix multiplication in JAGS)
for (t in 1:(n.occasions-2)){
   U[(t-1)*ns+(1:ns), (t-1)*ns+(1:ns)] <- ones
   for (j in (t+1):(n.occasions-1)){
      U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% psi[,t,] %*% dq[,t,]
      } # j
   } # t
U[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- ones
# Diagonal
for (t in 1:(n.occasions-2)){
   pr[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] %*% psi[,t,] %*% dp[,t,]
   # Above main diagonal
   for (j in (t+1):(n.occasions-1)){
      pr[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] %*% psi[,j,] %*% dp[,j,]
      } # j
   } # t
pr[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- psi[,n.occasions-1,] %*% dp[,n.occasions-1,]

# Below main diagonal
for (t in 2:(n.occasions-1)){
   for (j in 1:(t-1)){
      pr[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero
      } #j
   } #t

# Last column: probability of non-recapture
for (t in 1:((n.occasions-1)*ns)){
   pr[t,(n.occasions*ns-(ns-1))] <- 1-sum(pr[t,1:((n.occasions-1)*ns)])
   } #t
}
")


# bundle data
jags.data <- list(ice = ice.data, 
                  marr = ms.arr, n.occasions = ncol(ch), rel = rowSums(ms.arr), 
                  ns = ns, zero = matrix(0, ncol = ns, nrow = ns), ones = diag(ns), 
                  N = counts,
                  nest.obs.a = fecund.param$ns.obs.alpha[1:23],
                  nest.obs.b = fecund.param$ns.obs.beta[1:23])

# initial values
inits <- function(){list(
  mean.phi0 = runif(1, 0.2, 0.3), 
  mean.phiA = runif(1, 0.8, 0.9), mean.p = runif(1, 0.5, 0.6))}

# parameters monitored
parameters <- c("phiA", "phi0", "betaN", "beta", "beta.nest", 
                "mean.logit.alpha", "mean.phiA", "mean.phi0", "sigma.phi0", "sigma.phiA")

# MCMC settings
ni <- 11000; nt <- 1; nb <- 5000; nc <- 3

# Call JAGS from R (jagsUI), use autojags to run to convergence
out <- jags(jags.data, inits, parameters, "CMR.jags", 
                         n.chains = nc, n.burnin=nb, n.iter = ni,  
                         parallel = TRUE, n.adapt = 1000)

#saveRDS(YKD.constant.4.5, file = "YKD.constant.4.5.CMR.GAM.rds")
