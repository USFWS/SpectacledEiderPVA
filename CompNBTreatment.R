library(jagsUI)

# looking at the effects of changing the theta used in the NB draws from supplied as data (base)
# to supplied as prior (alt) (~N) using the estimate and SE from the NB regression of observed
# data against model estimates for RCP4.5 and RCP8.5 (see NegBinReg.R code)

# comparison using Lconstant scenario for lead, RCP4.5 scenario for ice.

##### Data
# years of predictions, to 2060
K <- 41

# Count Data, 1988 - 2019
counts <- read.csv("input_data/YKD_SPEI.csv", header = T)
counts$Nibb[counts$Year==2015] <- NA
counts$seNibb[counts$Year==2015] <- NA
observer <- c(1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
              5, 5, 5, 5, 5)

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

# Survival Covariates, 1988 - 2100
# lead exposure rates calculated in the model
#lead.exposure <- read.csv("input_data/lead.exposure.csv", header = T)
#lead.exposure <- lead.exposure[1:113,]
#lead.exposure$delta.const <- c(diff(lead.exposure$const.exp), NA)
#lead.exposure$delta.2008 <- c(diff(lead.exposure$exp.2008), NA)
#lead.exposure$delta.2026 <- c(diff(lead.exposure$exp.2026), NA)
ext.sea.ice <- read.csv("input_data/extreme.sea.ice.csv", header = T)
# year is winter year (Nov - Apr)
ice.data <- ext.sea.ice[ext.sea.ice$year >= 1988 & ext.sea.ice$year <= 2100,]
ice.data4.5 <- ifelse(is.na(ice.data$ice.obs), ice.data$ice.RCP4.5, 
                      ice.data$ice.obs)
ice.data8.5 <- ifelse(is.na(ice.data$ice.obs), ice.data$ice.RCP8.5, 
                      ice.data$ice.obs)
ice.datacurrent <- ifelse(is.na(ice.data$ice.obs), mean(ice.data$ice.obs, na.rm = T), 
                          ice.data$ice.obs)

# Fecundity Data, 1992 - 2015
fecund.param <- read.csv("input_data/fecundity.csv", header = T)
# convert mean and se of NS estimates to alpha and beta (Beta distribution) for
# input as covariates to the detection model (mark-recap)
fecund.param$ns.obs.alpha <- 
  round((fecund.param$ns.obs.gs-fecund.param$ns.obs.gs^2-fecund.param$sigma.ns.obs^2)*
          fecund.param$ns.obs.gs/fecund.param$sigma.ns.obs^2, 0)
fecund.param$ns.obs.beta <- round(fecund.param$ns.obs.alpha*(1-fecund.param$ns.obs.gs)
                                  /fecund.param$ns.obs.gs, 0)

#### RCP4.5 base

### Lead Constant, RCP4.5
cat(file = "YKD_IPM.jags", "
model {

## Priors
# state-space
sigma.obs[24] ~ dunif(800, 1400) # bounds consistent with s.e. observed since 2005
sigma.obs[28] ~ dunif(800, 1400)

# adding noise to the ice predictions s.t. interannual variation in the 
# future is consistent with the observed time series (1980 - 2018)
for (i in 1:(K+1)){
  nb.prob[i] <- nb.size/(obs.ice[n.occasions + BEFORE + AFTER - 1 + i] + nb.size)
  ice.draw[i] ~ dnegbin(nb.prob[i], nb.size)
} # i
    
ice.new <- c(obs.ice[1:(n.occasions - 1 + BEFORE + AFTER)], ice.draw[])
ice <- (ice.new - mean(ice.new))/sd(ice.new)

# lead exposure and decay rate
theta_0 ~ dbeta(6, 45) # mean of 0.12, with 90% between 0.5 and 0.2
decay ~ dnorm(10, 0.1) T(1,) # assumes highly likely decay rate is between 5-15 years

### COMMENT OUT LEAD SCENARIOS DEPENDING ON RUN (DECLINE IN 2008, DECLINE IN
### 2026, OR CONSTANT LEAD)

# constant lead
for (i in 1: (BEFORE+n.occasions-1+AFTER+K)){
  lead[i] <- theta_0
}

# decline in lead beginning in 2008 (21st year)
#for (i in 1:20){
#  lead[i] <- theta_0
#} # i
#for (i in 21:(BEFORE+n.occasions-1+AFTER+K)){
#  lead[i] <- 0.5^(1/decay)*lead[i-1]
#} # i

# decline in lead beginning in 2026 (39th year)
#for (i in 1:38){
#  lead[i] <- theta_0
#} # i
#for (i in 39:(BEFORE+n.occasions-1+AFTER+K)){
#  lead[i] <- 0.5^(1/decay)*lead[i-1]
#} # i

# productivity
mean.log.F ~ dnorm(-0.47,100) 
mean.logit.BP ~ dnorm(1.73, 1000) # from EE, 1.73 in RCP4.5, 1.54 in RCP8.5
betaF.ice ~ dnorm(2.57, 1000) # from EE, 2.57 in RCP4.5, 2.62 in RCP8.5
betaF.ice2 ~ dnorm(0.12, 1000) # from EE, 0.12 in RCP4.5, 0.06 in RCP8.5
betaF.ice3 ~ dnorm(-0.6, 1000) # from EE, -0.6 in RCP4.5, -0.66 in RCP8.5
tau.BP <- pow(sigma.BP, -2)
sigma.BP ~ dunif(0,1)
for (i in 1:2){
betaN.inv[i] ~ dgamma(100000, 1) # from EE
betaN[i] <- 1/betaN.inv[i] # DD for phi and F
} #i

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

omegaA ~ dbeta(1, 35) # from EE
omegaJ ~ dbeta(1, 9) # from EE

mean.p ~ dbeta(1,1)
mean.logit.p <- logit(mean.p)
tau.p <- pow(sigma.p, -2)
sigma.p ~ dunif(0, 1)
    
for (i in 1:4){
  beta[i] ~ dnorm(0, 100)
} # i

beta.nest ~ dnorm(0, 100)

kappa ~ dbeta(20, 20) #relate phiA to phi2 by lead rate

## State-space model for count data
# Model for the initial population size
N[1,1] ~ dunif(100, 300) # n1, 1YO 
N[2,1] ~ dunif(100, 350) # n2nB, 2YO non-breeding
N[3,1] ~ dunif(100, 250) # n2B, 2YO breeding
N[4,1] <- N[5,1] + N[6,1] # n3+, 3+YO breeding
N[5,1] ~ dunif(100, 350) # n3+ from n2nB
N[6,1] ~ dunif(100, 2500) # n3+ from n2B and n3+

# Process model
for (t in 1:(BEFORE+n.occasions-1+AFTER+K)){ # extended loop here
    N[1,t+1] <- F[t] * (N[3,t] + N[4,t]) 
    N[2,t+1] <- (phi1[t] + (1-phi1[t])*omegaJ)*(1-alpha[t])*N[1,t]
    N[3,t+1] <- (phi1[t] + (1-phi1[t])*omegaJ)*alpha[t]*N[1,t]
    N[4,t+1] <- N[5,t+1] + N[6,t+1]
    N[5,t+1] <- (phi2[t] + (1-phi2[t])*omegaA)*N[2,t]
    N[6,t+1] <- (phiA[t] + (1-phiA[t])*omegaA)*(N[3,t] + N[4,t])
} # t

# Observation model
tau.o <- pow(sigma.o, -2)
sigma.o ~ dgamma(1, 10)
tau.d <- pow(sigma.d, -2)
sigma.d ~ dgamma(2, 13) # from EE
for(i in 1:5){
  o[i] ~ dnorm(0, tau.o)
} # i
for (t in 1:(n.occasions+BEFORE+AFTER)){
  d[t] ~ dnorm(0, tau.d)
  count[t] ~ dnorm(2*exp(log(N[3,t] + N[4,t]) - o[obs[t]] - d[t]), tau.obs[t]) # count provided as data
  tau.obs[t] <- pow(sigma.obs[t], -2) # sigma.obs provided as data
} # t

## Fecundity Model
# Process model
for (t in 1:(n.occasions+BEFORE+AFTER+K-1)){ # extended loop here
  logit.BP[t] <- mean.logit.BP + betaF.ice*ice[t] + betaF.ice2*ice[t]*ice[t]
  + betaF.ice3*ice[t]*ice[t]*ice[t] + eps.BP[t]
  eps.BP[t] ~ dnorm(0, tau.BP)
  BP[t] <- ilogit(logit.BP[t])
  F[t] <- exp(log(phi0[t] + (1-phi0[t])*omegaJ) + mean.log.F 
          - betaN[1]*(N[3,t] + N[4,t]) + log(BP[t]))
} # t

## Multistate survival model
# process model
for (t in 1:(n.occasions-1)){
nest[t] ~ dbeta(nest.obs.a[t], nest.obs.b[t])
}

nest.s <- (nest - mean(nest))/sd(nest)

for (t in (1+BEFORE):(BEFORE+n.occasions-1)){
   logit.p[t] <- mean.logit.p + beta.nest*nest.s[t-BEFORE] + eps.p[t]
   eps.p[t] ~ dnorm(0, tau.p)
   p[t] <- ilogit(logit.p[t])
} # t

for (t in 1:(n.occasions - 1 + BEFORE + AFTER + K)){
    logit.phi0[t] <- mean.logit.phi0 + beta[3]*ice[t+1] 
                     + beta[4]*ice[t+1]*ice[t+1] + eps.phi0[t] 
    eps.phi0[t] ~ dnorm(0, tau.phi0)
    phi0[t] <- ilogit(logit.phi0[t])
    
    logit.phiA[t] <- logit(mean.phiA*(1-lead[t]*(1-kappa))) + beta[1]*ice[t+1] 
                     + beta[2]*ice[t+1]*ice[t+1] 
                     - betaN[2]*(N[3,t] + N[4,t]) + eps.phiA[t] 
    logit.phi2[t] <- logit(mean.phiA) + beta[1]*ice[t+1] 
                     + beta[2]*ice[t+1]*ice[t+1] 
                     - betaN[2]*(N[3,t] + N[4,t]) + eps.phiA[t]
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
  psi[1,t,2] <- phi0[t+BEFORE]
  psi[1,t,3] <- 0
  psi[1,t,4] <- 0
  psi[1,t,5] <- 0
  psi[2,t,1] <- 0
  psi[2,t,2] <- 0
  psi[2,t,3] <- phi1[t+BEFORE]*(1-alpha[t+BEFORE])
  psi[2,t,4] <- phi1[t+BEFORE]*alpha[t+BEFORE]
  psi[2,t,5] <- 0
  psi[3,t,1] <- 0
  psi[3,t,2] <- 0
  psi[3,t,3] <- 0
  psi[3,t,4] <- 0
  psi[3,t,5] <- phi2[t+BEFORE]
  psi[4,t,1] <- 0
  psi[4,t,2] <- 0
  psi[4,t,3] <- 0
  psi[4,t,4] <- 0
  psi[4,t,5] <- phiA[t+BEFORE]
  psi[5,t,1] <- 0
  psi[5,t,2] <- 0
  psi[5,t,3] <- 0
  psi[5,t,4] <- 0
  psi[5,t,5] <- phiA[t+BEFORE]
  po[1,t] <- 0
  po[2,t] <- 0
  po[3,t] <- 0
  po[4,t] <- p[t+BEFORE]
  po[5,t] <- p[t+BEFORE]
  
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

## Derived parameters
# Total population size and breeding population size
for (t in 1:(n.occasions + BEFORE + AFTER + K)){ # extended loop here
   Ntot[t] <- 2*(N[1,t] + N[2,t] + N[3,t] + N[4,t])
   Nb[t] <- 2*(N[3,t] + N[4,t])
}

# Check whether the population goes extinct in future
for (t in 1:K){ # extended loop here
   extinct[t] <- 1- step(Nb[n.occasions + BEFORE + AFTER + t] - 251) # quasiextinction of 250 
}
}
")


# bundle data
## nb.size for RCP4.5 and RCP8.5 = 2.8
jags.data <- list(obs.ice = ice.data4.5, nb.size = 2.8,
                  marr = ms.arr, n.occasions = ncol(ch), rel = rowSums(ms.arr), 
                  ns = ns, zero = matrix(0, ncol = ns, nrow = ns), ones = diag(ns), 
                  count = counts$Nibb, obs = observer, sigma.obs = counts$seNibb,
                  nest.obs.a = fecund.param$ns.obs.alpha[1:23],
                  nest.obs.b = fecund.param$ns.obs.beta[1:23],
                  K = K, BEFORE = 4, AFTER = 4)

# initial values
inits <- function(){list(
  mean.phi0 = runif(1, 0.2, 0.3), 
  mean.phiA = runif(1, 0.8, 0.9), mean.p = runif(1, 0.5, 0.6))}

# parameters monitored
parameters <- c("Nb", "phiA", "phi0", "mean.log.F", "betaN", "beta", "o", "sigma.d", 
                "sigma.o", "betaF.ice", "betaF.ice2", "betaF.ice3", "mean.logit.BP")

# MCMC settings
ni <- 100000; nt <- 1; nb <- 10000; nc <- 3

# Call JAGS from R (jagsUI), use autojags to run to convergence
NBbase.RCP4.5 <- jags(jags.data, inits, parameters, "YKD_IPM.jags", 
                         n.chains = nc, n.burnin=nb, n.iter = ni,  
                         parallel = TRUE, n.adapt = 1000)

saveRDS(NBbase.RCP4.5, file = "output/NBbase.RCP4.5.rds")


#### RCP4.5 alt

### Lead Constant, RCP4.5
cat(file = "YKD_IPM.jags", "
model {

## Priors
# state-space
sigma.obs[24] ~ dunif(800, 1400) # bounds consistent with s.e. observed since 2005
sigma.obs[28] ~ dunif(800, 1400)

# adding noise to the ice predictions s.t. interannual variation in the 
# future is consistent with the observed time series (1980 - 2018)
nb.size ~ dnorm(2.8, 2.2)
for (i in 1:(K+1)){
  nb.prob[i] <- nb.size/(obs.ice[n.occasions + BEFORE + AFTER - 1 + i] + nb.size)
  ice.draw[i] ~ dnegbin(nb.prob[i], nb.size)
} # i
    
ice.new <- c(obs.ice[1:(n.occasions - 1 + BEFORE + AFTER)], ice.draw[])
ice <- (ice.new - mean(ice.new))/sd(ice.new)

# lead exposure and decay rate
theta_0 ~ dbeta(6, 45) # mean of 0.12, with 90% between 0.5 and 0.2
decay ~ dnorm(10, 0.1) T(1,) # assumes highly likely decay rate is between 5-15 years

### COMMENT OUT LEAD SCENARIOS DEPENDING ON RUN (DECLINE IN 2008, DECLINE IN
### 2026, OR CONSTANT LEAD)

# constant lead
for (i in 1: (BEFORE+n.occasions-1+AFTER+K)){
  lead[i] <- theta_0
}

# decline in lead beginning in 2008 (21st year)
#for (i in 1:20){
#  lead[i] <- theta_0
#} # i
#for (i in 21:(BEFORE+n.occasions-1+AFTER+K)){
#  lead[i] <- 0.5^(1/decay)*lead[i-1]
#} # i

# decline in lead beginning in 2026 (39th year)
#for (i in 1:38){
#  lead[i] <- theta_0
#} # i
#for (i in 39:(BEFORE+n.occasions-1+AFTER+K)){
#  lead[i] <- 0.5^(1/decay)*lead[i-1]
#} # i

# productivity
mean.log.F ~ dnorm(-0.47,100) 
mean.logit.BP ~ dnorm(1.73, 1000) # from EE, 1.73 in RCP4.5, 1.54 in RCP8.5
betaF.ice ~ dnorm(2.57, 1000) # from EE, 2.57 in RCP4.5, 2.62 in RCP8.5
betaF.ice2 ~ dnorm(0.12, 1000) # from EE, 0.12 in RCP4.5, 0.06 in RCP8.5
betaF.ice3 ~ dnorm(-0.6, 1000) # from EE, -0.6 in RCP4.5, -0.66 in RCP8.5
tau.BP <- pow(sigma.BP, -2)
sigma.BP ~ dunif(0,1)
for (i in 1:2){
betaN.inv[i] ~ dgamma(100000, 1) # from EE
betaN[i] <- 1/betaN.inv[i] # DD for phi and F
} #i

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

omegaA ~ dbeta(1, 35) # from EE
omegaJ ~ dbeta(1, 9) # from EE

mean.p ~ dbeta(1,1)
mean.logit.p <- logit(mean.p)
tau.p <- pow(sigma.p, -2)
sigma.p ~ dunif(0, 1)
    
for (i in 1:4){
  beta[i] ~ dnorm(0, 100)
} # i

beta.nest ~ dnorm(0, 100)

kappa ~ dbeta(20, 20) #relate phiA to phi2 by lead rate

## State-space model for count data
# Model for the initial population size
N[1,1] ~ dunif(100, 300) # n1, 1YO 
N[2,1] ~ dunif(100, 350) # n2nB, 2YO non-breeding
N[3,1] ~ dunif(100, 250) # n2B, 2YO breeding
N[4,1] <- N[5,1] + N[6,1] # n3+, 3+YO breeding
N[5,1] ~ dunif(100, 350) # n3+ from n2nB
N[6,1] ~ dunif(100, 2500) # n3+ from n2B and n3+

# Process model
for (t in 1:(BEFORE+n.occasions-1+AFTER+K)){ # extended loop here
    N[1,t+1] <- F[t] * (N[3,t] + N[4,t]) 
    N[2,t+1] <- (phi1[t] + (1-phi1[t])*omegaJ)*(1-alpha[t])*N[1,t]
    N[3,t+1] <- (phi1[t] + (1-phi1[t])*omegaJ)*alpha[t]*N[1,t]
    N[4,t+1] <- N[5,t+1] + N[6,t+1]
    N[5,t+1] <- (phi2[t] + (1-phi2[t])*omegaA)*N[2,t]
    N[6,t+1] <- (phiA[t] + (1-phiA[t])*omegaA)*(N[3,t] + N[4,t])
} # t

# Observation model
tau.o <- pow(sigma.o, -2)
sigma.o ~ dgamma(1, 10)
tau.d <- pow(sigma.d, -2)
sigma.d ~ dgamma(2, 13) # from EE
for(i in 1:5){
  o[i] ~ dnorm(0, tau.o)
} # i
for (t in 1:(n.occasions+BEFORE+AFTER)){
  d[t] ~ dnorm(0, tau.d)
  count[t] ~ dnorm(2*exp(log(N[3,t] + N[4,t]) - o[obs[t]] - d[t]), tau.obs[t]) # count provided as data
  tau.obs[t] <- pow(sigma.obs[t], -2) # sigma.obs provided as data
} # t

## Fecundity Model
# Process model
for (t in 1:(n.occasions+BEFORE+AFTER+K-1)){ # extended loop here
  logit.BP[t] <- mean.logit.BP + betaF.ice*ice[t] + betaF.ice2*ice[t]*ice[t]
  + betaF.ice3*ice[t]*ice[t]*ice[t] + eps.BP[t]
  eps.BP[t] ~ dnorm(0, tau.BP)
  BP[t] <- ilogit(logit.BP[t])
  F[t] <- exp(log(phi0[t] + (1-phi0[t])*omegaJ) + mean.log.F 
          - betaN[1]*(N[3,t] + N[4,t]) + log(BP[t]))
} # t

## Multistate survival model
# process model
for (t in 1:(n.occasions-1)){
nest[t] ~ dbeta(nest.obs.a[t], nest.obs.b[t])
}

nest.s <- (nest - mean(nest))/sd(nest)

for (t in (1+BEFORE):(BEFORE+n.occasions-1)){
   logit.p[t] <- mean.logit.p + beta.nest*nest.s[t-BEFORE] + eps.p[t]
   eps.p[t] ~ dnorm(0, tau.p)
   p[t] <- ilogit(logit.p[t])
} # t

for (t in 1:(n.occasions - 1 + BEFORE + AFTER + K)){
    logit.phi0[t] <- mean.logit.phi0 + beta[3]*ice[t+1] 
                     + beta[4]*ice[t+1]*ice[t+1] + eps.phi0[t] 
    eps.phi0[t] ~ dnorm(0, tau.phi0)
    phi0[t] <- ilogit(logit.phi0[t])
    
    logit.phiA[t] <- logit(mean.phiA*(1-lead[t]*(1-kappa))) + beta[1]*ice[t+1] 
                     + beta[2]*ice[t+1]*ice[t+1] 
                     - betaN[2]*(N[3,t] + N[4,t]) + eps.phiA[t] 
    logit.phi2[t] <- logit(mean.phiA) + beta[1]*ice[t+1] 
                     + beta[2]*ice[t+1]*ice[t+1] 
                     - betaN[2]*(N[3,t] + N[4,t]) + eps.phiA[t]
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
  psi[1,t,2] <- phi0[t+BEFORE]
  psi[1,t,3] <- 0
  psi[1,t,4] <- 0
  psi[1,t,5] <- 0
  psi[2,t,1] <- 0
  psi[2,t,2] <- 0
  psi[2,t,3] <- phi1[t+BEFORE]*(1-alpha[t+BEFORE])
  psi[2,t,4] <- phi1[t+BEFORE]*alpha[t+BEFORE]
  psi[2,t,5] <- 0
  psi[3,t,1] <- 0
  psi[3,t,2] <- 0
  psi[3,t,3] <- 0
  psi[3,t,4] <- 0
  psi[3,t,5] <- phi2[t+BEFORE]
  psi[4,t,1] <- 0
  psi[4,t,2] <- 0
  psi[4,t,3] <- 0
  psi[4,t,4] <- 0
  psi[4,t,5] <- phiA[t+BEFORE]
  psi[5,t,1] <- 0
  psi[5,t,2] <- 0
  psi[5,t,3] <- 0
  psi[5,t,4] <- 0
  psi[5,t,5] <- phiA[t+BEFORE]
  po[1,t] <- 0
  po[2,t] <- 0
  po[3,t] <- 0
  po[4,t] <- p[t+BEFORE]
  po[5,t] <- p[t+BEFORE]
  
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

## Derived parameters
# Total population size and breeding population size
for (t in 1:(n.occasions + BEFORE + AFTER + K)){ # extended loop here
   Ntot[t] <- 2*(N[1,t] + N[2,t] + N[3,t] + N[4,t])
   Nb[t] <- 2*(N[3,t] + N[4,t])
}

# Check whether the population goes extinct in future
for (t in 1:K){ # extended loop here
   extinct[t] <- 1- step(Nb[n.occasions + BEFORE + AFTER + t] - 251) # quasiextinction of 250 
}
}
")


# bundle data
jags.data <- list(obs.ice = ice.data4.5,
                  marr = ms.arr, n.occasions = ncol(ch), rel = rowSums(ms.arr), 
                  ns = ns, zero = matrix(0, ncol = ns, nrow = ns), ones = diag(ns), 
                  count = counts$Nibb, obs = observer, sigma.obs = counts$seNibb,
                  nest.obs.a = fecund.param$ns.obs.alpha[1:23],
                  nest.obs.b = fecund.param$ns.obs.beta[1:23],
                  K = K, BEFORE = 4, AFTER = 4)

# initial values
inits <- function(){list(
  mean.phi0 = runif(1, 0.2, 0.3), 
  mean.phiA = runif(1, 0.8, 0.9), mean.p = runif(1, 0.5, 0.6))}

# parameters monitored
parameters <- c("Nb", "phiA", "phi0", "beta", "betaF.ice", "betaF.ice2", 
                "betaF.ice3", "mean.logit.BP", "nb.size")

# MCMC settings
ni <- 100000; nt <- 1; nb <- 10000; nc <- 3

# Call JAGS from R (jagsUI), use autojags to run to convergence
NBalt.RCP4.5 <- jags(jags.data, inits, parameters, "YKD_IPM.jags", 
                      n.chains = nc, n.burnin=nb, n.iter = ni,  
                      parallel = TRUE, n.adapt = 1000)

saveRDS(NBalt.RCP4.5, file = "output/NBalt.RCP4.5.rds")

NBbase.Nb.quant <- apply(NBbase.RCP4.5$sims.list$Nb, 2, quantile, 
                         probs=c(0.025,0.25,0.5,0.75,0.975))

NBalt.Nb.quant <- apply(NBalt.RCP4.5$sims.list$Nb, 2, quantile, 
                        probs=c(0.025,0.25,0.5,0.75,0.975))

plot(x=c(1988:2060), y=c(counts$Nibb, rep(NA,K)), type="p", pch=21, bg='black', 
     main="", xlab="Year", ylab="Breeding Birds", ylim=c(0,max(counts$Nibb, NBalt.Nb.quant, na.rm = TRUE)))
grid(col="lightgray")
polygon(x=c(1988:2060,rev(1988:2060)), # 95% CI
        y=c(NBalt.Nb.quant[1,],rev(NBalt.Nb.quant[5,])), col=rgb(1,0,0, alpha=0.25), border=FALSE)
polygon(x=c(1988:2060,rev(1988:2060)), #50% CI
        y=c(NBalt.Nb.quant[2,],rev(NBalt.Nb.quant[4,])), col=rgb(1,0,0, alpha=0.25), border=FALSE)
lines(x=c(1988:2060), y=NBalt.Nb.quant[3,], col=rgb(1,0,0, alpha=0.25), lwd=3)
polygon(x=c(1988:2060,rev(1988:2060)), # 95% CI
        y=c(NBbase.Nb.quant[1,],rev(NBbase.Nb.quant[5,])), col=rgb(0,0,1, alpha=0.25), border=FALSE)
polygon(x=c(1988:2060,rev(1988:2060)), #50% CI
        y=c(NBbase.Nb.quant[2,],rev(NBbase.Nb.quant[4,])), col=rgb(0,0,1, alpha=0.25), border=FALSE)
lines(x=c(1988:2060), y=NBbase.Nb.quant[3,], col=rgb(0,0,1, alpha=0.25), lwd=3)
