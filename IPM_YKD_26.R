# This file is based on IPM_MS.R, which was based on IPM_YKD.R, see below.
# Updates from IPM_MS.R are meant to inform the 2026 SSA
#   An important update is the CJS code, where a bug was found in the indexing and 
#   model was revised to allow for missing resight and banding years, see files under 
#   CMR_component folder and especially the CMR_simulation sub-folder. 


# Origin: Originally based on IPM_YKD.R, updating model code for manuscript (IPM_MS) 
# following work/discussions in winter 2021-22

# 1) Considering only 4 scenarios (lead 2008, RCP 4.5; lead constant, RCP 4.5;
# lead 2008, RCP 8.5; lead constant, RCP 8.5)

# 2) revising shape parameter in nb pulls to 4 (see nb_regress.R) s.t. the average
# variance of a simulated draw from model predicted ice from 1980 - 2019 is ~=
# the variance in the observed ice from 1980 - 2019.

# 3) revising EE priors for ice effects on breeding propensity,
# density dependence, and annual variation in aerial detection. See EE write up for
# all other EE priors. 1.27.2022 

# 4) update input data, ice (from Dan R.) and counts (from Chuck F.) 
################################################################################
# 2026 updates: (1) ice data, (2) CMR data and code for missing years; 
#               (3) count data, fecundity data, and observer data to 2026;
#               (4) revise neg. binomial size parameter for ice projections to 2.1,
#                   see code below. cannot reproduce Cat's original value of 4 
#                   (used in IPM_MS, see above). For each ice scenario, used a specific
#                   theta parameter: 4.5 -> 3.46, 8.5 -> 3.40, current (mean) -> 3.33. 
#                   Originally, 3.5 was used in the IPM for the SSA, IPM_YKD),
#               (5) simple ice effect to make shared between ducklings and >1 year old,
#                   cause correlation in survival between adults and ducklings, intercepts are not shared;
#               (6) used informed prior for alpha and p from CJS model because simulations show it is identifiable
#               (7) revised prior for VCF variation (EE question # 9), d in below, see ObsPriors.R
library(jagsUI)
library(tidyverse)
library(AKaerial)
################################################################################
#Define functions
#### Plot results function
plot_results <- function(out = out){
  #plot adults
  df <- data.frame(PhiA = apply(out$sims.list$phiA, 2, mean), 
                   upper = apply(out$sims.list$phiA, 2, quantile, probs = 0.975),
                   lower = apply(out$sims.list$phiA, 2, quantile, probs = 0.025), 
                   Year = 1992:2024)
  p1 <- ggplot(data = df, aes(x = Year, y = PhiA)) + 
    geom_pointrange(aes(ymin = lower, ymax = upper)) +
    scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1))
  # Not good!
  #plot ducklings
  df0 <- data.frame(Phi0 = apply(out$sims.list$phi0, 2, mean), 
                    upper = apply(out$sims.list$phi0, 2, quantile, probs = 0.975),
                    lower = apply(out$sims.list$phi0, 2, quantile, probs = 0.025), 
                    Year = 1992:2024)
  p2 <- ggplot(data = df0, aes(x = Year, y = Phi0)) + 
    geom_pointrange(aes(ymin = lower, ymax = upper)) +
    scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1))
  # similarly bad! 
  #plot p
  dfp <- data.frame(p = apply(out$sims.list$p, 2, mean), 
                    upper = apply(out$sims.list$p, 2, quantile, probs = 0.975),
                    lower = apply(out$sims.list$p, 2, quantile, probs = 0.025), 
                    Year = 1993:2025)
  p3 <- ggplot(data = dfp, aes(x = Year, y = p)) + 
    geom_pointrange(aes(ymin = lower, ymax = upper)) +
    scale_x_continuous(breaks = seq(1993, 2025, by = 4)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1))
  return(list(PhiA=p1, Phi0=p2, p=p3))
}
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
################################################################################
#Define JAGS model
# NOTE!!! COMMENTS REQUIRED IN MODEL AND INPUT DATA TO CHANGE FOR EACH SCENARIO
# RUN

### Lead 2008
cat(file = "YKD_IPM_L2008.jags", "
model {

## Priors
# state-space
sigma.obs[24] ~ dunif(800, 1400) # 2011, bounds consistent with s.e. observed since 2005
sigma.obs[28] ~ dunif(800, 1400) # 2015 removed in 2019 SSA because we have strong 
                                 # prior that population count was low due to observation process
                                 # for 2026 update, leave it in?
sigma.obs[33] ~ dunif(800, 1400) #2020 is missing
# adding noise to the ice predictions s.t. interannual variation in the 
# future is consistent with the observed time series (1980 - 2018)
for (i in 1:(K+1)){
  nb.prob[i] <- nb.size/(obs.ice[n.occasions + BEFORE + AFTER - 1 + i] + nb.size)
  ice.draw[i] ~ dnegbin(nb.prob[i], nb.size)
} # i
    
ice.new <- c(obs.ice[1:(n.occasions - 1 + BEFORE + AFTER)], ice.draw[])
ice.BP <- (ice.new - 60)/25 #standardize to mean and sd of the expert elicitation
ice <- (ice.new - 57.6)/25.3 #standardize to the ice used in the stand alone CJS model for priors

# lead exposure and decay rate
theta_0 ~ dbeta(6, 45) # mean of 0.12, with 90% between 0.5 and 0.2
decay ~ dnorm(10, 0.1) T(1,) # assumes highly likely decay rate is between 5-15 years

### COMMENT OUT LEAD SCENARIOS DEPENDING ON RUN (DECLINE IN 2008, DECLINE IN
### 2026, OR CONSTANT LEAD)

# constant lead
#for (i in 1: (BEFORE+n.occasions-1+AFTER+K)){
#  lead[i] <- theta_0
#}

# decline in lead beginning in 2008 (21st year)
for (i in 1:20){
  lead[i] <- theta_0
} # i
for (i in 21:(BEFORE+n.occasions-1+AFTER+K)){
  lead[i] <- 0.5^(1/decay)*lead[i-1]
} # i

# decline in lead beginning in 2026 (39th year)
#for (i in 1:38){
#  lead[i] <- theta_0
#} # i
#for (i in 39:(BEFORE+n.occasions-1+AFTER+K)){
#  lead[i] <- 0.5^(1/decay)*lead[i-1]
#} # i

# productivity
mean.log.F ~ dnorm(-0.47,100) 
betaBP[1:3] ~ dmnorm.vcov(muBP[1:3], SigmaBP[1:3, 1:3])
tau.BP <- pow(sigma.BP, -2)
sigma.BP ~ dunif(0,1)
for (i in 1:2){
betaN[i] ~ dgamma(3.1, 251.648) # from EE, see DDpriorsEO.R
} #i

# survival and breeding propensity
# mean.alpha.inv ~ dgamma(7.1, 1.95) T(1,) # from EE
# mean.logit.alpha <- logit(1/mean.alpha.inv) #from CJS model, this is identifiable
mean.alpha ~ dbeta(13 , 21) #from CJS model

mean.phi0 ~ dbeta(20, 60) # from cjs model 
mean.logit.phi0 <- logit(mean.phi0)
tau.phi0 <- pow(sigma.phi0, -2)
sigma.phi0 ~ dunif(0, 1)

mean.phiA ~ dbeta(84,16) # from cjs model
tau.phiA <- pow(sigma.phiA, -2)
sigma.phiA ~ dunif(0, 1)

omegaA ~ dbeta(1, 35) # from EE
omegaJ ~ dbeta(1, 9) # from EE

mean.p ~ dbeta(47, 53) #from CJS model 
mean.logit.p <- logit(mean.p)
tau.p <- pow(sigma.p, -2)
sigma.p ~ dunif(0, 1) #should consider this to vary more
#prior for covariates    
for (i in 1:2){
  beta[i] ~ dnorm(0, 100)
} # i

beta.nest ~ dnorm(0, 100) #prior for nest success effect on p

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
for(i in 1:6){
  o[i] ~ dnorm(0, tau.o)
} # i
for (t in 1:(n.occasions+BEFORE+AFTER)){
  x[t] ~ dbeta(4, 4) # from EE, see ObsPriors.R
  d[t] <- log( 1 + (x[t] * (0.8 - 0.2) + 0.2 - 0.5)/0.5 ) #standarize to (0.2, 0.8) and make multiplicative on log scale 
  count[t] ~ dnorm(2*exp(log(N[3,t] + N[4,t]) - o[obs[t]] - d[t]), tau.obs[t]) # count provided as data
  tau.obs[t] <- pow(sigma.obs[t], -2) # sigma.obs provided as data
} # t

## Fecundity Model
# Process model
for (t in 1:(n.occasions+BEFORE+AFTER+K-1)){ # extended loop here
  logit.BP[t] <- betaBP[1] + betaBP[2]*ice.BP[t] + betaBP[3]*ice.BP[t]*ice.BP[t]+ eps.BP[t]
  eps.BP[t] ~ dnorm(0, tau.BP)
  BP[t] <- ilogit(logit.BP[t])
  F[t] <- exp(log(phi0[t] + (1-phi0[t])*omegaJ) + mean.log.F 
          - betaN[1]*(N[3,t] + N[4,t])/1000 + log(BP[t]))
} # t

## Multistate survival model
# process model
for (t in 1:(n.occasions-1)){
nest[t] ~ dbeta(nest.obs.a[t], nest.obs.b[t]) #nest sucess observation from priors
}

nest.s <- (nest - mean(nest))/sd(nest)

for (t in (1+BEFORE):(BEFORE+n.occasions-1)){
   logit.p[t] <- mean.logit.p + beta.nest*nest.s[t-BEFORE] + eps.p[t]
   eps.p[t] ~ dnorm(0, tau.p)
   p[t] <- ilogit(logit.p[t])
} # t

for (t in 1:(n.occasions - 1 + BEFORE + AFTER + K)){ 
    logit.phi0[t] <- mean.logit.phi0 + beta[1]*ice[t+1]         #assume ice effect is shared with adults
                     + beta[2]*ice[t+1]*ice[t+1] + eps.phi0[t]  #t + 1 because 1987 is in ice data for fecundity
    eps.phi0[t] ~ dnorm(0, tau.phi0)
    phi0[t] <- ilogit(logit.phi0[t])
    
    logit.phiA[t] <- logit(mean.phiA*(1-lead[t]*(1-kappa))) + beta[1]*ice[t+1] 
                     + beta[2]*ice[t+1]*ice[t+1] 
                     - betaN[2]*(N[3,t] + N[4,t])/1000 + eps.phiA[t]
    logit.phi2[t] <- logit(mean.phiA) + beta[1]*ice[t+1] 
                     + beta[2]*ice[t+1]*ice[t+1] 
                     - betaN[2]*(N[3,t] + N[4,t])/1000 + eps.phiA[t] #assume phi2 = phiA
    eps.phiA[t] ~ dnorm(0, tau.phiA)
    phiA[t] <- ilogit(logit.phiA[t])
    phi2[t] <- ilogit(logit.phi2[t])
    # 
    phi1[t] <- phi2[t] #assume phi2 = phi1
    
    #logit.alpha[t] <- mean.logit.alpha
    #alpha[t] <- ilogit(logit.alpha[t])
    alpha[t] <- mean.alpha
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
  po[4,t] <- p[t+BEFORE] * resight.index[t + BEFORE]
  po[5,t] <- p[t+BEFORE] * resight.index[t + BEFORE]
  
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
      U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% psi[,j-1,] %*% dq[,j-1,]
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
}
")
################################################################################
### Lead 2026
cat(file = "YKD_IPM_L2026.jags", "
model {

## Priors
# state-space
sigma.obs[24] ~ dunif(800, 1400) # 2011, bounds consistent with s.e. observed since 2005
sigma.obs[28] ~ dunif(800, 1400) # 2015 removed in 2019 SSA because we have strong 
                                 # prior population count was low due to observation process
                                 # for 2026 update, leave it in?
sigma.obs[33] ~ dunif(800, 1400) #2020 is missing
# adding noise to the ice predictions s.t. interannual variation in the 
# future is consistent with the observed time series (1980 - 2018)
for (i in 1:(K+1)){
  nb.prob[i] <- nb.size/(obs.ice[n.occasions + BEFORE + AFTER - 1 + i] + nb.size)
  ice.draw[i] ~ dnegbin(nb.prob[i], nb.size)
} # i
    
ice.new <- c(obs.ice[1:(n.occasions - 1 + BEFORE + AFTER)], ice.draw[])
ice.BP <- (ice.new - 60)/25 #standardize to mean and sd of the expert elicitation
ice <- (ice.new - 57.6)/25.3 #standardize to the ice used in the stand salone CJS model for priors

# lead exposure and decay rate
theta_0 ~ dbeta(6, 45) # mean of 0.12, with 90% between 0.5 and 0.2
decay ~ dnorm(10, 0.1) T(1,) # assumes highly likely decay rate is between 5-15 years

### COMMENT OUT LEAD SCENARIOS DEPENDING ON RUN (DECLINE IN 2008, DECLINE IN
### 2026, OR CONSTANT LEAD)

# constant lead
#for (i in 1: (BEFORE+n.occasions-1+AFTER+K)){
#  lead[i] <- theta_0
#}

# # decline in lead beginning in 2008 (21st year)
# for (i in 1:20){
#   lead[i] <- theta_0
# } # i
# for (i in 21:(BEFORE+n.occasions-1+AFTER+K)){
#   lead[i] <- 0.5^(1/decay)*lead[i-1]
# } # i

# decline in lead beginning in 2026 (39th year)
for (i in 1:38){
 lead[i] <- theta_0
} # i
for (i in 39:(BEFORE+n.occasions-1+AFTER+K)){
 lead[i] <- 0.5^(1/decay)*lead[i-1]
} # i

# productivity
mean.log.F ~ dnorm(-0.47,100) 
betaBP[1:3] ~ dmnorm.vcov(muBP[1:3], SigmaBP[1:3, 1:3])
tau.BP <- pow(sigma.BP, -2)
sigma.BP ~ dunif(0,1)
for (i in 1:2){
betaN[i] ~ dgamma(3.1, 251.648) # from EE, see DDpriorsEO.R
} #i

# survival and breeding propensity
# mean.alpha.inv ~ dgamma(7.1, 1.95) T(1,) # from EE
# mean.logit.alpha <- logit(1/mean.alpha.inv) #from CJS model, this is identifiable
mean.alpha ~ dbeta(13 , 21) #from CJS model

mean.phi0 ~ dbeta(20, 60) # from cjs model 
mean.logit.phi0 <- logit(mean.phi0)
tau.phi0 <- pow(sigma.phi0, -2)
sigma.phi0 ~ dunif(0, 1)

mean.phiA ~ dbeta(84,16) # from cjs model
tau.phiA <- pow(sigma.phiA, -2)
sigma.phiA ~ dunif(0, 1)

omegaA ~ dbeta(1, 35) # from EE
omegaJ ~ dbeta(1, 9) # from EE

mean.p ~ dbeta(47, 53) #from CJS model 
mean.logit.p <- logit(mean.p)
tau.p <- pow(sigma.p, -2)
sigma.p ~ dunif(0, 1) #should consider this to vary more
#prior for covariates    
for (i in 1:2){
  beta[i] ~ dnorm(0, 100)
} # i

beta.nest ~ dnorm(0, 100) #prior for nest success effect on p

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
for(i in 1:6){
  o[i] ~ dnorm(0, tau.o)
} # i
for (t in 1:(n.occasions+BEFORE+AFTER)){
  x[t] ~ dbeta(4, 4) # from EE, see ObsPriors.R
  d[t] <- log( 1 + (x[t] * (0.8 - 0.2) + 0.2 - 0.5)/0.5 ) #standarize to (0.2, 0.8) and make multiplicative on log scale 
  count[t] ~ dnorm(2*exp(log(N[3,t] + N[4,t]) - o[obs[t]] - d[t]), tau.obs[t]) # count provided as data
  tau.obs[t] <- pow(sigma.obs[t], -2) # sigma.obs provided as data
} # t

## Fecundity Model
# Process model
for (t in 1:(n.occasions+BEFORE+AFTER+K-1)){ # extended loop here
  logit.BP[t] <- betaBP[1] + betaBP[2]*ice.BP[t] + betaBP[3]*ice.BP[t]*ice.BP[t]+ eps.BP[t]
  eps.BP[t] ~ dnorm(0, tau.BP)
  BP[t] <- ilogit(logit.BP[t])
  F[t] <- exp(log(phi0[t] + (1-phi0[t])*omegaJ) + mean.log.F 
          - betaN[1]*(N[3,t] + N[4,t])/1000 + log(BP[t]))
} # t

## Multistate survival model
# process model
for (t in 1:(n.occasions-1)){
nest[t] ~ dbeta(nest.obs.a[t], nest.obs.b[t]) #nest sucess observation from priors
}

nest.s <- (nest - mean(nest))/sd(nest)

for (t in (1+BEFORE):(BEFORE+n.occasions-1)){
   logit.p[t] <- mean.logit.p + beta.nest*nest.s[t-BEFORE] + eps.p[t]
   eps.p[t] ~ dnorm(0, tau.p)
   p[t] <- ilogit(logit.p[t])
} # t

for (t in 1:(n.occasions - 1 + BEFORE + AFTER + K)){ 
    logit.phi0[t] <- mean.logit.phi0 + beta[1]*ice[t+1]         #assume ice effect is shared with adults
                     + beta[2]*ice[t+1]*ice[t+1] + eps.phi0[t]  #t + 1 because 1987 is in ice data for fecundity
    eps.phi0[t] ~ dnorm(0, tau.phi0)
    phi0[t] <- ilogit(logit.phi0[t])
    
    logit.phiA[t] <- logit(mean.phiA*(1-lead[t]*(1-kappa))) + beta[1]*ice[t+1] 
                     + beta[2]*ice[t+1]*ice[t+1] 
                     - betaN[2]*(N[3,t] + N[4,t])/1000 + eps.phiA[t]
    logit.phi2[t] <- logit(mean.phiA) + beta[1]*ice[t+1] 
                     + beta[2]*ice[t+1]*ice[t+1] 
                     - betaN[2]*(N[3,t] + N[4,t])/1000 + eps.phiA[t] #assume phi2 = phiA
    eps.phiA[t] ~ dnorm(0, tau.phiA)
    phiA[t] <- ilogit(logit.phiA[t])
    phi2[t] <- ilogit(logit.phi2[t])
    # 
    phi1[t] <- phi2[t] #assume phi2 = phi1
    
    #logit.alpha[t] <- mean.logit.alpha
    #alpha[t] <- ilogit(logit.alpha[t])
    alpha[t] <- mean.alpha
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
  po[4,t] <- p[t+BEFORE] * resight.index[t + BEFORE]
  po[5,t] <- p[t+BEFORE] * resight.index[t + BEFORE]
  
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
      U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% psi[,j-1,] %*% dq[,j-1,]
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
}
")
################################################################################
### Lead Constant
cat(file = "YKD_IPM_Lconst.jags", "
model {

## Priors
# state-space
sigma.obs[24] ~ dunif(800, 1400) # 2011, bounds consistent with s.e. observed since 2005
sigma.obs[28] ~ dunif(800, 1400) # 2015 removed in 2019 SSA because we have strong 
                                 # prior population count was low due to observation process
                                 # for 2026 update, leave it in?
sigma.obs[33] ~ dunif(800, 1400) #2020 is missing
# adding noise to the ice predictions s.t. interannual variation in the 
# future is consistent with the observed time series (1980 - 2018)
for (i in 1:(K+1)){
  nb.prob[i] <- nb.size/(obs.ice[n.occasions + BEFORE + AFTER - 1 + i] + nb.size)
  ice.draw[i] ~ dnegbin(nb.prob[i], nb.size)
} # i
    
ice.new <- c(obs.ice[1:(n.occasions - 1 + BEFORE + AFTER)], ice.draw[])
ice.BP <- (ice.new - 60)/25 #standardize to mean and sd of the expert elicitation
ice <- (ice.new - 57.6)/25.3 #standardize to the ice used in the stand salone CJS model for priors

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
# for (i in 1:20){
#   lead[i] <- theta_0
# } # i
# for (i in 21:(BEFORE+n.occasions-1+AFTER+K)){
#   lead[i] <- 0.5^(1/decay)*lead[i-1]
# } # i

# decline in lead beginning in 2026 (39th year)
#for (i in 1:38){
#  lead[i] <- theta_0
#} # i
#for (i in 39:(BEFORE+n.occasions-1+AFTER+K)){
#  lead[i] <- 0.5^(1/decay)*lead[i-1]
#} # i

# productivity
mean.log.F ~ dnorm(-0.47,100) 
betaBP[1:3] ~ dmnorm.vcov(muBP[1:3], SigmaBP[1:3, 1:3])
tau.BP <- pow(sigma.BP, -2)
sigma.BP ~ dunif(0,1)
for (i in 1:2){
betaN[i] ~ dgamma(3.1, 251.648) # from EE, see DDpriorsEO.R
} #i

# survival and breeding propensity
# mean.alpha.inv ~ dgamma(7.1, 1.95) T(1,) # from EE
# mean.logit.alpha <- logit(1/mean.alpha.inv) #from CJS model, this is identifiable
mean.alpha ~ dbeta(13 , 21) #from CJS model

mean.phi0 ~ dbeta(20, 60) # from cjs model 
mean.logit.phi0 <- logit(mean.phi0)
tau.phi0 <- pow(sigma.phi0, -2)
sigma.phi0 ~ dunif(0, 1)

mean.phiA ~ dbeta(84,16) # from cjs model
tau.phiA <- pow(sigma.phiA, -2)
sigma.phiA ~ dunif(0, 1)

omegaA ~ dbeta(1, 35) # from EE
omegaJ ~ dbeta(1, 9) # from EE

mean.p ~ dbeta(47, 53) #from CJS model 
mean.logit.p <- logit(mean.p)
tau.p <- pow(sigma.p, -2)
sigma.p ~ dunif(0, 1) #should consider this to vary more
#prior for covariates    
for (i in 1:2){
  beta[i] ~ dnorm(0, 100)
} # i

beta.nest ~ dnorm(0, 100) #prior for nest success effect on p

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
for(i in 1:6){
  o[i] ~ dnorm(0, tau.o)
} # i
for (t in 1:(n.occasions+BEFORE+AFTER)){
  x[t] ~ dbeta(4, 4) # from EE, see ObsPriors.R
  d[t] <- log( 1 + (x[t] * (0.8 - 0.2) + 0.2 - 0.5)/0.5 ) #standarize to (0.2, 0.8) and make multiplicative on log scale 
  count[t] ~ dnorm(2*exp(log(N[3,t] + N[4,t]) - o[obs[t]] - d[t]), tau.obs[t]) # count provided as data
  tau.obs[t] <- pow(sigma.obs[t], -2) # sigma.obs provided as data
} # t

## Fecundity Model
# Process model
for (t in 1:(n.occasions+BEFORE+AFTER+K-1)){ # extended loop here
  logit.BP[t] <- betaBP[1] + betaBP[2]*ice.BP[t] + betaBP[3]*ice.BP[t]*ice.BP[t]+ eps.BP[t]
  eps.BP[t] ~ dnorm(0, tau.BP)
  BP[t] <- ilogit(logit.BP[t])
  F[t] <- exp(log(phi0[t] + (1-phi0[t])*omegaJ) + mean.log.F 
          - betaN[1]*(N[3,t] + N[4,t])/1000 + log(BP[t]))
} # t

## Multistate survival model
# process model
for (t in 1:(n.occasions-1)){
nest[t] ~ dbeta(nest.obs.a[t], nest.obs.b[t]) #nest sucess observation from priors
}

nest.s <- (nest - mean(nest))/sd(nest)

for (t in (1+BEFORE):(BEFORE+n.occasions-1)){
   logit.p[t] <- mean.logit.p + beta.nest*nest.s[t-BEFORE] + eps.p[t]
   eps.p[t] ~ dnorm(0, tau.p)
   p[t] <- ilogit(logit.p[t])
} # t

for (t in 1:(n.occasions - 1 + BEFORE + AFTER + K)){ 
    logit.phi0[t] <- mean.logit.phi0 + beta[1]*ice[t+1]         #assume ice effect is shared with adults
                     + beta[2]*ice[t+1]*ice[t+1] + eps.phi0[t]  #t + 1 because 1987 is in ice data for fecundity
    eps.phi0[t] ~ dnorm(0, tau.phi0)
    phi0[t] <- ilogit(logit.phi0[t])
    
    logit.phiA[t] <- logit(mean.phiA*(1-lead[t]*(1-kappa))) + beta[1]*ice[t+1] 
                     + beta[2]*ice[t+1]*ice[t+1] 
                     - betaN[2]*(N[3,t] + N[4,t])/1000 + eps.phiA[t]
    logit.phi2[t] <- logit(mean.phiA) + beta[1]*ice[t+1] 
                     + beta[2]*ice[t+1]*ice[t+1] 
                     - betaN[2]*(N[3,t] + N[4,t])/1000 + eps.phiA[t] #assume phi2 = phiA
    eps.phiA[t] ~ dnorm(0, tau.phiA)
    phiA[t] <- ilogit(logit.phiA[t])
    phi2[t] <- ilogit(logit.phi2[t])
    # 
    phi1[t] <- phi2[t] #assume phi2 = phi1
    
    #logit.alpha[t] <- mean.logit.alpha
    #alpha[t] <- ilogit(logit.alpha[t])
    alpha[t] <- mean.alpha
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
  po[4,t] <- p[t+BEFORE] * resight.index[t + BEFORE]
  po[5,t] <- p[t+BEFORE] * resight.index[t + BEFORE]
  
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
      U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% psi[,j-1,] %*% dq[,j-1,]
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
}
")
################################################################################
### Current conditions: lead decline starts in 2008 to 2026, then constant; 
#     climate is observed mean plus random annual deviation
cat(file = "YKD_IPM_current.jags", "
model {

## Priors
# state-space
sigma.obs[24] ~ dunif(800, 1400) # 2011, bounds consistent with s.e. observed since 2005
sigma.obs[28] ~ dunif(800, 1400) # 2015 removed in 2019 SSA because we have strong 
                                 # prior population count was low due to observation process
                                 # for 2026 update, leave it in?
sigma.obs[33] ~ dunif(800, 1400) #2020 is missing
# adding noise to the ice predictions s.t. interannual variation in the 
# future is consistent with the observed time series (1980 - 2018)
for (i in 1:(K+1)){
  nb.prob[i] <- nb.size/(obs.ice[n.occasions + BEFORE + AFTER - 1 + i] + nb.size)
  ice.draw[i] ~ dnegbin(nb.prob[i], nb.size)
} # i
    
ice.new <- c(obs.ice[1:(n.occasions - 1 + BEFORE + AFTER)], ice.draw[])
ice.BP <- (ice.new - 60)/25 #standardize to mean and sd of the expert elicitation
ice <- (ice.new - 57.6)/25.3 #standardize to the ice used in the stand salone CJS model for priors

# lead exposure and decay rate
theta_0 ~ dbeta(6, 45) # mean of 0.12, with 90% between 0.5 and 0.2
decay ~ dnorm(10, 0.1) T(1,) # assumes highly likely decay rate is between 5-15 years

### COMMENT OUT LEAD SCENARIOS DEPENDING ON RUN (DECLINE IN 2008, DECLINE IN
### 2026, OR CONSTANT LEAD)

# constant lead
#for (i in 1: (BEFORE+n.occasions-1+AFTER+K)){
#  lead[i] <- theta_0
#}

# decline in lead beginning in 2008 (21st year)
for (i in 1:20){
  lead[i] <- theta_0
} # i
for (i in 21:39){
  lead[i] <- 0.5^(1/decay)*lead[i-1]
} # i

for (i in 40:(BEFORE+n.occasions-1+AFTER+K)){
 lead[i] <- lead[39]
} # i

# productivity
mean.log.F ~ dnorm(-0.47,100) 
betaBP[1:3] ~ dmnorm.vcov(muBP[1:3], SigmaBP[1:3, 1:3])
tau.BP <- pow(sigma.BP, -2)
sigma.BP ~ dunif(0,1)
for (i in 1:2){
betaN[i] ~ dgamma(3.1, 251.648) # from EE, see DDpriorsEO.R
} #i

# survival and breeding propensity
# mean.alpha.inv ~ dgamma(7.1, 1.95) T(1,) # from EE
# mean.logit.alpha <- logit(1/mean.alpha.inv) #from CJS model, this is identifiable
mean.alpha ~ dbeta(13 , 21) #from CJS model

mean.phi0 ~ dbeta(20, 60) # from cjs model 
mean.logit.phi0 <- logit(mean.phi0)
tau.phi0 <- pow(sigma.phi0, -2)
sigma.phi0 ~ dunif(0, 1)

mean.phiA ~ dbeta(84,16) # from cjs model
tau.phiA <- pow(sigma.phiA, -2)
sigma.phiA ~ dunif(0, 1)

omegaA ~ dbeta(1, 35) # from EE
omegaJ ~ dbeta(1, 9) # from EE

mean.p ~ dbeta(47, 53) #from CJS model 
mean.logit.p <- logit(mean.p)
tau.p <- pow(sigma.p, -2)
sigma.p ~ dunif(0, 1) #should consider this to vary more
#prior for covariates    
for (i in 1:2){
  beta[i] ~ dnorm(0, 100)
} # i

beta.nest ~ dnorm(0, 100) #prior for nest success effect on p

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
for(i in 1:6){
  o[i] ~ dnorm(0, tau.o)
} # i
for (t in 1:(n.occasions+BEFORE+AFTER)){
  x[t] ~ dbeta(4, 4) # from EE, see ObsPriors.R
  d[t] <- log( 1 + (x[t] * (0.8 - 0.2) + 0.2 - 0.5)/0.5 ) #standarize to (0.2, 0.8) and make multiplicative on log scale 
  count[t] ~ dnorm(2*exp(log(N[3,t] + N[4,t]) - o[obs[t]] - d[t]), tau.obs[t]) # count provided as data
  tau.obs[t] <- pow(sigma.obs[t], -2) # sigma.obs provided as data
} # t

## Fecundity Model
# Process model
for (t in 1:(n.occasions+BEFORE+AFTER+K-1)){ # extended loop here
  logit.BP[t] <- betaBP[1] + betaBP[2]*ice.BP[t] + betaBP[3]*ice.BP[t]*ice.BP[t]+ eps.BP[t]
  eps.BP[t] ~ dnorm(0, tau.BP)
  BP[t] <- ilogit(logit.BP[t])
  F[t] <- exp(log(phi0[t] + (1-phi0[t])*omegaJ) + mean.log.F 
          - betaN[1]*(N[3,t] + N[4,t])/1000 + log(BP[t]))
} # t

## Multistate survival model
# process model
for (t in 1:(n.occasions-1)){
nest[t] ~ dbeta(nest.obs.a[t], nest.obs.b[t]) #nest sucess observation from priors
}

nest.s <- (nest - mean(nest))/sd(nest)

for (t in (1+BEFORE):(BEFORE+n.occasions-1)){
   logit.p[t] <- mean.logit.p + beta.nest*nest.s[t-BEFORE] + eps.p[t]
   eps.p[t] ~ dnorm(0, tau.p)
   p[t] <- ilogit(logit.p[t])
} # t

for (t in 1:(n.occasions - 1 + BEFORE + AFTER + K)){ 
    logit.phi0[t] <- mean.logit.phi0 + beta[1]*ice[t+1]         #assume ice effect is shared with adults
                     + beta[2]*ice[t+1]*ice[t+1] + eps.phi0[t]  #t + 1 because 1987 is in ice data for fecundity
    eps.phi0[t] ~ dnorm(0, tau.phi0)
    phi0[t] <- ilogit(logit.phi0[t])
    
    logit.phiA[t] <- logit(mean.phiA*(1-lead[t]*(1-kappa))) + beta[1]*ice[t+1] 
                     + beta[2]*ice[t+1]*ice[t+1] 
                     - betaN[2]*(N[3,t] + N[4,t])/1000 + eps.phiA[t]
    logit.phi2[t] <- logit(mean.phiA) + beta[1]*ice[t+1] 
                     + beta[2]*ice[t+1]*ice[t+1] 
                     - betaN[2]*(N[3,t] + N[4,t])/1000 + eps.phiA[t] #assume phi2 = phiA
    eps.phiA[t] ~ dnorm(0, tau.phiA)
    phiA[t] <- ilogit(logit.phiA[t])
    phi2[t] <- ilogit(logit.phi2[t])
    # 
    phi1[t] <- phi2[t] #assume phi2 = phi1
    
    #logit.alpha[t] <- mean.logit.alpha
    #alpha[t] <- ilogit(logit.alpha[t])
    alpha[t] <- mean.alpha
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
  po[4,t] <- p[t+BEFORE] * resight.index[t + BEFORE]
  po[5,t] <- p[t+BEFORE] * resight.index[t + BEFORE]
  
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
      U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% psi[,j-1,] %*% dq[,j-1,]
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
}
")
################################################################################
# Set up data
# years of predictions, to 2100
K <- 74
# Count Data, 1988 - 2019
eiders <- AKaerial::YKDVEst()
eiders[[2]]
counts <- eiders[[1]]  |> 
  select(Year, Nibb, seNibb) |>
  complete(Year = full_seq(Year, 1)) |>
  mutate(Nibb = ifelse(row_number() == 28, NA, Nibb),     #remove 2015, see comment above in JAGS model
         seNibb = ifelse(row_number() == 28, NA, seNibb))
write.csv(counts, file = "input_data/YKD_counts_2026.csv")
observer <- zoo::na.locf(as.numeric(factor(AKaerial::YKDVHistoric$output.table$Observer))) #give missing years the previous observer
####################################
# Mark-Recap Data, 1992 - 2025
mcr.kig <- read.csv("input_data/eh_na_1992-2025_adults_ducklings.csv", header = T) |>
  arrange(desc(is_duckling)) |>
  pivot_longer(cols = starts_with("s"), names_to = "Occasion", 
               values_to = "State") |>
  group_by(id_metalBand) |>
  mutate(State = if_else(State != 0 & is_duckling != 1, 5, State)) |> #change adult to state 5
  mutate(temp = which(State != 0)[2] - which(State != 0)[1]) |> #find time between first and second observations
  mutate(State = if_else(is_duckling == 1 & 
                           !is.na(temp) & #at least 2 observations
                           row_number() == which(State != 0)[2] & 
                           temp == 2, 4, State)) |> #replace to second year breeder state
  mutate(State = if_else(is_duckling == 1 & 
                           !is.na(temp) &
                           temp == 2 & #time between first and second observation = 2
                           row_number() %in% which(State != 0) & 
                           row_number() > which(State != 0)[2], 5, State)) |> #rows after second obs non-zero states
  mutate(State = if_else(is_duckling == 1 & 
                           !is.na(temp) &
                           temp > 2 & #time between first and second observation > 2
                           row_number() %in% which(State != 0) & 
                           row_number() > which(State != 0)[1], 5, State)) #rows of non-zero states
# Seems to have worked based on filtering and examining data
# There has to be a better way!
# Now pivot back to wide form
mcr.kig <- select(mcr.kig, -is_duckling, -temp) |>
  pivot_wider(names_from = "Occasion", values_from = "State")

# drop band id
ch <- as.matrix(mcr.kig[, -1]) 
##### 
# make m-array
ms.arr <- marray(ch)
####################################
# Survival Covariates, 1988 - 2100
# lead exposure rates calculated in the model
# sea ice
sea.ice <- readxl::read_xlsx("input_data/sea_ice_vars_1979_2025.xlsx") |>
  select(year_winter, ws_count_sic_ge95) |>
  filter(year_winter >= 1987) |> #starting ice data in 1987 because fecundity 
                                 # depends on ice experienced in t - 1, 
                                 # makes survival depend on ice at t+1, 
                                 # see code above
  rename(year = year_winter, ice.obs = ws_count_sic_ge95)
  
  # year is winter year starts in Nov ends in Apr, unlike in original data, the year for survival
  #   and ice match up. I think!
  # ice data does not seem to match that used for 2021 SSA, need to ask
# ice.data <- read.csv("input_data/extreme.sea.ice.csv", header = T) |>
#   drop_na()
# ggplot(data=ice.data, aes(x=ice.RCP4.5, y = ice.obs)) + geom_point()
# ggplot(data=ice.data, aes(x=ice.RCP8.5, y = ice.obs)) + geom_point()
# ggplot(data=ice.data, aes(x=year, y = ice.obs)) + geom_point()
# ggplot(data=ice.data, aes(x=year, y = ice.RCP4.5)) + geom_point()
# ggplot(data=ice.data, aes(x=year, y = ice.RCP8.5)) + geom_point()

ice.data <- read.csv("input_data/extreme.sea.ice.csv", header = T) |>
  filter(year >= 1988) |>
  mutate(year = year -1)
  # year is winter year (Nov - Apr) I think this data frame year is +1 relative to sea.ice
  # so 1987 (year that winter starts) ice in sea.ice is 1988 (year that winter ends) in ice.data
#swap in new ice data 1988 to 2024
ice.data$ice.obs[ice.data$year %in% 1987:2024] <- sea.ice$ice.obs
ice.data <- ice.data |> 
  mutate(ice_4.5 = ifelse(is.na(ice.obs), ice.RCP4.5, 
                      ice.obs),
         ice_8.5 = ifelse(is.na(ice.obs), ice.RCP8.5, 
                      ice.obs),
         ice_current = ifelse(is.na(ice.obs), mean(ice.obs, na.rm = T), 
                          ice.obs))
#estimate NB theta parameter
library(mgcv)
fit <- gam(formula=ice.obs~ice.RCP4.5, data=ice.data, family = "nb")
summary(fit)
theta4.5 <- fit$family$getTheta(TRUE)
fit <- gam(formula=ice.obs~ice.RCP8.5, data=ice.data, family = "nb")
summary(fit)
theta8.5 <- fit$family$getTheta(TRUE)
fit <- gam(formula=ice.obs~1, data=ice.data, family = "nb")
summary(fit)
theta.current <- fit$family$getTheta(TRUE)
####################################
# Fecundity Data, 1992 - 2015
fecund.param <- read.csv("input_data/fecundity.csv", header = T)
# convert mean and se of NS estimates to alpha and beta (Beta distribution) for
# input as covariates to the detection model (mark-recap)
fecund.param$ns.obs.alpha <- 
  round((fecund.param$ns.obs.gs-fecund.param$ns.obs.gs^2-fecund.param$sigma.ns.obs^2)*
          fecund.param$ns.obs.gs/fecund.param$sigma.ns.obs^2, 0)
fecund.param$ns.obs.beta <- round(fecund.param$ns.obs.alpha*(1-fecund.param$ns.obs.gs)
                                  /fecund.param$ns.obs.gs, 0)
#Need to add 2016 to 2025 fecundity data
# Can get 2019 and 2021 from Friendly: https://doi.org/10.1093/ornithology/ukaf008
#  Use apparent nest success, number of successful and failed nest as beta distribution parameters. 
#  From supplemental mat. in Friendly: 2019 = 118 successful, 47 failed; 2021 = 85 succ., 17 failed
#using the prior Beta(39, 11) as placeholder for now, looking at fig. 3 in Friendly and reported mean of 0.78
ns.a <- c(39, 39, 39, 118, 39, 85, 39, 39, 39)
ns.b <- c(11, 11, 11, 47, 11, 17, 11, 11, 11)
#Need to discuss with Dan which nest success data to use? See 2015 estimate of 
#  Friendly vs what is in the data here (used in 2021)
####################################
# from BP_prior_EO_revised, mean and VCV for BP ice covariates
muBP <- c(1.3956, -0.1632, -0.1048)
SigmaBP <- matrix(c(0.3496, -0.0516, -0.0716, -0.0516, 0.0335, 0.0159, -0.0716,
                  0.0159, 0.0190), nrow = 3, ncol = 3)
################################################################################
################################################################################
# set up JAGS data
# initial values
inits <- function(){list(
  mean.phi0 = 0.25, 
  mean.phiA = runif(1, 0.8, 0.9), 
  mean.p = runif(1, 0.5, 0.6), 
  sigma.phi0 = 0.5, 
  sigma.phiA = 0.5, 
  sigma.p = 0.5)}
# ---- resight.index construction (BEFORE-driven, readable) --------------------
# Which calendar years had NO resighting effort:
missing_resight_years <- c(2016, 2017, 2018, 2020)

# Timeline constants (single source of truth):
BEFORE      <- 4            # CJS occasions offset relative to count timeline
first_band  <- 1992         # calendar year of banding column 1
n.occ       <- ncol(ch)     # number of banding occasions/columns (= n.occasions)

# Map a missing RESIGHT calendar year to its resight.index position:
#   banding column of that year : (year - first_band + 1)
#   m-array recapture occasion  : (banding column - 1)   # standard m-array shift:
#                                                         # first-recapture is offset one occasion
#   resight.index position      : (recapture occasion + BEFORE)
# Composing: position = (year - first_band + 1) - 1 + BEFORE = year - first_band + BEFORE
gap_positions <- missing_resight_years - first_band + BEFORE

# Length must match the model's use: positions 1..((n.occ - 1) + BEFORE)
resight.index <- rep(1, (n.occ - 1) + BEFORE)
resight.index[gap_positions] <- 0
# ---- runtime alignment / sanity checks (run once, not needed in production) --
# 1. Positions that got zeroed, and the calendar years they map back to.
#    Back-map: position P -> recapture occasion (P - BEFORE) -> banding col (+1)
#              -> calendar year (first_band + col - 1).
zeroed_pos   <- which(resight.index == 0)
zeroed_years <- first_band + ((zeroed_pos - BEFORE) + 1) - 1
cat("Zeroed resight.index positions:", zeroed_pos, "\n")
cat("...map back to calendar years :", zeroed_years, "\n")
stopifnot(setequal(zeroed_years, missing_resight_years))

# 2. Length matches what the JAGS model indexes: max position = (n.occ-1)+BEFORE
cat("length(resight.index):", length(resight.index),
    " expected:", (n.occ - 1) + BEFORE, "\n")
stopifnot(length(resight.index) == (n.occ - 1) + BEFORE)

# 3. m-array columns for the gap recapture occasions.
#    Each recapture occasion k occupies m-array columns ((k-1)*ns + 1):(k*ns).
#    With a true effort gap these should be all zero (no first-recaptures possible).
ns <- 5
gap_occ <- zeroed_pos - BEFORE
for (k in gap_occ) {
  cols <- ((k - 1) * ns + 1):(k * ns)
  cat(sprintf("recapture occ %d (cols %d-%d) colSums: %s\n",
              k, min(cols), max(cols),
              paste(colSums(ms.arr)[cols], collapse = ", ")))
}
# ------------------------------------------------------------------------------
# ---------- bundle data
jags.data8.5 <- list(obs.ice = ice.data$ice_8.5, nb.size = theta8.5, muBP = muBP, SigmaBP=SigmaBP,
                     marr = ms.arr, n.occasions = ncol(ch), rel = rowSums(ms.arr), 
                     resight.index = resight.index,  
                     ns = 5,  
                     zero = matrix(0, ncol = 5, nrow = 5), 
                     ones = diag(5), 
                     count = counts$Nibb, obs = observer, sigma.obs = counts$seNibb,
                     nest.obs.a = c(fecund.param$ns.obs.alpha, ns.a),
                     nest.obs.b = c(fecund.param$ns.obs.beta, ns.b),
                     K = K, BEFORE = BEFORE, AFTER = 1)

# parameters monitored
parameters <- c("Nb", "phiA", "phi0", "F", "mean.phi0", "mean.phiA", "mean.alpha", 
                "mean.log.F", "betaN", "beta", "sigma.o", "d", "mean.p", "p",
                "sigma.phi0", "sigma.phiA", "betaBP")

# MCMC settings
ni <- 40000; nt <- 1; nb <- 20000; nc <- 4

# Call JAGS from R (jagsUI)
# lead 2008, ice 8.5
start <- Sys.time()
out <- jags(jags.data8.5, inits, parameters, "YKD_IPM_L2008.jags", 
                      n.chains = nc, n.burnin=nb, n.iter = ni,  
                      parallel = TRUE, n.adapt = 1000)
Sys.time() - start
saveRDS(out, file = "MS_Scenarios/YKD.L2008.8.5.rds")
################################################################################
# lead 2026, ice 8.5
start <- Sys.time()
out <- jags(jags.data8.5, inits, parameters, "YKD_IPM_L2026.jags", 
                      n.chains = nc, n.burnin=nb, n.iter = ni,  
                      parallel = TRUE, n.adapt = 1000)
Sys.time() - start
saveRDS(out, file = "MS_Scenarios/YKD.L2026.8.5.rds")
################################################################################
# lead constant, ice 8.5
start <- Sys.time()
out <- jags(jags.data8.5, inits, parameters, "YKD_IPM_Lconst.jags", 
                      n.chains = nc, n.burnin=nb, n.iter = ni,  
                      parallel = TRUE, n.adapt = 1000)
Sys.time() - start
saveRDS(out, file = "MS_Scenarios/YKD.Lconst.8.5.rds")
################################################################################
#Lead 2008, ice 4.5
jags.data4.5 <- list(obs.ice = ice.data$ice_4.5, nb.size = theta4.5, muBP = muBP, SigmaBP=SigmaBP,
                     marr = ms.arr, n.occasions = ncol(ch), rel = rowSums(ms.arr), 
                     resight.index = resight.index,  
                     ns = 5,  
                     zero = matrix(0, ncol = 5, nrow = 5), 
                     ones = diag(5), 
                     count = counts$Nibb, obs = observer, sigma.obs = counts$seNibb,
                     nest.obs.a = c(fecund.param$ns.obs.alpha, ns.a),
                     nest.obs.b = c(fecund.param$ns.obs.beta, ns.b),
                     K = K, BEFORE = BEFORE, AFTER = 1)
saveRDS(jags.data4.5, file = "MS_Scenarios/jags.data4.5.rds")
start <- Sys.time()
out <- jags(jags.data4.5, inits, parameters, "YKD_IPM_L2008.jags", 
                      n.chains = nc, n.burnin=nb, n.iter = ni,  
                      parallel = TRUE, n.adapt = 1000)
Sys.time() - start
saveRDS(out, file = "MS_Scenarios/YKD.L2008.4.5.rds")
################################################################################
#Lead 2026, ice 4.5
start <- Sys.time()
out <- jags(jags.data4.5, inits, parameters, "YKD_IPM_L2026.jags", 
                      n.chains = nc, n.burnin=nb, n.iter = ni,  
                      parallel = TRUE, n.adapt = 1000)
Sys.time() - start
saveRDS(out, file = "MS_Scenarios/YKD.L2026.4.5.rds")
################################################################################
#Lead Constant, ice 4.5
start <- Sys.time()
out <- jags(jags.data4.5, inits, parameters, "YKD_IPM_Lconst.jags", 
                      n.chains = nc, n.burnin=nb, n.iter = ni,  
                      parallel = TRUE, n.adapt = 1000)
Sys.time() - start
saveRDS(out, file = "MS_Scenarios/YKD.Lconst.4.5.rds")
################################################################################
# Current conditions
#   defined as: ice is sampled from mean ice with nb variance as observed, 
#   lead declines ar 2008 to 2019 (or should it be 2026?)
jags.data.current <- list(obs.ice = ice.data$ice_current, nb.size = theta.current, 
                          muBP = muBP, SigmaBP=SigmaBP,
                          marr = ms.arr, n.occasions = ncol(ch), rel = rowSums(ms.arr), 
                          resight.index = resight.index,  
                          ns = 5,  
                          zero = matrix(0, ncol = 5, nrow = 5), 
                          ones = diag(5), 
                          count = counts$Nibb, obs = observer, sigma.obs = counts$seNibb,
                          nest.obs.a = c(fecund.param$ns.obs.alpha, ns.a),
                          nest.obs.b = c(fecund.param$ns.obs.beta, ns.b),
                          K = K, BEFORE = BEFORE, AFTER = 1)
start <- Sys.time()
out <- jags(jags.data.current, inits, parameters, "YKD_IPM_current.jags", 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
Sys.time() - start
saveRDS(out, file = "MS_Scenarios/YKD.current.rds")
