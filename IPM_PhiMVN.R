# from  IPM_YKD.R, updating model code for manuscript (IPM_MS) following work/discussions
# in winter 2021-22

# 1) Considering only 4 scenarios (lead 2008, RCP 4.5; lead constant, RCP 4.5;
# lead 2008, RCP 8.5; lead constant, RCP 8.5)

# 2) revising shape parameter in nb pulls to 4 (see nb_regress.R) s.t. the average
# variance of a simulated draw from model predicted ice from 1980 - 2019 is ~=
# the variance in the observed ice from 1980 - 2019.

# 3) revising EE priors for ice effects on breeding propensity,
# density dependence, and annual variation in aerial detection. See EE write up for
# all other EE priors. 1.27.2022 

# 4) added a linear effect for minimal (<15%) ice on survival for phiA and phi0;
# priors for this and the betas for extreme ice (and its square) obtained from the CJS-only
# model and added as a MVN

# 5) update input data, ice (from Dan R.) and counts (from Chuck F.) 

# 6) revising indexing to reflect addition years of input data

library(jagsUI)

# Data
# years of predictions, to 2100
K <- 79 # (count) data observed through 2021

# Count Data, 1988 - 2021
counts <- read.csv("input_data/YKD_SPEI.csv", header = T)
observer <- c(1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
              5, 5, 5, 5, 5, 5, 5)

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
sea.ice <- read.csv("input_data/sea.ice.data.csv", header = T)
# year is winter year (Nov - Apr)

# subsetting years 1988 to 2100
ice.data <- sea.ice[sea.ice$year >= 1988 & sea.ice$year <= 2100,]
ext.ice.data4.5 <- ifelse(is.na(ice.data$ext.ice.obs), ice.data$ext.ice.RCP4.5, 
                      ice.data$ext.ice.obs)
ext.ice.data8.5 <- ifelse(is.na(ice.data$ext.ice.obs), ice.data$ext.ice.RCP8.5, 
                      ice.data$ext.ice.obs)
min.ice.data4.5 <- ifelse(is.na(ice.data$min.ice.obs), ice.data$min.ice.RCP4.5, 
                          ice.data$min.ice.obs)
min.ice.data8.5 <- ifelse(is.na(ice.data$min.ice.obs), ice.data$min.ice.RCP8.5, 
                          ice.data$min.ice.obs)


# Fecundity Data, 1992 - 2015
fecund.param <- read.csv("input_data/fecundity.csv", header = T)
# convert mean and se of NS estimates to alpha and beta (Beta distribution) for
# input as covariates to the detection model (mark-recap)
fecund.param$ns.obs.alpha <- 
  round((fecund.param$ns.obs.gs-fecund.param$ns.obs.gs^2-fecund.param$sigma.ns.obs^2)*
          fecund.param$ns.obs.gs/fecund.param$sigma.ns.obs^2, 0)
fecund.param$ns.obs.beta <- round(fecund.param$ns.obs.alpha*(1-fecund.param$ns.obs.gs)
                                  /fecund.param$ns.obs.gs, 0)

# from BP_prior_EO_revised, mean and VCV for BP ice covariates
muBP <- c(1.2811, -0.1444, -0.0724)
SigmaBP <- matrix(c(0.3121, -0.0465, -0.0671, -0.0465, 0.0400, 0.0186, -0.0671,
                    0.0186, 0.0197), nrow = 3, ncol = 3)

# from CJS only model, mean and VCV for phi ice covariates (order is Intercept (Agephi0), 
# AgePhiA additive effect, ice.min.s, ice.ext.s, and ice.ext.s^2)
muPhi <- c(-0.9621, 2.7400, -0.1585, -0.2699, -0.2384)
SigmaPhi <- matrix(c(0.0299, -0.0201, 0.0113, 0.0039, -0.0100, -0.0201, 0.0402, 
                     0, 0, 0, 0.0113, 0, 0.0544, 0.0404, -0.0116, 0.0039, 0, 0.0404,
                     0.0430, -0.0040, -0.0100, 0, -0.0116, -0.0040, 0.0103),
                   nrow = 5, ncol = 5)

 
#### JAGS set up
jags.data4.5 <- list(ext.obs.ice = ext.ice.data4.5, min.obs.ice = min.ice.data4.5, 
                     nb.size = 4, muBP = muBP, SigmaBP=SigmaBP, muPhi = muPhi,
                     SigmaPhi = SigmaPhi, marr = ms.arr, n.occasions = ncol(ch), rel = rowSums(ms.arr), 
                     ns = ns, zero = matrix(0, ncol = ns, nrow = ns), ones = diag(ns), 
                     count = counts$Nibb, obs = observer, sigma.obs = counts$seNibb, vcf = counts$mvcf,
                     nest.obs.a = fecund.param$ns.obs.alpha[1:23],
                     nest.obs.b = fecund.param$ns.obs.beta[1:23],
                     K = K, BEFORE = 4, AFTER = 6) 
                     # BEFORE = no. of count years before survival data (4)
                     # AFTER = no. of count years after survival data

jags.data8.5 <- list(ext.obs.ice = ext.ice.data8.5, min.obs.ice = min.ice.data8.5,
                     nb.size = 4, muBP = muBP, SigmaBP=SigmaBP, muPhi = muPhi,
                     SigmaPhi = SigmaPhi, marr = ms.arr, n.occasions = ncol(ch), rel = rowSums(ms.arr), 
                     ns = ns, zero = matrix(0, ncol = ns, nrow = ns), ones = diag(ns), 
                     count = counts$Nibb, obs = observer, sigma.obs = counts$seNibb, vcf = counts$mvcf,
                     nest.obs.a = fecund.param$ns.obs.alpha[1:23],
                     nest.obs.b = fecund.param$ns.obs.beta[1:23],
                     K = K, BEFORE = 4, AFTER = 4) 
                     # BEFORE = no. of count years before survival data (4)
                     # AFTER = no. of count years after survival data

# initial values
inits <- function(){list(
  mean.phi0 = runif(1, 0.2, 0.3), 
  mean.phiA = runif(1, 0.8, 0.9), mean.p = runif(1, 0.5, 0.6))}

# parameters monitored
parameters <- c("Nb", "phiA", "phi0", "F", "mean.phi0", "mean.phiA", "alpha", 
                "mean.log.F", "betaN", "betaPhi", "sigma.o", "d", 
                "betaBP")

# MCMC settings
ni <- 20000; nt <- 1; nb <- 10000; nc <- 3


### Lead Constant, RCP4.5
cat(file = "YKD_IPM.jags", "
model {

## Priors
# state-space
sigma.obs[24] ~ dunif(500, 1400) # 2011; bounds consistent with s.e. observed since 2005
sigma.obs[33] ~ dunif(500, 1400) # 2020; bounds consistent with s.e. observed since 2005

vcf[24] ~ dunif(2.1, 2.4) # 2011; bounds consistent with s.e. observed since 2005
vcf[33] ~ dunif(2.1, 2.4) # 2020; bounds consistent with s.e. observed since 2005

# adding noise to the ice predictions s.t. interannual variation in the 
# future is consistent with the observed time series (1980 - 2018)
for (i in 1:(K+1)){
  ext.nb.prob[i] <- nb.size/(ext.obs.ice[n.occasions + BEFORE + AFTER - 1 + i] + nb.size) 
  ext.ice.draw[i] ~ dnegbin(ext.nb.prob[i], nb.size)
  min.nb.prob[i] <- nb.size/(min.obs.ice[n.occasions + BEFORE + AFTER - 1 + i] + nb.size) 
  min.ice.draw[i] ~ dnegbin(min.nb.prob[i], nb.size)
} # i
    
ext.ice.new <- c(ext.obs.ice[1:(n.occasions - 1 + BEFORE + AFTER)], ext.ice.draw[])
ext.ice <- (ext.ice.new - mean(ext.ice.new))/sd(ext.ice.new)
min.ice.new <- c(min.obs.ice[1:(n.occasions - 1 + BEFORE + AFTER)], min.ice.draw[])
min.ice <- (min.ice.new - mean(min.ice.new))/sd(min.ice.new)

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
betaBP[1:3] ~ dmnorm.vcov(muBP[1:3], SigmaBP[1:3, 1:3])
tau.BP <- pow(sigma.BP, -2)
sigma.BP ~ dunif(0,1)
for (i in 1:2){
betaN[i] ~ dgamma(3.1, 251.648) # from EE, see DDpriorsEO.R
} #i

# survival and breeding propensity
mean.alpha.inv ~ dgamma(7.1, 1.95) T(1,) # from EE
mean.logit.alpha <- logit(1/mean.alpha.inv)

betaPhi[1:5] ~ dmnorm.vcov(muPhi[1:5], SigmaPhi[1:5, 1:5])

tau.phi0 <- pow(sigma.phi0, -2)
sigma.phi0 ~ dunif(0, 1)
tau.phiA <- pow(sigma.phiA, -2)
sigma.phiA ~ dunif(0, 1)

omegaA ~ dbeta(1, 35) # from EE
omegaJ ~ dbeta(1, 9) # from EE

mean.p ~ dbeta(1,1)
mean.logit.p <- logit(mean.p)
tau.p <- pow(sigma.p, -2)
sigma.p ~ dunif(0, 1)
    
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
for (t in 1:(BEFORE+n.occasions-1+AFTER+K)){ 
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
for(i in 1:5){
  o[i] ~ dnorm(0, tau.o)
} # i
for (t in 1:(n.occasions+BEFORE+AFTER)){
  dev[t] ~ dbeta(2.16, 11.56) T(,0.41) # from EE, see ObsPriors.R
  sign[t] ~ dbin(0.5, 1)
  d[t] <- ifelse(sign[t] == 1, dev[t], -dev[t])
  count[t] ~ dnorm(2*exp(log(N[3,t] + N[4,t]) - o[obs[t]] - log(1 + vcf[t]*d[t])), tau.obs[t]) # count and vcf provided as data
  tau.obs[t] <- pow(sigma.obs[t], -2) # sigma.obs provided as data
} # t

####START HERE####

## Fecundity Model
# Process model
for (t in 1:(n.occasions+BEFORE+AFTER+K-1)){ # extended loop here # should go to end - 1
  logit.BP[t] <- betaBP[1] + betaBP[2]*ext.ice[t] + betaBP[3]*ext.ice[t]*ext.ice[t] + eps.BP[t]
  eps.BP[t] ~ dnorm(0, tau.BP)
  BP[t] <- ilogit(logit.BP[t])
  F[t] <- exp(log(phi0[t] + (1-phi0[t])*omegaJ) + mean.log.F 
          - betaN[1]*(N[3,t] + N[4,t])/1000 + log(BP[t]))
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
    logit.phi0[t] <- betaPhi[1] + betaPhi[3]*min.ice[t+1] + betaPhi[4]*ext.ice[t+1]
                     + betaPhi[5]*ext.ice[t+1]*ext.ice[t+1] + eps.phi0[t] 
    eps.phi0[t] ~ dnorm(0, tau.phi0)
    phi0[t] <- ilogit(logit.phi0[t])
    
    logit.phiA[t] <- logit(ilogit(betaPhi[1] + betaPhi[2])*(1-lead[t]*(1-kappa))) 
                     + betaPhi[3]*min.ice[t+1] + betaPhi[4]*ext.ice[t+1]
                     + betaPhi[5]*ext.ice[t+1]*ext.ice[t+1]
                     - betaN[2]*(N[3,t] + N[4,t])/1000 + eps.phiA[t] 
    logit.phi2[t] <- betaPhi[1] + betaPhi[2] + betaPhi[3]*min.ice[t+1] 
                     + betaPhi[4]*ext.ice[t+1] + betaPhi[5]*ext.ice[t+1]*ext.ice[t+1] 
                     - betaN[2]*(N[3,t] + N[4,t])/1000 + eps.phiA[t]
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

# Call JAGS from R (jagsUI), use autojags to run to convergence
YKD.constant.4.5.PhiMVN <- jags(jags.data4.5, inits, parameters, "YKD_IPM.jags", 
                         n.chains = nc, n.burnin=nb, n.iter = ni,  
                         parallel = TRUE, n.adapt = 1000)

saveRDS(YKD.constant.4.5.PhiMVN, file = "MS_Scenarios/YKD.constant.4.5.PhiMVN.rds")
