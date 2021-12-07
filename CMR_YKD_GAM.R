###This is the CMR component pulled out of IPM_YKD.R in order to add a GAM to the ice functional response

library(jagsUI)
library(mgcv)
library(tidyverse)

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
max.ice <- (181 - mean(ice.data)) / sd(ice.data)  #max.ice is the upper bound on the ice covariate
min.ice <- (0 - mean(ice.data)) / sd(ice.data)  #min.ice is the lower bound on the ice covariate
ice.data <- (ice.data-mean(ice.data)) / sd(ice.data)

#make data for GAM, transform ice covariate to cubic spline regression basis functions
#get example file for IPM
#ignore random effect
#change data set to match IPM
mu <- -0.1*ice.data^2 + rnorm(length(ice.data), 0, 0.1)
y <- rbinom(length(ice.data), 1, prob=plogis(mu))
df=data.frame(y=y, ice=ice.data)
#augment data to include min and max value of ice covariate, this "sorta worked, hard to converge, 
#strange behavior beyond range, lambda seen to reture prior
#df=data.frame(y=c(y, 0, 0), ice=c(ice.data, min.ice, max.ice)) 
#so basis functions span this range.
gam.data <- mgcv::jagam(y~s(ice, k=5, bs='cr')-1, family=binomial, data = df, file="test_gam.txt")

### Lead Constant
cat(file = "CMR_GAM.jags", "
model {

## Priors

# survival and breeding propensity
mean.alpha.inv ~ dgamma(7.1, 1.95) T(1,) # from EE
mean.logit.alpha <- logit(1/mean.alpha.inv)

mean.phi0 ~ dbeta(15,45) # from cjs model
mean.logit.phi0 <- logit(mean.phi0)
tau.phi0 <- pow(sigma.phi0, -2)
sigma.phi0 ~ dgamma(1,2) #dunif(0, 1)

mean.phiA ~ dbeta(90,10) # from cjs model
tau.phiA <- pow(sigma.phiA, -2)
sigma.phiA ~ dgamma(1, 2) #dunif(0, 1)

mean.p ~ dbeta(1,1)
mean.logit.p <- logit(mean.p)
tau.p <- pow(sigma.p, -2)
sigma.p ~ dgamma(1, 2) #dunif(0, 1)
    
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
###################################################
##set up priors for GAM, basis dimention is K-1 = 4 # also tried 9 with slow mixing, trying 4 because the efd ~ 2 in fitted model 
## prior for s(ice)... 
K1 <- S1[1:4,1:4] * lambda[1]  + S1[1:4,5:8] * lambda[2]
b[1:4] ~ dmnorm(bzero[1:4],K1) 
## smoothing parameter priors CHECK...
for (i in 1:2) {
  lambda[i] ~ dgamma(0.05,0.005) #experimented with dgamma(0.05, 0.005), dgamma(0.01, 0.001), and dgamma(2, 0.01)
  # latter seems necessary for quick converge when predicting (extrapolating) far beyond range of data (to 150 days ice free), 
  # but is very informative and causes a very strong quadratic. Without extrapolating, the other prior seem fine and give
  # less of a prefect quadratic. 
}
#GAM component of linear predictor
#X is the ice design matrix without intercept
eta <- X%*%b
for (t in 1:(n.occasions - 1)){
    logit.phi0[t] <- mean.logit.phi0 + eta[t+1] + eps.phi0[t] 
    eps.phi0[t] ~ dnorm(0, tau.phi0)
    phi0[t] <- ilogit(logit.phi0[t])
    
    logit.phiA[t] <- logit(mean.phiA) + eta[t+1] - betaN*(N[t]) + eps.phiA[t] #WHY IS ice at t+1 but N at t?
    logit.phi2[t] <- logit(mean.phiA) + eta[t+1] - betaN*(N[t]) + eps.phiA[t]
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
                  X = gam.data$jags.data$X[1:25,],
                  S1 = gam.data$jags.data$S1,
                  bzero = gam.data$jags.data$zero,
                  marr = ms.arr, n.occasions = ncol(ch), rel = rowSums(ms.arr), 
                  ns = ns, zero = matrix(0, ncol = ns, nrow = ns), ones = diag(ns), 
                  N = counts,
                  nest.obs.a = fecund.param$ns.obs.alpha[1:23],
                  nest.obs.b = fecund.param$ns.obs.beta[1:23])

# initial values
inits <- function(){list(
  b = gam.data$jags.ini$b,
  lambda = gam.data$jags.ini$lambda,
  mean.phi0 = runif(1, 0.2, 0.3), 
  mean.phiA = runif(1, 0.8, 0.9), mean.p = runif(1, 0.5, 0.6))}

# parameters monitored
parameters <- c("phiA", "phi0", "betaN", "b", "lambda", "beta.nest", 
                "mean.logit.alpha", "mean.phiA", "mean.phi0", "sigma.phi0", "sigma.phiA", "sigma.p")

# MCMC settings
ni <- 10000; nt <- 1; nb <- 1000; nc <- 4

# Call JAGS from R (jagsUI)
out <- jags(jags.data, inits, parameters, "CMR_GAM.jags", 
                         n.chains = nc, n.burnin=nb, n.iter = ni,  
                         parallel = TRUE, n.adapt = 1000)

saveRDS(out, file = "CMR.GAM.rds")
out <- readRDS(file = "CMR.GAM.rds")

#can't seem to figure out how to use mgcv::sim2jam to plot smooth
# smooth.post <- list(b=cbind(qlogis(out$sims.list$mean.phiA), out$sims.list$b))
# jam <- sim2jam(smooth.post, gam.data$pregam)
#do it myself
lp <-  cbind(rep(1, dim(gam.data$pregam$X)[1]), gam.data$pregam$X) %*% 
  t(cbind(qlogis(out$sims.list$mean.phiA), out$sims.list$b))

# df <- data.frame(ice=c(ice.data, min.ice, max.ice), mlp = apply(lp, 1, mean),
#                  upperlp = apply(lp, 1, quantile, probs=0.9),
#                  lowerlp = apply(lp, 1, quantile, probs=0.1))
df <- data.frame(ice=ice.data, mlp = apply(lp, 1, mean),
                 upperlp = apply(lp, 1, quantile, probs=0.9),
                 lowerlp = apply(lp, 1, quantile, probs=0.1))
df <- arrange(df, ice)
gplot <- ggplot(data=df, aes(x=ice, y=mlp))+
  geom_ribbon(aes(x=ice, ymin=lowerlp, ymax=upperlp))+
  geom_line()+
  labs(x="Standardized ice", y="linear predictor for phiA")
print(gplot)
#now do this on the real scale
lp <- plogis(lp)
df <- data.frame(ice=jags.data$ice, mlp = apply(lp, 1, mean), 
                 upperlp = apply(lp, 1, quantile, probs=0.9),
                 lowerlp = apply(lp, 1, quantile, probs=0.1))
df <- arrange(df, ice)
gplot <- ggplot(data=df, aes(x=ice, y=mlp))+
  geom_ribbon(aes(x=ice, ymin=lowerlp, ymax=upperlp))+
  geom_line()+
  labs(x="Standardized ice", y="phiA", title="GAM for eider survival fit using JAGS")
print(gplot)

