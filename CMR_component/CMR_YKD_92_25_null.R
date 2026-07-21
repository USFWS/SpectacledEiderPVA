###This is the CMR component pulled out of IPM_YKD.R in order to test on own
# 2026-05-22 modify to add 1992 - 2025 data
# Additions:
#  (1) new data file is coded as separate capture histories for adults and juveniles, but not coded as state as before.
#      Original state codes used in 2021 were: 1 = duckling; 
#      2,3 = unobserved states; 2 = first year non-breeder; 3 = second year non-breeder;
#      4 = known second year breeder;
#      5 = adult.
#      To use Cat's original multiple state m-array format, need to add states to capture history.
#  (2) Original data 1992 - 2015 had no missing years. New data has no banding or 
#      re-sighting data for years 2016, 2017, 2018, and 2020. Need to deal with that.
#      (2.1) Proposed solution is to set re-sighting probability to zero for 
#            these years and estimate survival as a latent variable in Bayesian model. 
# The goal of this file (CMR_YKD_92_25_null.R) is to compare the individual sate-space formulation 
# to an m-array state-space formulation
#   --the individual based model is in CMR_IND_SS_null.R
#   --this file builds off of (CMR_YKD_92_25.R) where I was trying to add the 
#       missing year to the m-array model used in the IMP. 
#   --the idea is to use random effects for phi and p and get the two models as close as possible
#   --also mean to use the same priors and few or no covariates
#   -- the main problem is to be sure that the m-array model is structure 
#      correctly for the missing release and resight years.
# changed priors and structure the match Dan model as close as possible

library(tidyverse)
library(jagsUI)

# Data
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


### Lead Constant
cat(file = "CMR_null.jags", "
model {

## Priors

# survival and breeding propensity
# mean.alpha.inv ~ dgamma(7.1, 1.95) T(1,) # from EE
# mean.logit.alpha <- logit(1/mean.alpha.inv)
mean.logit.alpha ~ dnorm(0, 0.001) #make this equiv. to Dan's age-specific p effect
alpha <- ilogit(mean.logit.alpha)

#mean.phi0 ~ dbeta(15,45) # from cjs model
mean.logit.phi0 ~ dnorm(-2.8, (1/0.26)^2) #just to make it the same as the Ind. SS model
beta.phi0 <- mean.logit.phi0 - logit(mean.phiA)
tau.phi0 <- pow(sigma.phi0, -2)
sigma.phi0 ~ dunif(0, 1)

mean.phiA ~ dbeta(90,10) # from cjs model
tau.phiA <- pow(sigma.phiA, -2)
sigma.phiA ~ dunif(0, 1)

mean.p ~ dbeta(1,1)
mean.logit.p <- logit(mean.p)
tau.p <- pow(sigma.p, -2)
sigma.p ~ dunif(0, 1)
    
## Multistate survival model
# process model

for (t in 1:(n.occasions - 1)){ 
   logit.p[t] <- mean.logit.p + eps.p[t]
   eps.p[t] ~ dnorm(0, tau.p)
   p[t] <- ilogit(logit.p[t])
} 

for (t in 1:(n.occasions - 1)){
    # logit.phi0[t] <- mean.logit.phi0 + eps.phi0[t]
    # eps.phi0[t] ~ dnorm(0, tau.phi0)
    phi0[t] <- ilogit(logit.phiA[t] + beta.phi0) #make simple to match Dan, force correlation between adults and duckling, makes sense
    
    logit.phiA[t] <- logit(mean.phiA) + eps.phiA[t] 
    eps.phiA[t] ~ dnorm(0, tau.phiA)
    phiA[t] <- ilogit(logit.phiA[t])
    phi2[t] <- ilogit(logit.phiA[t])
    phi1[t] <- phi2[t]
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
  psi[2,t,3] <- phi1[t]*(1-alpha) #there is no breeding propensity in Dan's model, he has age-specific p
  psi[2,t,4] <- phi1[t]*alpha
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
# skip missing years: 2016, 2017, 2018, 2020
# marr indexes 121:135, 141:145
# for (t in 1:((n.occasions-1)*ns)){
for (t in marr.index){
   marr[t,1:(n.occasions*ns-(ns-1))] ~ dmulti(pr2[t,], rel[t])
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

# Now set probabilities to zero for no-effort year: 2016-2018, 2020   
for (i in 1:((n.occasions-1)*ns)){
   for (j in not.miss.index){
     pr2[i,j] <- pr[i,j]
   }
   for (j in miss.index){
     pr2[i,j] <- 0
   }
   #set probability of never recapture, need to recalculate because of new zeros
   pr2[i,(n.occasions*ns-(ns-1))] <- 1 - sum(pr[i,1:((n.occasions-1)*ns)])
}
   
}
")

# number of states for the multistate mark-recapture analysis
ns <- 5
# bundle data
jags.data <- list(marr = ms.arr, n.occasions = ncol(ch), rel = rowSums(ms.arr), 
                  marr.index = c(1:120, 136:140, 146:165), 
                  miss.index = c(116:130, 136:140),
                  not.miss.index = c(1:115, 131:135, 141:165),
                  ns = ns,  
                  zero = matrix(0, ncol = ns, nrow = ns), 
                  ones = diag(ns))

# initial values
inits <- function(){list(
  mean.phi0 = 0.25, 
  mean.phiA = runif(1, 0.8, 0.9), 
  mean.p = runif(1, 0.5, 0.6), 
  sigma.phi0 = 0.5, 
  sigma.phiA = 0.5, 
  sigma.p = 0.5)}

# parameters monitored
parameters <- c("phiA", "phi0", "p", 
                "alpha", "mean.phiA", "mean.phi0", "mean.logit.p", "mean.logit.alpha",
                "sigma.phi0", "sigma.phiA", "sigma.p")

# MCMC settings
ni <- 11000; nt <- 1; nb <- 5000; nc <- 3

# Call JAGS from R (jagsUI), use autojags to run to convergence
time <- Sys.time()
out <- jags(jags.data, inits, parameters, "CMR_null.jags", 
                         n.chains = nc, n.burnin=nb, n.iter = ni,  
                         parallel = TRUE, n.adapt = 1000)
Sys.time() - time
saveRDS(out, file = "test_null.rds")
#
jagsUI::traceplot(out, parameters = "phiA")

plot(1992:2024, out$mean$phiA)
plot(1992:2024, out$mean$p)
plot(out$mean$phiA, out$mean$p)
plot(out$sims.list$phiA[,1], out$sims.list$p[,1])
plot(out$sims.list$phiA[,33], out$sims.list$p[,33])
#plot adults
df <- data.frame(PhiA = apply(out$sims.list$phiA, 2, mean), 
                 upper = apply(out$sims.list$phiA, 2, quantile, probs = 0.975),
                 lower = apply(out$sims.list$phiA, 2, quantile, probs = 0.025), 
                 Year = 1992:2024)
ggplot(data = df, aes(x = Year, y = PhiA)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))
#plot ducklings
df0 <- data.frame(Phi0 = apply(out$sims.list$phi0, 2, mean), 
                 upper = apply(out$sims.list$phi0, 2, quantile, probs = 0.975),
                 lower = apply(out$sims.list$phi0, 2, quantile, probs = 0.025), 
                 Year = 1992:2024)
ggplot(data = df0, aes(x = Year, y = Phi0)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.01))
#plot p
dfp <- data.frame(p = apply(out$sims.list$p, 2, mean), 
                  upper = apply(out$sims.list$p, 2, quantile, probs = 0.975),
                  lower = apply(out$sims.list$p, 2, quantile, probs = 0.025), 
                  Year = 1992:2024)
ggplot(data = dfp, aes(x = Year, y = p)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.01))

AKaerial::YKDVHistoric$combined |>
  filter(Year >= 1992) |>
  select(Year, ibb, ibb.se) |>
  ggplot(aes(x = Year, y = ibb)) + geom_point() + 
  scale_x_continuous(breaks = seq(1992, 2024, by = 4))

