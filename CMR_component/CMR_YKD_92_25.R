###This is the CMR component pulled out of IPM_YKD.R in order to test on own
# 2026-07-20 modify to add 1992 - 2025 data after long debugg detailed in folder CMR_simulations
# Additions:
#  (1) new data file is coded as separate capture histories for adults and juveniles, 
#      but not coded as states as before.
#      Original state codes used in 2021 were: 
#      1 = duckling; 
#      2,3 = unobserved states; 2 = first year; 3 = second-year non-breeder;
#      4 = known second year breeder;
#      5 = adult.
#      To use Cat's original multiple state m-array format, need to add states to capture history.
#  (2) Original data 1992 - 2015 had no missing years. New data has no banding or 
#      re-sighting data for years 2016, 2017, 2018, and 2020. Need to deal with that.
#      (2.1) see code in CMR_simulations where missing data model is developed, 
#            debugged, tested against simulated data, and compared to individual capture history model
################################################################################
## define functions
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
library(tidyverse)
library(jagsUI)

# Data
# Count Data, 1988 - 2019
counts <- AKaerial::YKDVHistoric$combined |>
  filter(Year >= 1992) |>
  select(Year, ibb, ibb.se)
counts <- counts$ibb
counts[20] <- (counts[19]+counts[21])/2 #imput missing 2011 count
counts[29] <- (counts[28]+counts[30])/2 #imput missing 2020 count, assume mean

# Fecundity Data, 1992 - 2015
fecund.param <- read.csv("../input_data/fecundity.csv", header = T)
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

# Mark-Recap Data, 1992 - 2025
mcr.kig <- read.csv("../input_data/eh_na_1992-2025_adults_ducklings.csv", header = T) |>
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
# number of states for the multistate mark-recapture analysis
ns <- 5

# Survival Covariates, 1992 - 2025

sea.ice <- read.csv("../input_data/sea_ice_vars_1992_2025.csv", header = T) |>
# year is winter year (Nov - Apr), unlike in original data, the year for survival 
#   and ice match up. I think!
# ice data does not seem to match that used for 2021 SSA, need to ask
  select(year_winter, winter_index) |>
  mutate(ice_data = (winter_index - mean(winter_index))/sd(winter_index) )
################################################################################
## Start with simple model, no covariates, non-informative priors. 
## Get indexing and data handling correct!
cat(file = "CMR.missing.jags", "
model {

## Priors

# survival and breeding propensity
# mean.alpha.inv ~ dgamma(7.1, 1.95) T(1,) # from EE
# mean.logit.alpha <- logit(1/mean.alpha.inv)
# mean.logit.alpha ~ dnorm(0, 0.001) #make this equiv. to Dan's age-specific p effect
# alpha <- ilogit(mean.logit.alpha)
alpha ~ dbeta(1, 1)
#alpha <- 0.33

mean.phi0 ~ dbeta(1, 1) #dbeta(15,45) # from cjs model
#beta.phi0 ~ dnorm(0, 0.001)
tau.phi0 <- pow(sigma.phi0, -2)
sigma.phi0 ~ dunif(0, 1)

mean.phiA ~ dbeta(1, 1)
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
    logit.phi0[t] <- logit(mean.phi0) + eps.phi0[t]
    eps.phi0[t] ~ dnorm(0, tau.phi0)
    phi0[t] <- ilogit(logit.phi0[t]) # + beta.phi0) #make simple, force correlation between adults and duckling, makes sense
    
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
  po[4,t] <- p[t] * resight.index[t]
  po[5,t] <- p[t] * resight.index[t]
  
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
# banding.year indexes: !banding.year %in% 121:135, 141:145
for (t in 1:((n.occasions-1)*ns)){
# for (t in banding.year){
   marr[t,1:(n.occasions*ns-(ns-1))] ~ dmulti(pr[t,], rel[t])
   row.sum[t] <- sum(pr[t, 1:(n.occasions*ns-(ns-1))]) # just to check-debug
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
}
")
# # initial values
inits <- function(){list(
  mean.phi0 = 0.25, 
  mean.phiA = runif(1, 0.8, 0.9), 
  mean.p = runif(1, 0.5, 0.6), 
  sigma.phi0 = 0.5, 
  sigma.phiA = 0.5, 
  sigma.p = 0.5)}
# MCMC settings
ns <- 5
# make gap years
gap_calendar_years <- c(25:27, 29) #missing re-sight in calendar years 2016, 2017, 2018, and 2020
# resight.index is indexed by m-array recapture occasion = calendar year - 1
gap_occasions <- gap_calendar_years - 1
resight.index <- rep(1, ncol(ch) - 1)
resight.index[gap_occasions] <- 0
params <- c("phiA", "phi0", "p", "beta.phi0", "alpha",
                "mean.phiA", "mean.logit.p", "mean.phi0", 
                "sigma.phiA", "sigma.p", "row.sum")
jags.data <- list(marr = ms.arr, n.occasions = ncol(ch), rel = rowSums(ms.arr), 
                  resight.index = resight.index,  
                  ns = ns,  
                  zero = matrix(0, ncol = ns, nrow = ns), 
                  ones = diag(ns))
time <- Sys.time()
out <- jags(
  data = jags.data, inits = inits, parameters.to.save = params,
  model.file = "CMR.missing.jags", 
  n.chains = 3, n.adapt = 1000, n.burnin = 2000, n.iter = 6000,
  parallel = TRUE
)
Sys.time() - time
out$mean
plot_results(out)
#looks reasonable. 
################################################################################
## Add covariates
cat(file = "CMR.missing.covs.jags", "
model {

## Priors

# survival and breeding propensity
# mean.alpha.inv ~ dgamma(7.1, 1.95) T(1,) # from EE
# mean.logit.alpha <- logit(1/mean.alpha.inv)
# alpha <- ilogit(mean.logit.alpha)
alpha ~ dbeta(1, 1)
#alpha <- 0.33

mean.phi0 ~ dbeta(1, 1) #dbeta(15,45) # from cjs model
#beta.phi0 ~ dnorm(0, 0.001)
tau.phi0 <- pow(sigma.phi0, -2)
sigma.phi0 ~ dunif(0, 1)

mean.phiA ~ dbeta(1, 1)
tau.phiA <- pow(sigma.phiA, -2)
sigma.phiA ~ dunif(0, 1)

mean.p ~ dbeta(1,1)
mean.logit.p <- logit(mean.p)
tau.p <- pow(sigma.p, -2)
sigma.p ~ dunif(0, 1)

#prior for covariates
for (i in 1:2){
  beta[i] ~ dnorm(0, 100) #ice effect 
}
betaN.inv ~ dgamma(100000, 1) # from EE
betaN <- 1/betaN.inv # DD for phi and F

beta.nest ~ dnorm(0, 100) #prior for nest success effect on p

for (t in 1:(n.occasions-1)){
nest[t] ~ dbeta(nest.obs.a[t], nest.obs.b[t]) #nest sucess observation from priors
}

nest.s <- (nest - mean(nest))/sd(nest)
    
## Multistate survival model
# process model

for (t in 1:(n.occasions - 1)){
   logit.p[t] <- mean.logit.p  + beta.nest*nest.s[t] + eps.p[t]
   eps.p[t] ~ dnorm(0, tau.p)
   p[t] <- ilogit(logit.p[t])
} 

for (t in 1:(n.occasions - 1)){
    eps.phi0[t] ~ dnorm(0, tau.phi0)
    logit.phi0[t] <- logit(mean.phi0) 
                     + beta[1]*ice[t] #for now assume shared ice and density effect on ducklings
                     + beta[2]*ice[t]*ice[t] 
                     - betaN*(N[t]) + eps.phi0[t]  
                     
    eps.phiA[t] ~ dnorm(0, tau.phiA)
    logit.phiA[t] <- logit(mean.phiA) 
                     + beta[1]*ice[t] 
                     + beta[2]*ice[t]*ice[t] 
                     - betaN*(N[t]) + eps.phiA[t] 
    
    phi0[t] <- ilogit(logit.phi0[t])
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
  po[4,t] <- p[t] * resight.index[t]
  po[5,t] <- p[t] * resight.index[t]
  
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
# banding.year indexes: !banding.year %in% 121:135, 141:145
for (t in 1:((n.occasions-1)*ns)){
# for (t in banding.year){
   marr[t,1:(n.occasions*ns-(ns-1))] ~ dmulti(pr[t,], rel[t])
   row.sum[t] <- sum(pr[t, 1:(n.occasions*ns-(ns-1))]) # just to check-debug
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
}
")
params <- c("phiA", "phi0", "p", "beta.phi0", "alpha",
            "mean.phiA", "mean.logit.p", "mean.phi0", 
            "sigma.phiA", "sigma.p", "beta", "betaN", "beta.nest")
jags.data <- list(marr = ms.arr, n.occasions = ncol(ch), rel = rowSums(ms.arr), 
                  resight.index = resight.index,  
                  ns = ns,  
                  zero = matrix(0, ncol = ns, nrow = ns), 
                  ones = diag(ns), 
                  ice = sea.ice$ice_data,
                  N = counts,
                  nest.obs.a = c(fecund.param$ns.obs.alpha, ns.a),
                  nest.obs.b = c(fecund.param$ns.obs.beta, ns.b))
time <- Sys.time()
out <- jags(
  data = jags.data, inits = inits, parameters.to.save = params,
  model.file = "CMR.missing.covs.jags", 
  n.chains = 3, n.adapt = 1000, n.burnin = 2000, n.iter = 6000,
  parallel = TRUE
)
Sys.time() - time
out$mean
plot_results(out)
hist(out$sims.list$beta[,2])
hist(out$sims.list$betaN)
hist(out$sims.list$beta.nest)
################################################################################
