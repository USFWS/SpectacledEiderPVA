# DRAFT code started 2026-04-28
# modifying the year-specific phi model to ice_days + ice_days2 + age
# ice data are in: output/sic_beringSea_4pixels/sea_ice_vars_1995_2025.csv
# relevant variable is ws_index
# with year_winter being the year at the start of the winter season Nov year t through April year t+1
# runs and converges 2026-04-29

# CJS model phi ~ icedays + icedays^2 + age first year, p ~ random year, fixed age-3year old (phi_ice2_a_p_rT_a)

# model runs with data "data/eh_na_1992-2025_adults_ducklings.csv"

# Fit Hierarchical Bayesian Cormack-Joly-Seber model to mark-resight encounter histories from spectacled eiders at Kigigak 1992-2025
# no resight effort in 2015-2018 and 2020, with these years set to NA in the encounter history and indexed out of the likelihood for p
# code is based on Kery and Schaub 2012, Chpt 7 (but indexing of p is changed to not offset p from t)

# load data (encounter history created with file "create_encounterHistory_1992-2025_20250731.R")

# Erik Osnas, 2026-06-01, added random time effect on adult survival. Model originally from Dan R., fit_cjs_phi_ice2_a-p_rT_a_20260428.R
#    --changed hyperpriors on p and phi to dunif(0, 1)
#    --fitting this just to examine difference in comp. time and compare estimates to m-array format
#       --m-array format takes 20-30minute of comp time; this ind. format in Dan's original model takes about 8 hours. 
# load packages
#library(rjags)
library(jagsUI) #package to bridge R to JAGS and work with JAGS output
library(ggplot2)

set.seed(12) # make random results reproducible

# load encounter history data
CH <- read.csv("input_data/eh_na_1992-2025_adults_ducklings.csv") 
CH <- CH[order(CH$id_metalBand), ] # sort descending by id_metalBand
# remove column names and id_metalBand for JAGS input
ch <- CH # new data frame to modify
colnames(ch) <- NULL # remove col names
ch <- ch[ ,2:35] # remove first column for id_metalBand and last indicator variable column is_duckling

# create duckling indicator vector for creating an age matrix
is_duckling <- CH$is_duckling

# create a vector to index years with resight effort and exclude years without resight effort, for indexing p in model likelihood
is_resight <- c(1:24, 28, 30:34) # in time series, no resight effort during 2016, 2017, 2018, 2020 (years 25, 26, 27, 29)

# set index for the loop
n_resight <- length(is_resight)

# create vector indicating occasion of first capture/marking for each individual
get.first<- function(x) min(which(x != 0)) # return the index of the first non-zero value of each row
#f_old <- apply(ch, 1, get.first) # apply get.first to rows in ch
f <- apply(ch, 1, function(x) min(which(x == 1)))
hist(f,
     breaks = seq(1,ncol(ch),1),
     freq = FALSE)

# create matrix indicating third summer (earliest a female would return to breed) for individuals marked as ducklings and 0 otherwise
# for individs marked as ducklings, p_3rd summer will be estimated (2nd summer no resight is dealt with separately with "avail" matrix)
nind <- nrow(ch) # number of marked individuals in encounter history 1992-2025
n.occasions <- ncol(ch) # number of occasions 1992-2025 including years with no resight effort (2016-2018, 2020)

age_return <- matrix(0, nind, n.occasions) # create template matrix of all zeros

for (i in 1:nind) { # for each individual
  if (is_duckling[i] == 1) { # if it is indicated to be a duckling
    t2 <- f[i] + 1 # assign indicator for 2nd summer when individ will not return to breeding area
    t3 <- f[i] + 2  # assign indicator for third summer (summer of they turn 2 years old) as 2 years after marking occasion (f)
    if(t2 <= n.occasions) { # if second summer is within the range of encounter occasions
      age_return[i, t2] <- 0 # assign NA as individ should not have estimated resight for 2nd summer (actually asign 0 to just not estimate)
    }
    if (t3 <= n.occasions) { # if third summer is within the range of encounter occasions (i.e., don't project past range of data)
      age_return[i, t3] <- 1 # create indicator value of 1 for the third summer of individs marked as ducklings
    }
  }
  # adults stay 0
}

# deal with individs marked as ducklings not being available to resight in their 2nd summer when they remain at sea
# create indicator matrix with all 1s except for the 2nd summer of ducklings
# use this as a filter on mu2
avail <- matrix(1, nind, n.occasions)
for (i in 1:nind) { # for each individual
  if (is_duckling[i] == 1) { # if it is indicated to be a duckling
    t2 <- f[i] + 1 # assign indicator for 2nd summer when individ will not return to breeding area
    if(t2 <= n.occasions) { # if second summer is within the range of encounter occasions
      avail[i, t2] <- 0 # assign 0 to not estimate)
    }
  }
}

# setup age effect (first year) on phi
# create a first year indicator matrix for the first year of individs marked as ducklings for a first year survival effect on phi (derived value based on assumption that 2nd year survival = adult survival)
first_yr_duckling <- matrix(0, nind, n.occasions-1) # start with matrix of zeros, indexed to n.occasions-1 to match phi
for (i in 1:nind) { # for each individual
  if (is_duckling[i] == 1) { # if it is indicated to be a duckling
    first_yr <- f[i]  # assign indicator for 1st t (year) after a duckling is marked
    if(first_yr <= n.occasions) { # if first year is within the range of encounter occasions
      first_yr_duckling[i, first_yr] <- 1 # assign 1
    }
  }
}
  
# load sea ice data
ice_all <- read.csv("input_data/sea_ice_vars_1992_2025.csv")
ice_days <- ice_all$ws_count_sic_ge95 
# center ice_days
ice_days_c <- as.numeric(scale(ice_days, center = TRUE, scale = FALSE))

# set capture history data to y
y <- as.matrix(ch)

#################################################################################################
# CJS model: time-specific fixed phi with fixed age effect and RANDOM time p and fixed age p ####
#################################################################################################
# write to JAGS
writeLines("
    model {
  # Priors for Phi: quadratic ice-day effect + age effect for year 1 ducklings
  beta_ice0 ~ dnorm(0, 0.001) # prior for phi intercept
  beta_ice1 ~ dnorm(0, 0.001) # prior for phi ice_day linear effect
  beta_ice2 ~ dnorm(0, 0.001) # prior for phi ice_day^2 effect
  beta_first ~ dnorm(0, 0.001) # prior for the age effect on phi

  # Priors for p: time-varying random effect and fixed age effect
  for (t in 1:(n.occasions)) { # index to occasions p[2]=t[2], p[1] not in likelihood because t > f[i] and p[1] = f; NB this is different from Kery and Schaub indexing, they offset t-1 so p[1] = t[2]
    alpha_p[t] ~ dnorm(0, tau_p)  # draw prior for random effect for time from normal distribution with mean 0 and precision tau_p; each alpha_p[t] is drawn from the same distirbition
  }

  tau_p <- pow(sigma_p, -2)       # assign precision (inverse variance) for time random effect using sigma_p where 1/sigma_p^2 is the precision; sigma_p^2 is the variance 
  mu_p ~ dnorm(0, 0.001)          # draw prior for global intercept for p from normal distribution with mean 0 and precision 0.001
  beta_age ~ dnorm(0, 0.001)      # draw priors for fixed effect of age from normal distribution with mean 0 and precision 0.001
  sigma_p ~ dunif(0, 1)           # draw SD hyperparameter (prior for a prior) of time random effect from uniform distribution between 0 and 5; influences how much p varies over t

# priors: beta_phi[t], alpha_p[t], tau_p, mu_p, beta_age, sigma_p (where t is 34 occasions)
# prior for year random effect on phi
tau_phi <- pow(sigma_phi, -2)
sigma_phi ~ dunif(0, 1)
for (t in 1:(n.occasions-1)){
  eps.phi[t] ~ dnorm(0, tau_phi)
}
# Logit survival model with year effect + first-year duckling effect
  for (i in 1:nind) {
    for (t in 1:(n.occasions-1)) {
      logit(phi[i,t]) <- beta_ice0 + beta_ice1 * ice_days_c[t] + beta_ice2 * pow(ice_days_c[t], 2) + 
        beta_first * first_yr_duckling[i,t] + eps.phi[t]
  }
}
  # Likelihood
  for (i in 1:nind) {
    for (t in 1:(f[i]-1)){ # index of occasions before first capture for each individual i (needed because is_resight[k] (years with resight effort) in detection would otherwise be evaluated at k < f)
      z[i,t] <- 0 # set z to 0 pre-capture to ensure z exists everywhere and indicates an individ was not avail for detetection at t < f
    }
      z[i, f[i]] <- 1               # Latent state conditional on first capture (f[i]) 
     # State process 
    for (t in (f[i] + 1):n.occasions) { # indexed to n.occasions because individs are observed at each t
      mu1[i,t] <- phi[i, t - 1] * z[i, t - 1] # assign phi for preceding occasion conditional on z of preceding occasion (dead individs stay dead)
      z[i,t] ~ dbern(mu1[i,t])             # draw state z from Bernoulli distrib with probability mu1 for each individ and occasion
       } # close time
     } # close individuals
      # detection model
      for (i in 1:nind) {
        for (t in 1:n.occasions) {
      logit(p[i, t]) <- mu_p + alpha_p[t] + beta_age * age_return[i, t] # covariates on detection: random year, fixed age-class 3rd summer
     } # close time
      } # close individuals
      # Observation model
    for (i in 1:nind){ # loop over individs
      for (k in 1:n_resight){ # loop over occasions with resight effort (index vector n_resight)
      mu2[i, is_resight[k]] <-  # assign effective detection prob (i.e., prob that individ i is observed at time t), conditional on:
      p[i, is_resight[k]] * # individual detection probability
      z[i, is_resight[k]] * # latent state to ensure only individs that not already dead are resighted
      step(is_resight[k] - f[i]) * # indicator function to prevent modeling detection at or before marking
      avail[i,is_resight[k]] # indicator matrix to filter the 2nd summer of individs marked as ducklings
      y[i, is_resight[k]] ~ dbern(mu2[i, is_resight[k]]) # observation model with detection drawn from Bernoulli dist with prob mu2 if alive and post capture, otherwise 0 (detection not possble)
     } # close time
      } # close individuals
      
  # derived values
  for (t in 1:(n.occasions - 1)) {

  logit(phi_adult[t]) <- 
    beta_ice0 + 
    beta_ice1 * ice_days_c[t] + 
    beta_ice2 * pow(ice_days_c[t], 2) + 
    eps.phi[t]

  logit(phi_first[t]) <- 
    beta_ice0 + 
    beta_ice1 * ice_days_c[t] + 
    beta_ice2 * pow(ice_days_c[t], 2) + 
    beta_first + 
    eps.phi[t]

}

    } # close model statement
    ", con ="model.txt") # save model script as model.txt file (referenced in jags command later) and close writeLines function
# parameters estimated: mu_p, alpha_p[t-1], beta_age...

phi_ice2_a_p_rT_a <- 'model.txt' # save model as text file

# create matrix of known z states to include as data for the model (z = 1 after first capture to last encounter, NA otherwise, including NA for f)
z_known <- function(ch){ # encounter history ch as input
  ch_temp <- ch # change object name for tweaking here
  state <- as.matrix(ch_temp) # create state from ch as starting matrix to modify with known z information
  state[is.na(state)] <- 0 # turn columns of NAs in years of no resight effort to 0 for this function because all zeros will eventually be turned to NA unless they occur between two 1's in a row
  for (i in 1:dim(ch_temp)[1]){ # for each individual (row in state)
    n1 <- min(which(state[i, ] == 1)) # assign n1 as the minimum occasion with a 1 
    n2 <- max(which(state[i, ] == 1)) # assign n2 as the maximum occasion with a 1
    state[i, n1:n2] <- 1 # assign 1s between the min and max occasions with 1 (individual known alive even if not seen)
    state[i, n1] <- NA # assign NAs from the first occasion to n1
  }
  state[state == 0] <- NA # all occasions with 0 in ch assign NA in state
  return(state)
}

# Initial values ####
# z initial values based on ch, all known z changed to NA, keep ch values otherwise (i.e., 0, or consider changing to 1's?)
z_init <- function(ch) { # capture history as input
  ch[is.na(ch)] <- 0  # Replace NA with 0 (occasions with no resight effort)
  ch <- as.matrix(ch) # convert to matrix to avoid errors caused by data frame subsetting
  for (i in 1:nrow(ch)) { # for each individual, row
    if (sum(ch[i, ]) == 1) { # for individuals without any resights, skip
      next
    }
    n1 <- min(which(ch[i, ] == 1)) # assign n1 for min occasion with 1 in the row
    n2 <- max(which(ch[i, ] == 1)) # assign n2 for max occasion with 1 in the row
    ch[i, n1:n2] <- NA # assign NA between n1 and n2
    ch[i, n1] <- NA # assign NA for first capture
  }
  for(i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA # assign NA for occasions before first capture
  }
  return(ch) # output the modified ch_temp matrix as z_init
}

# all required initial values for estimated parameters as starting points for the estimation
inits <- function() {
  list(
    z = z_init(ch),
    #z = z.init,
    beta_ice0 = rnorm(0, 0.5),   # logit-scale survival effects, random init values centered on 0
    beta_ice1 = rnorm(0, 0.1),
    beta_ice2 = rnorm(0, 0.1),
    mu_p = rnorm(1, 0, 1), # logit-scale intercept for detection model, random init values centered on 0
    beta_age = rnorm(1, 0, 1), # age effect for detection model
    alpha_p = rnorm(n.occasions, 0, 1), # random time effect on detection
    sigma_p = runif(1, 0.1, 1), # SD of random time effect on detection, init must be positive, set as small positive value
    beta_first = rnorm(1, 0, 0.5) # age effect on phi
  )
}

# parameters to save
params <- c("beta_ice0", "beta_ice1", "beta_ice2", "beta_first", "phi_adult", 
            "phi_first", "mu_p", "beta_age", "sigma_p", "alpha_p") # "deviance", "z"
# run the model ####
start_time <- Sys.time()   # record start to measure model run time
# set MCMC parameters for test run and run model in JAGS with jags function
results <- jags(data = list(y = ch,
                            age_return = age_return,
                            n.occasions = n.occasions, 
                            nind = nind,
                            f = f,
                            z = z_known(ch),
                            is_resight = is_resight,
                            n_resight = n_resight,
                            avail = avail, 
                            first_yr_duckling = first_yr_duckling,
                            ice_days_c = ice_days_c),
                inits = inits,
                parameters.to.save = params,# what parameters to monitor for each iteration 
                model.file = "model.txt", # where to find the model code
                n.chains = 3, # 2 independent chains to permit assessing R hat (1)
                n.iter = 50000, # run each chain for 10k iterations (100)
                n.burnin = 10000, # discard the first 2k iterations (20)
                n.thin = 10, # keep every 5th value (1)
                parallel = TRUE) # run chains in parallel with multiple processor cores

end_time <- Sys.time()     # record end time
# calc total model run time and display
runtime <- end_time - start_time
print(runtime)
saveRDS(results, "results.RDS")
# results summaries
summary(results)
results
traceplot(results)

# create plot

# extract posterior samples
beta0 <- results$sims.list$beta_ice0
beta1 <- results$sims.list$beta_ice1
beta2 <- results$sims.list$beta_ice2
beta_first <- results$sims.list$beta_first
# create sequence of ice_days_c values from min to max
ice_seq <- seq(min(ice_days_c), max(ice_days_c), length.out = 200)
# define inverse logit function
ilogit <- function(x) exp(x) / (1 + exp(x))
# calc phi for adults
phi_adult_mat <- sapply(ice_seq, function(x) {
  ilogit(beta0 + beta1 * x + beta2 * x^2)
})
# calc phi for first year birds
phi_first_mat <- sapply(ice_seq, function(x) {
  ilogit(beta0 + beta1 * x + beta2 * x^2 + beta_first)
})
# get 95% credible intervals
phi_adult_df <- data.frame(
  ice_days_c = ice_seq,
  mean = apply(phi_adult_mat, 2, mean),
  lower = apply(phi_adult_mat, 2, quantile, 0.025),
  upper = apply(phi_adult_mat, 2, quantile, 0.9750))

ggplot(data = phi_adult_df, aes(x=ice_days_c, y = mean)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper))


df <- data.frame(Year = 1992:2024, 
                 PhiA = results$mean$phi_adult,
                 upper_PhiA = results$q97.5$phi_adult,
                 lower_PhiA = results$q2.5$phi_adult)
p_post <- plogis(matrix(results$sims.list$mu_p, nrow = length(results$sims.list$mu_p), 
                        ncol = 34, byrow = FALSE) + 
                   results$sims.list$alpha_p[,]) #beta_age = 0 is adults, see above
df$p <- apply(p_post, 2, mean)[-1] #I think this is correct due to comment above, ask Dan
df$upper_p <- apply(p_post, 2, quantile, probs = 0.975)[-1]
df$lower_p <- apply(p_post, 2, quantile, probs = 0.025)[-1]

ggplot(data = df, aes(x=Year, y = PhiA)) + 
  geom_pointrange(aes(ymin = lower_PhiA, ymax = upper_PhiA)) +
  scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))
ggplot(data = df, aes(x=Year, y = p)) + 
  geom_pointrange(aes(ymin = lower_p, ymax = upper_p)) +
  scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))

#Compare to m-array formatted model
out <- readRDS("test.rds")
df_marray <- data.frame(PhiA = apply(out$sims.list$phiA, 2, mean), 
                 upper_PhiA = apply(out$sims.list$phiA, 2, quantile, probs = 0.975),
                 lower_PhiA = apply(out$sims.list$phiA, 2, quantile, probs = 0.025), 
                 Year = 1992:2024 + 0.2, 
                 p = out$mean$p, 
                 upper_p = out$q97.5$p, 
                 lower_p = out$q2.5$p,
                 Type = "M-array")
ggplot(data = df_marray, aes(x = Year, y = PhiA)) + 
  geom_pointrange(aes(ymin = lower_PhiA, ymax = upper_PhiA)) +
  scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))

df$Type <- "Ind."
df2 <- rbind(df, df_marray)
ggplot(data = df2, aes(x = Year, y = PhiA, col = Type)) + 
  geom_pointrange(aes(ymin = lower_PhiA, ymax = upper_PhiA)) +
  scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))
ggplot(data = df2, aes(x = Year, y = p, col = Type)) + 
  geom_pointrange(aes(ymin = lower_p, ymax = upper_p)) +
  scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))
