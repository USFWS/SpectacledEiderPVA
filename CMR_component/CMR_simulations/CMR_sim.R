# simulate a capture-recapture history data to test model
# main issue to to make sure model structure with missing release and resight years gives correct estimates
# year varying phi and p
# age and stage structure as in IPM for eiders
# 1 = ducklings
# 2 = second year, not observable
# 3 = third year non-breeders, not observable
# 4 = third year breeders, observable
# 5 = adult breeders
library(tidyverse)
library(jagsUI)

################################################################################
# function to simulate data for multi-state capture history
# input is number of occasions and true parameter values
# states described above and below
ch_data_sim <- function(
    n_occasions = 34, #re-sight occasions by year
    
    d_marked = rep(1000, each = n_occasions), #number of ducklings each year
    a_marked = rep(1000, each = n_occasions), #number of adults
    
    # Probabilities (Adjust these to fit your specific study)
    phi0 = rep(0.3, n_occasions),   # Survival probability of ducklings marked in nest to 1-year old
    phi1 = rep(0.80, n_occasions),   # Survival probability of second year 
    phi2 = rep(0.80, n_occasions),  # Survival probability of third year
    phiA = rep(0.80, n_occasions),  # Survival probability of adults
    alpha = 0.33,                   # Breeding probability of second year
    p_adult = rep(0.4, n_occasions)# Re-sight probability of adults
){
  #set up ch matrix
  ch_d <- ch_a <- data.frame(NULL)
  
  # Simulate individual histories
  # set.seed(1234) # Set seed for reproducibility
  for(j in 1:(n_occasions - 1)){
    # Initialize empty list and matrix for histories
    ch_alist <- ch_dlist <- list()
    #ducklings first
    for(i in 1:d_marked[j]) {
      # State tracking: 1 = Alive (Young, first year), 2 = Alive (Second year), 
      # 3 = Alive (3rd year but not breeding), 4 = Alive (3rd year, breeding), 
      # 5 = Alive Adults, 0 = Dead/Not in sight
      state <- 1 
      history <- numeric(n_occasions - j + 1)
      
      # Occasion 1: All marked as young and released
      history[1] <- 1
      
      # Occasion 2: Can't be re-captured/re-sighted. Survival is the only process.
      if(runif(1) > phi0[j]) {
        state <- 0 # Died in first year
      } else {
        state <- 2 # Survived and transitioned to second year
      }
      history[2] <- 0 # Set to 0 because they can't be seen in Year 2
      # Occasion 3 onwards: Adult stage
      if(length(history) > 2){
        for(t in 3:(n_occasions - j + 1)) {
          if(state == 0) {
            history[t] <- 0
          } else {
            if(state == 2){
              # Survival step
              if(runif(1) > phi1[j + t - 2]) {
                state <- 0
                history[t] <- 0
              } else {
                if(runif(1) < alpha){
                  state <- 4
                  if(runif(1) < p_adult[j + t - 2]) {
                    history[t] <- 1
                  } else {
                    history[t] <- 0
                  }
                } else {
                  state <- 3
                  history[t] <- 0
                }
              }
            } else if(state == 3){
              if(runif(1) > phi2[j + t - 2]) {
                state <- 0
                history[t] <- 0
              } else {
                state <- 5
                if(runif(1) < p_adult[j + t - 2]) {
                  history[t] <- 1
                } else {
                  history[t] <- 0
                }
              }
            } else if(state == 4){
              if(runif(1) > phiA[j + t - 2]) {
                state <- 0
                history[t] <- 0
              } else {
                state <- 5
                if(runif(1) < p_adult[j + t - 2]) {
                  history[t] <- 1
                } else {
                  history[t] <- 0
                }
              }
            } else if(state == 5){
              if(runif(1) > phiA[j + t - 2]) {
                state <- 0
                history[t] <- 0
              } else {
                if(runif(1) < p_adult[j + t - 2]) {
                  history[t] <- 1
                } else {
                  history[t] <- 0
                }
              }
            }
          } # closes the `else {` from `if(state == 0)`
        } # closes the `for(t ...)` loop
      }
      
      ch_dlist[[i]] <- c(rep(0, j - 1), history)
    }
    #now adults
    for(i in 1:a_marked[j]) {
      # State tracking: 1 = Alive (Young, first year), 2 = Alive (Second year), 
      # 3 = Alive (3rd year but not breeding), 4 = Alive (3rd year, breeding), 
      # 5 = Alive Adults, 0 = Dead/Not in sight
      state <- 5 # assume all marked not as ducklings are adults
      history <- numeric(n_occasions - j + 1)
      
      # Occasion 1: All marked as adult and released
      history[1] <- 1
      
      # Occasion 2 onwards: Adult stage
      for(t in 2:(n_occasions - j +1)) {
        if(state == 0) {
          history[t] <- 0
        } else {
          # Survival step
          if(runif(1) > phiA[j + t - 1]) {
            state <- 0 # Died
            history[t] <- 0
          } else {
            # Detection step
            if(runif(1) < p_adult[j + t - 1]) {
              history[t] <- 1
            } else {
              history[t] <- 0
            }
          }
        }
      }
      ch_alist[[i]] <- c(rep(0, j - 1), history)
    }
    
    #ch_list <- c(ch_dlist, ch_alist)
    # Convert list to a dataframe
    ch_d_matrix <- do.call(rbind, ch_dlist)
    colnames(ch_d_matrix) <- paste0("Year", 1:n_occasions)
    ch_a_matrix <- do.call(rbind, ch_alist)
    colnames(ch_a_matrix) <- paste0("Year", 1:n_occasions)
    ch_d <- rbind(ch_d, as.data.frame(ch_d_matrix))
    ch_a <- rbind(ch_a, as.data.frame(ch_a_matrix))
  }
  
  ch <- rbind(data.frame(ch_d, is_duckling = 1), data.frame(ch_a, is_duckling = 0)) |>
    rownames_to_column(var = "id_metalBand")
  return(ch)
}
#set some years (columns) to zeros and remove those release years
ch <- ch_data_sim()
#prep data for JAGS model
#make data formatting a function
mcr <- function(ch, starts.with = "Y"){
  x <- ch |>
    arrange(desc(is_duckling)) |>
    pivot_longer(cols = starts_with(starts.with), names_to = "Occasion", 
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
  # There has to be a better way!
  # Now pivot back to wide form
  x <- select(x, -is_duckling, -temp) |>
    pivot_wider(names_from = "Occasion", values_from = "State")
  
  x <- as.matrix(x[, -1])   # drop band id
  return(x)
}
#apply it
ch.sim <- mcr(ch=ch)
#saveRDS(ch.sim, "ch.sim_1.RDS")
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
####
ms.arr <- marray(ch.sim)
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
    scale_y_continuous(breaks = seq(0, 1, by = 0.01))
  # Not good!
  #plot ducklings
  df0 <- data.frame(Phi0 = apply(out$sims.list$phi0, 2, mean), 
                    upper = apply(out$sims.list$phi0, 2, quantile, probs = 0.975),
                    lower = apply(out$sims.list$phi0, 2, quantile, probs = 0.025), 
                    Year = 1992:2024)
  p2 <- ggplot(data = df0, aes(x = Year, y = Phi0)) + 
    geom_pointrange(aes(ymin = lower, ymax = upper)) +
    scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.01))
  # similarly bad! 
  #plot p
  dfp <- data.frame(p = apply(out$sims.list$p, 2, mean), 
                    upper = apply(out$sims.list$p, 2, quantile, probs = 0.975),
                    lower = apply(out$sims.list$p, 2, quantile, probs = 0.025), 
                    Year = 1993:2025)
  p3 <- ggplot(data = dfp, aes(x = Year, y = p)) + 
    geom_pointrange(aes(ymin = lower, ymax = upper)) +
    scale_x_continuous(breaks = seq(1993, 2025, by = 4)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.01))
  return(list(PhiA=p1, Phi0=p2, p=p3))
}

### modify model for complete data (no missing years)
cat(file = "CMR.jags", "
model {

## Priors

# survival and breeding propensity
# mean.alpha.inv ~ dgamma(7.1, 1.95) T(1,) # from EE
# mean.logit.alpha <- logit(1/mean.alpha.inv)
# mean.logit.alpha ~ dnorm(0, 0.001) #make this equiv. to Dan's age-specific p effect
# alpha <- ilogit(mean.logit.alpha)
alpha ~ dbeta(1, 1)

#mean.phi0 ~ dbeta(15,45) # from cjs model
beta.phi0 ~ dnorm(0, 0.001)
# tau.phi0 <- pow(sigma.phi0, -2)
# sigma.phi0 ~ dunif(0, 1)

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
    # logit.phi0[t] <- mean.logit.phi0 + eps.phi0[t]
    # eps.phi0[t] ~ dnorm(0, tau.phi0)
    phi0[t] <- ilogit(logit.phiA[t] + beta.phi0) #make simple, force correlation between adults and duckling, makes sense
    
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
for (t in 1:((n.occasions-1)*ns)){
# for (t in marr.index){
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

# # Now set probabilities to zero for no-effort year: 2016-2018, 2020   
# for (i in 1:((n.occasions-1)*ns)){
#    for (j in not.miss.index){
#      pr2[i,j] <- pr[i,j]
#    }
#    for (j in miss.index){
#      pr2[i,j] <- 0
#    }
#    #set probability of never recapture, need to recalculate because of new zeros
#    pr2[i,(n.occasions*ns-(ns-1))] <- 1 - sum(pr[i,1:((n.occasions-1)*ns)])
# }
   
}
")
# number of states for the multistate mark-recapture analysis
ns <- 5
# bundle data
jags.data <- list(marr = ms.arr, n.occasions = ncol(ch.sim), rel = rowSums(ms.arr), 
                  # marr.index = c(1:120, 136:140, 146:165), 
                  # miss.index = c(116:130, 136:140),
                  # not.miss.index = c(1:115, 131:135, 141:165),
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
# MCMC settings
ni <- 11000; nt <- 1; nb <- 5000; nc <- 3

parameters <- c("phiA", "phi0", "p", "beta.phi0", "alpha",
                "mean.phiA", "mean.logit.p", "mean.phi0", 
                "sigma.phiA", "sigma.p")
time <- Sys.time()
out <- jags(jags.data, inits, parameters, "CMR.jags", 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
Sys.time() - time
#plot adults
df <- data.frame(PhiA = apply(out$sims.list$phiA, 2, mean), 
                 upper = apply(out$sims.list$phiA, 2, quantile, probs = 0.975),
                 lower = apply(out$sims.list$phiA, 2, quantile, probs = 0.025), 
                 Year = 1992:2024)
ggplot(data = df, aes(x = Year, y = PhiA)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.01))
#looks good
#plot ducklings
df0 <- data.frame(Phi0 = apply(out$sims.list$phi0, 2, mean), 
                  upper = apply(out$sims.list$phi0, 2, quantile, probs = 0.975),
                  lower = apply(out$sims.list$phi0, 2, quantile, probs = 0.025), 
                  Year = 1992:2024)
ggplot(data = df0, aes(x = Year, y = Phi0)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.01))
#about 0.065 too low
#plot p
dfp <- data.frame(p = apply(out$sims.list$p, 2, mean), 
                  upper = apply(out$sims.list$p, 2, quantile, probs = 0.975),
                  lower = apply(out$sims.list$p, 2, quantile, probs = 0.025), 
                  Year = 1992:2024)
ggplot(data = dfp, aes(x = Year, y = p)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.01))
#close, about 0.01 -0.02 too low
hist(out$sims.list$beta.phi0)
hist(out$sims.list$mean.logit.p)
hist(out$sims.list$phi0)
mean(out$mean$phi0)
#so phi0 is about 0.1 too low. 
# given that phi0 is the only free parameter, it should adjust to phi0*phi1 = 0.3*0.8 ~= 0.23
out$mean$mean.logit.p
mean(out$mean$p)
#just about 0.01-0.02 too low; fairly good
hist(out$sims.list$sigma.p)
hist(out$sims.list$sigma.phiA)
plot(out$sims.list$mean.logit.p, out$sims.list$beta.phi0)
traceplot(out, parameters = "beta.phi0")
traceplot(out, parameters = "mean.logit.p")
plot(out$sims.list$beta.phi0, out$sims.list$alpha)
#calculate bias
mean(out$mean$phiA - 0.8)
mean(out$mean$phi0 - 0.3)
mean(out$mean$p - 0.4)
mean(apply(out$sims.list$phiA*out$sims.list$phi0, 2, mean) - 0.3*0.8)
################################################################################
## test above with missing years
# How to handle the missing years? (1) set column of pr to 0 manualy [did not work!] or 
#   (2) set p to zero in code. Same?
cat(file = "CMR.missing.jags", "
model {

## Priors

# survival and breeding propensity
# mean.alpha.inv ~ dgamma(7.1, 1.95) T(1,) # from EE
# mean.logit.alpha <- logit(1/mean.alpha.inv)
# mean.logit.alpha ~ dnorm(0, 0.001) #make this equiv. to Dan's age-specific p effect
# alpha <- ilogit(mean.logit.alpha)
alpha ~ dbeta(1, 1)

#mean.phi0 ~ dbeta(15,45) # from cjs model
beta.phi0 ~ dnorm(0, 0.001)
# tau.phi0 <- pow(sigma.phi0, -2)
# sigma.phi0 ~ dunif(0, 1)

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
    # logit.phi0[t] <- mean.logit.phi0 + eps.phi0[t]
    # eps.phi0[t] ~ dnorm(0, tau.phi0)
    phi0[t] <- ilogit(logit.phiA[t] + beta.phi0) #make simple, force correlation between adults and duckling, makes sense
    
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
# for (t in 1:((n.occasions-1)*ns)){
for (t in banding.year){
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

# # Now set probabilities to zero for no-effort year: 2016-2018, 2020
# for (i in 1:((n.occasions-1)*ns)){
#    for (j in resight.year){
#      pr2[i,j] <- pr[i,j]
#    }
#    for (j in no.resight.year){
#      pr2[i,j] <- 0
#    }
#    #set probability of never recapture, need to recalculate because of new zeros
#    pr2[i,(n.occasions*ns-(ns-1))] <- 1 - sum(pr[i,1:((n.occasions-1)*ns)])
# }
   
}
")

#Now set up data
no.resight.year <- c(116:130, 136:140)
ms.arr.missing <- ms.arr
ms.arr.missing[,no.resight.year] <- 0 #missing resight years
ms.arr.missing[c(121:135, 141:145),] <- 0 #missing banding years
sum(ms.arr.missing[,no.resight.year])
sum(ms.arr[,no.resight.year])
sum(ms.arr.missing[c(121:135, 141:145),])
sum(ms.arr[c(121:135, 141:145),])
#make jags data
jags.data <- list(marr = ms.arr.missing, n.occasions = ncol(ch.sim), rel = rowSums(ms.arr.missing), 
                  banding.year = c(1:120, 136:140, 146:165),
                  no.resight.year = no.resight.year,
                  resight.index = as.numeric( ! 1:165 %in% no.resight.year ),
                  resight.year = c(1:115, 131:135, 141:165),
                  ns = ns,  
                  zero = matrix(0, ncol = ns, nrow = ns), 
                  ones = diag(ns))
time <- Sys.time()
out <- jags(jags.data, inits, parameters, "CMR.missing.jags", 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
Sys.time() - time
#plot adults
df <- data.frame(PhiA = apply(out$sims.list$phiA, 2, mean), 
                 upper = apply(out$sims.list$phiA, 2, quantile, probs = 0.975),
                 lower = apply(out$sims.list$phiA, 2, quantile, probs = 0.025), 
                 Year = 1992:2024)
ggplot(data = df, aes(x = Year, y = PhiA)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))
# Not good!
#plot ducklings
df0 <- data.frame(Phi0 = apply(out$sims.list$phi0, 2, mean), 
                  upper = apply(out$sims.list$phi0, 2, quantile, probs = 0.975),
                  lower = apply(out$sims.list$phi0, 2, quantile, probs = 0.025), 
                  Year = 1992:2024)
ggplot(data = df0, aes(x = Year, y = Phi0)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))
# similarly bad! 
#plot p
dfp <- data.frame(p = apply(out$sims.list$p, 2, mean), 
                  upper = apply(out$sims.list$p, 2, quantile, probs = 0.975),
                  lower = apply(out$sims.list$p, 2, quantile, probs = 0.025), 
                  Year = 1992:2024)
ggplot(data = dfp, aes(x = Year, y = p)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  scale_x_continuous(breaks = seq(1992, 2024, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))
# Not working!
##After much thought and using Claude Opus 4.8 to help debug, there is a proposed problem in the indexing of U
##  The proposed fix os to change the U line to:
#   U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% psi[,j-1,] %*% dq[,j-1,]
#   Note the j-1 index instead of t
#   We will implement this and test it against the simulated data above. 
#   Will also add a diagnostic to check that the multinomial probability sum to 1
# for (t in banding.year){
#   row.sum[t] <- sum(pr[t, 1:(n.occasions*ns-(ns-1))])
# }
# We also found a problem in the data simulation function that led to too many duckling being resight, given parameters. 
#   fixed in above on 2026-06-26
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
#saveRDS(ch.sim, "ch.sim_1.RDS")
# ch.sim  <- readRDS("ch.sim_1.RDS") 
# ms.arr <- marray(ch.sim)

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
ni <- 11000; nt <- 1; nb <- 5000; nc <- 3
#parameters to monitor, add row.sum
parameters <- c("phiA", "phi0", "p", "beta.phi0", "alpha",
                "mean.phiA", "mean.logit.p", "mean.phi0", 
                "sigma.phiA", "sigma.p", "row.sum")
#######################
#first run with no missing years (gaps)
jags.data <- list(marr = ms.arr, n.occasions = ncol(ch.sim), rel = rowSums(ms.arr), 
                  banding.year = c(1:165),
                  #no.resight.year = no.resight.year,
                  resight.index = rep(1, 33), 
                  #resight.year = c(1:115, 131:135, 141:165),
                  ns = ns,  
                  zero = matrix(0, ncol = ns, nrow = ns), 
                  ones = diag(ns))
time <- Sys.time()
out <- jags(jags.data, inits, parameters, "CMR.missing.jags", 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
Sys.time() - time
out$mean
plot_results(out)
plot(out$sims.list$alpha, out$sims.list$mean.phi0)
# Code used during debugging:
# #calculate observed resight frequency of ducklings from sim data
# ch.duckling <- ch.sim[1:33000,]
# #remove releases that can never be seen 
# ch.duckling <- ch.duckling[ch.duckling[,"Year33"] != 1,]
# get.first <- function(x) min(which(x!=0))
# first.duck <- apply(ch.duckling, 1, get.first)
# index <- cbind(c(1:3200), first.duck+2)
# second.duck <- ch.duckling[index]
# sum(second.duck > 0)/length(second.duck) #seems to low
# #try working on m-array
# ma <- ms.arr[seq(1, 165, by=5),]
# inx <- cbind(1:32, seq(9, 165, by = 5))
# mean(ma[inx]/1000)
# sum(ch[1:1000, 4] > 0)/1000
# sum(ch[1001:2000, 5] > 0)/1000
# sum(ch[2001:3000, 6] > 0)/1000
# x <- c()
# for(i in 1:31){x[i] <- sum(ch[(1:1000) + (i*1000), (4 + i)] > 0)/1000}
# mean(x)
#######################
## Now with gap years
# #Now set up data
no.resight.year <- c(116:130, 136:140)
ms.arr.missing <- ms.arr
ms.arr.missing[,no.resight.year] <- 0 #missing resight years
#ms.arr.missing[c(121:135, 141:145),] <- 0 #missing banding years
# or is it:
ms.arr.missing[c(116:130, 136:140),] <- 0 #missing banding years
# sum(ms.arr.missing[,no.resight.year])
# sum(ms.arr[,no.resight.year])
# sum(ms.arr.missing[c(121:135, 141:145),])
# sum(ms.arr[c(121:135, 141:145),])
#make jags data
jags.data <- list(marr = ms.arr.missing, n.occasions = ncol(ch.sim), rel = rowSums(ms.arr.missing), 
                  banding.year = c(1:115, 131:135, 141:165),
                  resight.index = as.numeric( 1:33 %in% c(1:23, 27, 29:33) ),
                  ns = ns,  
                  zero = matrix(0, ncol = ns, nrow = ns), 
                  ones = diag(ns))
time <- Sys.time()
out <- jags(jags.data, inits, parameters, "CMR.missing.jags", 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
Sys.time() - time
out$mean
plot_results(out)
plot(out$sims.list$alpha, out$sims.list$mean.phi0)
#try with different banding year and resight years missing:
no.resight.year <- c(116:130, 136:140)
ms.arr.missing <- ms.arr
ms.arr.missing[,no.resight.year] <- 0 #missing resight years
ms.arr.missing[c(121:135, 141:145),] <- 0 #missing banding years
#make jags data
jags.data <- list(marr = ms.arr.missing, n.occasions = ncol(ch.sim), rel = rowSums(ms.arr.missing), 
                  banding.year = c(1:120, 136:140, 146:165),
                  resight.index = as.numeric( 1:33 %in% c(1:23, 27, 29:33) ),
                  ns = ns,  
                  zero = matrix(0, ncol = ns, nrow = ns), 
                  ones = diag(ns))
time <- Sys.time()
out <- jags(jags.data, inits, parameters, "CMR.missing.jags", 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
Sys.time() - time
out$mean
plot_results(out)
#OK, strange results where p and phiA go high just before the missing years, but 
# the p estimate during missing years seem to be as expected, where the estimates 
# return the posterior mean and r.e. variance. Not the case for phiA where the estimate 
# is much too low.
# Try simulating data with no missing year, but variation in survival and p across years
x = rnorm(34, qlogis(0.8), 0.5)
y = rnorm(34, qlogis(0.4), 0.5)
ch <- ch_data_sim(phiA = plogis(x), phi0 = plogis(x - 2), p_adult = plogis(y))
ch.sim <- mcr(ch=ch)
ms.arr <- marray(ch.sim)
jags.data <- list(marr = ms.arr, n.occasions = ncol(ch.sim), rel = rowSums(ms.arr), 
                  banding.year = c(1:165),
                  resight.index = rep(1, 33), 
                  ns = ns,  
                  zero = matrix(0, ncol = ns, nrow = ns), 
                  ones = diag(ns))
time <- Sys.time()
out <- jags(jags.data, inits, parameters, "CMR.missing.jags", 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
Sys.time() - time
out$mean
plots <- plot_results(out)
plots[["PhiA"]]+geom_point(data = data.frame(Year = 1992:2024, PhiA = plogis(x)[-34]), aes(x=Year, y=PhiA, col="red"))
plots[["Phi0"]]+geom_point(data = data.frame(Year = 1992:2024, PhiA = plogis(x - 2)[-34]), aes(x=Year, y=PhiA, col="red"))
plots[["p"]]+geom_point(data = data.frame(Year = 1992:2024, PhiA = plogis(y)[-34]), aes(x=Year, y=PhiA, col="red"))
#That seems fairly good, now add in a block of missing data:
no.resight.year <- 10
ms.arr.missing <- ms.arr
ms.arr.missing[,(no.resight.year*ns - 4):(no.resight.year*ns)] <- 0 #missing resight years
ms.arr.missing[(no.resight.year*ns - 4):(no.resight.year*ns)+5,] <- 0 #missing banding years
jags.data <- list(marr = ms.arr, n.occasions = ncol(ch.sim), rel = rowSums(ms.arr.missing), 
                  banding.year = c(1:165),
                  resight.index = c(rep(1, 9), 0, rep(1, 23)), 
                  ns = ns,  
                  zero = matrix(0, ncol = ns, nrow = ns), 
                  ones = diag(ns))
time <- Sys.time()
out <- jags(jags.data, inits, parameters, "CMR.missing.jags", 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
Sys.time() - time
################################################################################
# fit model with a gap and monitor the product of phiA across the gap
ch <- ch_data_sim(n_occasions = 30, phi0 = rep(0.8, 30), d_marked = rep(1000, 30), a_marked = rep(1000, 30))
ch.sim <- mcr(ch=ch)
# make gap years
gap_calendar_years <- c(5, 6, 7) #missing re-sight in calendar years 5, 6, 7
#remove banding years
get.first <- function(x) min(which(x!=0)) #need indicator of banding year
first <- apply(ch.sim, 1, get.first)
missing <- first %in% gap_calendar_years
ch.sim <- ch.sim[!missing,] #drop missing banding years
ch.sim[, gap_calendar_years] <- 0
# resight.index is indexed by m-array recapture occasion = calendar year - 1
gap_occasions <- gap_calendar_years - 1
resight.index <- rep(1, ncol(ch.sim) - 1)
resight.index[gap_occasions] <- 0
#now make m-array
ms.arr <- marray(ch.sim) 
rowSums(ms.arr)[21:35] #check missing banding years sum to 0
which(colSums(ms.arr) == 0)    # should now equal columns 16:30 (occasions 4,5,6 x5 blocking)
which(resight.index == 0)      # should equal 4,5,6 (Option A)
parameters <- c("phiA", "phi0", "p", "beta.phi0", "alpha",
                "mean.phiA", "mean.logit.p", "mean.phi0", 
                "sigma.phiA", "sigma.p", "row.sum")
jags.data <- list(marr = ms.arr, n.occasions = ncol(ch.sim), rel = rowSums(ms.arr), 
                  resight.index = resight.index,  
                  ns = ns,  
                  zero = matrix(0, ncol = ns, nrow = ns), 
                  ones = diag(ns))
time <- Sys.time()
out <- jags(jags.data, inits, parameters, "CMR.missing.jags", 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
Sys.time() - time
out$mean
plot(out$mean$phiA, pch = 1, ylim = c(0, 1))
points(out$mean$phi0, pch = 2)
points(out$mean$p, pch = 3)
mean(apply(out$sims.list$phiA[,5:7], 1, prod))
#truth
0.8^3 #works!!!!!!!!
###############################
## Now try again with random year variation and gaps like the eider, as above:
x = rnorm(34, qlogis(0.8), 0.5)
y = rnorm(34, qlogis(0.4), 0.5)
ch <- ch_data_sim(phiA = plogis(x), phi0 = plogis(x - 2), p_adult = plogis(y))
ch.sim <- mcr(ch=ch)
# make gap years
gap_calendar_years <- c(24:26, 28) #missing re-sight in calendar years 24:26, 28
#remove banding years
get.first <- function(x) min(which(x!=0)) #need indicator of banding year
first <- apply(ch.sim, 1, get.first)
missing <- first %in% gap_calendar_years
ch.sim <- ch.sim[!missing,] #drop missing banding years
ch.sim[, gap_calendar_years] <- 0
# resight.index is indexed by m-array recapture occasion = calendar year - 1
gap_occasions <- gap_calendar_years - 1
resight.index <- rep(1, ncol(ch.sim) - 1)
resight.index[gap_occasions] <- 0
#now make m-array
ms.arr <- marray(ch.sim) 

parameters <- c("phiA", "phi0", "p", "beta.phi0", "alpha",
                "mean.phiA", "mean.logit.p", "mean.phi0", 
                "sigma.phiA", "sigma.p", "row.sum")
jags.data <- list(marr = ms.arr, n.occasions = ncol(ch.sim), rel = rowSums(ms.arr), 
                  resight.index = resight.index,  
                  ns = ns,  
                  zero = matrix(0, ncol = ns, nrow = ns), 
                  ones = diag(ns))
time <- Sys.time()
out <- jags(jags.data, inits, parameters, "CMR.missing.jags", 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
Sys.time() - time
plots <- plot_results(out)
plots[["PhiA"]]+geom_point(data = data.frame(Year = 1992:2024, PhiA = plogis(x)[-34]), aes(x=Year, y=PhiA, col="red"))
plots[["Phi0"]]+geom_point(data = data.frame(Year = 1992:2024, PhiA = plogis(x - 2)[-34]), aes(x=Year, y=PhiA, col="red"))
plots[["p"]]+geom_point(data = data.frame(Year = 1992:2024, PhiA = plogis(y)[-34]), aes(x=Year, y=PhiA, col="red"))
#WORKS!!!
# Now test with low variance in truth
x = rnorm(34, qlogis(0.8), 0.01)
y = rnorm(34, qlogis(0.4), 0.01)
y[c(10,15)] <- qlogis(rep(0.05, 2))
ch <- ch_data_sim(phiA = plogis(x), phi0 = plogis(x - 2), p_adult = plogis(y))
ch.sim <- mcr(ch=ch)
# make gap years
gap_calendar_years <- c(24:26, 28) #missing re-sight in calendar years 24:26, 28
#remove banding years
get.first <- function(x) min(which(x!=0)) #need indicator of banding year
first <- apply(ch.sim, 1, get.first)
missing <- first %in% gap_calendar_years
ch.sim <- ch.sim[!missing,] #drop missing banding years
ch.sim[, gap_calendar_years] <- 0
# resight.index is indexed by m-array recapture occasion = calendar year - 1
gap_occasions <- gap_calendar_years - 1
resight.index <- rep(1, ncol(ch.sim) - 1)
resight.index[gap_occasions] <- 0
#now make m-array
ms.arr <- marray(ch.sim) 
ns = 5
parameters <- c("phiA", "phi0", "p", "beta.phi0", "alpha",
                "mean.phiA", "mean.logit.p", "mean.phi0", 
                "sigma.phiA", "sigma.p")
# initial values
inits <- function(){list(
  mean.phi0 = 0.25, 
  mean.phiA = runif(1, 0.8, 0.9), 
  mean.p = runif(1, 0.5, 0.6), 
  sigma.phi0 = 0.5, 
  sigma.phiA = 0.5, 
  sigma.p = 0.5)}
# MCMC settings
ni <- 11000; nt <- 1; nb <- 5000; nc <- 3

jags.data <- list(marr = ms.arr, n.occasions = ncol(ch.sim), rel = rowSums(ms.arr), 
                  resight.index = resight.index,  
                  ns = ns,  
                  zero = matrix(0, ncol = ns, nrow = ns), 
                  ones = diag(ns))
time <- Sys.time()
out <- jags(jags.data, inits, parameters, "CMR.missing.jags", 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
Sys.time() - time
out$mean
plots <- plot_results(out)
plots[["PhiA"]]+geom_point(data = data.frame(Year = 1992:2024, PhiA = plogis(x)[-34]), aes(x=Year, y=PhiA, col="red"))
plots[["Phi0"]]+geom_point(data = data.frame(Year = 1992:2024, PhiA = plogis(x - 2)[-34]), aes(x=Year, y=PhiA, col="red"))
plots[["p"]]+geom_point(data = data.frame(Year = 1992:2025, PhiA = plogis(y)), aes(x=Year, y=PhiA, col="red"))
#Interest one-year lag dependancy in p after a very large low estimate? Is this sampling covariation?
plot(out$sims.list$p[,9], out$sims.list$p[,10]) #not sampling correlation. This is strange!
#small bias (< 0.01) in PhiA (low) and Phi0 (high), is this sampling correlation?
plot(out$sims.list$mean.phiA, out$sims.list$mean.phi0) #small negative correlation
cor(out$sims.list$mean.phiA, out$sims.list$mean.phi0)
# [1] -0.1814083
plot(out$sims.list$phiA[,9], out$sims.list$p[,10])
################################################################################
#Claude may be stumped; It seems this is likely statisitical confounding and not a bug; 
#  it suggests that my idea to vary p to high is a good test. Doing that here:
y[c(10,15)] <- qlogis(rep(0.95, 2))
ch <- ch_data_sim(phiA = plogis(x), phi0 = plogis(x - 2), p_adult = plogis(y))
ch.sim <- mcr(ch=ch)
# make gap years
gap_calendar_years <- c(24:26, 28) #missing re-sight in calendar years 24:26, 28
#remove banding years
get.first <- function(x) min(which(x!=0)) #need indicator of banding year
first <- apply(ch.sim, 1, get.first)
missing <- first %in% gap_calendar_years
ch.sim <- ch.sim[!missing,] #drop missing banding years
ch.sim[, gap_calendar_years] <- 0
# resight.index is indexed by m-array recapture occasion = calendar year - 1
gap_occasions <- gap_calendar_years - 1
resight.index <- rep(1, ncol(ch.sim) - 1)
resight.index[gap_occasions] <- 0
#now make m-array
ms.arr <- marray(ch.sim) 
jags.data <- list(marr = ms.arr, n.occasions = ncol(ch.sim), rel = rowSums(ms.arr), 
                  resight.index = resight.index,  
                  ns = ns,  
                  zero = matrix(0, ncol = ns, nrow = ns), 
                  ones = diag(ns))
time <- Sys.time()
out <- jags(jags.data, inits, parameters, "CMR.missing.jags", 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
Sys.time() - time
out$mean
plots <- plot_results(out)
plots[["PhiA"]]+geom_point(data = data.frame(Year = 1992:2024, PhiA = plogis(x)[-34]), aes(x=Year, y=PhiA, col="red"))
plots[["Phi0"]]+geom_point(data = data.frame(Year = 1992:2024, PhiA = plogis(x - 2)[-34]), aes(x=Year, y=PhiA, col="red"))
plots[["p"]]+geom_point(data = data.frame(Year = 1992:2025, PhiA = plogis(y)), aes(x=Year, y=PhiA, col="red"))
# it seem to be a really confounding, not a bug. At least that is our current hypothesis, which we feel is likely.
################################################################################
# To examine the practical effect of the confounding, I will set the detection to 
#  a less extreme values to see the effects
y[c(10,15)] <- qlogis(rep(0.3, 2))
ch <- ch_data_sim(phiA = plogis(x), phi0 = plogis(x - 2), p_adult = plogis(y))
ch.sim <- mcr(ch=ch)
# make gap years
gap_calendar_years <- c(24:26, 28) #missing re-sight in calendar years 24:26, 28
#remove banding years
get.first <- function(x) min(which(x!=0)) #need indicator of banding year
first <- apply(ch.sim, 1, get.first)
missing <- first %in% gap_calendar_years
ch.sim <- ch.sim[!missing,] #drop missing banding years
ch.sim[, gap_calendar_years] <- 0
# resight.index is indexed by m-array recapture occasion = calendar year - 1
gap_occasions <- gap_calendar_years - 1
resight.index <- rep(1, ncol(ch.sim) - 1)
resight.index[gap_occasions] <- 0
#now make m-array
ms.arr <- marray(ch.sim) 
jags.data <- list(marr = ms.arr, n.occasions = ncol(ch.sim), rel = rowSums(ms.arr), 
                  resight.index = resight.index,  
                  ns = ns,  
                  zero = matrix(0, ncol = ns, nrow = ns), 
                  ones = diag(ns))
time <- Sys.time()
out <- jags(jags.data, inits, parameters, "CMR.missing.jags", 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
Sys.time() - time
out$mean
plots <- plot_results(out)
plots[["PhiA"]]+geom_point(data = data.frame(Year = 1992:2024, PhiA = plogis(x)[-34]), aes(x=Year, y=PhiA, col="red"))
plots[["Phi0"]]+geom_point(data = data.frame(Year = 1992:2024, PhiA = plogis(x - 2)[-34]), aes(x=Year, y=PhiA, col="red"))
plots[["p"]]+geom_point(data = data.frame(Year = 1992:2025, PhiA = plogis(y)), aes(x=Year, y=PhiA, col="red"))
################################################################################
#  Fit an individual capture history model to test against for the structural "low p" issue. 
#  Reduced sample size and occasions to make comp time reasonable. 
## ============================================================================
## Individual-history (state-space) multistate CJS model
## Built to EXACTLY mirror the m-array model's biology and parameters,
## so it can serve as an independent cross-check.
##
## Only the likelihood differs: here we sample a latent true-state z[i,t]
## for every individual-occasion (Kery & Schaub 2012, ch. 9), instead of
## using summary m-array counts. Same states, transitions, survival and
## detection structure, and same priors as the m-array version.
##
## States (must match the m-array model):
##   1 = duckling (released, never resighted again as a duckling)
##   2 = 1-year-old            (unobservable)
##   3 = 2-year-old non-breeder(unobservable)
##   4 = 2-year-old breeder    (observable)
##   5 = adult (3+)            (observable)
##   6 = dead                  (added explicitly here; in the m-array it was
##                               the complement 1 - sum(survival))
## ============================================================================

library(jagsUI)   # or R2jags / rjags

## ----------------------------------------------------------------------------
## 1.  JAGS MODEL
## ----------------------------------------------------------------------------
## NOTE ON PARAMETER EQUIVALENCE WITH THE M-ARRAY MODEL:
##   phi1 = phi2 = phiA  (all share logit.phiA[t]), exactly as in the m-array.
##   phi0 is free here (its own mean + random effect), matching the model
##     version you used when you found the low-p phenomenon. To use the
##     offset link instead, replace the phi0 block with:
##       phi0[t] <- ilogit(logit.phiA[t] + beta.phi0)
##   p[t] carries resight.index[t] so structural gap years force p = 0,
##     exactly as po[.,t] <- p[t]*resight.index[t] did in the m-array.

cat("
model {

  # ---- Priors (identical to the m-array model) ----
  alpha ~ dbeta(1, 1)                    # or supply as data to fix it

  mean.phiA ~ dbeta(1, 1)
  sigma.phiA ~ dunif(0, 1)
  tau.phiA <- pow(sigma.phiA, -2)

  mean.phi0 ~ dbeta(1, 1)
  sigma.phi0 ~ dunif(0, 1)
  tau.phi0 <- pow(sigma.phi0, -2)

  mean.p ~ dbeta(1, 1)
  mean.logit.p <- logit(mean.p)
  sigma.p ~ dunif(0, 1)
  tau.p <- pow(sigma.p, -2)

  for (t in 1:(n.occasions - 1)) {
    logit.phiA[t] <- logit(mean.phiA) + eps.phiA[t]
    eps.phiA[t] ~ dnorm(0, tau.phiA)
    phiA[t] <- ilogit(logit.phiA[t])
    phi2[t] <- phiA[t]
    phi1[t] <- phiA[t]

    logit.phi0[t] <- logit(mean.phi0) + eps.phi0[t]
    eps.phi0[t] ~ dnorm(0, tau.phi0)
    phi0[t] <- ilogit(logit.phi0[t])
  }

  for (t in 1:n.occasions) {
    logit.p[t] <- mean.logit.p + eps.p[t]
    eps.p[t] ~ dnorm(0, tau.p)
    p[t] <- ilogit(logit.p[t]) * resight.index[t]
  }

  # ---- State-transition matrix PSI[state_t, t, state_t+1] ----
  # Rows = state at t; columns = state at t+1; includes dead state 6.
  # Complement of each alive row is death (1 - survival), matching the
  # m-array's implicit mortality.
  for (t in 1:(n.occasions - 1)) {
    PSI[1,t,1] <- 0
    PSI[1,t,2] <- phi0[t]
    PSI[1,t,3] <- 0
    PSI[1,t,4] <- 0
    PSI[1,t,5] <- 0
    PSI[1,t,6] <- 1 - phi0[t]

    PSI[2,t,1] <- 0
    PSI[2,t,2] <- 0
    PSI[2,t,3] <- phi1[t] * (1 - alpha)
    PSI[2,t,4] <- phi1[t] * alpha
    PSI[2,t,5] <- 0
    PSI[2,t,6] <- 1 - phi1[t]

    PSI[3,t,1] <- 0
    PSI[3,t,2] <- 0
    PSI[3,t,3] <- 0
    PSI[3,t,4] <- 0
    PSI[3,t,5] <- phi2[t]
    PSI[3,t,6] <- 1 - phi2[t]

    PSI[4,t,1] <- 0
    PSI[4,t,2] <- 0
    PSI[4,t,3] <- 0
    PSI[4,t,4] <- 0
    PSI[4,t,5] <- phiA[t]
    PSI[4,t,6] <- 1 - phiA[t]

    PSI[5,t,1] <- 0
    PSI[5,t,2] <- 0
    PSI[5,t,3] <- 0
    PSI[5,t,4] <- 0
    PSI[5,t,5] <- phiA[t]
    PSI[5,t,6] <- 1 - phiA[t]

    PSI[6,t,1] <- 0
    PSI[6,t,2] <- 0
    PSI[6,t,3] <- 0
    PSI[6,t,4] <- 0
    PSI[6,t,5] <- 0
    PSI[6,t,6] <- 1
  }

  # ---- Observation matrix PO[true_state, t, observed_event] ----
  # Observed events: 1..5 = seen in that state; 6 = not seen.
  # Only states 4 and 5 are observable (with prob p[t]); all else -> not seen.
  # (Observed state is known because age is known from banding year:
  #  a resighted duckling-banded bird at age 2 is a breeder (state 4);
  #  at age >= 3 it is an adult (state 5). See data-prep below.)
  for (t in 1:n.occasions) {
    PO[1,t,1] <- 0
    PO[1,t,2] <- 0
    PO[1,t,3] <- 0
    PO[1,t,4] <- 0
    PO[1,t,5] <- 0
    PO[1,t,6] <- 1

    PO[2,t,1] <- 0
    PO[2,t,2] <- 0
    PO[2,t,3] <- 0
    PO[2,t,4] <- 0
    PO[2,t,5] <- 0
    PO[2,t,6] <- 1

    PO[3,t,1] <- 0
    PO[3,t,2] <- 0
    PO[3,t,3] <- 0
    PO[3,t,4] <- 0
    PO[3,t,5] <- 0
    PO[3,t,6] <- 1

    PO[4,t,1] <- 0
    PO[4,t,2] <- 0
    PO[4,t,3] <- 0
    PO[4,t,4] <- p[t]
    PO[4,t,5] <- 0
    PO[4,t,6] <- 1 - p[t]

    PO[5,t,1] <- 0
    PO[5,t,2] <- 0
    PO[5,t,3] <- 0
    PO[5,t,4] <- 0
    PO[5,t,5] <- p[t]
    PO[5,t,6] <- 1 - p[t]

    PO[6,t,1] <- 0
    PO[6,t,2] <- 0
    PO[6,t,3] <- 0
    PO[6,t,4] <- 0
    PO[6,t,5] <- 0
    PO[6,t,6] <- 1
  }

  # ---- Likelihood: state process + observation process ----
  # z[i, first[i]] is supplied as DATA (known release state).
  for (i in 1:nind) {
    for (t in (first[i] + 1):n.occasions) {
      z[i, t] ~ dcat(PSI[z[i, t-1], t-1, ])
      y[i, t] ~ dcat(PO[z[i, t], t, ])
    }
  }

  # ---- Derived quantities for easy comparison to truth ----
  mean.phiA.out <- mean.phiA
  mean.phi0.out <- mean.phi0
}
", file = "CMR_individual.jags")


## ----------------------------------------------------------------------------
## 2.  SIMULATE DATA (small + fast, for the cross-check)
## ----------------------------------------------------------------------------
## Key point: the low-p confounding is driven by DESIGN STRUCTURE
## (a near-invisible year between normal ones), NOT by sample size. So shrink
## everything to make the slow individual model tractable in minutes:
##   - few occasions (enough to contain: normal years, one low-p year, and
##     a couple of years after it)
##   - modest numbers marked
## You only need to reproduce the PATTERN and its DIRECTION, not measure it
## precisely.

n.occasions <- 12

# --- low-p test: set true detection very low in years 5 and 6 ---
p_true <- rep(0.4, n.occasions)
p_true[c(5, 8)] <- 0.95           # <- the extreme low-detection years

# (For the REVERSE test, instead do: p_true[c(5,6)] <- 0.95)

phi0_true <- 0.30
phiA_true <- 0.80
alpha_true <- 0.33

# Use your existing simulator, but feed it the per-year p vector and small
# release sizes. (ch_data_sim already accepts p_adult as a vector.)
ch <- ch_data_sim(
  n_occasions = n.occasions,
  phi0        = rep(phi0_true, n.occasions),
  phi1        = rep(phiA_true, n.occasions),
  phi2        = rep(phiA_true, n.occasions),
  phiA        = rep(phiA_true, n.occasions),
  alpha       = alpha_true,
  p_adult     = p_true,
  d_marked    = rep(150, n.occasions),   # small!
  a_marked    = rep(150, n.occasions)    # small!
)

## ----------------------------------------------------------------------------
## 3.  CONVERT 0/1 HISTORIES  ->  OBSERVED-STATE MATRIX y  (and known z)
## ----------------------------------------------------------------------------
## Because age is known from the banding year, the observed state is
## deterministic given a sighting:
##   duckling-banded, seen at age 2      -> state 4 (breeder)
##   duckling-banded, seen at age >= 3   -> state 5 (adult)
##   adult-banded, seen                  -> state 5
## Age 0 = release (t == first), age 1 = unobservable.

CH <- as.matrix(ch[, grep("^Year", names(ch))])
is.duck <- ch$is_duckling
nind <- nrow(CH)

first <- apply(CH, 1, function(x) min(which(x != 0)))
init.state <- ifelse(is.duck == 1, 1, 5)

# Observed-event matrix: default 6 (= not seen)
y <- matrix(6, nrow = nind, ncol = n.occasions)
for (i in 1:nind) {
  f <- first[i]
  for (t in 1:n.occasions) {
    if (CH[i, t] == 1 && t > f) {          # a resighting (exclude the release)
      if (is.duck[i] == 1) {
        age <- t - f
        if (age == 2)      y[i, t] <- 4
        else if (age >= 3) y[i, t] <- 5
        # age 1 is never observable; if it appears, the sim is wrong
      } else {
        y[i, t] <- 5
      }
    }
  }
}

# Known latent states z (supplied as DATA): release state + observed states.
# NA everywhere else -> sampled by JAGS.
z.known <- matrix(NA, nrow = nind, ncol = n.occasions)
for (i in 1:nind) {
  z.known[i, first[i]] <- init.state[i]
  seen4 <- which(y[i, ] == 4); if (length(seen4)) z.known[i, seen4] <- 4
  seen5 <- which(y[i, ] == 5); if (length(seen5)) z.known[i, seen5] <- 5
}

## ----------------------------------------------------------------------------
## 4.  INITIAL VALUES FOR LATENT z  (must be a legal state path)
## ----------------------------------------------------------------------------
## Provide inits ONLY for latent (NA) entries; NA for the known ones.
## Use the age-implied alive trajectory (all-alive is a valid dcat path since
## every alive->alive transition used here has positive probability).
z.init.fun <- function() {
  zi <- matrix(NA, nrow = nind, ncol = n.occasions)
  for (i in 1:nind) {
    f <- first[i]
    if (f < n.occasions) {
      for (t in (f + 1):n.occasions) {
        if (is.na(z.known[i, t])) {
          if (is.duck[i] == 1) {
            age <- t - f
            zi[i, t] <- if (age == 1) 2 else if (age == 2) 4 else 5
          } else {
            zi[i, t] <- 5
          }
        }
      }
    }
  }
  zi
}

## ----------------------------------------------------------------------------
## 5.  BUNDLE + RUN
## ----------------------------------------------------------------------------
# resight.index: 1 where there IS resighting effort. For this test there are
# no STRUCTURAL gaps (the low-p years still have effort, just poor detection),
# so it is all ones. Set entries to 0 only to impose true structural gaps.
resight.index <- rep(1, n.occasions)

jags.data <- list(
  y             = y,
  z             = z.known,       # known states as data; NA = latent
  first         = first,
  nind          = nind,
  n.occasions   = n.occasions,
  resight.index = resight.index
  # To FIX alpha (as in earlier tests), add:  alpha = alpha_true
  # and delete the alpha ~ dbeta(1,1) prior line from the model.
)

inits <- function() list(
  z          = z.init.fun(),
  mean.phiA  = runif(1, 0.7, 0.9),
  mean.phi0  = runif(1, 0.2, 0.4),
  mean.p     = runif(1, 0.3, 0.5),
  sigma.phiA = runif(1, 0.05, 0.3),
  sigma.phi0 = runif(1, 0.05, 0.3),
  sigma.p    = runif(1, 0.05, 0.3),
  alpha      = runif(1, 0.25, 0.45)
)

params <- c("phiA", "phi0", "p", "alpha",
            "mean.phiA", "mean.phi0", "mean.p",
            "sigma.phiA", "sigma.phi0", "sigma.p")
start.time <- Sys.time()
out.ind <- jags(
  data = jags.data, inits = inits, parameters.to.save = params,
  model.file = "CMR_individual.jags",
  n.chains = 3, n.adapt = 1000, n.burnin = 2000, n.iter = 6000,
  parallel = TRUE
)
Sys.time() - start.time
#take a little more than 1 hr with the above samples sizes
plot(out.ind$mean$p)
#fairly clear there is no year-after effect. Although, the year after estimate both 
#  are slightly below the mean but well within other estimates. 
#Now run the reverse test and report results here in comments:
#  the after year are reversed and both are also slightly above the mean of other years, 
#  but within the range of normal years. Not super clear but suggests the statistical 
#  issue hypothesis is shared between models. 
#
#Simulate many (30) data sets and compare p and phiA across ind-based and m-array 
#  estimates across reps. Will save the p and phiA mean vectors, sigma.p and sigma.phiA
#  for each model type. Also need to save design metadata, code might be enough. 
#  Possible option thought is to repeat this over two samples size to examine 
#  shrinkage effects. Will hold off on this for now and use the moderate sample 
#  size (12 years, 100 individuals). 
nreps <- 40
band.size <- 150
n.occasions <- 12
p_true <- rep(0.4, n.occasions)
phi0_true <- 0.30
phiA_true <- 0.80
alpha_true <- 0.33
resight.index <- rep(1, n.occasions)

sim.results <- list( #list to hold results
  design.data = list("nreps" = nreps, "band.size" = band.size, 
                     "n.occasions" = n.occasions, 
                     phi0_true = phi0_true, 
                     phiA_true = phiA_true,
                     alpha_true = alpha_true, 
                     resight.index = resight.index),
  ind.model = vector("list", nreps), 
  ma.model =  vector("list", nreps)
  ) 
for (j in 1:nreps) {
  if ( j %% 2 != 0 ) p_true[c(6)] <- 0.05           # <- the extreme low-detection years
  if ( j %% 2 == 0 ) p_true[c(6)] <- 0.95           # <- the extreme high-detection years
  
  #sim data
  ch <- ch_data_sim(
    n_occasions = n.occasions,
    phi0        = rep(phi0_true, n.occasions),
    phi1        = rep(phiA_true, n.occasions),
    phi2        = rep(phiA_true, n.occasions),
    phiA        = rep(phiA_true, n.occasions),
    alpha       = alpha_true,
    p_adult     = p_true,
    d_marked    = rep(band.size, n.occasions),   # small!
    a_marked    = rep(band.size, n.occasions)    # small!
  )
  CH <- as.matrix(ch[, grep("^Year", names(ch))])
  is.duck <- ch$is_duckling
  nind <- nrow(CH)
  
  first <- apply(CH, 1, function(x) min(which(x != 0)))
  init.state <- ifelse(is.duck == 1, 1, 5)
  
  # Observed-event matrix: default 6 (= not seen)
  y <- matrix(6, nrow = nind, ncol = n.occasions)
  for (i in 1:nind) {
    f <- first[i]
    for (t in 1:n.occasions) {
      if (CH[i, t] == 1 && t > f) {          # a resighting (exclude the release)
        if (is.duck[i] == 1) {
          age <- t - f
          if (age == 2)      y[i, t] <- 4
          else if (age >= 3) y[i, t] <- 5
          # age 1 is never observable; if it appears, the sim is wrong
        } else {
          y[i, t] <- 5
        }
      }
    }
  }
  
  # Known latent states z (supplied as DATA): release state + observed states.
  # NA everywhere else -> sampled by JAGS.
  z.known <- matrix(NA, nrow = nind, ncol = n.occasions)
  for (i in 1:nind) {
    z.known[i, first[i]] <- init.state[i]
    seen4 <- which(y[i, ] == 4); if (length(seen4)) z.known[i, seen4] <- 4
    seen5 <- which(y[i, ] == 5); if (length(seen5)) z.known[i, seen5] <- 5
  }
  
  jags.data <- list(
    y             = y,
    z             = z.known,       # known states as data; NA = latent
    first         = first,
    nind          = nind,
    n.occasions   = n.occasions,
    resight.index = resight.index
    # To FIX alpha (as in earlier tests), add:  alpha = alpha_true
    # and delete the alpha ~ dbeta(1,1) prior line from the model.
  )
  
  inits <- function() list(
    z          = z.init.fun(),
    mean.phiA  = runif(1, 0.7, 0.9),
    mean.phi0  = runif(1, 0.2, 0.4),
    mean.p     = runif(1, 0.3, 0.5),
    sigma.phiA = runif(1, 0.05, 0.3),
    sigma.phi0 = runif(1, 0.05, 0.3),
    sigma.p    = runif(1, 0.05, 0.3),
    alpha      = runif(1, 0.25, 0.45)
  )
  
  params <- c("phiA", "phi0", "p", "mean.phiA", "mean.p", "sigma.phiA", "sigma.p")
  print(paste0("Starting run ", j))
  #fit individual model
  out.ind <- jags(
    data = jags.data, inits = inits, parameters.to.save = params,
    model.file = "CMR_individual.jags",
    n.chains = 3, n.adapt = 1000, n.burnin = 2000, n.iter = 6000,
    parallel = TRUE
  )
  #store results
  sim.results$ind.model[[j]] = list("p_true" = p_true, 
                                    "phiA" = out.ind$mean$phiA, 
                                    "p" = out.ind$mean$p, 
                                    "mean.phiA" = out.ind$mean$mean.phiA, 
                                    "mean.p" = out.ind$mean$mean.p, 
                                    "sigma.phiA" = out.ind$mean$sigma.phiA, 
                                    "sigma.p" = out.ind$mean$sigma.p)
  #prep data for m-array
  ns = 5
  ch.sim <- mcr(ch=ch)
  ms.arr <- marray(ch.sim) 
  jags.data <- list(marr = ms.arr, n.occasions = ncol(ch.sim), rel = rowSums(ms.arr), 
                    resight.index = resight.index,  
                    ns = ns,  
                    zero = matrix(0, ncol = ns, nrow = ns), 
                    ones = diag(ns))
  #fit m-array model
  out.ma <- jags(
    data = jags.data, inits = inits, parameters.to.save = params,
    model.file = "CMR.missing.jags", 
    n.chains = 3, n.adapt = 1000, n.burnin = 2000, n.iter = 6000,
    parallel = TRUE
  )
  #store results
  sim.results$ma.model[[j]] = list("p_true" = p_true, 
                                    "phiA" = out.ma$mean$phiA, 
                                    "p" = out.ma$mean$p, 
                                    "mean.phiA" = out.ma$mean$mean.phiA, 
                                    "mean.p" = out.ma$mean$mean.p, 
                                    "sigma.phiA" = out.ma$mean$sigma.phiA, 
                                    "sigma.p" = out.ma$mean$sigma.p)
} #loop over reps
saveRDS(sim.results, "sim.results.RDS")
sim.results <- readRDS("sim.results.RDS")
# summarize and plot results
df <- data.frame(NULL)
for (i in 1:nreps) {
  tmp <- cbind(data.frame(Model = "Individual", Dataset = i, Type = ifelse( i %% 2 != 0, "Low", "High")), 
               matrix(sim.results$ind.model[[i]]$p, 1, 12))
  df <- rbind(df, tmp)
  tmp <- cbind(data.frame(Model = "M-array", Dataset = i, Type = ifelse( i %% 2 != 0, "Low", "High")), 
               matrix(c(NA, sim.results$ma.model[[i]]$p), 1, 12))
  df <- rbind(df, tmp)
}
names(df)[4:15] <- paste0("p", 1:12)
df <- pivot_longer(df, 
                   cols = starts_with("p"), 
                   values_to = "Probability", 
                   names_to = "index",
                   names_pattern = "p(\\d+)")
df <- drop_na(df) |> mutate(index = as.numeric(index)) |>
  mutate(index = index - 6)
ggplot(data = df, aes(x = index, y = Probability, col = Model, group = Model)) + 
  geom_point(position = position_dodge(width = 0.5),) + 
  scale_x_continuous(breaks = -5:6, minor_breaks = NULL) + 
  scale_y_continuous(breaks = c(1:10)/10) + 
  geom_line(aes(col = Type, group = Dataset)) + 
  geom_segment(aes(x=-5, y=0.4, xend = -0.5, yend = 0.4), col = "black") + 
  geom_segment(aes(x=0.5, y=0.4, xend = 6, yend = 0.4), col = "black") +
  geom_segment(aes(x=-0.5, y=0.95, xend = 0.5, yend = 0.95), col = "black") + 
  geom_segment(aes(x=-0.5, y=0.05, xend = 0.5, yend = 0.05), col = "black")
#compute the paired difference between models
df2 <- df |> filter(index != 1) |>
  arrange(Dataset, Type, index, Model) |>
  group_by(Dataset, Type, index) |>
  summarise(Diff = first(Probability) - last(Probability)) |>
  ungroup() |>
  group_by(Type, index) |>
  summarise(meanDiff = mean(Diff), sdDiff = sd(Diff)) |>
  ungroup()
ggplot(data = df2, aes(x= index, y = meanDiff, col = Type)) + 
  geom_ribbon(aes(ymin = meanDiff - sdDiff, ymax = meanDiff + sdDiff, 
                  fill = Type, col = Type), alpha = 0.5) + 
  geom_line() + geom_hline(yintercept = 0) + 
  scale_x_continuous(breaks = -5:6, minor_breaks = NULL) + 
  labs(y="Mean difference between models")
#look at deviation from truth by type
mdf <- data.frame(index = rep(-5:6, times = 2), 
                  Type = c(rep("Low", 12), rep("High", 12)), 
                  Truth = c(sim.results$ind.model[[1]]$p_true, sim.results$ind.model[[2]]$p_true))
df2 <- df |> 
  left_join(mdf) |>
  mutate(Diff = Probability - Truth) |>
  group_by(index, Model, Type) |> 
  summarise(meanDiff = mean(Diff), sdDiff = sd(Diff), n = n(), seDiff = sdDiff/sqrt(n))
ggplot(data = df2, aes(x = index, y = meanDiff, col = Type)) + geom_line() + 
  geom_ribbon(aes(ymin = meanDiff - 2*seDiff, ymax = meanDiff + 2*seDiff, 
                  fill = Type, col = NULL), alpha = 0.3) + 
  geom_line(lwd = 1) + geom_hline(yintercept = 0) + 
  scale_x_continuous(breaks = -6:6, minor_breaks = NULL) + 
  labs(y="Mean difference from truth")
#### Survival -----------------------------
# summarize and plot results
df <- data.frame(NULL)
for (i in 1:nreps) {
  tmp <- cbind(data.frame(Model = "Individual", Dataset = i, Type = ifelse( i %% 2 != 0, "Low", "High")), 
               matrix(sim.results$ind.model[[i]]$phiA, 1, 11))
  df <- rbind(df, tmp)
  tmp <- cbind(data.frame(Model = "M-array", Dataset = i, Type = ifelse( i %% 2 != 0, "Low", "High")), 
               matrix(sim.results$ma.model[[i]]$phiA, 1, 11))
  df <- rbind(df, tmp)
}
names(df)[4:14] <- paste0("phiA", 1:11)
df <- pivot_longer(df, 
                   cols = starts_with("phiA"), 
                   values_to = "Survival", 
                   names_to = "index",
                   names_pattern = "phiA(\\d+)")
df <- drop_na(df) |> mutate(index = as.numeric(index)) |>
  mutate(index = index - 6)
ggplot(data = df, aes(x = index, y = Survival, col = Model, group = Model)) + 
  geom_point(position = position_dodge(width = 0.5),) + 
  scale_x_continuous(breaks = -5:6, minor_breaks = NULL) + 
  scale_y_continuous(breaks = c(1:10)/10) + 
  geom_line(aes(col = Type, group = Dataset)) + 
  geom_hline(yintercept = 0.8)
#compute the paired difference between models
df2 <- df |>
  arrange(Dataset, Type, index, Model) |>
  group_by(Dataset, Type, index) |>
  summarise(Diff = first(Survival) - last(Survival)) |>
  ungroup() |>
  group_by(Type, index) |>
  summarise(meanDiff = mean(Diff), sdDiff = sd(Diff)) |>
  ungroup()
ggplot(data = df2, aes(x= index, y = meanDiff, col = Type)) + 
  geom_ribbon(aes(ymin = meanDiff - sdDiff, ymax = meanDiff + sdDiff, 
                  fill = Type, col = Type), alpha = 0.5) + 
  geom_line() + geom_hline(yintercept = 0) + 
  labs(y="Mean difference in survival between models")
#look at deviation from truth by type
df2 <- df |> 
  mutate(Diff = Survival - 0.8) |>
  group_by(index, Model, Type) |> 
  summarise(meanDiff = mean(Diff), sdDiff = sd(Diff), n = n(), seDiff = sdDiff/sqrt(n))
ggplot(data = df2, aes(x = index, y = meanDiff, col = Type)) + geom_line() + 
  geom_ribbon(aes(ymin = meanDiff - 2*seDiff, ymax = meanDiff + 2*seDiff, 
                  fill = Type, col = NULL), alpha = 0.3) + 
  geom_line(lwd = 1) + geom_hline(yintercept = 0) + 
  scale_x_continuous(breaks = -6:6, minor_breaks = NULL) + 
  labs(y="Mean difference from truth for Survival")
## ----------------------------------------------------------------------------
## 6.  HOW TO COMPARE TO THE M-ARRAY MODEL  (read this!)
## ----------------------------------------------------------------------------
## INDEXING CORRESPONDENCE FOR p  (avoid a false off-by-one alarm):
##   Individual model: p[t] = detection at CALENDAR occasion t, t = 1..n.
##                     p[1] is unused (occasion 1 is release-only).
##   M-array model:    p[k] = detection at RECAPTURE occasion k = calendar
##                     year k+1,  k = 1..(n-1).
##   => m-array p[k]  corresponds to  individual p[k+1].
##   So the low-p years here (calendar 5,6) are m-array p[4], p[5].
##
## WHAT TO CHECK:
##   1. Does the individual model show the SAME pattern -- detection in the
##      year AFTER a low-p year biased low, with compensating survival?
##   2. Does the REVERSE (p_true[c(5,6)] <- 0.95) flip the sign, as it did
##      for the m-array?
##   If both models agree, the phenomenon lives in the shared statistics,
##   not in the m-array-specific code (U-recursion, cell probs, data hand-off).
##   That is exactly the code the individual model does NOT share, so agreement
##   rules it out as the source.
##
## Plot both for a quick visual:
##   plot(2:n.occasions, out.ind$mean$p[2:n.occasions], type="b",
##        ylim=c(0,1), xlab="calendar occasion", ylab="p")
##   abline(v = c(5,6), lty = 3)   # the low-p years; watch year 7