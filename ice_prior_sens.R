# Do prior sensitivity anaysis on ice survival coefficient
# Swap out the quadriatic effect prior with something that forces the relationship to be concave down
# use existing model, and gsub to to the change
# Currently in all models the prior is defined as:
#  #prior for covariates    
#  for (i in 1:2){
#    beta[i] ~ dnorm(0, 100) 
#  } # i
# This is SD = 0.1, which is not hugely informative for the quadratic
# Switch to beta[2] ~ dnorm(0, 100) T(, 0)
library(jagsUI)
# set up JAGS data
# initial values
inits <- function(){list(
  mean.phi0 = 0.25, 
  mean.phiA = runif(1, 0.8, 0.9), 
  mean.p = runif(1, 0.5, 0.6), 
  sigma.phi0 = 0.5, 
  sigma.phiA = 0.5, 
  sigma.p = 0.5)}
# parameters monitored
parameters <- c("Nb", "phiA", "phi0", "F", "mean.phi0", "mean.phiA", "mean.alpha", 
                "mean.log.F", "betaN", "beta", "sigma.o", "d", "mean.p", "p",
                "sigma.phi0", "sigma.phiA", "betaBP")
# MCMC settings
ni <- 40000; nt <- 1; nb <- 20000; nc <- 4
#-------------------------------------------------------------------------------
# First current conditions
jags.data <- readRDS("MS_Scenarios/jags.data.YKD.current.rds")
#read in model string
model_string <- paste(readLines("YKD_IPM_current.jags"), collapse = "\n")
old_prior <- "for (i in 1:2){\n  beta[i] ~ dnorm(0, 100)\n} # i"
new_prior <- "beta[1] ~ dnorm(0, 100)\nbeta[2] ~ dnorm(0, 100) T(,0)"
updated_model_string <- gsub(pattern = old_prior, 
                             replacement = new_prior, 
                             x = model_string, 
                             fixed = TRUE)
#check model
cat(updated_model_string)
#write to temp file
tmp_model_file <- tempfile(fileext = ".txt")
writeLines(updated_model_string, con = tmp_model_file)
#fit model
out <- jags(jags.data, inits, parameters, tmp_model_file, 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
saveRDS(out, file = "MS_Scenarios/YKD.current.ice.prior2.rds")
#-------------------------------------------------------------------------------
#Lead 2008, ice 4.5
jags.data <- readRDS("MS_Scenarios/jags.data.YKD.4.5.rds")
#read in model string
model_string <- paste(readLines("YKD_IPM_L2008.jags"), collapse = "\n")
old_prior <- "for (i in 1:2){\n  beta[i] ~ dnorm(0, 100)\n} # i"
new_prior <- "beta[1] ~ dnorm(0, 100)\nbeta[2] ~ dnorm(0, 100) T(,0)"
updated_model_string <- gsub(pattern = old_prior, 
                             replacement = new_prior, 
                             x = model_string, 
                             fixed = TRUE)
#check model
cat(updated_model_string)
#write to temp file
tmp_model_file <- tempfile(fileext = ".txt")
writeLines(updated_model_string, con = tmp_model_file)
#fit model
out <- jags(jags.data, inits, parameters, tmp_model_file, 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
saveRDS(out, file = "MS_Scenarios/YKD.L2008.4.5.ice.prior2.rds")
#-------------------------------------------------------------------------------
#Lead 2026, ice 4.5
jags.data <- readRDS("MS_Scenarios/jags.data.YKD.4.5.rds")
#read in model string
model_string <- paste(readLines("YKD_IPM_L2026.jags"), collapse = "\n")
old_prior <- "for (i in 1:2){\n  beta[i] ~ dnorm(0, 100)\n} # i"
new_prior <- "beta[1] ~ dnorm(0, 100)\nbeta[2] ~ dnorm(0, 100) T(,0)"
updated_model_string <- gsub(pattern = old_prior, 
                             replacement = new_prior, 
                             x = model_string, 
                             fixed = TRUE)
#check model
cat(updated_model_string)
#write to temp file
tmp_model_file <- tempfile(fileext = ".txt")
writeLines(updated_model_string, con = tmp_model_file)
#fit model
out <- jags(jags.data, inits, parameters, tmp_model_file, 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
saveRDS(out, file = "MS_Scenarios/YKD.L2026.4.5.ice.prior2.rds")
#-------------------------------------------------------------------------------
# lead constant, ice 4.5
jags.data <- readRDS("MS_Scenarios/jags.data.YKD.4.5.rds")
#read in model string
model_string <- paste(readLines("YKD_IPM_Lconst.jags"), collapse = "\n")
old_prior <- "for (i in 1:2){\n  beta[i] ~ dnorm(0, 100)\n} # i"
new_prior <- "beta[1] ~ dnorm(0, 100)\nbeta[2] ~ dnorm(0, 100) T(,0)"
updated_model_string <- gsub(pattern = old_prior, 
                             replacement = new_prior, 
                             x = model_string, 
                             fixed = TRUE)
#check model
cat(updated_model_string)
#write to temp file
tmp_model_file <- tempfile(fileext = ".txt")
writeLines(updated_model_string, con = tmp_model_file)
#fit model
out <- jags(jags.data, inits, parameters, tmp_model_file, 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
saveRDS(out, file = "MS_Scenarios/YKD.Lconst.4.5.ice.prior2.rds")
#-------------------------------------------------------------------------------
#Lead 2008, ice 8.5
jags.data <- readRDS("MS_Scenarios/jags.data.YKD.8.5.rds")
#read in model string
model_string <- paste(readLines("YKD_IPM_L2008.jags"), collapse = "\n")
old_prior <- "for (i in 1:2){\n  beta[i] ~ dnorm(0, 100)\n} # i"
new_prior <- "beta[1] ~ dnorm(0, 100)\nbeta[2] ~ dnorm(0, 100) T(,0)"
updated_model_string <- gsub(pattern = old_prior, 
                             replacement = new_prior, 
                             x = model_string, 
                             fixed = TRUE)
#check model
cat(updated_model_string)
#write to temp file
tmp_model_file <- tempfile(fileext = ".txt")
writeLines(updated_model_string, con = tmp_model_file)
#fit model
out <- jags(jags.data, inits, parameters, tmp_model_file, 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
saveRDS(out, file = "MS_Scenarios/YKD.L2008.8.5.ice.prior2.rds")
#-------------------------------------------------------------------------------
#Lead 2026, ice 8.5
jags.data <- readRDS("MS_Scenarios/jags.data.YKD.8.5.rds")
#read in model string
model_string <- paste(readLines("YKD_IPM_L2026.jags"), collapse = "\n")
old_prior <- "for (i in 1:2){\n  beta[i] ~ dnorm(0, 100)\n} # i"
new_prior <- "beta[1] ~ dnorm(0, 100)\nbeta[2] ~ dnorm(0, 100) T(,0)"
updated_model_string <- gsub(pattern = old_prior, 
                             replacement = new_prior, 
                             x = model_string, 
                             fixed = TRUE)
#check model
cat(updated_model_string)
#write to temp file
tmp_model_file <- tempfile(fileext = ".txt")
writeLines(updated_model_string, con = tmp_model_file)
#fit model
out <- jags(jags.data, inits, parameters, tmp_model_file, 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
saveRDS(out, file = "MS_Scenarios/YKD.L2026.8.5.ice.prior2.rds")
#-------------------------------------------------------------------------------
#Lead constant, ice 8.5
jags.data <- readRDS("MS_Scenarios/jags.data.YKD.8.5.rds")
#read in model string
model_string <- paste(readLines("YKD_IPM_Lconst.jags"), collapse = "\n")
old_prior <- "for (i in 1:2){\n  beta[i] ~ dnorm(0, 100)\n} # i"
new_prior <- "beta[1] ~ dnorm(0, 100)\nbeta[2] ~ dnorm(0, 100) T(,0)"
updated_model_string <- gsub(pattern = old_prior, 
                             replacement = new_prior, 
                             x = model_string, 
                             fixed = TRUE)
#check model
cat(updated_model_string)
#write to temp file
tmp_model_file <- tempfile(fileext = ".txt")
writeLines(updated_model_string, con = tmp_model_file)
#fit model
out <- jags(jags.data, inits, parameters, tmp_model_file, 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
saveRDS(out, file = "MS_Scenarios/YKD.Lconst.8.5.ice.prior2.rds")
#-------------------------------------------------------------------------------
# Now do ACP
#-------------------------------------------------------------------------------
#current conditions
jags.data <- readRDS("MS_Scenarios/jags.data.ACP.current.rds")
#read in model string
model_string <- paste(readLines("ACP_IPM.jags"), collapse = "\n")
old_prior <- "for (i in 1:2){\n  beta[i] ~ dnorm(0, 100)\n} # i"
new_prior <- "beta[1] ~ dnorm(0, 100)\nbeta[2] ~ dnorm(0, 100) T(,0)"
updated_model_string <- gsub(pattern = old_prior, 
                             replacement = new_prior, 
                             x = model_string, 
                             fixed = TRUE)
#check model
cat(updated_model_string)
#write to temp file
tmp_model_file <- tempfile(fileext = ".txt")
writeLines(updated_model_string, con = tmp_model_file)
#fit model
out <- jags(jags.data, inits, parameters, tmp_model_file, 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
saveRDS(out, file = "MS_Scenarios/ACP.current.ice.prior2.rds")
#-------------------------------------------------------------------------------
# RCP 4.5
jags.data <- readRDS("MS_Scenarios/jags.data.ACP.4.5.rds")
#read in model string
model_string <- paste(readLines("ACP_IPM.jags"), collapse = "\n")
old_prior <- "for (i in 1:2){\n  beta[i] ~ dnorm(0, 100)\n} # i"
new_prior <- "beta[1] ~ dnorm(0, 100)\nbeta[2] ~ dnorm(0, 100) T(,0)"
updated_model_string <- gsub(pattern = old_prior, 
                             replacement = new_prior, 
                             x = model_string, 
                             fixed = TRUE)
#check model
cat(updated_model_string)
#write to temp file
tmp_model_file <- tempfile(fileext = ".txt")
writeLines(updated_model_string, con = tmp_model_file)
#fit model
out <- jags(jags.data, inits, parameters, tmp_model_file, 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
saveRDS(out, file = "MS_Scenarios/ACP.4.5.ice.prior2.rds")
#-------------------------------------------------------------------------------
#RCP 8.5
jags.data <- readRDS("MS_Scenarios/jags.data.ACP.8.5.rds")
#read in model string
model_string <- paste(readLines("ACP_IPM.jags"), collapse = "\n")
old_prior <- "for (i in 1:2){\n  beta[i] ~ dnorm(0, 100)\n} # i"
new_prior <- "beta[1] ~ dnorm(0, 100)\nbeta[2] ~ dnorm(0, 100) T(,0)"
updated_model_string <- gsub(pattern = old_prior, 
                             replacement = new_prior, 
                             x = model_string, 
                             fixed = TRUE)
#check model
cat(updated_model_string)
#write to temp file
tmp_model_file <- tempfile(fileext = ".txt")
writeLines(updated_model_string, con = tmp_model_file)
#fit model
out <- jags(jags.data, inits, parameters, tmp_model_file, 
            n.chains = nc, n.burnin=nb, n.iter = ni,  
            parallel = TRUE, n.adapt = 1000)
saveRDS(out, file = "MS_Scenarios/ACP.8.5.ice.prior2.rds")
