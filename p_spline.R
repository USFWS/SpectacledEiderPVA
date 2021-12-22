#write R and Jags code to implement p-splines based on description in Wood 2017,
#  pp. 204-206
library(mgcv)

#calculate basis, p. 205 in Wood.
bspline <- function(x, k, i, m=2)
{
  if(m==-1) {
    res <- as.numeric(x<k[i+1]&x>=k[i])
  } else {
    z0 <- (x-k[i])/(k[i+m+1]-k[i])
    z1 <- (k[i+m+2]-x)/(k[i+m+2]-k[i+1])
    res <- z0*bspline(x,k,i,m-1)+z1*bspline(x,k,i+1,m-1)
  }
  res
}

#try plotting the p-spline basis
# knot number, must be k + m + 2 with
# m + 1 the order of the basis, e.g., m = 2 for a cubic spline
# must define m+1 first and last knots outside the range of prediction
# see Wood, 2017, p. 204.
# ice range is 0 to 150, so we need to define knots < 0 and > 150
# if m = 2, then we need 3 knots < 0 and 3 > 150. 
# they should be equal spacing as the 'interior knots' in the range we care 
# about (i think this is needed)
k = 10 #number of knots in range of data we care about
m = 2 #order of basis
kp <- seq(0, 150, length=k)
kp <- c(-((m+1):1)*kp[2]+kp[1], kp, (1:(m+1))*kp[2]+kp[k]) #augment knots

x <- seq(0, 150, by=1)
f<- matrix(0, length(x), k)
for( i in 1:k){
  f[,i] <- bspline(x,k=kp, i=i, m=2)
}

beta <- diag(k) #matrix(c(0,0,0,0,1,0,0,0,0,0), k, 1) #rep(1, k), k, 1) #rnorm(n=k)
for(i in 1:k){
  if(i==1) plot(x, f%*%beta[,i], type="l", ylim=c(0, 1))
  lines(x, f%*%beta[,i])
}
#try random betas
for(i in 1:60){
beta <- matrix(rnorm(n=k), k, 1)
plot(x, f%*%beta, type="l", ylim=c(-3, 3))
Sys.sleep(1)
}

################################################################################
################################################################################
# now can we implement it in JAGS?
ext.sea.ice <- read.csv("input_data/extreme.sea.ice.csv", header = T)
# year is winter year (Nov - Apr)
ice.data <- ext.sea.ice[ext.sea.ice$year >= 1992 & ext.sea.ice$year <= 2100,] #replace to 2100 for projections, 2016 for data
ice.data4.5 <- ifelse(is.na(ice.data$ice.obs), ice.data$ice.RCP4.5, 
                      ice.data$ice.obs)
ice.data8.5 <- ifelse(is.na(ice.data$ice.obs), ice.data$ice.RCP8.5, 
                      ice.data$ice.obs)
ice.datacurrent <- ifelse(is.na(ice.data$ice.obs), mean(ice.data$ice.obs, na.rm = T), 
                          ice.data$ice.obs)
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
#strange behavior beyond range, lambda seen to return prior
#df=data.frame(y=c(y, 0, 0), ice=c(ice.data, min.ice, max.ice)) 
#so basis functions span this range.
gam.data <- mgcv::jagam(y~s(ice, k=5, bs='ps', m=c(2,2))-1, family=binomial, data = df, file="test_gam.txt")

cat(file = "p_spline.jags", "
model {
  eta <- X %*% b ## linear predictor
  for (i in 1:n) { mu[i] <-  ilogit(eta[i]) } ## expected response
  for (i in 1:n) { y[i] ~ dbin(mu[i],w[i]) } ## response 
  ## prior for s(ice)... 
  K1 <- S1[1:4,1:4] * lambda[1]  + S1[1:4,5:8] * lambda[2]
  b[1:4] ~ dmnorm(zero[1:4],K1) 
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
  #predict new ice
  for (i in 1:(K+1)){
    nb.prob[i] <- nb.size/(obs.ice[i] + nb.size)
    ice.draw[i] ~ dnegbin(nb.prob[i], nb.size)
  }
}
")

# bundle data
jags.data <- list(ice = ice.data, 
                  X = gam.data$jags.data$X[1:25,],
                  S1 = gam.data$jags.data$S1,
                  bzero = gam.data$jags.data$zero,
                  y=y,
                  nb.size=3.5)

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