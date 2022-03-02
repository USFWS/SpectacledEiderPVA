#Code to explore priors on latent detection/observation process, d_t in manuscript
# Question 9 in EE: annual variation in YKD detection rate
# aggregate responses were low = 0.05, median = 0.14, high = 0.29 for 10th, 50th, and 90th percentiles
library(rriskDistributions)
ee <- c(0.05, 0.14, 0.29) #annual variation in detection probability (dp), 
# interpreted as quantiles of the difference in detection from the average
# abs(dp) = abs(p_t - p_hat) or dropping the abs() 
# 1/dp = 1/(p_t - p_hat)
# log(dp) = - log(p_t - p_hat)

par <- get.beta.par(p=c(0.1, 0.5, 0.9), q=ee)
hist(rbeta(10000, par[[1]], par[[2]])) 
#this is the distribution for absolute value of differences in detection probability, (p_t - p_hat)
#Now translate that to log-scale d for IPM observation model:
# because the observed estimates are centered on mean detection (VCF, see Lewis 2019), 
# we need to translate the above to something consistent with deviations on the log-scale,
# Y ~ Normal(2b, sigma_obs), log(b) = log(mu) - o - d
# Y = breeding birds (males + females) from survey, sigma_obs is estimated from survey design, 
# b = breeding females from population model, o = observer effect estimated in IPM,
# d = annual multiplicative deviation in VCF or 1/detection (unobserved, detection = 1/VCF)
# o = average observer effect
# mu = model expected breeding females
# Note that Y has already been adjusted for average VCF, Y = Y_obs*VCF_bar; therefore
# in a year where the true VCF is > average (detection was low), Y_obs will be < average 
# and the VCF should have been adjusted up and the deviation is > 1. In other words, 
# we should have inflated Y_obs to a greater degree, causing Y to be greater than we measured. 
# But this is not observed, so we add the log offset to the linear predictor to adjust the expected 
# response for any given model-predicted expectation mu. This causes lower predictions of mu to be 
# favored in the draw of the MCMC chain. The observer effect, however, is estimated as an average 
# effect of an observer over the multiple years. 
#
# One the log scale, we add - log(1 + vcf[t]*d[t]) to the linear predictor, where vcf[t] is 
#  the average vcf for year t (averaged across the 4 eider strata, see Lewis 2019). 
#  vcf[t] is calculated as Nibb[t]/Y[t], where Nibb is from Lewis and Y is the observed index.  
#
#Below is just exploratory and assumes d = 0.5
x <- rbeta(100000, par[[1]], par[[2]])
x <- x[x<0.5]
hist(x) #this is what we elicited in question 9
x2 <- c(0.5+x,0.5-x)
hist(x2, breaks=100) #this is the associated detection probability, assuming and average of 0.5
hist(1/x2, xlim=c(0, 10), breaks=10000)  #this is the VCF = 1/detection, note the dip at 2
x3 <- 1/x2 - 2 #this is the deviation in the VCF: VCF - average VCF = 1/p - 1/0.5 = VCF - 2
hist(x3, xlim=c(-2, 10), breaks=10000) #note the dip at 0
x4 <- 1/(3*x2) #vcf = 2
hist(x4, xlim=c(0, 4), breaks=8000) #this is the multiplicative factor for the VCF
#log scale
hist(log(x4), xlim=c(0, 4), breaks=1000)

#find distributions, x4 > 1
quants <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
par2 <- get.gamma.par(p=quants, q=quantile(log(x4[x4>1]), probs=quants))
#x4 <= 1
hist(log(x4[x4<=1])) #range is log(0.5) < x4 <= 0
x5 <- -log(x4[x4<=1])/-log(0.5) #standardize to 0,1 interval
hist(x5)
par3 <- get.beta.par(p=quants, q=quantile(x5, probs=quants))

#now simulate
a <- rgamma(100000, par2[[1]], par2[[2]])
b <- rbeta(100000, par3[[1]], par3[[2]])*log(0.5)
x <- ifelse(rbinom(100000, 1, 0.5), a, b)
hist(x, xlim=c(-1, 4), breaks=100)
#there must be a smarter way!
#ignore the dip at zero? This is just the prior, after all
par4 <- get.gamma.par(p=quants, q=quantile(x-log(0.5), probs=quants))
hist(rgamma(10000, par4[[1]], par4[[2]]) + log(0.5) )

#Explore the effects of d on the expected population size and counts:
counts <- read.csv("input_data/YKD_SPEI.csv", header = T)
Nsim <- 1000
d <- rbeta(length(counts$Nibb)*Nsim, par[[1]], par[[2]])
d <- ifelse(runif(length(counts$Nibb)*Nsim)<0.5, 1, -1)*d
d <- matrix(d, Nsim, length(counts$Nibb))
n <- matrix(counts$Nibb, Nsim, length(counts$Nibb), byrow=TRUE)
jit <- 0.2
vcf <- matrix(counts$mvcf, Nsim, length(counts$Nibb), byrow=TRUE)
nibb <- exp(log(n) - log(1 + vcf*d))
mnibb <- apply(nibb, 2, mean, na.rm=TRUE)
pnibb <- apply(nibb, 2, quantile, probs=c(0.025, 0.975), na.rm=TRUE)
plot(1:34-jit, counts$Nibb, pch=16, ylim=c(0, 30000))
points(1:34+jit, mnibb, pch=1)
arrows(x0=1:34+jit, y0=pnibb[1,], y1=pnibb[2,], length=0, col="lightgray")
arrows(x0=1:34-jit, y0=counts$Nibb-2*counts$seNibb, y1=counts$Nibb+2*counts$seNibb, length=0, col="lightgray")
#Jensen's gap = exp(mu)*(exp(vd/2) - 1)
vd <- apply( - log(1 + vcf*d), 2, var, na.rm=TRUE)
jg <- counts$Nibb*(exp(vd/2) - 1)
points(1:34+jit, mnibb - jg, pch=2)
