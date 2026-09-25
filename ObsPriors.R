#Code to explore priors on latent detection/observation process, d_t in manuscript
# Question 9 in EE: annual variation in YKD detection rate
# aggregate responses were low = 0.05, median = 0.14, high = 0.29 for 10th, 50th, and 90th percentiles
library(rriskDistributions)
ee <- c(0.05, 0.14, 0.29) #annual variation in detection probability abs(dp), 
# interpreted as quantiles of the difference in detection from the average on abs scale
par <- get.beta.par(p=c(0.1, 0.5, 0.9), q=ee)
hist(rbeta(10000, par[[1]], par[[2]])) 
#this is the distribution for absolute value of differences in detection probability, (p_hat + d)
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
# N_t = y_t/(p + e_t) = y_t/[p*(1 + e_t/p)], so (1 + e_t/p) is the annual multiplicative effect 
#   as the deviation from the average p.
#   or on the log scale log(N_t) = log(y_t/p) - log(1 - e_t/p), 
#   so in the IPM we will use d_t = log(1 + e_t/p)
#
# One issue is extreme values of e_t or d_t, especially values near e_t = -p
#   This would cause detection --> 0, which is unrealistic and cause the population estimate to --> Inf
#   Therefore, we will set an upper bound on the deviation, 0.3 is consistent with the EE.
#
# Below is just exploratory and assumes p = 0.5 on average. This is close to the average across all years 
#   and is what we state in the EE.
#
x <- rbeta(100000, par[[1]], par[[2]])
x <- x[x<0.5]
hist(x) #this is what we elicited in question 9
qbeta(ee, par[[1]], par[[2]])
x2 <- c(0.5+x,0.5-x)
hist(x2, breaks=100) #this is the associated detection probability, assuming and average of 0.5
hist(1/x2, xlim=c(0, 10), breaks=10000)  #this is the VCF = 1/detection, note the dip at 2
#add a constraint < 0.3
x <- x[x < 0.3]
hist(x)
x2 <- c(0.5+x,0.5-x)
hist(x2, breaks=100) # p + e_t is between 0.2 and 0.8, and is almost never average
################################################################################
#try something simpler
#try a beta ranging between 0.2 and 0.8
qq <- quantile(x2, probs=c(0.1, 0.25, 0.5, 0.75, 0.9))
par2 <- get.beta.par(p=c(0.1, 0.25, 0.5, 0.75, 0.9), q=qq)
hist(rbeta(10000, par2[[1]], par2[[2]]))
hist(rbeta(10000, par2[[1]], par2[[2]]) * (0.8 - 0.2) + 0.2)
x <- rbeta(10000, par2[[1]], par2[[2]]) * (0.8 - 0.2) + 0.2
e_t <- x - 0.5
hist(e_t)
d <- 1 + e_t/0.5
hist(d)
hist(log(d))
############
#Explore the effects of d on the expected population size and counts:
counts <- read.csv("input_data/YKD_SPEI.csv", header = T)
Nsim <- 1000
dp <- (rbeta(length(counts$Nibb)*Nsim, par2[[1]], par2[[2]]) * (0.8 - 0.2) + 0.2) - 0.5
d <- log(1 + dp/0.5)
d <- matrix(d, Nsim, length(counts$Nibb))
n <- matrix(counts$Nibb, Nsim, length(counts$Nibb), byrow=TRUE)
jit <- 0.2
vcf <- matrix(counts$mvcf, Nsim, length(counts$Nibb), byrow=TRUE)
nibb <- exp(log(n) - d)
mnibb <- apply(nibb, 2, mean, na.rm=TRUE)
pnibb <- apply(nibb, 2, quantile, probs=c(0.025, 0.975), na.rm=TRUE)
plot(1:34-jit, counts$Nibb, pch=16, ylim=c(0, 30000))
points(1:34+jit, mnibb, pch=1)
arrows(x0=1:34+jit, y0=pnibb[1,], y1=pnibb[2,], length=0, col="lightgray")
arrows(x0=1:34-jit, y0=counts$Nibb-2*counts$seNibb, y1=counts$Nibb+2*counts$seNibb, length=0, col="lightgray")
# That looks reasonable, par2 ~= beta( 4, 4) * (0.8 - 0.2) + 0.2)
################################################################################
# for ACP
ee <- c(0.02, 0.08, 0.15) #annual variation in detection probability abs(dp), 
# interpreted as quantiles of the difference in detection from the average on abs scale
par <- get.beta.par(p=c(0.1, 0.5, 0.9), q=ee)
hist(rbeta(10000, par[[1]], par[[2]])) 
#for ACP assume mean detection is 0.75
x <- rbeta(100000, par[[1]], par[[2]])
x <- x[x<0.25]
hist(x) #this is what we elicited in question 9
qbeta(ee, par[[1]], par[[2]])
x2 <- c(0.75+x,0.75-x)
hist(x2, breaks=100) #this is the associated detection probability, assuming and average of 0.75
#try a beta ranging between 0.5 and 1.0
x3 <-  (x2 - 0.5) / (1 - 0.5) 
qq <- quantile(x3, probs=c(0.1, 0.25, 0.5, 0.75, 0.9))
par2 <- get.beta.par(p=c(0.1, 0.25, 0.5, 0.75, 0.9), q=qq)
hist(rbeta(10000, par2[[1]], par2[[2]]))
hist(rbeta(10000, par2[[1]], par2[[2]]) * (1 - 0.5) + 0.5)
x <- rbeta(10000, par2[[1]], par2[[2]]) * (1 - 0.5) + 0.5
e_t <- x - 0.75
hist(e_t)
d <- 1 + e_t/0.75
hist(d)
hist(log(d))
par2
################################################################################
# # other explorations below:
# ee2 <- c(0.05, 0.14, 0.29)*2 #annual variation in detection probability abs(dp), standardized to 0, 1 range 
# par2 <- get.beta.par(p=c(0.1, 0.5, 0.9), q=ee2)
# x <- rbeta(10000, par2[[1]], par2[[2]])/2 #set back to 0, 0.5 range
# hist(x)
# quantile(x, probs = c(0.1, 0.5, 0.9))
# # that looks good, but these are deviations on the probability scale, need them on the log scale
# # (1 + d/p) is the multiplicative factor 
# hist(1 + x/0.5) #positive deviations
# hist(1 - x/0.5) #negative deviations
# hist(c(1 + x/0.5, 1 - x/0.5), breaks = 100)
# hist(log(c(1 + x/0.5, 1 - x/0.5)), breaks = 100)
# #can this be approximated by a beta?
# qq <- quantile(c(1 + x/0.5, 1 - x/0.5)/2, probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
# par3 <- get.beta.par(p=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), q=qq)
# hist(rbeta(10000, par3[[1]], par3[[2]])*2)
# 
# #Explore the effects of d on the expected population size and counts:
# counts <- read.csv("input_data/YKD_SPEI.csv", header = T)
# Nsim <- 1000
# dp <- rbeta(length(counts$Nibb)*Nsim, par3[[1]], par3[[2]])
# d <- rnorm(length(counts$Nibb)*Nsim, 0, dp)
# d <- matrix(d, Nsim, length(counts$Nibb))
# n <- matrix(counts$Nibb, Nsim, length(counts$Nibb), byrow=TRUE)
# jit <- 0.2
# vcf <- matrix(counts$mvcf, Nsim, length(counts$Nibb), byrow=TRUE)
# nibb <- exp(log(n) - d)
# mnibb <- apply(nibb, 2, mean, na.rm=TRUE)
# pnibb <- apply(nibb, 2, quantile, probs=c(0.025, 0.975), na.rm=TRUE)
# plot(1:34-jit, counts$Nibb, pch=16, ylim=c(0, 30000))
# points(1:34+jit, mnibb, pch=1)
# arrows(x0=1:34+jit, y0=pnibb[1,], y1=pnibb[2,], length=0, col="lightgray")
# arrows(x0=1:34-jit, y0=counts$Nibb-2*counts$seNibb, y1=counts$Nibb+2*counts$seNibb, length=0, col="lightgray")
# #Now try with par2
# dp <- rbeta(length(counts$Nibb)*Nsim, par2[[1]], par2[[2]])
# dp <- ifelse(rbinom(length(counts$Nibb)*Nsim, size = 1, prob = 0.5) == 1, 1 + dp, 1 - dp)
# hist(dp, breaks = 100)
# d <- rnorm(length(counts$Nibb)*Nsim, 0, dp)
# d <- matrix(d, Nsim, length(counts$Nibb))
# n <- matrix(counts$Nibb, Nsim, length(counts$Nibb), byrow=TRUE)
# jit <- 0.2
# vcf <- matrix(counts$mvcf, Nsim, length(counts$Nibb), byrow=TRUE)
# nibb <- exp(log(n) - d)
# mnibb <- apply(nibb, 2, mean, na.rm=TRUE)
# pnibb <- apply(nibb, 2, quantile, probs=c(0.025, 0.975), na.rm=TRUE)
# plot(1:34-jit, counts$Nibb, pch=16, ylim=c(0, 30000))
# points(1:34+jit, mnibb, pch=1)
# arrows(x0=1:34+jit, y0=pnibb[1,], y1=pnibb[2,], length=0, col="lightgray")
# arrows(x0=1:34-jit, y0=counts$Nibb-2*counts$seNibb, y1=counts$Nibb+2*counts$seNibb, length=0, col="lightgray")
# #Will use the beta(par3[[1]], par3[[2]])*2 as the EE value seem way too much
# # even these this beta seems too much, maybe a beta 4,4?
# hist(rbeta(10000, 3, 3)*2)
