# find prior for sigma d from EE
# log scale, d is the log deviation in the 1/detection or the 'VCF' from the average detection or VCF. 
ee.dat <- data.frame(Question=rep("Q10", 4), Lowest=c(0.01, 0.1, 0.025, 0.03), Best=c(0.05, 0.125, 0.05, 0.1), Highest=c(0.12, 0.25, 0.1, 0.17), Confidence=c(0.80, 0.50, 0.85, 0.60))
head(ee.dat)
ee.dat$Question <- as.factor(ee.dat$Question)

#standardize to average of two observer, make Expert response proportional to 1/sqrt(2)
ee.dat[,2:4] <- ee.dat[,2:4]/sqrt(2)
head(ee.dat)

# create standardized 80% bounds using equations in Hemmings et al. (2018)
s.CI <- 0.8
ee.dat$Low80 <- ee.dat$Best - ((ee.dat$Best - ee.dat$Lowest)*(s.CI/ee.dat$Confidence))
ee.dat$Hi80 <- ee.dat$Best + ((ee.dat$Highest - ee.dat$Best)*(s.CI/ee.dat$Confidence))
# truncating where necessary
ee.dat$Low80 <- ifelse(ee.dat$Low80 < 0, 0, ee.dat$Low80)
ee.dat$Hi80 <- ifelse(ee.dat$Hi80 > 0.25, 0.25, ee.dat$Hi80)

# Q10
#Best
round(mean(ee.dat$Best, na.rm = T), 2)
#[1] 0.06
#Low
round(mean(ee.dat$Low80, na.rm = T), 2)
#[1] 0.02
#High
round(mean(ee.dat$Hi80, na.rm = T), 2)
#[1] 0.13
#find prior from EE response
# log scale!
a=3
b=10
#x <- abs(exp(rnorm(10000, 0, rgamma(10000, a,b))) - 1)
#sigma_d <- runif(10000, 0.02, 0.13)
#sigma_d <- rbeta(10000, 6, 94)*(0.13-0.02) + 0.02
sigma_d <- rgamma(10000, 1, 13) # or use rgamma(10000, 1, 13)? results not sensitive
hist(sigma_d)
mean(sigma_d)
quantile(sigma_d, probs=c(0.1, 0.2, 0.5, 0.8, 0.9))

x <- rnorm(10000, 0, sigma_d)
hist(x)
hist(0.75*exp(x) - 0.75)
quantile(0.75*exp(x) - 0.75, probs=c(0.1, 0.2, 0.5, 0.8, 0.9))
#for YKD
sigma_d <- rgamma(10000, 2, 13) # or use rgamma(10000, 1, 13)? results not sensitive
hist(sigma_d)
mean(sigma_d)
quantile(sigma_d, probs=c(0.1, 0.2, 0.5, 0.8, 0.9))
x <- rnorm(10000, 0, sigma_d)
hist(x)
hist(0.5*exp(x) - 0.5)
sd(0.5*exp(x) - 0.5)
quantile(0.5*exp(x) - 0.5, probs=c(0.1, 0.2, 0.5, 0.8, 0.9))

#for sigma.o
a=1
b=10
x <- abs(exp(rnorm(10000, 0, rgamma(10000, a,b))) - 1)
mean(x)
hist(x)
quantile(x, probs=c(0.1, 0.5, 0.9))
