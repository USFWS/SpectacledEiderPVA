# using NB regression to model observed ice days as a function of model predicted
# ice days


library(MASS)
library(tidyverse)

set.seed(89)

# note that this data was updated from D.Rizzolo to include recent years (19-20) and 
# revised numbers from 1980-2018
ext.sea.ice <- read.csv("input_data/extreme.sea.ice.csv", header = T)
head(ext.sea.ice)

ice.data <- ext.sea.ice[is.na(ext.sea.ice$ice.obs) == FALSE,]
head(ice.data)

# looking at the variances in the different time series
var(ice.data$ice.obs) # 627.61
var(ice.data$ice.RCP4.5) # 150.82
var(ice.data$ice.RCP8.5) # 166.94
# issue is that RCP modelled ice data has much less noise than the realized 
# observations (to some degree a consequence of the RCP time series each being
# a weighted average of 8 climate models)

# looking to add uncertainty to the RCP data series to reflect the uncertainty
# in the ice covariates in the future

ice.data$max.day <- 181

### RCP4.5
# relationship between observed data and RCP4.5 data
plot(ice.obs~ ice.RCP4.5, data = ice.data) # there is no relationship at all
abline(a = 0, b = 1)

summary(m1 <- glm.nb(ice.obs~ offset(log(max.day)) + ice.RCP4.5, data = ice.data))
# theta = 2.159, s.e. = 0.504
# looks like I used the ice data before the most recent years data was appended 
# (up to 2016) to develop the theta numbers in the previous runs

# how do simulated pulls look when allowing theta to vary?
niter <- 1000
RCP4.5.sim <- matrix(NA, nrow = niter, ncol = nrow(ext.sea.ice))
for (i in 1: nrow(ext.sea.ice)){
 for (j in 1:niter){
   # theta <- rnorm(1, mean = m1$theta, sd = m1$SE.theta) # mean and sd from nbregression
   RCP4.5.sim[j,i] <- rnbinom(1, size = 4, mu = ext.sea.ice$ice.RCP4.5[i])
   RCP4.5.sim[j,i] <- ifelse(RCP4.5.sim[j,i] > 180, 180, RCP4.5.sim[j,i])
   # note that draws > 180 can't be used...
 }
}

# variance of the draws over the observed time period (1980 - 2020)
RCP4.5.var.sim <- apply(RCP4.5.sim[,c(30:70)], 1, var)
hist(RCP4.5.var.sim)
sum(RCP4.5.var.sim>var(ice.data$ice.obs))/1000 # 57%

plot(c(ice.data$ice.obs, rep(NA,80)) ~ c(1980:2100), ylim = c(0,181))
for (i in 1:50){
  lines(RCP4.5.sim[i,c(30:150)] ~ c(1980:2100), col = i)
}
points(ice.data$ice.obs~ c(1980:2020), pch = 16)
points(ext.sea.ice$ice.RCP4.5[30:150] ~ c(1980:2100), pch = 16, col = "red")

RCP4.5.sim <- data.frame(RCP4.5.sim)
RCP4.5.sim.obs <- RCP4.5.sim[,c(30:70)]
RCP4.5.sim.obs <- RCP4.5.sim.obs %>%
  pivot_longer(cols = starts_with("X"),
               names_to = "year_count", names_prefix = "wk", values_to = "ice.days")
  
RCP4.5.sim.obs$year <- as.factor(rep(c(1980:2020), 1000))

ggplot(RCP4.5.sim.obs, aes(x = year, y = ice.days)) +
  geom_violin() + 
  stat_summary(fun=mean, geom="point", shape=23, size=2, col = "red") +
  geom_point(data = ice.data, aes(x = as.factor(year), y = ice.RCP4.5), col = "red") +
  geom_point(data = ice.data, aes(x = as.factor(year), y = ice.obs), col = "blue") +
  theme(axis.text.x = element_blank())

### RCP8.5
# relationship between observed data and RCP8.5 data
plot(ice.obs~ ice.RCP8.5, data = ice.data) # there is no relationship at all
abline(a = 0, b = 1)

summary(m1 <- glm.nb(ice.obs~ offset(log(max.day)) + ice.RCP8.5, data = ice.data))
# theta = 2.153, s.e. = 0.502
# very close to the RCP4.5 results. should be able to use the same theta for both?

# how do simulated pulls look when allowing theta to vary?
niter <- 1000
RCP8.5.sim <- matrix(NA, nrow = niter, ncol = nrow(ext.sea.ice))
for (i in 1: nrow(ext.sea.ice)){
  for (j in 1:niter){
    # theta <- rnorm(1, mean = m1$theta, sd = m1$SE.theta) # mean and sd from nbregression
    RCP8.5.sim[j,i] <- rnbinom(1, size = 4, mu = ext.sea.ice$ice.RCP8.5[i])
    RCP8.5.sim[j,i] <- ifelse(RCP8.5.sim[j,i] > 180, 180, RCP8.5.sim[j,i])
    # note that draws > 180 can't be used...
  }
}

# variance of the draws over the observed time period (1980 - 2020)
RCP8.5.var.sim <- apply(RCP8.5.sim[,c(30:70)], 1, var)
hist(RCP8.5.var.sim)
sum(RCP8.5.var.sim>var(ice.data$ice.obs))/1000 # 54%

plot(c(ice.data$ice.obs, rep(NA,80)) ~ c(1980:2100), ylim = c(0,181))
for (i in 1:50){
  lines(RCP8.5.sim[i,c(30:150)] ~ c(1980:2100), col = i)
}
points(ice.data$ice.obs~ c(1980:2020), pch = 16)
points(ext.sea.ice$ice.RCP8.5[30:150] ~ c(1980:2100), pch = 16, col = "red")

RCP8.5.sim <- data.frame(RCP8.5.sim)
RCP8.5.sim.obs <- RCP8.5.sim[,c(30:70)]
RCP8.5.sim.obs <- RCP8.5.sim.obs %>%
  pivot_longer(cols = starts_with("X"),
               names_to = "year_count", names_prefix = "wk", values_to = "ice.days")

RCP8.5.sim.obs$year <- as.factor(rep(c(1980:2020), 1000))

ggplot(RCP8.5.sim.obs, aes(x = year, y = ice.days)) +
  geom_violin() + 
  stat_summary(fun=mean, geom="point", shape=23, size=2, col = "red") +
  geom_point(data = ice.data, aes(x = as.factor(year), y = ice.RCP8.5), col = "red") +
  geom_point(data = ice.data, aes(x = as.factor(year), y = ice.obs), col = "blue") +
  theme(axis.text.x = element_blank())