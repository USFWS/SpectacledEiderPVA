# using NB regression to model observed ice days as a function of model predicted
# ice days

# NEED TO UPDATE THIS WITH RIZZOLO'S data if we decide to switch over.
library(MASS)
library(tidyverse)


ext.sea.ice <- read.csv("input_data/extreme.sea.ice.csv", header = T)
head(ext.sea.ice)

ice.data <- ext.sea.ice[is.na(ext.sea.ice$ice.obs) == FALSE,]
head(ice.data)
var(ice.data$ice.obs) # 725.76
var(ice.data$ice.RCP4.5) # 141.03
var(ice.data$ice.RCP8.5) # 161.42

ice.data$max.day <- 181

### RCP4.5
# relationship between observed data and RCP4.5 data
plot(ice.obs~ ice.RCP4.5, data = ice.data) # there is no relationship at all
summary(m1 <- glm.nb(ice.obs~ offset(log(max.day)) + ice.RCP4.5, data = ice.data))
# theta = 2.768, s.e. = 0.672
# looks like I used the ice data before the most recent years data was appended to develop
# the theta numbers in the previous runs

#how do simulated pulls look when allowing theta to vary?
niter <- 1000
RCP4.5.sim <- matrix(NA, nrow = niter, ncol = nrow(ice.data))
for (i in 1: nrow(ice.data)){
 for (j in 1:niter){
   theta <- rnorm(1, mean = m1$theta, sd = m1$SE.theta) # mean and sd from nbregression
   RCP4.5.sim[j,i] <- rnbinom(1, size = theta, mu = ice.data$ice.RCP4.5[i])
   RCP4.5.sim[j,i] <- ifelse(RCP4.5.sim[j,i] > 180, 180, RCP4.5.sim[j,i])
   # note that draws > 180 can't be used...
 }
}


RCP4.5.sim <- data.frame(RCP4.5.sim)
RCP4.5.sim <- pivot_longer(RCP4.5.sim, cols = c(1:39), names_to = "case", values_to = "ice.days")
RCP4.5.sim$year <- as.factor(rep(c(1980:2018), 1000))

ggplot(RCP4.5.sim, aes(x = year, y = ice.days)) +
  geom_violin() + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=2, col = "red") +
  geom_point(data = ice.data, aes(x = as.factor(year), y = ice.RCP4.5), col = "red") +
  geom_point(data = ice.data, aes(x = as.factor(year), y = ice.obs), col = "blue") +
  theme(axis.text.x = element_blank())

### RCP8.5
# relationship between observed data and RCP8.5 data
plot(ice.obs~ ice.RCP8.5, data = ice.data) # there is no relationship at all
summary(m2 <- glm.nb(ice.obs~ offset(log(max.day)) + ice.RCP8.5, data = ice.data))
# theta = 2.777, s.e. = 0.675
# looks like I used the ice data before we added the most recent years info to develop
# the theta numbers in the previous runs

#how do simulated pulls look when allowing theta to vary?
niter <- 1000
RCP8.5.sim <- matrix(NA, nrow = niter, ncol = nrow(ice.data))
for (i in 1: nrow(ice.data)){
  for (j in 1:niter){
    theta <- rnorm(1, mean = m2$theta, sd = m2$SE.theta) # mean and sd from nbregression
    RCP8.5.sim[j,i] <- rnbinom(1, size = theta, mu = ice.data$ice.RCP8.5[i])
    RCP8.5.sim[j,i] <- ifelse(RCP8.5.sim[j,i] > 180, 180, RCP8.5.sim[j,i])
    # note that draws > 180 can't be used...
  }
}


RCP8.5.sim <- data.frame(RCP8.5.sim)
RCP8.5.sim <- pivot_longer(RCP8.5.sim, cols = c(1:39), names_to = "case", values_to = "ice.days")
RCP8.5.sim$year <- as.factor(rep(c(1980:2018), 1000))

ggplot(RCP8.5.sim, aes(x = year, y = ice.days)) +
  geom_violin() + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=2, col = "red") +
  geom_point(data = ice.data, aes(x = as.factor(year), y = ice.RCP4.5), col = "red") +
  geom_point(data = ice.data, aes(x = as.factor(year), y = ice.obs), col = "blue") +
  theme(axis.text.x = element_blank())