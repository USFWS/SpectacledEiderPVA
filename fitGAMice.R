# examine the effect of ice on survival using a GAM
# GAM based on annual estimate from Bayesian CMR model, CMR.jags/CMR_YKD.R
library(mgcv)
library(jagsUI)
library(dplyr)
out <- readRDS(file="CMR.GAM.rds")

plot(1:23, out$mean$phiA, ylim=c(0.5, 1))

#load ice data
ext.sea.ice <- read.csv("input_data/extreme.sea.ice.csv", header = T)
# year is winter year (Nov - Apr)
ice.data <- ext.sea.ice[ext.sea.ice$year >= 1992 & ext.sea.ice$year <= 2016,]
ice.data <- ice.data$ice.obs
# ice.data <- (ice.data-mean(ice.data)) / sd(ice.data)
# ice.data <- ice.data[-c(1,25)]
#plot survival verses ice
plot(ice.data[-c(1,25)], out$mean$phiA)

fit <- gam(phiA~s(ice), data=data.frame(phiA=out$mean$phiA, ice=ice.data[-c(1,25)]))
plot(fit, residuals = TRUE, cex=5)

#now add minimum ice
min.ice <- read.csv("input_data/minimal.sea.ice.csv", header = T) %>%
  filter(year >= 1992 & year <= 2016) %>%
  select(ice.obs)
fit2 <-  gam(phiA~s(ice.min), data=data.frame(phiA=out$mean$phiA, ice.min=min.ice$ice.obs[-c(1,25)]))
plot(fit2, residuals = TRUE, cex=5)

#add both to same model
plot(ice.data, min.ice$ice.obs)
# looks too correlated
cor(ice.data, min.ice$ice.obs)
#[1] -0.8837191
cov(data.frame(ice.data, min.ice$ice.obs))
eigen(cov(data.frame(ice.data, min.ice$ice.obs)))
df = data.frame(phiA=out$mean$phiA, 
                ice.min=min.ice$ice.obs[-c(1,25)],
                ice.ext=ice.data[-c(1,25)])
fit3 <- gam(phiA~s(ice.min, ice.ext, k=5), data=df)
vis.gam(fit3, plot.type = "contour", labcex=1)
points(df$ice.min,df$ice.ext, pch=16)

#fit separately
fit4 <- gam(phiA~s(ice.min, k=5)+s(ice.ext, k=5), data=df)
plot(fit4)

fitlm <- lm(phiA~ice.min+ice.ext+I(ice.ext^2), data=df)

AIC(fit, fit2, fit3, fit4)

#Now look at phi0
df <- data.frame(df, phi0=out$mean$phi0)
plot(df$ice.ext, df$phi0, pch=16)
plot(df$ice.min, df$phi0, pch=16)
fit <- gam(phi0~s(ice.ext), data=df)
fit2 <-  gam(phi0~s(ice.min), data=df)
fit3 <- gam(phi0~s(ice.min, ice.ext, k=5), data=df)
fit4 <- gam(phi0~s(ice.min, k=5)+s(ice.ext, k=5), data=df)
AIC(fit, fit2, fit3, fit4)
vis.gam(fit3, plot.type = "contour", labcex=1)
points(df$ice.min,df$ice.ext, pch=16)

#fit them together
library(tidyr)
df <- df %>% pivot_longer(cols=c(1,4), names_to="Age", values_to="Phi")
fit <- gam(Phi~factor(Age)+s(ice.ext), data=df)
fit2 <-  gam(Phi~factor(Age)+s(ice.min), data=df)
fit3 <- gam(Phi~factor(Age)+s(ice.min, ice.ext, k=5), data=df)
fit4 <- gam(Phi~factor(Age)+s(ice.min, k=5)+s(ice.ext, k=5), data=df)
fitlm <- lm(Phi~Age+ice.min+ice.ext+I(ice.ext^2), data=df)
AIC(fit, fit2, fit3, fit4, fitlm)
plot(fit2)
vis.gam(fit3, plot.type = "contour", labcex=1, view=c("ice.min", "ice.ext"))
points(df$ice.min,df$ice.ext, pch=16)
plot(fit4)
vcov(fitlm)

#put it on logit scale and center
df <- df %>% 
  mutate(logitPhi=qlogis(Phi), ice.min.s=scale(ice.min), ice.ext.s=scale(ice.ext))
fitlm <- lm(logitPhi~Age+ice.min.s+ice.ext.s+I(ice.ext.s^2), data=df)
summary(fitlm)
vcov(fitlm)

#Use PCA to make orthogonal covariates between min and ext ice
# first find PCA based on all ice data
ext.sea.ice <- read.csv("input_data/extreme.sea.ice.csv", header = T)
min.sea.ice <- read.csv("input_data/minimal.sea.ice.csv", header = T)
plot(ext.sea.ice$ice.obs, min.ice$ice.obs, pch=16, ylim=c(0, 200))
points(ext.sea.ice$ice.RCP4.5, min.ice$ice.RCP4.5, pch=1)
points(ext.sea.ice$ice.RCP8.5, min.ice$ice.RCP8.5, pch=2)
df2 <- ext.sea.ice %>% pivot_longer(cols=2:4, names_to = "type", values_to="ext") %>%
  filter(!is.na(ext))
df3 <- min.sea.ice %>% pivot_longer(cols=2:4, names_to = "type", values_to="min") %>%
  filter(!is.na(min)) %>%
  left_join(df2, by=c("year", "type")) 
df2 <- df3
rm(df3)
pca <- princomp(formula=~ext+min, data=df2)
plot(df2$ext, df2$min, pch=16)
ext1 <- c(0, 80)
min1 <- (pca$loadings[1,1]/pca$loadings[2,1])*ext1 + pca$center[2]
lines(ext1, min1)
