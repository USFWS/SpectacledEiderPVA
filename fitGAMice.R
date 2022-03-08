# examine the effect of ice on survival using a GAM
# GAM based on annual estimate from Bayesian CMR model, CMR.jags/CMR_YKD.R
library(mgcv)
library(jagsUI)
library(tidyverse)
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
plot(ext.sea.ice$ice.obs, min.sea.ice$ice.obs, pch=16, ylim=c(0, 200))
points(ext.sea.ice$ice.RCP4.5, min.sea.ice$ice.RCP4.5, pch=1)
points(ext.sea.ice$ice.RCP8.5, min.sea.ice$ice.RCP8.5, pch=2)
df2 <- ext.sea.ice %>% pivot_longer(cols=2:4, names_to = "type", values_to="ext") %>%
  filter(!is.na(ext))
df3 <- min.sea.ice %>% pivot_longer(cols=2:4, names_to = "type", values_to="min") %>%
  filter(!is.na(min)) %>%
  left_join(df2, by=c("year", "type")) 
df2 <- df3
rm(df3)
pca <- princomp(formula=~ext+min, data=df2)
plot(df2$ext, df2$min, pch=1)
x=100
arrows(x0=pca$loadings[1,1]*(-x) + pca$center[1], x1=pca$loadings[1,1]*x + pca$center[1], 
       y0=pca$loadings[2,1]*(-x) + pca$center[2], y1=pca$loadings[2,1]*x + pca$center[2],
       length=0)
x=10
arrows(x0=pca$loadings[1,2]*(-x) + pca$center[1], x1=pca$loadings[1,2]*x + pca$center[1], 
       y0=pca$loadings[2,2]*(-x) + pca$center[2], y1=pca$loadings[2,2]*x + pca$center[2],
       length=0)

points(pca$center[1], pca$center[2], pch=16, col="red")

plot(pca$scores[,1], pca$scores[,2])
pca$loadings
#compute pca scores from scratch
scores <- scale(as.matrix(df2[,c("ext", "min")]), scale=FALSE)%*%pca$loadings


#fit linear model to pca scores
df <- df2 %>% 
  bind_cols(data.frame(pc1=pca$scores[,1], pc2=pca$scores[,2])) %>%
  filter(year > 1992 & year < 2016, type == "ice.obs") %>%
  bind_cols(data.frame(PhiA=out$mean$phiA), data.frame(Phi0=out$mean$phi0)) %>%
  pivot_longer(cols=c(7,8), names_to="Age", values_to="Phi") %>% 
  mutate(logitPhi=qlogis(Phi), ice.min.s=scale(min), ice.ext.s=scale(ext)) 
ggplot(data=df) + geom_point(aes(x=ext, y=min))
ggplot(data=df) + geom_point(aes(x=pc1, y=pc2, col="red"))
ggplot(data = df) + 
  geom_point(aes(x=pc1, y=Phi, col=Age)) + 
  geom_smooth(aes(x=pc1, y=Phi, col=Age), method="gam")
ggplot(data = df) + 
  geom_point(aes(x=pc2, y=Phi, col=Age)) + 
  geom_smooth(aes(x=pc2, y=Phi, col=Age), method="gam")

fitgam <- gam(logitPhi~Age+s(pc1, pc2, k=10), data=df)
summary(fitgam)
vis.gam(fitgam, view=c("pc1", "pc2"), plot.type = "contour")
points(df$pc1, df$pc2)

fitlm1 <- lm(logitPhi~Age+pc1+I(pc1^2), data=df)
fitlm2 <- lm(logitPhi~Age+pc1+pc2+I(pc1^2), data=df)
fitlm3 <- lm(logitPhi~Age+pc1+pc2+I(pc1^2) + I(pc2^2), data=df)
fitlm4 <- lm(logitPhi~Age+pc1+pc2+I(pc1^2) + I(pc2^2) + I(pc1^2):I(pc2^2), data=df)

lapply(list(fitlm1, fitlm2, fitlm3, fitlm4), summary)
AIC(fitlm1, fitlm2, fitlm3, fitlm4)
vcov(fitlm1)

#Is the pca any better than just ext ice?
fitlm5 <- lm(logitPhi~Age+ext+I(ext^2), data=df)
fitlm6 <- lm(logitPhi~Age+min+I(min^2), data=df)
fitlm7 <- lm(logitPhi~Age+ext+min+I(ext^2) + I(min^2), data=df)
AIC(fitlm1, fitlm5, fitlm6, fitlm7)

#ggplots
ggplot(data=df2) + 
  geom_point(data=df, aes(x=ext, y=min, shape="open circle", size=1)) + 
  geom_point(aes(x=ext, y=min, col=type)) + 
  geom_text(data = df2, aes(x=ext, y=min, label=year), 
    nudge_x = 0.25, nudge_y = 0.25, 
    check_overlap = T) + 
  geom_segment(aes(x=0, y=180, xend=180, yend=0))

#all obs ice
df3 <- filter(df2, type == "ice.obs")
ggplot(data=df) + 
  geom_point(aes(x=ext, y=min, col=type)) + 
  # geom_text(data = df2, aes(x=ext, y=min, label=year), 
  #           nudge_x = 0.25, nudge_y = 0.25, 
  #           check_overlap = T) + 
  geom_segment(aes(x=0, y=180, xend=180, yend=0))
pca <- princomp(formula=~ext+min, data=df3)

df3 <- df2 %>% pivot_wider(names_from = type, values_from = c("min", "ext")) %>%
  drop_na()
p1 <- ggplot(data = df2) +
  geom_point(aes(x=year, y=min, col=type))
p2 <- ggplot(data = df2) +
  geom_point(aes(x=year, y=ext, col=type))
ggpubr::ggarrange(p1, p2, ncol=1)
#use pca of ext and min ice projections in linear predictor.