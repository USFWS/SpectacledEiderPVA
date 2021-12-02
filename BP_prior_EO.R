#Explore EE responses for extreme sea ice effect on breeding
library(ggplot2)
low=c(0.5, 0.78, 1.0) #35 days extreme ice
high=c(0.46, 0.72, 0.96) #85 days of extreme ice

mean = 0.9 #60 extreme ice days

df <- data.frame(ice=c(35, 60, 85), bp=c(0.78, 0.9, 0.72), upper=c(1.0, 1.0, 0.96), lower=c(0.5, 0.8, 0.46))

gplot <- ggplot(data=df) +
  geom_ribbon(aes(x=ice, ymin=lower, ymax=upper)) +
  geom_line(aes(x=ice, y=bp))
print(gplot)

#constraint: sec. derivative <=0 ?
#intercept ice = 0, b0 ~dbeta(90, 10)

#try to simulate some data consistent with the EE response
#needs to be improved, use rriskDistributions::get.XXX.par functions
#average ice, logit(BP)
x1 <- rnorm(100, 2.24, 0.34)
#high ice
x2 <- rbeta(10000, 4.5, 1.725)
quantile(x2, probs=c(0.1, 0.5, 0.9))
hist(x2)
mean(qlogis(x2))
sd(qlogis(x2))
x2 <- rnorm(100, 1.16, 1.03)
#low ice
x3 <- rbeta(10000, 5.2, 1.725)
quantile(x3, probs=c(0.1, 0.5, 0.9))
hist(x3)
mean(qlogis(x3))
sd(qlogis(x3))
x3 <- rnorm(100, 1.34, 1)

#plot each rep, unfortunate that I named these x and especially in the order I did
plot(c(-1, 0, 1), plogis(c(x3[1], x1[1], x2[1])), type='l', ylim=c(0, 1), xlim=c(-3, 3))
for(i in 2:100){
  lines(c(-1, 0, 1), plogis(c(x3[i], x1[i], x2[i])))
}
#find and store best polynomial fit
betas <- matrix(0, 100, 3)
objfun <- function(beta, y){
  x <- c(-1, 0, 1)
  X <- matrix(c(1,1,1,x,x^2), 3, 3)
  sum((y-X%*%beta)^2)
}
for(i in 1:100){
  fit <- optim(par=c(0,1,0), fn=objfun, y=c(x3[i], x1[i], x2[i]))
  betas[i,] <- fit$par
}
#plot fits
x <- seq(-3, 3, by=0.01)
plot(x, matrix(c(rep(1, length(x)),x,x^2),length(x),3)%*%betas[1,], type='l', ylim=c(-3, 4), xlim=c(-3, 3))
for(i in 2:100){
  lines(x, matrix(c(rep(1, length(x)),x,x^2),length(x),3)%*%betas[i,])
}
#on real scale
plot(x, plogis(matrix(c(rep(1, length(x)),x,x^2),length(x),3)%*%betas[1,]), type='l', ylim=c(0, 1), xlim=c(-3, 3))
for(i in 2:100){
  lines(x, plogis(matrix(c(rep(1, length(x)),x,x^2),length(x),3)%*%betas[i,]))
}
#find mean and cov of betas
mbetas <- apply(betas, 2, mean)
vcv <- cov(betas)
#now simulate from MVnormal and plot to see if this looks right
preds <- mgcv::rmvn(n=100, mu=mbetas, V=vcv)
plot(x, plogis(matrix(c(rep(1, length(x)),x,x^2),length(x),3)%*%preds[1,]), type='l', ylim=c(0, 1), xlim=c(-3, 3))
for(i in 2:100){
  lines(x, plogis(matrix(c(rep(1, length(x)),x,x^2),length(x),3)%*%preds[i,]))
}


#below is exploratory stuff, no long applicable
# df <- data.frame(ice=rep(c(-1, 0, 1), each=100), bp=c(x3, x1, x2))
# plot(df$ice, df$bp)
# 
# library(mgcv)
# fit <- gam(bp~s(ice, k=3, bs="cr"), family=gaussian, data = df)
# plot(fit)
# df2 <- data.frame(ice=seq(-1, 1, by=0.1))
# preds <- predict(fit, newdata=df2, type="response", se.fit=TRUE)
# plot(df2$ice, plogis(preds$fit))
# 
# #why not just fit a GLM?
# df$ice2 <- df$ice*df$ice
# fit2 <- glm(bp~ice + ice2, data=df, family=gaussian)
# summary(fit2)
# summary(fit2)$cov.scaled
# 
# 
