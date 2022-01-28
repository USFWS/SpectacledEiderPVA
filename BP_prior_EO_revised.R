library(ggplot2)
library(rriskDistributions)

#Explore EE responses for extreme sea ice effect on breeding
vlow=c(0.7, 0.78, 0.96) # 0 days of extreme ice
low=c(0.44, 0.69, 0.9) #35 days extreme ice
high=c(0.41, 0.64, 0.83) #85 days of extreme ice
vhigh=c(0.51, 0.68, 0.81) # 110 days of extreme ice

mean = 0.9 #60 extreme ice days

df <- data.frame(ice=c(0, 35, 60, 85, 110), bp=c(0.78, 0.69, 0.9, 0.64, 0.68), 
                 upper=c(0.96, 0.9, 1.0, 0.83, 0.81), lower=c(0.7, 0.44, 0.8, 0.41, 0.51))

gplot <- ggplot(data=df) +
  geom_ribbon(aes(x=ice, ymin=lower, ymax=upper)) +
  geom_line(aes(x=ice, y=bp))
print(gplot)

#constraint: sec. derivative <=0 ?
#intercept ice = 0, b0 ~dbeta(90, 10)

#try to simulate some data consistent with the EE response
#revised from BP_prior_EO to use rriskDistributions
ee35 <- low
par35 <- get.beta.par(p=c(0.1, 0.5, 0.9), q=ee35)
quantile(rbeta(10000, par35[[1]], par35[[2]]), c(0.1, 0.5, 0.9))
# 10%       50%       90% 
# 0.4272811 0.6923437 0.8921043 
sim35 <- rbeta(100, par35[[1]], par35[[2]])

base60 <- c(0.8,0.9,1.0)
par60 <- get.beta.par(p=c(0.1, 0.5, 0.9), q=base60, tol = 1e8)
quantile(rbeta(10000, par60[[1]], par60[[2]]), c(0.1, 0.5, 0.9))
# 10%       50%       90% 
# 0.7987395 0.9002246 0.9608277
sim60 <- rbeta(100, par60[[1]], par60[[2]])

ee85 <- high
par85 <- get.beta.par(p=c(0.1, 0.5, 0.9), q=ee85)
quantile(rbeta(10000, par85[[1]], par85[[2]]), c(0.1, 0.5, 0.9))
# 10%       50%       90% 
# 0.4115011 0.6410260 0.8313530 
sim85 <- rbeta(100, par85[[1]], par85[[2]])

plot(c(-1, 0, 1), c(sim35[1], sim60[1], sim85[1]), type='l', ylim=c(0, 1), xlim=c(-3, 3))
for(i in 2:100){
  lines(c(-1, 0, 1), c(sim35[i], sim60[i], sim85[i]))
}

#find and store best polynomial fit
betas <- matrix(0, 100, 3)
objfun <- function(beta, y){
  x <- c(-1, 0, 1)
  X <- matrix(c(1,1,1,x,x^2), 3, 3)
  sum((y-X%*%beta)^2)
}
for(i in 1:100){
  fit <- optim(par=c(0,1,0), fn=objfun, y=qlogis(c(sim35[i], sim60[i], sim85[i])))
  betas[i,] <- fit$par
}
#plot fits
x <- seq(-3, 3, by=0.01)
plot(x, matrix(c(rep(1, length(x)),x,x^2),length(x),3)%*%betas[1,], type='l', ylim=c(-3, 4), xlim=c(-3, 3))
for(i in 2:100){
  lines(x, matrix(c(rep(1, length(x)),x,x^2),length(x),3)%*%betas[i,])
}
for(i in 1:100){
  lines(c(-1, 0, 1), qlogis(c(sim35[i], sim60[i], sim85[i])), col = "red")
}

#on real scale
plot(x, plogis(matrix(c(rep(1, length(x)),x,x^2),length(x),3)%*%betas[1,]), type='l', ylim=c(0, 1), xlim=c(-3, 3))
for(i in 2:100){
  lines(x, plogis(matrix(c(rep(1, length(x)),x,x^2),length(x),3)%*%betas[i,]))
}
for(i in 1:100){
  lines(c(-1, 0, 1), qlogis(c(sim35[i], sim60[i], sim85[i]), col = "red")
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
