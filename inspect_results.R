library(ggplot2)
library(jagsUI)
#Inspect and summarize output from a JAGS model
out <- readRDS(file = "MS_Scenarios/YKD.L2008.4.5.rds")
dat <- readRDS(file = "MS_Scenarios/jags.data4.5.rds")
out$mean
out$Rhat
lapply(out$Rhat, max, na.rm = TRUE)
lapply(out$n.eff, min, na.rm = TRUE)

#plot a sample of Population trajectories
nreps = 100
s <- sample(1:length(out$sims.list$Nb[,1]), nreps)
out2 <- out$sims.list$Nb[s,]
x <- 1:length(out$sims.list$Nb[1,])
plot(x, out2[1,], type = "l", ylim = c(0, 15000), 
     col = "lightgray")
for( i in 2:nreps){
  lines(x, out2[i,], col = "lightgray")
}
lines(x[-c(1:39)], out$mean$Nb[-c(1:39)], pch = 16)
points(1:39, out$mean$Nb[1:39])
#points(1:39, dat$count, pch = 16, cex = 0.5) #for YKD
points(20:39, dat$count[-c(1:19)], pch = 16, cex = 0.5) #for ACP
abline(h=250)
hist(out$sims.list$Nb[,40])
## plot just observed data years
plot(1:39, out2[1,1:39], type = "l", ylim = c(0, 20000), 
     col = "lightgray")
for( i in 2:nreps){
  lines(1:39, out2[i,1:39], col = "lightgray")
}
lines(x[-c(1:39)], out$mean$Nb[-c(1:39)], pch = 16)
points(1:39, out$mean$Nb[1:39])
#points(1:39, dat$count, pch = 16, cex = 0.5) #for YKD
points(20:39, dat$count[-c(1:19)], pch = 16, cex = 0.5) #for ACP
abline(h=250)
################################################################################
#calculate extinction probability
q_extinct <- function (x, threshold = 250) {
  s <- which(x <= threshold)
  t <- ifelse(length(s) == 0, length(x) + 1, min(s))
  return(t)
}
nyears <- length(out$mean$Nb)
qe <- apply(out$sims.list$Nb, 1, q_extinct)
#probability dist. extinction time
hist(qe[qe != nyears + 1], breaks = 1:nyears, probability = TRUE)
#cummulative dist extinction time
cum_qe <- c()
nsims <- length(qe)
for(i in 1:nyears){
  cum_qe[i] <- sum(qe <= i)/nsims
}
max(cum_qe)
cum_qe[39 + 20]
cum_qe[39 + 50]
cum_qe[39 + 74]
plot(1988:2100, cum_qe, col = c(rep("lightgray", 39), rep("black", 113-39)))
plot(2027:2100, cum_qe[-c(1:39)], xlim = c(2030, 2100))
################################################################################
#survival
out2 <- out$sims.list$phiA[,5:40]
df <- data.frame(PhiA = apply(out2, 2, mean), 
                 upper = apply(out2, 2, quantile, probs = 0.975),
                 lower = apply(out2, 2, quantile, probs = 0.025), 
                 Year = 1992:2027)
ggplot(data = df, aes(x = Year, y = PhiA)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  scale_x_continuous(breaks = seq(1992, 2027, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))
################################################################################
hist(out$sims.list$beta[,2])
traceplot(out, c("beta", "betaN"))
traceplot(out, c("betaBP"))
################################################################################
df <- data.frame(Year = 1988:2026, d = out$mean$d, upper = out$q97.5$d, lower = out$q2.5$d)
ggplot(data = df, aes( x = Year, y = d)) + geom_pointrange(aes(ymin = lower, ymax = upper)) + 
  scale_x_continuous(breaks = seq(1988, 2026, by = 2))
#interesting 2004 looks account for by d, maybe should add 2015 back in?
