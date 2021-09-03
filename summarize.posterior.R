#summarize IPM results
library(jagsUI)

extinct <- function(x){
  mean(x <= 250)*100
}

# Count Data, 1988 - 2019
counts <- read.csv("input_data/YKD_SPEI.csv", header = T)
counts <- rbind(data.frame(Year=1988:2006, Nibb=NA, seNibb=NA, Index=NA, seIndex=NA), counts)
observer <- c(rep(1, 19), c(1, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5, 6, 7)+1)

post <- readRDS("YKD_IPM_Scenarios/YKD.L2026.4.5.rds")

summary(post)
print(post)
hist(post$summary[,"Rhat"])
post$summary[,"Rhat"][post$summary[,"Rhat"]>1.1]

post$mean

ext.prob <- apply(post$sims.list$Nb, 2, FUN = extinct)
max(ext.prob)

plot(ext.prob, xlab = "Year", ylab = "QE probability (%)")

#plot results 
#counts
jit=0.2
plot(counts$Year, counts$Nibb, xlim=c(1987, 2100), ylim=c(0, 20000), 
     xlab="Year", ylab="Count")
arrows(x0=counts$Year, 
       y0=counts$Nibb - 1.96*counts$seNibb, 
       y1=counts$Nibb + 1.96*counts$seNibb, length=0)
points(1:length(post$mean$Nb)+jit+1987, post$mean$Nb, pch=16)
arrows(x0=1:length(post$mean$Nb)+jit+1987,  
       y0=post$q2.5$Nb, 
       y1=post$q97.5$Nb, length=0, lwd=2)
sim.plot <- sample(1:nrow(post$sims.list$Nb), 20, replace = FALSE)
for(i in 1:20){lines(2019:2100, post$sims.list$Nb[sim.plot[i],32:113])}

hist(post$sims.list$betaN[,1], probability = TRUE, col="lightgray")
x=seq(0, 0.1, length=1000)
lines(x, dgamma(x, 3, 3300))
legend("topleft", legend=c("Posterior", "Prior"), pch=c(22, NA), 
       lty=c(NA, 1), lwd=c(NA, 1), col=c(1), pt.bg=c("lightgray", NA))

hist(post$sims.list$betaN[,2], probability = TRUE, col="lightgray")
x=seq(0, 0.1, length=1000)
lines(x, dgamma(x, 3, 3300))
legend("topleft", legend=c("Posterior", "Prior"), pch=c(22, NA), 
       lty=c(NA, 1), lwd=c(NA, 1), col=c(1), pt.bg=c("lightgray", NA))

hist(post$sims.list$mean.log.F, probability = TRUE, col="lightgray")
x=seq(-1, 0, length=1000)
lines(x, dnorm(x, -0.47, 0.1))
legend("topleft", legend=c("Posterior", "Prior"), pch=c(22, NA), 
       lty=c(NA, 1), lwd=c(NA, 1), col=c(1), pt.bg=c("lightgray", NA))

post$mean$o

hist(post$sims.list$sigma.o, probability = TRUE, col="lightgray")
x=seq(0, 1, length=1000)
lines(x, dgamma(x, 1, 10))
legend("topright", legend=c("Posterior", "Prior"), pch=c(22, NA), 
       lty=c(NA, 1), lwd=c(NA, 1), col=c(1), pt.bg=c("lightgray", NA))

hist(post$sims.list$sigma.d, probability = TRUE, col="lightgray")
x=seq(0, 1, length=1000)
lines(x, dgamma(x, 1, 13))
legend("topright", legend=c("Posterior", "Prior"), pch=c(22, NA), 
       lty=c(NA, 1), lwd=c(NA, 1), col=c(1), pt.bg=c("lightgray", NA))

