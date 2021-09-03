ext.sea.ice <- read.csv("input_data/extreme.sea.ice.csv", header = T)
ice.data <- ext.sea.ice[ext.sea.ice$year >= 1988 & ext.sea.ice$year <= 2100,]
ice.data4.5 <- ifelse(is.na(ice.data$ice.obs), ice.data$ice.RCP4.5, 
                      ice.data$ice.obs)
ice.data8.5 <- ifelse(is.na(ice.data$ice.obs), ice.data$ice.RCP8.5, 
                      ice.data$ice.obs)

ice8.5 <- scale(ice.data8.5)
ice4.5 <- scale(ice.data4.5)

#MIDPOINT
hold4.5 <- expand.grid(a = seq(0.25, 0.75, length.out = 10), b = seq(2.5,3, length.out = 10), c =  seq(-1.5,-0.5, length.out = 10))

test4.5 <- c((35-mean(ice.data4.5))/sd(ice.data4.5),
             (60-mean(ice.data4.5))/sd(ice.data4.5),
             (85-mean(ice.data4.5))/sd(ice.data4.5))

pred <- matrix(NA, nrow = length(test4.5), ncol = nrow(hold4.5))
error <- rep(NA, length = nrow(hold4.5))
for (i in 1:nrow(hold4.5)){
  pred[,i] <- plogis(hold4.5$a[i] + hold4.5$b[i]*test4.5 + hold4.5$c[i]*test4.5^2)
  error[i] <- sum(abs(pred[1,i] - 0.776) + abs(pred[2,i] - 0.9) + abs(pred[3,i] - .72))
}
which(error == min(error))
pred[,365]
hold4.5[365,]

BP4.5 <- plogis(0.472 + 2.83*ice4.5 + - 1.17*ice4.5^2) # middle point
plot(BP4.5 ~ ice.data4.5, ylim = c(0,1))

#LOW POINT
hold4.5 <- expand.grid(a = seq(-0.9,-0.6, length.out = 10), b = seq(2.6,3, length.out = 10), c =  seq(-1.3,-0.9, length.out = 10))

test4.5 <- c((35-mean(ice.data4.5))/sd(ice.data4.5),
             (60-mean(ice.data4.5))/sd(ice.data4.5),
             (85-mean(ice.data4.5))/sd(ice.data4.5))

pred <- matrix(NA, nrow = length(test4.5), ncol = nrow(hold4.5))
error <- rep(NA, length = nrow(hold4.5))
for (i in 1:nrow(hold4.5)){
  pred[,i] <- plogis(hold4.5$a[i] + hold4.5$b[i]*test4.5 + hold4.5$c[i]*test4.5^2)
  error[i] <- sum(abs(pred[1,i] - 0.495) + abs(pred[2,i] - 0.9) + abs(pred[3,i] - 0.46))
}
which(error == min(error))
pred[,295]
hold4.5[295,]

BP4.5low <- plogis(-0.77 + 3*ice4.5 + - 1.21*ice4.5^2)
lines(BP4.5low[order(ice.data4.5)] ~ ice.data4.5[order(ice.data4.5)])

# HIGH POINT
hold4.5 <- expand.grid(a = seq(-3,4, length.out = 10), b = seq(-3,3, length.out = 10), c =  seq(-3,3, length.out = 10))

test4.5 <- c((35-mean(ice.data4.5))/sd(ice.data4.5),
             (60-mean(ice.data4.5))/sd(ice.data4.5),
             (85-mean(ice.data4.5))/sd(ice.data4.5))

pred <- matrix(NA, nrow = length(test4.5), ncol = nrow(hold4.5))
error <- rep(NA, length = nrow(hold4.5))
for (i in 1:nrow(hold4.5)){
  pred[,i] <- plogis(hold4.5$a[i] + hold4.5$b[i]*test4.5 + hold4.5$c[i]*test4.5^2)
  error[i] <- sum(abs(pred[1,i] - 1) + abs(pred[2,i] - 0.9) + abs(pred[3,i] - 0.93))
}
which(error == min(error))
pred[,539]
hold4.5[539,]

BP4.5high <- plogis(3.22 -1*ice4.5 + 0.33*ice4.5^2) # middle point
lines(BP4.5high[order(ice.data4.5)] ~ ice.data4.5[order(ice.data4.5)])

# B0: center on 0.472, going -0.77 to 3.22
# B1: center on 2.83, going -1 to 3
# B2:center on -1.17, going -1.21 to 0.33


quantile(rnorm(100000, 0.472, 2.15), c(0.1, 0.9))
1/2.15^2 # 0.22
quantile(rnorm(100000, 2.83, 3), c(0.1, 0.9)) 
1/3^2 # 0.11
quantile(rnorm(100000, -1.17, 1.2), c(0.1, 0.9)) 
1/1.2^2 #0.69

b0 <- rnorm(100000, .472, 2.15)
b0 <- b0[b0 > -0.77 & b0 < 3.22]
b0 <- sample(b0, 10000, replace = FALSE)
mean()
hist(b0)
b1 <- rnorm(100000, 2.83, 3)
b1 <- b1[b1 > -1 & b1 < 3]
b1 <- sample(b1, 10000, replace = FALSE)
hist(b1)
b2 <- rnorm(100000, -1.17, 1.2)
b2 <- b2[b2 > -1.21 & b2 < 0.33]
b2 <- sample(b2, 10000, replace = FALSE)
hist(b2)

test4.5.hold35 <- rep(NA, 10000)

BP4.5.hold <- plogis(b0[1] + b1[1]*ice4.5 + b2[1]*ice4.5^2)
test4.5.hold[1] <- plogis(b0[1] + b1[1]*test4.5[1] + b2[1]*test4.5[1]^2)
plot(BP4.5.hold ~ ice.data4.5, ylim = c(0,1))

for (i in 2:100){
  BP4.5sim <- plogis(b0[i] + b1[i]*ice4.5 + b2[i]*ice4.5^2)
  lines(BP4.5sim[order(ice.data4.5)] ~ ice.data4.5[order(ice.data4.5)])
}

for (i in 2:10000){
  BP4.5sim <- plogis(b0[i] + b1[i]*ice4.5 + b2[i]*ice4.5^2)
  test4.5.hold[i] <- plogis(b0[i] + b1[i]*test4.5[1] + b2[i]*test4.5[1]^2)
}

hist(test4.5.hold)
quantile(test4.5.hold, c(.1, .9)) 

# looking at model posteriors
#YKD.constant.8.5
BP.hold <- matrix(NA, nrow = 1000, ncol = 113)
BP.post.sample <- sample(1:40000, 1000, replace = FALSE)
YKD.constant.8.5samp <- data.frame(YKD.constant.8.5$sims.list$mean.logit.BP[BP.post.sample],
                                   YKD.constant.8.5$sims.list$betaF.ice[BP.post.sample],
                                   YKD.constant.8.5$sims.list$betaF.ice2[BP.post.sample])
names(YKD.constant.8.5samp) <- c("b0", "b1", "b2")
for (i in 1:1000){
  for (j in 1:113){
    BP.hold[i,j] <- plogis(YKD.constant.8.5samp$b0[i] + 
                             YKD.constant.8.5samp$b1[i]*ice8.5[j] +
                             YKD.constant.8.5samp$b2[i]*ice8.5[j]*ice8.5[j])
  }
}

BP.post.mean <- apply(BP.hold, 2, mean)
plot(BP.post.mean)
