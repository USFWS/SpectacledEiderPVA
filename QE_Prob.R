library(jagsUI)

YKD.L2008.4.5 <- readRDS("YKD_IPM_Scenarios/YKD.L2008.4.5.rds")
YKD.L2026.4.5 <- readRDS("YKD_IPM_Scenarios/YKD.L2026.4.5.rds")
YKD.constant.4.5 <- readRDS("YKD_IPM_Scenarios/YKD.constant.4.5.rds")
YKD.L2008.8.5 <- readRDS("YKD_IPM_Scenarios/YKD.L2008.8.5.rds")
YKD.L2026.8.5 <- readRDS("YKD_IPM_Scenarios/YKD.L2026.8.5.rds")
YKD.constant.8.5 <- readRDS("YKD_IPM_Scenarios/YKD.constant.8.5.rds")
YKD.L2019.current <- readRDS("YKD_IPM_Scenarios/YKD.L2019.current.rds")
ACP.current <- readRDS("ACP_IPM_Scenarios/ACP.current.rds")
ACP.4.5 <- readRDS("ACP_IPM_Scenarios/ACP.4.5.rds")
ACP.8.5 <- readRDS("ACP_IPM_Scenarios/ACP.8.5.rds")


QE <- function (x) {
  min(which(x <= 250))
}

n_iter <- 40000

YKD.L2008.4.5QE <- apply(YKD.L2008.4.5$sims.list$Nb, 1, QE)
YKD.L2008.4.5QE <- ifelse(YKD.L2008.4.5QE == "Inf", 0, YKD.L2008.4.5QE)
hist(YKD.L2008.4.5QE[YKD.L2008.4.5QE != 0])
length(YKD.L2008.4.5QE[YKD.L2008.4.5QE > 33 & YKD.L2008.4.5QE < 74])/n_iter*100
length(YKD.L2008.4.5QE[YKD.L2008.4.5QE >73])/n_iter*100

YKD.L2026.4.5QE <- apply(YKD.L2026.4.5$sims.list$Nb, 1, QE)
YKD.L2026.4.5QE <- ifelse(YKD.L2026.4.5QE == "Inf", 0, YKD.L2026.4.5QE)
hist(YKD.L2026.4.5QE[YKD.L2026.4.5QE != 0])
length(YKD.L2026.4.5QE[YKD.L2026.4.5QE > 33 & YKD.L2026.4.5QE < 74])/n_iter*100
length(YKD.L2026.4.5QE[YKD.L2026.4.5QE >73])/n_iter*100

YKD.constant.4.5QE <- apply(YKD.constant.4.5$sims.list$Nb, 1, QE)
YKD.constant.4.5QE <- ifelse(YKD.constant.4.5QE == "Inf", 0, YKD.constant.4.5QE)
hist(YKD.constant.4.5QE[YKD.constant.4.5QE != 0])
length(YKD.constant.4.5QE[YKD.constant.4.5QE > 33 & YKD.constant.4.5QE < 74])/n_iter*100
length(YKD.constant.4.5QE[YKD.constant.4.5QE >73])/n_iter*100

YKD.L2008.8.5QE <- apply(YKD.L2008.8.5$sims.list$Nb, 1, QE)
YKD.L2008.8.5QE <- ifelse(YKD.L2008.8.5QE == "Inf", 0, YKD.L2008.8.5QE)
hist(YKD.L2008.8.5QE[YKD.L2008.8.5QE != 0])
length(YKD.L2008.8.5QE[YKD.L2008.8.5QE > 33 & YKD.L2008.8.5QE < 74])/n_iter*100
length(YKD.L2008.8.5QE[YKD.L2008.8.5QE >73])/n_iter*100

YKD.L2026.8.5QE <- apply(YKD.L2026.8.5$sims.list$Nb, 1, QE)
YKD.L2026.8.5QE <- ifelse(YKD.L2026.8.5QE == "Inf", 0, YKD.L2026.8.5QE)
hist(YKD.L2026.8.5QE[YKD.L2026.8.5QE != 0])
length(YKD.L2026.8.5QE[YKD.L2026.8.5QE > 33 & YKD.L2026.8.5QE < 74])/n_iter*100
length(YKD.L2026.8.5QE[YKD.L2026.8.5QE >73])/n_iter*100

YKD.constant.8.5QE <- apply(YKD.constant.8.5$sims.list$Nb, 1, QE)
YKD.constant.8.5QE <- ifelse(YKD.constant.8.5QE == "Inf", 0, YKD.constant.8.5QE)
hist(YKD.constant.8.5QE[YKD.constant.8.5QE != 0])
length(YKD.constant.8.5QE[YKD.constant.8.5QE > 33 & YKD.constant.8.5QE < 74])/n_iter*100
length(YKD.constant.8.5QE[YKD.constant.8.5QE >73])/n_iter*100

YKD.L2019.currentQE <- apply(YKD.L2019.current$sims.list$Nb, 1, QE)
YKD.L2019.currentQE <- ifelse(YKD.L2019.currentQE == "Inf", 0, YKD.L2019.currentQE)
hist(YKD.L2019.currentQE[YKD.L2019.currentQE != 0])
length(YKD.L2019.currentQE[YKD.L2019.currentQE > 33 & YKD.L2019.currentQE < 74])/n_iter*100
length(YKD.L2019.currentQE[YKD.L2019.currentQE >73])/n_iter*100

ACP.currentQE <- apply(ACP.current$sims.list$Nb, 1, QE)
ACP.currentQE <- ifelse(ACP.currentQE == "Inf", 0, ACP.currentQE)
hist(ACP.currentQE[ACP.currentQE != 0])
length(ACP.currentQE[ACP.currentQE > 33 & ACP.currentQE < 74])/n_iter*100
length(ACP.currentQE[ACP.currentQE >73])/n_iter*100

ACP.4.5QE <- apply(ACP.4.5$sims.list$Nb, 1, QE)
ACP.4.5QE <- ifelse(ACP.4.5QE == "Inf", 0, ACP.4.5QE)
hist(ACP.4.5QE[ACP.4.5QE != 0])
length(ACP.4.5QE[ACP.4.5QE > 33 & ACP.4.5QE < 74])/n_iter*100
length(ACP.4.5QE[ACP.4.5QE >73])/n_iter*100

ACP.8.5QE <- apply(ACP.8.5$sims.list$Nb, 1, QE)
ACP.8.5QE <- ifelse(ACP.8.5QE == "Inf", 0, ACP.8.5QE)
hist(ACP.8.5QE[ACP.8.5QE != 0])
length(ACP.8.5QE[ACP.8.5QE > 33 & ACP.8.5QE < 74])/n_iter*100
length(ACP.8.5QE[ACP.8.5QE >73])/n_iter*100

