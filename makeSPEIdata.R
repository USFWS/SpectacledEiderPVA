#Produce VCF corrected SPEI population estimates for input to IPM
#Erik Osnas 20220126

#load function that makes estimates
source("popest.R")
#first YKD
dat <- read.csv("raw_data/YKDV_TransLevel_SPEI_1988to2021.csv")

vcfData <- data.frame(strata = c("High", "Low", "LowE", "LowN", "LowS", "Medium", "Pooled"),
                      vcf = c(3.09395, 1.354825, 1.354825, 1.354825, 1.354825, 2.459412, 2.734084), 
                      se = c(0.1873677, 0.1502611, 0.1502611, 0.1502611, 0.1502611, 0.1757193, 0.1502795))
dat <- dat[dat$strata%in%vcfData$strata,]
#Need total number of possible transects for each strata. These vary by year. 
#Based on Current (as of 20190729) YKD 'tyler density' strata (original used for paper is different)
#from file YKDV_M_1988to2021.csv
Mdat <- read.csv("raw_data/YKDV_M_1988to2021.csv")
Mdat <- Mdat[Mdat$strata%in%vcfData$strata,]
Mdat$M <- round(Mdat$M)

table(dat$Obs_Type, dat$Num)

results <- popest(dat=dat, vcfData=vcfData, Mdat=Mdat)
print(results$plot)
ggsave("YKD_SPEI.png", plot=results$plot, device = "png",  width = 7, height = 7)
results$popest
results$popest$mvcf <- results$popest$Nibb/results$popest$Index
write.csv(results$popest, file="input_data/YKD_SPEI.csv", row.names=FALSE)
