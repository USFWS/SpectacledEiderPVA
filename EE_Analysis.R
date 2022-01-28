library(ggplot2)
library(tidyverse)

# using EERound3 which includes additional responses from PF and DR on extreme
# sea ice day effects on BP (0 days and 110 days)
ee.dat <- read.csv("input_data/EERound3.csv", header = T)
head(ee.dat)
ee.dat$Question <- as.factor(ee.dat$Question)

# create standardized 80% bounds using equations in Hemmings et al. (2018)
s.CI <- 0.8
ee.dat$Low80 <- ee.dat$Best - ((ee.dat$Best - ee.dat$Lowest)*(s.CI/ee.dat$Range))
ee.dat$Hi80 <- ee.dat$Best + ((ee.dat$Highest - ee.dat$Best)*(s.CI/ee.dat$Range))
# truncating where necessary
ee.dat$Low80 <- ifelse(ee.dat$Low80 < ee.dat$LowerBound, ee.dat$LowerBound, ee.dat$Low80)
ee.dat$Hi80 <- ifelse(ee.dat$Hi80 > ee.dat$UpperBound, ee.dat$UpperBound, ee.dat$Hi80)

#calculate logical to indicate out of range
ee.dat$out <- ifelse(ee.dat$Round2==1,0,1)

# Q1
ee.Q1 <- ee.dat[ee.dat$Question == 1 & ee.dat$out == 0,]
Q1.Best <- round(mean(ee.Q1$Best, na.rm = T), 2)
Q1.Low <- round(mean(ee.Q1$Low80, na.rm = T), 2)
Q1.Hi <- round(mean(ee.Q1$Hi80, na.rm = T), 2)
ggplot(ee.Q1, aes(x = CodeName, y = Best)) + 
  geom_pointrange(aes(ymin = Low80, ymax = Hi80)) +
  geom_hline(yintercept = Q1.Best, col = "red", alpha = 0.5) + 
  geom_hline(yintercept = Q1.Low, col = "red", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = Q1.Hi, col = "red", linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(breaks = c(0.8, 0.85, Q1.Best, Q1.Low, Q1.Hi)) +
  labs(x = "Code names") + 
  labs(y = "Adult Mean Survival (ACP)") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  coord_flip()

# Q2
ee.Q2 <- ee.dat[ee.dat$Question == 2 & ee.dat$out == 0,]
Q2.Best <- round(mean(ee.Q2$Best, na.rm = T), 2)
Q2.Low <- round(mean(ee.Q2$Low80, na.rm = T), 2)
Q2.Hi <- round(mean(ee.Q2$Hi80, na.rm = T), 2)
ggplot(ee.Q2, aes(x = CodeName, y = Best)) + 
  geom_pointrange(aes(ymin = Low80, ymax = Hi80)) +
  geom_hline(yintercept = Q2.Best, col = "red", alpha = 0.5) + 
  geom_hline(yintercept = Q2.Low, col = "red", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = Q2.Hi, col = "red", linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(breaks = c(0, 0.25, 0.75, Q2.Best, Q2.Low, Q2.Hi)) +
  labs(x = "Code names") + 
  labs(y = "Juvenile Mean Survival (ACP)") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  coord_flip()

# Q3
ee.Q3 <- ee.dat[ee.dat$Question == 3 & ee.dat$out == 0,]
Q3.Best <- round(mean(ee.Q3$Best, na.rm = T), 2)
Q3.Low <- round(mean(ee.Q3$Low80, na.rm = T), 2)
Q3.Hi <- round(mean(ee.Q3$Hi80, na.rm = T), 2)
ggplot(ee.Q3, aes(x = CodeName, y = Best)) + 
  geom_pointrange(aes(ymin = Low80, ymax = Hi80)) +
  geom_hline(yintercept = Q3.Best, col = "red", alpha = 0.5) + 
  geom_hline(yintercept = Q3.Low, col = "red", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = Q3.Hi, col = "red", linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(breaks = c(0.50, 0.75, Q3.Best, Q3.Low, Q3.Hi)) +
  labs(x = "Code names") + 
  labs(y = "Juvenile Dispersal Rate") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  coord_flip()


# Q4
ee.Q4 <- ee.dat[ee.dat$Question == 4  & ee.dat$out == 0,]
Q4.Best <- round(mean(ee.Q4$Best, na.rm = T), 2)
Q4.Low <- round(mean(ee.Q4$Low80, na.rm = T), 2)
Q4.Hi <- round(mean(ee.Q4$Hi80, na.rm = T), 2)
ggplot(ee.Q4, aes(x = CodeName, y = Best)) + 
  geom_pointrange(aes(ymin = Low80, ymax = Hi80)) +
  geom_hline(yintercept = Q4.Best, col = "red", alpha = 0.5) + 
  geom_hline(yintercept = Q4.Low, col = "red", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = Q4.Hi, col = "red", linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(breaks = c(0.25, 0.50, 0.75, Q4.Best, Q4.Low, Q4.Hi)) +
  labs(x = "Code names") + 
  labs(y = "Adult Dispersal Rate") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  coord_flip()

# Q5.1
ee.Q5.1 <- ee.dat[ee.dat$Question == 5.1  & ee.dat$out == 0,]
Q5.1.Best <- round(mean(ee.Q5.1$Best, na.rm = T), 2)
Q5.1.Low <- round(mean(ee.Q5.1$Low80, na.rm = T), 2)
Q5.1.Hi <- round(mean(ee.Q5.1$Hi80, na.rm = T), 2)
ggplot(ee.Q5.1, aes(x = CodeName, y = Best)) + 
  geom_pointrange(aes(ymin = Low80, ymax = Hi80)) +
  geom_hline(yintercept = Q5.1.Best, col = "red", alpha = 0.5) + 
  geom_hline(yintercept = Q5.1.Low, col = "red", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = Q5.1.Hi, col = "red", linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(breaks = c(.2, 2.5, Q5.1.Best, Q5.1.Low, Q5.1.Hi)) +
  labs(x = "Code names") + 
  labs(y = "Productivity after 85 Ext Sea Ice Days") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  coord_flip()

# Q5.2
ee.Q5.2 <- ee.dat[ee.dat$Question == 5.2 & ee.dat$out == 0,]
Q5.2.Best <- round(mean(ee.Q5.2$Best, na.rm = T), 2)
Q5.2.Low <- round(mean(ee.Q5.2$Low80, na.rm = T), 2)
Q5.2.Hi <- round(mean(ee.Q5.2$Hi80, na.rm = T), 2)
ggplot(ee.Q5.2, aes(x = CodeName, y = Best)) + 
  geom_pointrange(aes(ymin = Low80, ymax = Hi80)) +
  geom_hline(yintercept = Q5.2.Best, col = "red", alpha = 0.5) + 
  geom_hline(yintercept = Q5.2.Low, col = "red", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = Q5.2.Hi, col = "red", linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(breaks = c(.2, 2.5, Q5.2.Best, Q5.2.Low, Q5.2.Hi)) +
  labs(x = "Code names") + 
  labs(y = "Productivity after 35 Ext Sea Ice Days") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  coord_flip()

# Q5.3
ee.Q5.3 <- ee.dat[ee.dat$Question == 5.3 & ee.dat$out == 0,]
Q5.3.Best <- round(mean(ee.Q5.3$Best, na.rm = T), 2)
Q5.3.Low <- round(mean(ee.Q5.3$Low80, na.rm = T), 2)
Q5.3.Hi <- round(mean(ee.Q5.3$Hi80, na.rm = T), 2)
ggplot(ee.Q5.3, aes(x = CodeName, y = Best)) + 
  geom_pointrange(aes(ymin = Low80, ymax = Hi80)) +
  geom_hline(yintercept = Q5.3.Best, col = "red", alpha = 0.5) + 
  geom_hline(yintercept = Q5.3.Low, col = "red", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = Q5.3.Hi, col = "red", linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(breaks = c(.2, 2.5, Q5.3.Best, Q5.3.Low, Q5.3.Hi)) +
  labs(x = "Code names") + 
  labs(y = "Productivity after 0 Ext Sea Ice Days") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  coord_flip()

# Q5.4
ee.Q5.4 <- ee.dat[ee.dat$Question == 5.4 & ee.dat$out == 0,]
Q5.4.Best <- round(mean(ee.Q5.4$Best, na.rm = T), 2)
Q5.4.Low <- round(mean(ee.Q5.4$Low80, na.rm = T), 2)
Q5.4.Hi <- round(mean(ee.Q5.4$Hi80, na.rm = T), 2)
ggplot(ee.Q5.4, aes(x = CodeName, y = Best)) + 
  geom_pointrange(aes(ymin = Low80, ymax = Hi80)) +
  geom_hline(yintercept = Q5.4.Best, col = "red", alpha = 0.5) + 
  geom_hline(yintercept = Q5.4.Low, col = "red", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = Q5.4.Hi, col = "red", linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(breaks = c(.2, 2.5, Q5.4.Best, Q5.4.Low, Q5.4.Hi)) +
  labs(x = "Code names") + 
  labs(y = "Productivity after 110 Ext Sea Ice Days") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  coord_flip()

# Q5 by respondent
ee.Q5 <- ee.dat %>%
  filter(Question == 5.1|Question == 5.2|Question == 5.3| Question == 5.4)
ee.Q5$Day <- NA
ee.Q5$Day[ee.Q5$Question == 5.1] <- 85
ee.Q5$Day[ee.Q5$Question == 5.2] <- 35
ee.Q5$Day[ee.Q5$Question == 5.3] <- 0
ee.Q5$Day[ee.Q5$Question == 5.4] <- 110

ggplot(ee.Q5, aes(x = Day, y = Best)) +
  geom_pointrange(aes(ymin = Low80, ymax = Hi80)) + 
  geom_point(aes(x = 60, y = 0.9), col = "red") +
  labs(x = "Extreme Sea Ice Days") +
  labs(y = "Breeding Propensity") +
  facet_wrap(~ Respondent)

# Q6
ee.Q6 <- ee.dat[ee.dat$Question == 6  & ee.dat$out == 0,]
Q6.Best <- round(mean(ee.Q6$Best, na.rm = T), 0)
Q6.Low <- round(mean(ee.Q6$Low80, na.rm = T), 0)
Q6.Hi <- round(mean(ee.Q6$Hi80, na.rm = T), 0)
ggplot(ee.Q6, aes(x = CodeName, y = Best)) + 
  geom_pointrange(aes(ymin = Low80, ymax = Hi80)) +
  geom_hline(yintercept = Q6.Best, col = "red", alpha = 0.5) + 
  geom_hline(yintercept = Q6.Low, col = "red", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = Q6.Hi, col = "red", linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(breaks = c(100000, Q6.Best, Q6.Low, Q6.Hi)) +
  labs(x = "Code names") + 
  labs(y = "Max Equil. Population Size (YKD)") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  coord_flip()

# Q7
ee.Q7 <- ee.dat[ee.dat$Question == 7  & ee.dat$out == 0,]
Q7.Best <- round(mean(ee.Q7$Best, na.rm = T), 0)
Q7.Low <- round(mean(ee.Q7$Low80, na.rm = T), 0)
Q7.Hi <- round(mean(ee.Q7$Hi80, na.rm = T), 0)
ggplot(ee.Q7, aes(x = CodeName, y = Best)) + 
  geom_pointrange(aes(ymin = Low80, ymax = Hi80)) +
  geom_hline(yintercept = Q7.Best, col = "red", alpha = 0.5) + 
  geom_hline(yintercept = Q7.Low, col = "red", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = Q7.Hi, col = "red", linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(breaks = c(50000, 75000, 100000, Q7.Best, Q7.Low, Q7.Hi)) +
  labs(x = "Code names") + 
  labs(y = "Max Equil. Population Size (ACP)") +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 11)) +
  coord_flip()

# Q8
ee.Q8 <- ee.dat[ee.dat$Question == 8  & ee.dat$out == 0,]
Q8.Best <- round(median(ee.Q8$Best, na.rm = T), 2)
Q8.Low <- round(median(ee.Q8$Low80, na.rm = T), 2)
Q8.Hi <- round(median(ee.Q8$Hi80, na.rm = T), 2)
ggplot(ee.Q8, aes(x = CodeName, y = Best)) + 
  geom_pointrange(aes(ymin = Low80, ymax = Hi80)) +
  geom_hline(yintercept = Q8.Best, col = "red", alpha = 0.5) + 
  geom_hline(yintercept = Q8.Low, col = "red", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = Q8.Hi, col = "red", linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(breaks = c(0.75, 1, Q8.Best, Q8.Low, Q8.Hi)) +
  labs(x = "Code names") + 
  labs(y = "Breeding Propensity 2YO") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  coord_flip()

# Q9
ee.Q9 <- ee.dat[ee.dat$Question == 9 & ee.dat$out == 0,]
Q9.Best <- round(mean(ee.Q9$Best, na.rm = T), 2)
Q9.Low <- round(mean(ee.Q9$Low80, na.rm = T), 2)
Q9.Hi <- round(mean(ee.Q9$Hi80, na.rm = T), 2)
ggplot(ee.Q9, aes(x = CodeName, y = Best)) + 
  geom_pointrange(aes(ymin = Low80, ymax = Hi80)) +
  geom_hline(yintercept = Q9.Best, col = "red", alpha = 0.5) + 
  geom_hline(yintercept = Q9.Low, col = "red", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = Q9.Hi, col = "red", linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(breaks = c(0, 0.5, Q9.Best, Q9.Low, Q9.Hi)) +
  labs(x = "Code names") + 
  labs(y = "Variation in Aerial Detection (YKD)") +
  coord_flip() + 
  theme(legend.position = "none",
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11))

# Q10
ee.Q10 <- ee.dat[ee.dat$Question == 10  & ee.dat$out == 0,]
Q10.Best <- round(mean(ee.Q10$Best, na.rm = T), 2)
Q10.Low <- round(mean(ee.Q10$Low80, na.rm = T), 2)
Q10.Hi <- round(mean(ee.Q10$Hi80, na.rm = T), 2)
ggplot(ee.Q10, aes(x = CodeName, y = Best)) + 
  geom_pointrange(aes(ymin = Low80, ymax = Hi80)) +
  geom_hline(yintercept = Q10.Best, col = "red", alpha = 0.5) + 
  geom_hline(yintercept = Q10.Low, col = "red", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = Q10.Hi, col = "red", linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(breaks = c(0.25, Q10.Best, Q10.Low, Q10.Hi)) +
  labs(x = "Code names") + 
  labs(y = "Variation in Aerial Detection (ACP)") +
  coord_flip() + 
  theme(legend.position = "none",
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11))

