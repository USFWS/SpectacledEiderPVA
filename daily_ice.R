library(tidyverse)
ice <- read.table(file="raw_data/obs_ice_raw.txt", header=TRUE) %>%
  mutate(date = as.Date(date, format="%m/%d/%Y"), 
         winter = ifelse(month >= 11, year, year-1), 
         dow = ifelse(month >= 11, doy - 304, doy+60))


ggplot(data=ice, aes(x=dow, y=min4_bs3))+
  geom_path()
