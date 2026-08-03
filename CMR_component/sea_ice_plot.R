# plot 2022 ice data and compare to ice data from Dan

ice.new <- read.csv("../input_data/sea_ice_vars_1992_2025.csv", header = T) |>
  # year is winter year (Nov - Apr), unlike in original data, the year for survival 
  #   and ice match up. I think!
  # ice data does not seem to match that used for 2021 SSA, need to ask
  select(year_winter, winter_count_sic_ge95) |>
  mutate(ice_data = (winter_count_sic_ge95 - mean(winter_count_sic_ge95))/sd(winter_count_sic_ge95) ) |>
  filter(year_winter < 2021)

ice.old <- read.csv("../input_data/extreme.sea.ice.csv", header = T) |>
  filter(year > 1991) |>
  drop_na() |>
  mutate(ice.sd = (ice.obs-mean(ice.obs)) / sd(ice.obs))
# year is winter year (Nov - Apr)
#ice.old <- ext.sea.ice[ext.sea.ice$year >= 1992 & ext.sea.ice$year <= 2016,]
#ice.old <- ice.data$ice.obs
#ice.old <- (ice.data-mean(ice.data)) / sd(ice.data)

plot(ice.old$ice.sd, ice.new$ice_data)
cor(ice.old$ice.sd, ice.new$ice_data) #not a strong correlation at all!
View(cbind(ice.old, ice.new))
plot(ice.new$ice_data[-29], ice.old$ice.sd[-1])
