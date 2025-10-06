### Expand marine mammal dive depth data
### AVC Jan 2025


library(tidyverse)


dive_depth <- read.csv("Data/MM_dive_time.csv") %>% 
  mutate(max_depth = as.numeric(max_depth)) %>% 
  group_by(Species) %>% 
  complete(max_depth = 1:max(max_depth)) %>% 
  mutate(bin10m = cut_width(max_depth, width = 10, boundary = 1, labels = FALSE)) %>% 
  mutate(binEqual = cut_number(max_depth, n = 10, labels = FALSE)) %>% 
  fill(time_per_m, depth_bin_m, surrogate, percent_time, .direction ="up") %>% 
  group_by(Species, bin10m) %>% 
  mutate(time_10m = sum(time_per_m, na.rm = FALSE)) %>% 
  group_by(Species,binEqual) %>% 
  mutate(time_equalBin = sum(time_per_m, na.rm = FALSE)) %>% 
  rename("depth" = "max_depth", "time_in_bin" = "percent_time")
  
write.csv(dive_depth, file = "Data/MM_dive_time_expand.csv", row.names = FALSE)
