### Expand marine mammal dive depth data
### AVC Jan 2025


library(tidyverse)


dive_depth <- read.csv("Data/MM_dive_time.csv") %>% 
  mutate(max_depth = as.numeric(max_depth)) %>% 
  group_by(Species) %>% 
  complete(max_depth = 0:max(max_depth)) %>% 
  fill(time_per_m, depth_bin_m, surrogate, percent_time, .direction ="up") %>% 
  rename("depth" = "max_depth", "depth_bin_max" = "depth_bin_m", "time_in_bin" = "percent_time")
  
write.csv(dive_depth, file = "Data/MM_dive_time_expand.csv")
