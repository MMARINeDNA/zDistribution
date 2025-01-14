#### 3D Distribution
#### wrangling data
#### January 2025
#### EKJ

library(tidyverse)
library(mgcv)

load("./ProcessedData/detect_data_meta.RData")

# need to expand to 1/0 for each possible species+sample combo

all_species_samples <- detect_data_meta %>%
  expand(nesting(Sample_name, primer, dilution, techRep, station, depth, volume, 
                 bottom.depth.consensus, utm.lat, utm.lon), BestTaxon) 

nondetections <- all_species_samples %>% 
  anti_join(detect_data_meta) %>%
  mutate(Detected = 0)

detections <- detect_data_meta %>%
  select(names(all_species_samples)) %>%
  mutate(Detected = 1)

binary_data <- rbind.data.frame(detections, nondetections)

# now want to collapse across tech reps, so we have one record per sample



         
         