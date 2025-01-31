### 3D distribution
### Density distribution by depth
### January 2025
### AVC

#### Set up environment --------------------------------------------------------

library(tidyverse)
library(PNWColors)
library(ggridges)

load("./ProcessedData/detect_species_meta.RData")
metadata <- read.csv("./Data/Hake_2019_metadata.csv")

#### Add metadata --------------------------------------------------------------

detect_data_meta <- detect_data %>% 
  left_join(metadata, by = c("NWFSCsampleID" = "sampleID"))

#### Collapse by station/species -----------------------------------------------

detect_by_station <- detect_data_meta %>% 
  group_by(station, depth, techRep, BestTaxon) %>% 
  mutate(totReads = sum(nReads)) %>% 
  mutate(detect = case_when(totReads > 1 ~ 1,
                            TRUE ~ 0)) %>% 
  slice_head() %>% 
  ungroup() %>% 
  filter(detect == 1)

#### Bubbleplot ----------------------------------------------------------------
detectDepth_bubble <- ggplot(detect_by_station, aes(y = BestTaxon, x = depth, 
                              fill = BestTaxon, color = BestTaxon)) +
  geom_count() +
  theme_minimal() + 
  coord_flip() +
  scale_x_reverse() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = c(pnw_palette("Cascades",12, type = "continuous"),
                               pnw_palette("Sunset",12, type = "continuous")[1:11])) +
  scale_color_manual(values = c(pnw_palette("Cascades",12, type = "continuous"),
                                pnw_palette("Sunset",12, type = "continuous")[1:11])) +
  theme(legend.position = "none")

#### Ridgeplot -----------------------------------------------------------------

lowdetect_subset <- detect_by_station %>% 
  group_by(BestTaxon) %>% 
  mutate(nDetect = n()) %>% 
  filter(nDetect < 3)

detectDepth_ridge <- ggplot(detect_by_station, aes(y = BestTaxon, x = depth, 
                              fill = BestTaxon, color = BestTaxon)) +
  geom_density_ridges(scale = 1, 
                      bandwidth = 100, 
                      jittered_points = TRUE,
                      point_alpha = 1,
                      point_shape = 21,
                      alpha = 0.6) +
  geom_point(data = lowdetect_subset, aes()) +
  theme_minimal() + 
  coord_flip() +
  scale_x_reverse() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = c(pnw_palette("Cascades",12, type = "continuous"),
                               pnw_palette("Sunset",12, type = "continuous")[1:12])) +
  scale_color_manual(values = c(pnw_palette("Cascades",12, type = "continuous"),
                                pnw_palette("Sunset",12, type = "continuous")[1:12])) +
  theme(legend.position = "none")

#### Save figures -------------------------------------------------------------

save(detectDepth_bubble, detectDepth_ridge, file = "./ProcessedData/detectDepth_plots.Rdata")

pdf(file = "./Figures/detectDepth_bubble.pdf")
detectDepth_bubble
dev.off()

pdf(file = "./Figures/detectDepth_ridge.pdf")
detectDepth_ridge
dev.off()

