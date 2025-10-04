#### 3D Distribution
#### Figure 1 Sampling map
#### AVC October 2025

library(tidyverse)
library(ggOceanMaps)

## Get data --------------------------------------------------------------------

load("./ProcessedData/detect_data.Rdata")


## Get station locations -------------------------------------------------------

station_loc <- detect_data %>% 
  group_by(station) %>% 
  slice_head() %>% 
  select(station, lat, lon)

## Sample map ------------------------------------------------------------------

map <- basemap(limits = c(min(station_loc$lon)-0.2,
                          max(station_loc$lon)+0.2,
                          min(station_loc$lat)-0.1,
                          max(station_loc$lat)+0.1),
               bathy.style = "rcb") +
  ggspatial::geom_spatial_point(data = station_loc, 
                                aes(x = lon, y = lat),
                                size = 1.5, alpha = 0.8,
                                color = "grey20") +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  guides(fill = "none") +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme(plot.margin = margin(0, 0, 0, 0))

map

save(map, file = "./Figures/sampling_map.Rdata")
