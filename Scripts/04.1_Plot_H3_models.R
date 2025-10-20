#### 3D Distribution
#### Figures for Q3 models
#### October 2025
#### AVC

library(tidyverse)
library(PNWColors)
library(mgcv)
library(patchwork)
library(terra)
library(ggOceanMaps)
library(sf)

load("./ProcessedData/m3.0models_preds.Rdata")
load("./ProcessedData/detect_data.RData")
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")
metadata <- read.csv("./Data/Hake_2019_metadata.csv")

### m3.0c max POD depth map for three species ----------------------------------

# pull depth of max POD
maxPOD_depth <- m3.0c_sePreds %>% 
  group_by(BestTaxon, lat,lon) %>% #3816 groups
  arrange(desc(mu), .by_group = TRUE) %>% 
  slice_head() %>% 
  ungroup() %>% 
  mutate(ci95 = high-low)

# pull depth of detections
pos_detect <- detect_data %>% 
  filter(BestTaxon %in% c("Lagenorhynchus obliquidens",
                          "Megaptera novaeangliae",
                          "Berardius bairdii")) %>% 
  filter(Detected == 1) %>% 
  mutate(depth = case_when(depth %in% c(48,50)~50,
                           depth %in% c(467,485,495,500)~500,
                           TRUE~depth))

#create coastline shapefile
world <- st_read("Data/ne_10m_land/ne_10m_land.shp")

data_bbox <- st_as_sf(maxPOD_depth, coords = c("lon", "lat"), crs = 4326) %>%
  st_bbox() %>%
  st_as_sfc()  %>%
  st_buffer(dist = 2)

westcoast_land <- st_crop(world, data_bbox)

# depth of max POD map
depth_max_detect <- ggplot(westcoast_land) +
  geom_tile(data = maxPOD_depth, aes(x = lon, y = lat, 
                                          fill = depth)) +
  theme_minimal() +
  geom_sf(fill = "grey50", colour = NA) +
  scale_fill_viridis_c(name = "Depth of max POD (m)",
                       option = "mako",
                       trans = "reverse",
                       begin = 0.4, end = 0.9) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        legend.position = "right",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = margin(0, 0, 0, 0)) +
  ggspatial::geom_spatial_point(data = pos_detect, 
                                aes(x = lon, y = lat, 
                                    color = as.factor(depth)),
                                size = 1,
                                alpha = 0.8,
                                stroke = 1,
                                position = position_jitter(width = 0.05, 
                                                           height = 0.05)) +
  scale_color_viridis_d("Detection depth", option = "rocket",
                        direction = -1, begin = 0.3, end = 0.8) +
  facet_grid(~BestTaxon, labeller = labeller(BestTaxon = 
                                                     c("Berardius bairdii" = "Bbai",
                                                       "Lagenorhynchus obliquidens" = "Lobl",
                                                       "Megaptera novaeangliae" = "Mnov")))

depth_max_detect

# 95% CI of max POD

POD_depth_CI <- ggplot(westcoast_land) +
  geom_tile(data = maxPOD_depth, aes(x = lon, y = lat, 
                                          fill = ci95)) +
  theme_minimal() +
  geom_sf(fill = "grey50", colour = NA) +
  scale_fill_viridis_c(name = "95% CI width",
                       option = "rocket",
                       trans = "reverse") +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        legend.position = "right",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = margin(0, 0, 0, 0)) +
  facet_wrap(~BestTaxon, labeller = labeller(BestTaxon = 
                                               c("Berardius bairdii" = "Bbai",
                                                 "Lagenorhynchus obliquidens" = "Lobl",
                                                 "Megaptera novaeangliae" = "Mnov")))

POD_depth_CI

save(depth_max_detect, POD_depth_CI, file = "./Figures/H3.0c_map.Rdata")

### m3.0c dept-lat variability figure ------------------------------------------

m3.0c_predictions <- expand_grid(depth = 0:500, 
                                 lat = c(34.39501,41.47946,48.5639),
                                 lon = (min(metadata$lon) + max(metadata$lon))/2,
                                 BestTaxon = as.factor(c("Lagenorhynchus obliquidens",
                                                         "Megaptera novaeangliae",
                                                         "Berardius bairdii")))

m3.0cpreds <- predict(m3.0c, m3.0c_predictions, type = "response", se.fit = TRUE)
m3.0c_sePreds <- data.frame(m3.0c_predictions,
                           mu   = exp(m3.0cpreds$fit),
                           low  = exp(m3.0cpreds$fit - 0.674 * m3.0cpreds$se.fit), #50% CI
                           high = exp(m3.0cpreds$fit + 0.674 * m3.0cpreds$se.fit)) %>% 
  left_join(mmEcoEvo, by = c("BestTaxon" = "Species"))

###Try this 

#--- 1. Compute tighter facet-specific limits (based on central 90% of mu)
facet_lims <- m3.0c_sePreds %>%
  group_by(lat, abbrev) %>%
  summarise(
    y_min = quantile(mu, 0.05, na.rm = TRUE),
    y_max = quantile(mu, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

#--- 2. Define a vector of colors (one per plot)
facet_colors <- c("#88A2B9", "#88A2B9", "#88A2B9", 
                  "#41476B", "#41476B", "#41476B",
                  "#2D4030", "#2D4030", "#2D4030")

#--- 3. Split data by facet and plot each subset with its own zoom window and color
split_data <- m3.0c_sePreds %>%
  split(list(.$lat, .$abbrev), drop = TRUE)

plots <- vector("list", length(split_data))

for (i in seq_along(split_data)) {
  df <- split_data[[i]]
  lims <- facet_lims %>%
    filter(lat == unique(df$lat), abbrev == unique(df$abbrev))
  
  color_i <- facet_colors[(i - 1) %% length(facet_colors) + 1]  # cycle through 3 colors
  
  plots[[i]] <- ggplot(df, aes(x = depth)) +
    geom_smooth(
      aes(ymin = low, ymax = high, y = mu),
      stat = "identity", linewidth = 1,
      color = color_i, fill = scales::alpha(color_i, 0.3)
    ) +
    coord_cartesian(ylim = c(lims$y_min, lims$y_max), expand = FALSE) +
    labs(
      title = paste(unique(df$abbrev), unique(df$lat), sep = "  |  "),
      x = "Depth (m)",
      y = "POD"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      legend.position = "none"
    )
}

#--- 4. Combine all facets into a single layout
m3.0clat_POD <- wrap_plots(plots, axis_titles = "collect") &
  theme(legend.position = "none")

save(m3.0clat_POD, file = "./Figures/H3.0c_plot.Rdata")

png(file = "./Figures/H3.0c_plot.png")
m3.0clat_POD
dev.off()
