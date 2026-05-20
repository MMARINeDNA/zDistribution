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
library(cowplot)
library(ggspatial)
library(marmap)

load("ProcessedData/m3.0models_preds_0.05degree.Rdata")
# load("./ProcessedData/detect_data.Rdata")
# load("./ProcessedData/detect_data_clean.RData")
# mmEcoEvo <- read.csv("./Data/MM_metadata.csv")
metadata <- read.csv("./Data/Hake_2019_metadata.csv")

det_colors <- c("500" = "#03051AFF", "300" = "#611F53FF", "150" = "#CB1B4FFF","50" = "#F16445FF","0" = "#F69C73FF")

#"#751F58FF" "#AC185AFF" "#DB2946FF" 

### m3.0c max POD depth map for three species ----------------------------------

# pull depth of max POD
maxPOD_depth <- m3.0c_sePreds %>% 
  group_by(BestTaxon, utm.lat,utm.lon) %>% #3816 groups
  arrange(desc(mu), .by_group = TRUE) %>% 
  mutate(max_mu = first(mu), max_low50 = first(low50), 
         max_high50 = first(high50), max_depth = first(depth)) %>%
  filter(mu > max_low50 & mu < max_high50) %>%
  mutate(depth_min = min(depth), depth_max = max(depth)) %>% 
  slice_head() %>% 
  ungroup() %>% 
  mutate(depthWidth = depth_max - depth_min) %>% 
  mutate(ci95 = high-low) %>% 
  ungroup()
  
# create convex hull study area
study_area <- st_as_sf(detections_sf, coords = c("utm.lon", "utm.lat"), crs = 32610) %>%
  summarise(geometry = st_union(geometry)) %>%
  st_convex_hull() 

# convert POD to sf
maxPOD_depth_sf <- maxPOD_depth %>% 
  mutate(lon_plain = utm.lon, lat_plain = utm.lat) %>% 
  st_as_sf(coords = c("utm.lon", "utm.lat"), crs = 32610) 

# clip POD <0.005
maxPOD_depth_clipped <- st_join(study_area, maxPOD_depth_sf, left = TRUE) %>% 
  mutate(depth = case_when(max_mu < 0.005~NA,
                           TRUE~depth)) %>% 
  mutate(depthWidth = case_when(max_mu < 0.005~NA,
                           TRUE~depthWidth))

# pull positive detections
pos_detect <- detections_sf %>% 
  filter(Detected == 1) %>% 
  mutate(depth = case_when(depth %in% c(48,50)~50,
                           depth %in% c(467,485,495,500)~500,
                           TRUE~depth))
  

### Create coastline shapefile -------------------------------------------------
world <- st_read("Data/ne_10m_land/ne_10m_land.shp")

data_bbox <- st_transform(detections_sf, 4326) %>%
  st_bbox() %>%
  st_as_sfc() %>% 
  st_buffer(dist = 70000)

westcoast_land <- st_crop(world, data_bbox)
westcoast_land_utm <- st_transform(westcoast_land, 32610)

### Get bathymetry data --------------------------------------------------------

# lon‐range and lat‐range:
lon1 <- min(detections_sf$lon); lon2 <- max(detections_sf$lon) 
lat1 <- min(detections_sf$lat); lat2 <- max(detections_sf$lat)

# Download bathymetry
bath <- getNOAA.bathy(lon1 = lon1, lon2 = lon2,
                      lat1 = lat1, lat2 = lat2,
                      resolution = 2)  # “1” ~ 1-minute (~1.8 km) resolution

# Convert bathy to a data.frame 
bath_df <- fortify.bathy(bath) 
#sort x and y values
bath_grid <- bath_df %>%
  arrange(y, x)
#create a matrix
zmat <- bath_grid %>%
  pivot_wider(names_from = x, values_from = z) %>%
  column_to_rownames("y") %>%
  as.matrix()
zmat <- t(zmat)
#pull x and y values
xvals <- as.numeric(rownames(zmat))
yvals <- as.numeric(colnames(zmat))
#create contour lines
cl <- contourLines(x = xvals, y = yvals, z = zmat, levels = c(-500, -1000, -2000))
#convert to sf
contour_sf <- do.call(rbind,
  lapply(cl, function(x) {st_sf(level = x$level,
                                geometry = st_sfc(st_linestring(cbind(x$x, x$y)),
                                                  crs = 4326))}))
#convert to utm
contour_sf <- st_transform(contour_sf, 32610)

#save(maxPOD_depth, pos_detect, maxPOD_depth_clipped, file = "./ProcessedData/H3.0c_pred_flt.Rdata")

#### Max POD plot with matching color scale ------------------------------------

depth_max_detect <- ggplot() +
  geom_tile(data = maxPOD_depth_clipped, 
            aes(x=lon_plain, y = lat_plain, fill = depth)) +
  geom_sf(data = westcoast_land_utm, fill = "grey50", colour = NA) +
  scale_fill_viridis_c(name = "Depth (m)",
                       option = "mako",
                       trans = "reverse",
                       begin = 0.4, end = 0.9, na.value = "transparent") +
  ggspatial::geom_spatial_point(data = pos_detect,
                                aes(x = lon, y = lat,
                                    color = as.factor(depth)),
                                size = 1,
                                alpha = 0.8,
                                stroke = 1,
                                position = position_jitter(width = 6500,
                                                           height = 6500)) +
  scale_color_manual(values = det_colors) +
  facet_wrap(~BestTaxon, labeller = label_wrap_gen(width=10)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "none", #c(0.54, 0.45),    # <<-- Adjust to place legend over land
        legend.justification = c("left"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        strip.text = element_text(face = "italic")) +
  geom_sf(data = contour_sf, color = "grey70", linewidth = 0.3) +
  guides(color = guide_legend("Detection\ndepth (m)")) 

depth_max_detect

ci95_plot <- ggplot() +
  geom_tile(data = maxPOD_depth_clipped, 
            aes(x = lon_plain, y = lat_plain, 
                fill = depthWidth)) +
  geom_sf(data = westcoast_land_utm, fill = "grey50", colour = NA) +
  theme_minimal() +
  geom_sf(fill = "grey50", colour = NA) +
  scale_fill_viridis_c(name = "50% CI\ndepth range (m)",
                       option = "magma",
                       trans = "reverse",
                       begin = 0.15, end = 1, na.value = "transparent") +
  facet_wrap(~BestTaxon, labeller = label_wrap_gen(width=10), ncol = 1) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "none", #c(0.54, 0.45),    # <<-- Adjust to place legend over land
        legend.justification = c("left"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        strip.text = element_text(face = "italic")) +
  geom_sf(data = contour_sf, color = "grey70", linewidth = 0.3)

depth_max_detect | ci95_plot

save(depth_max_detect, ci95_plot, file = "./Figures/H3.0c_map.Rdata")

#### Plot all as one grid ####
sp_labels <- c("Berardius bairdii" = "northern giant\nbottlenose whale", "Lagenorhynchus obliquidens" = "Pacific\nwhite-sided dolphin", 
               "Megaptera novaeangliae" = "humpback whale")

depth_max_detect_long <- maxPOD_depth_clipped %>%
  #select(-geometry) %>%
  select(BestTaxon, lon_plain, lat_plain, geometry, depth, depthWidth) %>%
  pivot_longer(depth:depthWidth, names_to = "pred_type", values_to = "pred_value")

depth_max_detect_grid <- ggplot() +
  geom_tile(data = subset(depth_max_detect_long, pred_type == "depth"),
            aes(x = lon_plain, y = lat_plain, fill = pred_value)) +
  scale_fill_viridis_c(name = "Maximized detection\ndepth (m)",
                       option = "mako",
                       trans = "reverse",
                       begin = 0.4, end = 0.9, na.value = "grey90") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = subset(depth_max_detect_long, pred_type == "depthWidth"),
            aes(x = lon_plain, y = lat_plain, fill = pred_value)) +
  scale_fill_viridis_c(name = "Model uncertainty:\n50% CI strip width,\ndepth (m)",
                       option = "magma",
                       trans = "reverse",
                       begin = 0.15, end = 1,
                       na.value = "grey90") +
  geom_sf(data=westcoast_land_utm, fill = "grey50", colour = NA) +
  ggspatial::geom_spatial_point(data = pos_detect,
                                aes(x = lon, y = lat,
                                    color = as.factor(depth)),
                                size = 1,
                                alpha = 0.8,
                                stroke = 1,
                                position = position_jitter(width = 6500,
                                                           height = 6500)) +
  scale_color_manual(values = det_colors) +
  facet_grid(BestTaxon~pred_type,
             labeller = labeller(BestTaxon = as_labeller(sp_labels))) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "right", #c(0.54, 0.45),    # <<-- Adjust to place legend over land
        legend.justification = c("center"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        strip.text.x = element_blank()) +
  geom_sf(data = contour_sf, color = "grey70", linewidth = 0.3) +
  guides(color = guide_legend("Detection\ndepth (m)"))

depth_max_detect_grid

#### m3.1 Marginal effect of density -------------------------------------------

Dmarg_plot <- ggplot(D_grid) +
  geom_ribbon(aes(D, ymin = low, ymax = high, fill = species), alpha = 0.3) +
  geom_ribbon(aes(D, ymin = low50, ymax = high50, fill = species), stat = "identity", alpha = 0.3, color = NA) +
  geom_line(aes(D, fit, color = species), linewidth = 1) +
  geom_rug(data = detect_data_dens %>%
             mutate(D_jittered = D + runif(n(), 0.0001, 0.001)),
           aes(D_jittered), color = "grey70") +
  geom_rug(data = detect_data_dens %>% filter(Detected == 1) %>%
             mutate(D_jittered = D + runif(n(), 0.0001, 0.001)), 
           aes(D_jittered, color = species), sides = "t") +
  scale_fill_manual(values = c(pnw_palette("Cascades",5, type = "continuous")[4:5],
                               pnw_palette("Sunset",1, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Cascades",5, type = "continuous")[4:5],
                                pnw_palette("Sunset",1, type = "continuous"))) +
  facet_wrap(~species, scales = "free", 
             labeller = labeller(species = as_labeller(sp_labels))) +
  theme_minimal() +
  labs(y = "Predicted detection probability", x = "Predicted Density", color = "Species", fill = "Species") 

#### Save plots ----------------------------------------------------------------
save(depth_max_detect, ci95_plot, 
     depth_max_detect_grid, Dmarg_plot, 
     file = "./Figures/H3.0c_map_clean.Rdata")

