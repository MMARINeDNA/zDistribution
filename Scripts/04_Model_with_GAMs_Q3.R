#### 3D Distribution
#### Models testing Q3: Depth distribution across xy
#### Fall 2025
#### EKJ&AVC

library(mgcv)
library(tidyverse)
library(PNWColors)
library(marmap)
library(terra)
library(sf)

load("./ProcessedData/detect_data.Rdata")
load("./ProcessedData/detect_data_clean.RData")
detect_data <- detect_data %>% mutate(BestTaxon = as.factor(BestTaxon))
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")

### Get bathymetry for pred grid -----------------------------------------------

bathy <- getNOAA.bathy(lon1 = min(detect_data_clean$lon), 
                       lon2 = max(detect_data_clean$lon), 
                       lat1 = min(detect_data_clean$lat),  
                       lat2 = max(detect_data_clean$lat),
                       resolution = 1)

bathy_raster <- marmap::as.raster(bathy)
bathy_r <- rast(bathy_raster)

### Get density estimates from Becker et al. 2020 tech memo --------------------

density_files <- list.files(path = "./Data/CCE_model_run", pattern = "predgrid_yearly", full.names = TRUE)

density_data <- density_files %>%
  set_names(basename(.)) %>%
  map_dfr(read_csv, .id = "file") %>% 
  filter(year == 2014) %>% 
  separate(file, into = c(NA,NA,"species",NA), sep = "_") 

#subset detect data to target the three main species
detect_data_sub <- detect_data_clean %>% 
  filter(BestTaxon %in% c("Berardius bairdii", 
                          "Lagenorhynchus obliquidens",
                          "Megaptera novaeangliae"))

#join data using a spatial join
detections_sf <- detect_data_sub %>%
  mutate(lat_plain = lat, lon_plain = lon,
         utm.lat_plain = utm.lat, utm.lon_plain = utm.lon) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  rename("lat" = lat_plain, "lon" = lon_plain) %>% 
  select(-utm.lat_plain, -utm.lon_plain) %>% 
  st_transform(32610) 

coords <- st_coordinates(detections_sf)

detections_sf <- detections_sf %>%
  mutate(utm.lon = coords[,1], utm.lat = coords[,2])

density_sf <- density_data %>%
  st_as_sf(coords = c("mlon", "mlat"), crs = 4326) %>% 
  st_transform(32610) 

detect_data_dens <- detections_sf %>%
  group_by(BestTaxon) %>%
  group_modify(~ {density_sub <- density_sf %>% filter(species == unique(.y$BestTaxon))
    st_join(.x, density_sub, join = st_nearest_feature)}) %>%
  ungroup() %>% 
  select(-geometry)


# Q3: Does depth distribution of detections vary across xy spatial distribution or species?
# H3.0: Depth distribution of detections does not vary across xy spatial distribution or species.
  #H3.0a: xy variability + z variability
  #H3.0b: xyz covariability
  #H3.0c: xyz covariability by species
  #H3.1: z and density covariability by species

### H3.0a: Depth distribution plus xy space distribution -----------------------

m3.0a <- bam(Detected ~ s(depth) + s(utm.lat, utm.lon),
             family = "binomial",
             data = detect_data_dens,
             method = "fREML",
             discrete = FALSE,
             nthreads = 4)
summary(m3.0a)

# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(depth)            2.763  3.269  23.80 4.74e-05 ***
#   s(utm.lat,utm.lon) 14.600 18.878  80.77  < 2e-16 ***
#   ---
#   
# R-sq.(adj) =  0.0229   Deviance explained = 7.88%

gam.check(m3.0a)
# k'   edf k-index p-value
# s(depth)            9.00  2.76     0.9    0.25
# s(utm.lat,utm.lon) 29.00 14.60     0.9    0.38

AIC(m3.0a)
#1629

### H3.0b: Depth distribution smoothed by xy space distribution ----------------

m3.0b <- bam(Detected ~ te(depth, utm.lat, utm.lon),
             family = "binomial",
             data = detect_data_dens,
             method = "fREML",
             discrete = TRUE,
             nthreads = 4)
summary(m3.0b)
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# te(depth,utm.lat,utm.lon) 28.72  36.26  142.1  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0261   Deviance explained = 11.2%
AIC(m3.0b)
#1596

### m3.0b predictions ----------------------------------------------------------

m3.0b_pred_grid <- expand_grid(depth = seq(0,500, by = 100),
                               utm.lat = seq(min(detect_data$utm.lat, na.rm = TRUE),
                                         max(detect_data$utm.lat, na.rm = TRUE),
                                         by = 100),
                               utm.lon = seq(min(detect_data$utm.lon, na.rm = TRUE),
                                         max(detect_data$utm.lon, na.rm = TRUE),
                                         by = 100))


m3.0bpreds <- predict.gam(m3.0b, m3.0b_pred_grid, se.fit = TRUE)

m3.0b_sePreds <- data.frame(m3.0b_pred_grid,
                            mu   = binomial()$linkinv(m3.0bpreds$fit),
                            low  = binomial()$linkinv(m3.0bpreds$fit - 1.96 * m3.0bpreds$se.fit),
                            high = binomial()$linkinv(m3.0bpreds$fit + 1.96 * m3.0bpreds$se.fit))

### H3.0c: Depth smoothed over xy with shape and intercept variable by species -

m3.0c <-
  bam(Detected ~ 
        # main effects of space, depth, taxon
        ti(utm.lon, utm.lat,
           d=2,
           k=20,
           bs="tp")+
         ti(depth,
            k=5,
            bs="ts")+
        ti(BestTaxon,
           k=16,
           bs="re")+
        # interaction between *everything*
        ti(utm.lon, utm.lat, depth, BestTaxon,
           d=c(2,1,1),
           k=c(20, 5, 16),
           bs=c("tp","ts", "re"))+
        # space-taxon effect
        ti(utm.lon, utm.lat, BestTaxon,
           d=c(2,1),
           k=c(20,16),
           bs=c("tp","re"))+
        # depth-taxon effect
        ti(depth, BestTaxon,
           k=c(5,16),
           bs=c("ts","re")),
      family = "binomial",
      method = "fREML",
      data = detect_data_dens,
      discrete = TRUE)

summary(m3.0c)
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq  p-value    
# ti(utm.lon,utm.lat)                 6.867e+00    8.4  12.809 0.137625    
# ti(depth)                           1.583e-05    4.0   0.000 0.077249 .  
# ti(BestTaxon)                       1.624e+00    2.0   7.366 0.002203 ** 
#   ti(BestTaxon,depth,utm.lon,utm.lat) 3.763e+01  228.0 109.688  < 2e-16 ***
#   ti(BestTaxon,utm.lon,utm.lat)       8.200e+00   55.0  15.952 0.000601 ***
#   ti(depth,BestTaxon)                 7.326e+00   12.0  35.598  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.094   Deviance explained = 26.2%

### WITHOUT DEPTH OR LAT/LON

# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq  p-value    
# ti(BestTaxon)                        1.582      2   6.703  0.00331 ** 
#   ti(BestTaxon,depth,utm.lon,utm.lat) 37.062    228 107.975  < 2e-16 ***
#   ti(BestTaxon,utm.lon,utm.lat)       13.765     57  34.000 1.86e-06 ***
#   ti(depth,BestTaxon)                  7.371     12  36.356  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0915   Deviance explained = 25.7%
# fREML = 7855.4  Scale est. = 1         n = 7740

AIC(m3.0c)
# 1410 with all terms
# 1413 with non-significant terms (depth and lat/lon) removed

gam.check(m3.0c)
# k'    edf k-index p-value
# ti(BestTaxon)                         3.00   1.68      NA      NA
# ti(BestTaxon,depth,utm.lon,utm.lat) 228.00  40.89      NA      NA
# ti(BestTaxon,utm.lon,utm.lat)        57.00  17.08      NA      NA
# ti(BestTaxon,depth)                  12.00   6.95      NA      NA

#mean squared Pearson residual dispersion parameter
sum(residuals(m3.0c, type = "pearson")^2) / df.residual(m3.0c)

detect_data_dens$resid_response <- residuals(m3.0c,
                                             type = "response")

detect_data_dens$resid_pearson <- residuals(m3.0c, type = "pearson")

ggplot(detect_data_dens,
       aes(utm.lon, utm.lat,
           color = resid_pearson), size = 2) +
  geom_point(alpha = 0.7) +
  facet_wrap(~BestTaxon) +
  scale_color_gradient2(midpoint = 0) +
  coord_equal()

### m3.0c predictions ----------------------------------------------------------

m3.0c_pred_grid <- expand_grid(depth = seq(from = 0, to = 500, by = 10),
                               utm.lat = seq(min(detect_data_dens$utm.lat, na.rm = TRUE),
                                         max(detect_data_dens$utm.lat, na.rm = TRUE),
                                         by = 5000),
                               utm.lon = seq(min(detect_data_dens$utm.lon, na.rm = TRUE),
                                         max(detect_data_dens$utm.lon, na.rm = TRUE),
                                         by = 5000),
                               BestTaxon = as.factor(c("Lagenorhynchus obliquidens",
                                                        "Megaptera novaeangliae",
                                                        "Berardius bairdii")))

coords <- m3.0c_pred_grid %>%
  distinct(utm.lon, utm.lat)

coords_sf <- coords %>%
  st_as_sf(coords = c("utm.lon", "utm.lat"), crs = 32610) %>%
  st_transform(4326)

coords <- coords %>%
  mutate(lon = st_coordinates(coords_sf)[,1], 
         lat = st_coordinates(coords_sf)[,2])

m3.0c_pred_grid <- m3.0c_pred_grid %>%
  left_join(coords,
            by = c("utm.lon", "utm.lat"))

m3.0c_pred_grid <- st_as_sf(m3.0c_pred_grid,
                           coords = c("utm.lon", "utm.lat"), crs = 32610)

m3.0c_pred_grid <- m3.0c_pred_grid %>% 
  mutate(utm.lon = st_coordinates(m3.0c_pred_grid)[,1],
         utm.lat = st_coordinates(m3.0c_pred_grid)[,2])

###crop to MURI data extent
#convex hull from MURI detections
dat_hull <- detections_sf %>%
  st_union() %>%
  st_convex_hull()

#crop grid
m3.0c_pred_grid_crop <- m3.0c_pred_grid[dat_hull, ]

m3.0c_pred_grid <- m3.0c_pred_grid_crop %>%
  mutate(utm.lon = st_coordinates(m3.0c_pred_grid_crop)[,1],
         utm.lat = st_coordinates(m3.0c_pred_grid_crop)[,2])

# mask depths that aren't real
m3.0c_pred_grid$seafloor <- terra::extract(bathy_r,
                                                 m3.0c_pred_grid[, c("lon", "lat")])[,2]

m3.0c_pred_grid_trimmed <- m3.0c_pred_grid[is.na(m3.0c_pred_grid$seafloor) | 
                                          m3.0c_pred_grid$depth <= pmax(abs(m3.0c_pred_grid$seafloor), 0),]


# response predictions
m3.0cpreds <- predict.bam(m3.0c, m3.0c_pred_grid_trimmed,
                          se.fit = TRUE)

m3.0c_sePreds <- data.frame(m3.0c_pred_grid_trimmed,
                            mu   = binomial()$linkinv(m3.0cpreds$fit),
                            low  = binomial()$linkinv(m3.0cpreds$fit - 1.96 * m3.0cpreds$se.fit),
                            high = binomial()$linkinv(m3.0cpreds$fit + 1.96 * m3.0cpreds$se.fit),
                            low50  = binomial()$linkinv(m3.0cpreds$fit - 0.674 * m3.0cpreds$se.fit),
                            high50 = binomial()$linkinv(m3.0cpreds$fit + 0.674 * m3.0cpreds$se.fit))


#Marginal effect of depth
# ref_lon <- mean(detect_data_dens$utm.lon, na.rm = TRUE)
# ref_lat <- mean(detect_data_dens$utm.lat, na.rm = TRUE)
# 
# depth_grid <- expand.grid(
#   depth = seq(from = 0, to = 500, by = 10),
#   utm.lon = ref_lon,
#   utm.lat = ref_lat,
#   BestTaxon = as.factor(c("Lagenorhynchus obliquidens",
#                           "Megaptera novaeangliae",
#                           "Berardius bairdii")))
# 
# pred <- predict.bam(m3.0c, newdata = depth_grid, se.fit = TRUE)
# 
# depth_grid$fit <- binomial()$linkinv(pred$fit)
# 
# depth_grid$low <- binomial()$linkinv(pred$fit - 1.96 * pred$se.fit)
# 
# depth_grid$high <- binomial()$linkinv(pred$fit + 1.96 * pred$se.fit)
# 
# ggplot(depth_grid, aes(depth, fit, color = BestTaxon, fill = BestTaxon)) +
#   geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, color = NA) +
#   geom_line(linewidth = 1.2) +
#   labs(y = "Predicted detection probability",
#     x = "Depth") +
#   facet_wrap(~BestTaxon, scales = "free")

#### H3.1: Detection smoothed over depth + density by species -----------------

m3.1 <- bam(Detected ~ s(D, by = BestTaxon, bs="ts") + s(depth, bs = "ts") + BestTaxon,
             family = "binomial",
             method = "fREML",
             data = detect_data_dens,
             discrete = TRUE)

summary(m3.1)
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(D):BestTaxonBerardius bairdii          0.0003021      3   0.00    0.333    
# s(D):BestTaxonLagenorhynchus obliquidens 1.8203224      9  22.63 3.71e-06 ***
#   s(D):BestTaxonMegaptera novaeangliae     0.9846625      6  50.87  < 2e-16 ***
#   s(depth)                                 2.4472630      9  21.17 1.91e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0323   Deviance explained =   10%

AIC(m3.1)
#1570

#mean squared Pearson residual dispersion parameter
sum(residuals(m3.1, type = "pearson")^2) / df.residual(m3.1)

detect_data_dens$resid_pearson <- residuals(m3.1, type = "pearson")

ggplot(detect_data_dens,
       aes(utm.lon, utm.lat,
           color = resid_pearson), size = 2) +
  geom_point(alpha = 0.7) +
  facet_wrap(~BestTaxon) +
  scale_color_gradient2(midpoint = 0) +
  coord_equal()

### m3.1 predictions ----------------------------------------------------------

m3.1_pred_grid <- expand_grid(density_data %>%
                              rename("BestTaxon" = species) %>%
                              select(BestTaxon, mlat, mlon, D),
                              depth = c(0,300,500)) %>%
                  mutate(lat = mlat, lon=mlon)

m3.1_pred_grid <- st_as_sf(m3.1_pred_grid,
                           coords = c("mlon", "mlat"), crs = 4326) %>%
  st_transform(32610)

###crop to MURI data extent
#convex hull from MURI detections
dat_hull <- detections_sf %>%
  st_union() %>%
  st_convex_hull()

#crop grid
m3.1_pred_grid_crop <- m3.1_pred_grid[dat_hull, ]

m3.1_pred_grid <- m3.1_pred_grid_crop %>%
  mutate(utm.lon = st_coordinates(m3.1_pred_grid_crop)[,1],
         utm.lat = st_coordinates(m3.1_pred_grid_crop)[,2])


# mask depths that aren't real
m3.1_pred_grid$seafloor <- terra::extract(bathy_r,
                                           m3.1_pred_grid[, c("lon", "lat")])[,2]

m3.1_pred_grid_trimmed <- m3.1_pred_grid[is.na(m3.1_pred_grid$seafloor) |
                                             m3.1_pred_grid$depth <= pmax(abs(m3.1_pred_grid$seafloor), 0),]

# response predictions
m3.1preds <- predict.bam(m3.1, m3.1_pred_grid_trimmed,
                         se.fit = TRUE)

m3.1_sePreds <- data.frame(m3.1_pred_grid_trimmed,
                           mu   = binomial()$linkinv(m3.1preds$fit),
                           low  = binomial()$linkinv(m3.1preds$fit - 1.96 * m3.1preds$se.fit),
                           high = binomial()$linkinv(m3.1preds$fit + 1.96 * m3.1preds$se.fit),
                           low50  = binomial()$linkinv(m3.1preds$fit - 0.674 * m3.1preds$se.fit),
                           high50 = binomial()$linkinv(m3.1preds$fit + 0.674 * m3.1preds$se.fit))

#Marginal effect of density predictions ----------------------------------------
D_grid_lags <- expand.grid(depth = median(detect_data_dens$depth, na.rm = TRUE),
                           D = seq(min(detect_data_dens %>% 
                filter(BestTaxon == "Lagenorhynchus obliquidens") %>% 
                pull(D), 
              na.rm = TRUE),
          max(detect_data_dens %>% 
                filter(BestTaxon == "Lagenorhynchus obliquidens") %>% 
                pull(D), na.rm = TRUE),
          length.out = 100),
  BestTaxon = as.factor(c("Lagenorhynchus obliquidens"))) %>% 
  mutate(species = BestTaxon) %>%
  mutate(D_jittered = D + runif(n(), 0.001, 0.05))

D_grid_mnov <- expand.grid(depth = median(detect_data_dens$depth, na.rm = TRUE),
                           D = seq(min(detect_data_dens %>% 
                filter(BestTaxon == "Megaptera novaeangliae") %>% 
                pull(D), 
              na.rm = TRUE),
          max(detect_data_dens %>% 
                filter(BestTaxon == "Megaptera novaeangliae") %>% 
                pull(D), na.rm = TRUE),
          length.out = 100),
  BestTaxon = as.factor(c("Megaptera novaeangliae"))) %>% 
  mutate(species = BestTaxon) %>%
  mutate(D_jittered = D + runif(n(), 0.001, 0.005))

D_grid_bbai <- expand.grid(depth = median(detect_data_dens$depth, na.rm = TRUE),
                           D = seq(min(detect_data_dens %>% 
                filter(BestTaxon == "Berardius bairdii") %>% 
                pull(D), 
              na.rm = TRUE),
          max(detect_data_dens %>% 
                filter(BestTaxon == "Berardius bairdii") %>% 
                pull(D), na.rm = TRUE),
          length.out = 100),
  BestTaxon = as.factor(c("Berardius bairdii"))) %>% 
  mutate(species = BestTaxon) %>%
  mutate(D_jittered = D + runif(n(), 0.001, 0.002))

D_grid <- bind_rows(D_grid_mnov, D_grid_bbai, D_grid_lags)

pred_D <- predict.bam(m3.1, newdata = D_grid, se.fit = TRUE)

D_grid$fit <- binomial()$linkinv(pred_D$fit)
D_grid$low <- binomial()$linkinv(pred_D$fit - 1.96 * pred_D$se.fit)
D_grid$high <- binomial()$linkinv(pred_D$fit + 1.96 * pred_D$se.fit)
D_grid$low50  <- binomial()$linkinv(pred_D$fit - 0.674 * pred_D$se.fit)
D_grid$high50 <- binomial()$linkinv(pred_D$fit + 0.674 * pred_D$se.fit)

### Save -----------------------------------------------------------------------

save(detections_sf, m3.0a, m3.0b, m3.0c, m3.1, 
     m3.0b_sePreds, m3.0cpreds, m3.0c_sePreds, m3.1_sePreds, 
     D_grid, detect_data_dens,
     file = "./ProcessedData/m3.0models_preds_0.05degree.Rdata")

###Some plots for m3.1 ---------------------------------------------------------
#Detection rate by depth and density estimate (3D interactive)
# library(plotly)
# 
# plot_ly(
#   m3.1_sePreds,
#   x = ~depth,
#   y = ~D,
#   z = ~mu,
#   color = ~BestTaxon,
#   type = "scatter3d",
#   mode = "markers"
# )

#Marginal effect of depth
# depth_grid <- expand.grid(
#   depth = seq(0, 500, by = 5),
#   D = median(detect_data_dens$D, na.rm = TRUE),
#   BestTaxon = as.factor(c("Lagenorhynchus obliquidens",
#                           "Megaptera novaeangliae",
#                           "Berardius bairdii")))
# 
# pred_depth <- predict.bam(m3.1, newdata = depth_grid, se.fit = TRUE)
# 
# depth_grid$fit <- binomial()$linkinv(pred_depth$fit)
# depth_grid$low <- binomial()$linkinv(pred_depth$fit - 1.96 * pred_depth$se.fit)
# depth_grid$high <- binomial()$linkinv(pred_depth$fit + 1.96 * pred_depth$se.fit)
# 
# ggplot(depth_grid, aes(depth, fit)) +
#   geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2) +
#   geom_line(linewidth = 1) +
#   facet_wrap(~BestTaxon) +
#   labs(y = "Predicted detection probability")

### Big ole' plot --------------------------------------------------------------
# Create coastline shapefile 
# world <- st_read("Data/ne_10m_land/ne_10m_land.shp")
# 
# data_bbox <- st_as_sf(detect_data_dens, coords = c("lon", "lat"), 
#                       crs = 4326) %>%
#   st_bbox() %>%
#   st_as_sfc() %>% 
#   st_buffer(dist = 70000)
# 
# westcoast_land <- st_crop(world, data_bbox)
# 
# g <- m3.1_sePreds %>%
# #  arrange(BestTaxon, utm.lon, utm.lat, desc(mu)) %>%
#   mutate(depth = as.factor(depth)) %>% 
#   group_by(BestTaxon) %>% 
#   filter(abs(scale(mu)) < 3) %>% 
#   mutate(mu = case_when(mu < 0.005~NA,
#                         TRUE~mu)) %>% 
#   ungroup()
# 
# pp <- ggplot() +
#   geom_sf(data = westcoast_land, fill = "grey50", colour = NA) +
#   geom_tile(aes(x=lon, y=lat, fill=mu), data=g) +
# #  geom_sf(aes(colour = mu), data=st_as_sf(g, coords=c("lon", "lat")), geom="tile") +
#   facet_grid(BestTaxon~depth) +
#   theme_minimal() +
#   scale_fill_viridis_c(name = "Detection rate",
#                        option = "mako",
#                        #trans = "reverse",
#                        begin = 0.4, end = 0.9, na.value = "grey70")
# 
# ggsave(pp, file="thing1.pdf", width=10, height=16)
# 
# ###m3.0c
# g2 <- m3.0c_sePreds %>%
#   filter(depth %in% c(0, 300, 500)) %>% 
#   mutate(depth = as.factor(depth)) %>% 
#   group_by(BestTaxon) %>% 
#   filter(abs(scale(mu)) < 5) %>% 
#   mutate(mu = case_when(mu < 0.001~NA,
#                         TRUE~mu)) %>% 
#   ungroup() #%>% 
# 
# g2_sf <- g2 %>%
#   st_as_sf(coords = c("utm.lon", "utm.lat"), crs = 32610, remove = FALSE)
# 
# westcoast_land_utm <- st_transform(westcoast_land, 32610)
# 
# detect_pts <- detect_data_dens %>%
#   filter(Detected == 1) %>%
#   st_as_sf(coords = c("utm.lon", "utm.lat"), crs = 32610)
# 
# pp2 <- ggplot() +
#   #geom_sf(data = g2_sf, aes(color = mu)) +
#   geom_raster(aes(x=utm.lon, y=utm.lat, fill=mu), data=g2_sf, interpolate = FALSE) +
#    #geom_tile(aes(x=lon, y=lat, fill=mu), width = 0.4,
#    #          height = 0.35, data=g2) +
#   scale_fill_viridis_c(name = "Detection rate",
#                        option = "mako",
#                        #trans = "reverse",
#                        begin = 0.4, end = 0.9, na.value = "grey70") +
#   geom_sf(data = westcoast_land_utm, fill = "grey50", colour = NA, inherit.aes = FALSE) +
#   # geom_sf(data = detect_pts,
#   #         size = 1, alpha = 0.8, stroke = 1, inherit.aes = FALSE) +
#   facet_grid(BestTaxon~depth) +
#   theme_minimal() 
# 
# ggsave(pp2, file="thing2.pdf", width=10, height=16)

