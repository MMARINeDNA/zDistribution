#### 3D Distribution
#### Models testing Q3: Depth distribution across xy
#### Fall 2025
#### EKJ&AVC

library(mgcv)
library(tidyverse)
library(PNWColors)

load("./ProcessedData/detect_data.RData")
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")

# Q3: Does depth distribution of detections vary across xy spatial distribution?
# H3.0: Depth distribution of detections does not vary across xy spatial distribution.
  #H3.0a: xy variability + z variability
  #H3.0b: xyz covariability
  #H3.0c: xyz covariability by species
# H3.1: Depth distribution of detections does not correspond with oceanographic variability.
  #H3.1a: z and oceanography covariability
  #H3.1b: z and oceanography covariability by species

### H3.0a: Depth distribution plus xy space distribution -----------------------

m3.0a <- bam(Detected ~ s(depth) + s(utm.lat, utm.lon),
             family = "binomial",
             data = detect_data,
             method = "fREML",
             discrete = TRUE,
             nthreads = 40)
summary(m3.0a)
#both depth and xy are significant
AIC(m3.0a)
#4708

### H3.0b: Depth distribution smoothed by xy space distribution ----------------

m3.0b <- bam(Detected ~ te(depth, lat, lon),
             family = "binomial",
             data = detect_data,
             method = "fREML",
             discrete = TRUE,
             nthreads = 40)
summary(m3.0b)
#te(depth, utm.lat, utm.lon) is significant
AIC(m3.0b)
#4669

### m3.0b predictions ----------------------------------------------------------

m3.0b_pred_grid <- expand_grid(depth = seq(0,500, by = 100),
                               lat = seq(min(detect_data$lat, na.rm = TRUE),
                                         max(detect_data$lat, na.rm = TRUE),
                                         by = 0.1),
                               lon = seq(min(detect_data$lon, na.rm = TRUE),
                                         max(detect_data$lon, na.rm = TRUE),
                                         by = 0.1))


m3.0bpreds <- predict.gam(m3.0b, m3.0b_pred_grid, se.fit = TRUE)

m3.0b_sePreds <- data.frame(m3.0b_pred_grid,
                            mu   = exp(m3.0bpreds$fit),
                            low  = exp(m3.0bpreds$fit - 1.96 * m3.0bpreds$se.fit),
                            high = exp(m3.0bpreds$fit + 1.96 * m3.0bpreds$se.fit))

### H3.0c: Depth smoothed over xy with shape and intercept variable by species -

Sys.time()
m3.0c <- bam(Detected ~ te(depth, lat, lon, 
                           k=c(10,5,5),
                           by = as.factor(BestTaxon)),
             family = "binomial",
             method = "fREML",
             data = detect_data,
             discrete = TRUE,
             nthreads = 40)
Sys.time()

save(m3.0c, file = "./ProcessedData/m3.0c.RData")

### m3.0c predictions ----------------------------------------------------------

m3.0c_pred_grid <- expand_grid(depth = seq(0,500, by = 100),
                               lat = seq(min(detect_data$lat, na.rm = TRUE),
                                         max(detect_data$lat, na.rm = TRUE),
                                         by = 0.1),
                               lon = seq(min(detect_data$lon, na.rm = TRUE),
                                         max(detect_data$lon, na.rm = TRUE),
                                         by = 0.1),
                               BestTaxon = as.factor(unique(detect_data$BestTaxon)))


m3.0cpreds <- predict.bam(m3.0c, m3.0c_pred_grid, 
                          se.fit = TRUE,
                          discrete = TRUE,
                          n.threads = 40)

m3.0c_sePreds <- data.frame(m3.0c_pred_grid,
                            mu   = exp(m3.0cpreds$fit),
                            low  = exp(m3.0cpreds$fit - 1.96 * m3.0cpreds$se.fit),
                            high = exp(m3.0cpreds$fit + 1.96 * m3.0cpreds$se.fit))

save(m3.0cpreds,m3.0c_sePreds, file = "m3.0c_preds.Rdata")

### Halfway save ---------------------------------------------------------------

save(m3.0a, m3.0b, m3.0c, m3.0b_sePreds, m3.0c_sePreds, 
     file = "m3.0models_preds.Rdata")

### NOT RUN:

### H3.1a: Depth and oceanographic variables -----------------------------------

# m3.1a <- bam(Detected ~ te(depth, SST, bottom_depth,
#                           k = c(10,5,5)),
#                 family = "binomial",
#                 data = "REML",
#                 discrete = TRUE,
#                 nthreads = 40)
# 
# summary(m3.1a)
# #
# AIC(m3.1a)
# #
# 
# save(m3.1a, file = "./ProcessedData/m3.1a.RData")

### H3.1b: Depth and oceanographic variables by species ------------------------ 

# m3.1b <- bam(Detected ~ te(depth, SST, bottom_depth,     
#                            k=c(10,5,5),  
#                            by = as.factor(BestTaxon)), 
#              family = "binomial",   
#              data = detect_data,  
#              discrete = TRUE,
#              nthreads = 40)   
# 
# save(m3.1b, file = "m3.1b.RData")   
# 
# summary(m3.1b)
# #
# AIC(m3.1b)
# #
# 
# ### M3.1b predictions ----------------------------------------------------------
# 
# m3.1b_pred_grid <- expand_grid(depth = seq(0,500, by = 100),
#                                SST = seq(min(detect_data$SST, na.rm = TRUE),
#                                          max(detect_data$SST, na.rm = TRUE),
#                                          by = 0.1),
#                                lon = seq(min(detect_data$lon, na.rm = TRUE),
#                                          max(detect_data$lon, na.rm = TRUE),
#                                          by = 0.1),
#                                BestTaxon = as.factor(unique(detect_data$BestTaxon)))
# 
# 
# m3.1bpreds <- predict.bam(m3.1b, m3.1b_pred_grid,
#                           se.fit = TRUE,
#                           discrete = TRUE,
#                           n.threads = 40)
# 
# m3.1b_sePreds <- data.frame(m3.1b_pred_grid,
#                             mu   = exp(m3.1bpreds$fit),
#                             low  = exp(m3.1bpreds$fit - 1.96 * m3.1bpreds$se.fit),
#                             high = exp(m3.1bpreds$fit + 1.96 * m3.1bpreds$se.fit))
# 
# save(m3.1bpreds,m3.1b_sePreds, file = "m3.1b_preds.Rdata")

#TODO:
#1. Which oceanographic variables to use?
#2. Pull from sat data
#3.  Run 3.1 models
