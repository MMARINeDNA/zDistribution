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
# H3.1: Depth distribution of detection does vary across xy spatial distribution agnostic to oceanography (e.g. upwelling).
# H3.2: Depth distribution of detection does vary across xy spatial distribution according to oceanography (e.g. upwelling).

### H3.0a: Depth distribution plus xy space distribution -----------------------

m3.0a <- gam(Detected ~ s(depth) + s(utm.lat, utm.lon), 
          family = "binomial",
          data = detect_data)
summary(m3.0a)
#both depth and xy are significant
AIC(m3.0a)
#4708

### H3.0b: Depth distribution smoothed by xy space distribution ----------------

m3.0b <- gam(Detected ~ te(depth, utm.lat, utm.lon), 
             family = "binomial",
             data = detect_data)
summary(m3.0b)
#te(depth, utm.lat, utm.lon) is significant
AIC(m3.0b)
#4669

### H3.0c: Depth smoothed over xy with shape and intercept variable by species -

m3.0c <- gam(Detected ~ te(depth, utm.lat, utm.lon, by = as.factor(BestTaxon)) + 
               as.factor(BestTaxon),
          family = "binomial",
          method = "REML",
          data = detect_data)
          #nthreads = 40) #might want to start running these on Hyak

summary(m3.0c)
#
AIC(m3.0c)
#

save(m3.0c, file = "./ProcessedData/m3.0c.RData")
### H3.1a: Depth and oceanographic variables plus xy distribution --------------

m3.1a <- gam(Detected ~ te(depth, SST, bottom_depth) + s(utm.lat, utm.lon),
                family = "binomial",
                data = "REML")

summary(m3.1a)
#
AIC(m3.1a)
#

save(m3.1a, file = "./ProcessedData/m3.1a.RData")

### H3.1b: Depth and oceanographic variables by species ------------------------

m3.1b <- gam(Detected ~ te(depth, SST, bottom_depth, by = as.factor(BestTaxon)) + 
               as.factor(BestTaxon) + s(utm.lat, utm.lon), 
             #remove s(utm.lat, utm.lon) if it's not significant in m3.1a
             family = "binomial",
             data = "REML")
              #nthreads = 40) #might want to start running these on Hyak

summary(m3.1b)
#
AIC(m3.1b)
#

#TODO:
#1. Run on Hyak
#2. Which oceanographic variables to use?
