#### 3D Distribution
#### Models testing Q1: detection probability across depth
#### January2025
#### EKJ&AVC

library(mgcv)
library(tidyverse)
library(PNWColors)

load("./ProcessedData/detect_species_meta.RData")
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")

# Q1: Does the probability of a detecting cetaceans in eDNA samples vary with sample depth?
# H0: Probability of detection does not vary with depth.
# H1: Probability of detection varies across depth agnostic to species or functional group.
# H2: Probability of detection varies across depth according to species or functional group.

### H1: POD by depth alone -----------------------------------------------------
# basic model with no species-specific terms
m1.0 <- gam(Detected ~ s(depth), family = "binomial", data = detect_species_meta) 
# depth p-value = 0.085
# AIC 3694.004
m1.0_predictions <- data.frame(depth = 0:500)
m1.0_predictions$pred <- predict.gam(m1.0, m1.0_predictions, type = "response")
plot(m1.0_predictions$depth, m1.0_predictions$pred, 
     type = "l", xlab = "Depth", ylab = "P(Detection)")

m1.0preds <- predict(m1.0, m1.0_predictions, se.fit = TRUE)

m1.0_sePreds <- data.frame(m1.0_predictions,
                           mu   = exp(m1.0preds$fit),
                           low  = exp(m1.0preds$fit - 1.96 * m1.0preds$se.fit),
                           high = exp(m1.0preds$fit + 1.96 * m1.0preds$se.fit))

ggplot(m1.0_sePreds, aes(x = depth, y = mu)) +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", 
              alpha = 0.2, color = "#74677e", fill = "#74677e") +
  ylab("POD") +
  theme_minimal()

### H2: POD by depth across species --------------------------------------------
# this model will have a different intercept for each species, but spline will be same shape
m1.1 <- gam(Detected ~ s(depth) + BestTaxon, family = "binomial", data = detect_species_meta)
m1.1_predictions <- expand_grid(depth = 0:500, BestTaxon = as.factor(unique(detect_species_meta$BestTaxon)))
m1.1_predictions$pred <- predict.gam(m1.1, m1.1_predictions, type = "response")
# AIC = 3257

ggplot(m1.1_predictions) +
  geom_line(aes(x=depth, y = pred, group = BestTaxon)) +
  xlab("Depth")+
  ylab("P(Detection)")+
  theme_bw()

m1.1preds <- predict(m1.1, m1.1_predictions, se.fit = TRUE)
m1.1_sePreds <- data.frame(m1.1_predictions,
                           mu   = exp(m1.1preds$fit),
                           low  = exp(m1.1preds$fit - 1.96 * m1.1preds$se.fit),
                           high = exp(m1.1preds$fit + 1.96 * m1.1preds$se.fit)) %>% 
  left_join(mmEcoEvo, by = c("BestTaxon" = "Species"))

ggplot(m1.1_sePreds, aes(x = depth, y = mu, color = BestTaxon, fill = BestTaxon)) +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", 
              alpha = 0.2) +
  ylab("POD") +
  scale_fill_manual(values = c(pnw_palette("Cascades",12, type = "continuous"),
                               pnw_palette("Sunset",12, type = "continuous")[1:12])) +
  scale_color_manual(values = c(pnw_palette("Cascades",12, type = "continuous"),
                                pnw_palette("Sunset",12, type = "continuous")[1:12])) +
  facet_wrap(~Family, scales = "free_y") +
  coord_cartesian(ylim = c(0,0.1)) +
  theme_minimal()


# this model will have separate smooths for each species
m1.2 <- gam(Detected ~ s(depth, by = as.factor(BestTaxon)), family = "binomial", data = detect_species_meta)
# Depth significant for some taxa
# AIC 3355.708
m1.2_predictions <- expand_grid(depth = 0:500, BestTaxon = as.factor(unique(detect_species_meta$BestTaxon)))
m1.2_predictions$pred <- predict.gam(m1.2, m1.2_predictions, type = "response")

ggplot(m1.2_predictions) +
  geom_line(aes(x=depth, y = pred, group = BestTaxon, color = BestTaxon)) +
  xlab("Depth")+
  ylab("P(Detection)")+
  theme_bw()

m1.2preds <- predict(m1.2, m1.2_predictions, se.fit = TRUE)
m1.2_sePreds <- data.frame(m1.2_predictions,
                      mu   = exp(m1.2preds$fit),
                      low  = exp(m1.2preds$fit - 1.96 * m1.2preds$se.fit),
                      high = exp(m1.2preds$fit + 1.96 * m1.2preds$se.fit)) %>% 
  left_join(mmEcoEvo, by = c("BestTaxon" = "Species"))

ggplot(m1.2_sePreds, aes(x = depth, color = BestTaxon, fill = BestTaxon)) +
  geom_line(aes(y = mu)) +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity") +
  scale_fill_manual(values = c(pnw_palette("Cascades",9, type = "continuous"),
                               pnw_palette("Sunset",9, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Cascades",9, type = "continuous"),
                                pnw_palette("Sunset",9, type = "continuous"))) +
  facet_wrap(~BestTaxon, scales = "free_y") +
  geom_rug(data = detect_species_meta, aes(x=depth), color = "grey")+
  geom_rug(data = filter(detect_species_meta, Detected == 1), aes(x=depth))+
  
  #coord_cartesian(ylim = c(0,0.25)) +
  theme_minimal() +
  theme(legend.position = "none")

#save(m1.2, file = "./ProcessedData/m1.2.RData")

### H2a: POD by depth across taxonomic family ----------------------------------
m1.2a <- gam(Detected ~ s(depth, by = as.factor(Family)), 
            family = "binomial", data = detect_species_meta)
summary(m1.2a)
#significant for some families but not e.g. grey whales or bowhead whales (too few detections?)
AIC(m1.2a)
#AIC 4057.063: fit seems better by species than family

m1.2a_predictions <- expand_grid(depth = 0:500, Family = as.factor(unique(detect_species_meta$Family)))
m1.2a_predictions$pred <- predict.gam(m1.2a, m1.2a_predictions, type = "response")

ggplot(m1.2a_predictions) +
  geom_line(aes(x=depth, y = pred, group = Family, color = Family)) +
  xlab("Depth")+
  ylab("POD")+
  theme_bw()

m1.2apreds <- predict(m1.2a, m1.2a_predictions, se.fit = TRUE)
m1.2a_sePreds <- data.frame(m1.2a_predictions,
                           mu   = exp(m1.2apreds$fit),
                           low  = exp(m1.2apreds$fit - 1.96 * m1.2apreds$se.fit),
                           high = exp(m1.2apreds$fit + 1.96 * m1.2apreds$se.fit))

ggplot(m1.2a_sePreds, aes(x = depth, y = mu, color = Family, fill = Family)) +
  geom_point() +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity") +
  scale_fill_manual(values = c(pnw_palette("Bay",8, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Bay",8, type = "continuous"))) +
  facet_wrap(~Family, scales = "free_y") +
  #coord_cartesian(ylim = c(0,0.25)) +
  theme_minimal()

### H2b:POD by depth across prey category --------------------------------------

m1.2b <- gam(Detected ~ s(depth, by = as.factor(Prey.family)), 
             family = "binomial", data = detect_species_meta)
summary(m1.2b)
#significant for all three types
AIC(m1.2b)
#AIC 4075.112: even worse!

m1.2b_predictions <- expand_grid(depth = 0:500, Prey.family = as.factor(unique(detect_species_meta$Prey.family)))
m1.2b_predictions$pred <- predict.gam(m1.2b, m1.2b_predictions, type = "response")

ggplot(m1.2b_predictions) +
  geom_line(aes(x=depth, y = pred, group = Prey.family, color = Prey.family)) +
  xlab("Depth")+
  ylab("POD")+
  theme_bw()

m1.2bpreds <- predict(m1.2b, m1.2b_predictions, se.fit = TRUE)
m1.2b_sePreds <- data.frame(m1.2b_predictions,
                            mu   = exp(m1.2bpreds$fit),
                            low  = exp(m1.2bpreds$fit - 1.96 * m1.2bpreds$se.fit),
                            high = exp(m1.2bpreds$fit + 1.96 * m1.2bpreds$se.fit))

ggplot(m1.2b_sePreds, aes(x = depth, y = mu, color = Prey.family, fill = Prey.family)) +
  geom_line() +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity") +
  scale_fill_manual(values = c(pnw_palette("Moth",5, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Moth",5, type = "continuous"))) +
  facet_wrap(~Prey.family, scales = "free_y") +
  theme_minimal() +
  xlab("Depth")+
  ylab("POD")

### H2c: POD by time-at-depth --------------------------------------------------

m1.2c <- gam(Detected ~ s(time_per_m), 
             family = "binomial", data = detect_species_divetime)
summary(m1.2c)
#p<2e-16
AIC(m1.2c)
#AIC 3185.768 - still better across depth by individual species

m1.2c_predictions <- data.frame(time_per_m = seq(min(detect_species_divetime$time_per_m),max(detect_species_divetime$time_per_m), by = 0.1))
m1.2c_predictions$pred <- predict.gam(m1.2c, m1.2c_predictions, type = "response")

ggplot(m1.2c_predictions) +
  geom_line(aes(x=time_per_m, y = pred)) +
  xlab("Time at Depth")+
  ylab("POD")+
  theme_bw()

m1.2cpreds <- predict(m1.2c, m1.2c_predictions, se.fit = TRUE)

m1.2c_sePreds <- data.frame(m1.2c_predictions,
                           mu   = exp(m1.2cpreds$fit),
                           low  = exp(m1.2cpreds$fit - 1.96 * m1.2cpreds$se.fit),
                           high = exp(m1.2cpreds$fit + 1.96 * m1.2cpreds$se.fit))

ggplot(m1.2c_sePreds, aes(x = time_per_m, y = mu)) +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", 
              alpha = 0.2, color = "#74677e", fill = "#74677e") +
  ylab("POD") +
  xlab("Time at depth") +
  theme_minimal()

### H2d: POD by time-at-depth across Family ------------------------------------

m1.2d <- gam(Detected ~ s(time_per_m, by = as.factor(Family)), 
             family = "binomial", data = detect_species_divetime)
summary(m1.2d)
#significant for all but bowhead and grey whale
AIC(m1.2d)
#3732.61 - still higher than depth by species

m1.2d_predictions <- expand_grid(time_per_m = seq(min(detect_species_divetime$time_per_m),max(detect_species_divetime$time_per_m), by = 0.1),
                                 Family = as.factor(unique(detect_species_divetime$Family)))
m1.2d_predictions$pred <- predict.gam(m1.2d, m1.2d_predictions, type = "response")

ggplot(m1.2d_predictions) +
  geom_line(aes(x=time_per_m, y = pred, group = Family, color = Family)) +
  xlab("Time at depth")+
  ylab("POD")+
  theme_bw()

m1.2dpreds <- predict(m1.2d, m1.2d_predictions, se.fit = TRUE)

m1.2d_sePreds <- data.frame(m1.2d_predictions,
                            mu   = exp(m1.2dpreds$fit),
                            low  = exp(m1.2dpreds$fit - 1.96 * m1.2dpreds$se.fit),
                            high = exp(m1.2dpreds$fit + 1.96 * m1.2dpreds$se.fit))

ggplot(m1.2d_sePreds, aes(x = time_per_m, y = mu, color = Family, fill = Family)) +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", 
              alpha = 0.2) +
  
  ylab("POD") +
  xlab("Time at depth") +
  scale_fill_manual(values = c(pnw_palette("Bay",8, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Bay",8, type = "continuous"))) +
  facet_wrap(~Family, scales = "free_y") +
  theme_minimal()

### H2e: POD by time-at-depth across species -----------------------------------
## This one takes a really long time to run.
m1.2e <- gam(Detected ~ s(time_per_m, by = as.factor(BestTaxon)), 
             family = "binomial", data = detect_species_divetime)
summary(m1.2e)
#significant for some species but not all (e.g. lissos, grey whale, bowhead, minke, 
#blue whale, killer whale, all beaked whales, harbor seal)
AIC(m1.2e)
#3428.249 - better than by family but still not as good as by depth

m1.2e_predictions <- expand_grid(time_per_m = min(detect_species_divetime$time_per_m):max(detect_species_divetime$time_per_m),
                                 BestTaxon = as.factor(unique(detect_species_divetime$BestTaxon)))
m1.2e_predictions$pred <- predict.gam(m1.2e, m1.2e_predictions, type = "response")

ggplot(m1.2e_predictions) +
  geom_line(aes(x=time_per_m, y = pred, group = BestTaxon, color = BestTaxon)) +
  xlab("Time at depth")+
  ylab("POD")+
  theme_bw()

m1.2epreds <- predict(m1.2e, m1.2e_predictions, se.fit = TRUE)

m1.2e_sePreds <- data.frame(m1.2e_predictions,
                            mu   = exp(m1.2epreds$fit),
                            low  = exp(m1.2epreds$fit - 1.96 * m1.2epreds$se.fit),
                            high = exp(m1.2epreds$fit + 1.96 * m1.2epreds$se.fit))

ggplot(m1.2e_sePreds, aes(x = time_per_m, y = mu, color = BestTaxon, fill = BestTaxon)) +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", 
              alpha = 0.2) +
  
  ylab("POD") +
  xlab("Time at depth") +
  scale_fill_manual(values = c(pnw_palette("Bay",24, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Bay",24, type = "continuous"))) +
  facet_wrap(~BestTaxon, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none")

### H2f: POD by depth + time-at-depth across species ---------------------------
## This one takes a really long time to run.
m1.2f <- gam(Detected ~ s(time_per_m, by = as.factor(BestTaxon)) +
               s(depth, by = as.factor(BestTaxon)), 
             family = "binomial", data = detect_species_divetime)
summary(m1.2f)

AIC(m1.2f)
# 3477.885

m1.2f_predictions <- expand_grid(depth = 0:500, time_per_m = min(detect_species_divetime$time_per_m):max(detect_species_divetime$time_per_m),
                                 BestTaxon = as.factor(unique(detect_species_divetime$BestTaxon)))
m1.2f_predictions$pred <- predict.gam(m1.2f, m1.2f_predictions, type = "response")

ggplot(m1.2f_predictions) +
  geom_line(aes(x=time_per_m, y = pred, group = BestTaxon, color = BestTaxon)) +
  xlab("Time at depth")+
  ylab("POD")+
  theme_bw()

m1.2fpreds <- predict(m1.2f, m1.2f_predictions, se.fit = TRUE)

m1.2f_sePreds <- data.frame(m1.2f_predictions,
                            mu   = exp(m1.2fpreds$fit),
                            low  = exp(m1.2fpreds$fit - 1.96 * m1.2fpreds$se.fit),
                            high = exp(m1.2fpreds$fit + 1.96 * m1.2fpreds$se.fit)) %>% 
  left_join(mmEcoEvo, by = c("BestTaxon" = "Species"))

ggplot(m1.2f_sePreds, aes(x = depth, y = mu, color = BestTaxon, fill = BestTaxon)) +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", 
              alpha = 0.2) +
  
  ylab("POD") +
  xlab("Depth") +
  scale_fill_manual(values = c(pnw_palette("Bay",24, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Bay",24, type = "continuous"))) +
  facet_wrap(Broad_taxa~BestTaxon, scales = "free") +
  theme_minimal() +
  theme(legend.position = "none")

### H2g: POD by time-at-depth across depths ------------------------------------

bins <- c(0, 100, 300, 500)
detect_species_depthbin <- detect_species_divetime %>% 
  mutate(depth_bin = cut(depth, bins, include.lowest = TRUE))

m1.2g <- gam(Detected ~ s(time_per_m, by = as.factor(depth_bin)), 
             family = "binomial", data = detect_species_depthbin)
summary(m1.2g)

AIC(m1.2g)
# 3876.487

m1.2g2 <- gam(Detected ~ s(time_per_m, by = depth), 
             family = "binomial", data = detect_species_divetime)
summary(m1.2g2)

AIC(m1.2g2)
# 3991.547
m1.2g_predictions <- expand_grid(time_per_m = min(detect_species_depthbin$time_per_m):max(detect_species_depthbin$time_per_m),
                                 depth_bin = as.factor(unique(detect_species_depthbin$depth_bin)))
m1.2g_predictions$pred <- predict.gam(m1.2g, m1.2g_predictions, type = "response")


m1.2gpreds <- predict(m1.2g, m1.2g_predictions, se.fit = TRUE)

m1.2g_sePreds <- data.frame(m1.2g_predictions,
                            mu   = exp(m1.2gpreds$fit),
                            low  = exp(m1.2gpreds$fit - 1.96 * m1.2gpreds$se.fit),
                            high = exp(m1.2gpreds$fit + 1.96 * m1.2gpreds$se.fit))

#very large uncertainty bands
ggplot(m1.2g_sePreds, aes(x = time_per_m, y = mu, color = depth_bin, fill = depth_bin)) +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", 
              alpha = 0.2) +
  
  ylab("POD") +
  xlab("Time at depth") +
  scale_fill_manual(values = c(pnw_palette("Bay",5, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Bay",5, type = "continuous"))) +
  facet_wrap(~depth_bin, scales = "free_y") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal()

### Aggregate model AIC --------------------------------------------------------

modelAIC <- AIC(m1.0, m1.1, m1.2, m1.2a, m1.2b, m1.2c, m1.2d, m1.2e, m1.2f, m1.2g)

########

save(m1.0, m1.1, m1.2, m1.2a, m1.2b, m1.2c, m1.2d, m1.2e, m1.2f, m1.2g,
     file = "ProcessedData/H1models.RData")
