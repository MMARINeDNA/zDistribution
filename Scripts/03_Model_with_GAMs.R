
library(mgcv)
library(tidyverse)

load("./ProcessedData/detect_species_meta.RData")

# Q1: Does the probability of a detecting cetaceans in eDNA samples vary with sample depth?
# H0: Probability of detection does not vary with depth.
# H1: Probability of detection varies across depth agnostic to species or functional group.
# H2: Probability of detection varies across depth according to species or functional group.

# basic model with no species-specific terms
m1.0 <- gam(Detected ~ s(depth), family = "binomial", data = detect_species_meta) 
# depth p-value = 0.085
# AIC 3694.004
m1.0_predictions <- data.frame(depth = 0:500)
m1.0_predictions$pred <- predict.gam(m1.0, m1.0_predictions, type = "response")
plot(m1.0_predictions$depth, m1.0_predictions$pred, 
     type = "l", xlab = "Depth", ylab = "P(Detection)")

# this model will have a different intercept for each species, but spline will be same shape
m1.1 <- gam(Detected ~ s(depth) + BestTaxon, family = "binomial", data = detect_species_meta)
m1.1_predictions <- expand_grid(depth = 0:500, BestTaxon = as.factor(unique(detect_species_meta$BestTaxon)))
m1.1_predictions$pred <- predict.gam(m1.1, m1.1_predictions, type = "response")

ggplot(m1.1_predictions) +
  geom_line(aes(x=depth, y = pred, group = BestTaxon)) +
  xlab("Depth")+
  ylab("P(Detection)")+
  theme_bw()

# this model will have separate smooths for each species
m1.2 <- gam(Detected ~ s(depth, by = as.factor(BestTaxon)), family = "binomial", data = detect_species_meta)
# Depth significant for some taxa
# AIC 3355.708
m1.2_predictions <- expand_grid(depth = 0:500, BestTaxon = as.factor(unique(detect_species_meta$BestTaxon)))
m1.2_predictions$pred <- predict.gam(m1.2, m1.2_predictions, type = "response")

ggplot(m1.2_predictions) +
  geom_line(aes(x=depth, y = pred, group = BestTaxon)) +
  xlab("Depth")+
  ylab("P(Detection)")+
  theme_bw()

save(m1.2, file = "./ProcessedData/m1.2.RData")
