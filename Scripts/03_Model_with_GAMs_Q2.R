library(mgcv)
library(tidyverse)
library(PNWColors)

load("./ProcessedData/detect_species_meta.RData")

# Q2: Does the probability of a detecting cetaceans in eDNA samples change with the number of technical replicates?

# H0: Detection does not vary with # technical replicates.
# H1: Detection varies with # technical replicated agnostic to species/fuctional group or depth.

m2.0 <- gam(Detected ~ nTechReps, 
            family = "binomial", data = detect_species_meta)

# H2: Detection varies with # technical replicates according to species/functional 
# group, depth, or a combination of the two.

m2.1 <- gam(Detected ~ s(depth) + nTechReps, 
            family = "binomial", data = detect_species_meta)

m2.2 <- gam(Detected ~ s(depth) + BestTaxon + nTechReps,
            family = "binomial", data = detect_species_meta)

m2.3 <- gam(Detected ~ s(depth, by = as.factor(BestTaxon)) + nTechReps, 
            family = "binomial", data = detect_species_meta)

save(m2.0, m2.1, m2.2, m2.3,
     file = "./ProcessedData/m2models.RData")
