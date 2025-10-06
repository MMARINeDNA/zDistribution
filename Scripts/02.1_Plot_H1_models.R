#### 3D Distribution
#### Figures for Q1 models
#### October 2025
#### AVC

library(tidyverse)
library(PNWColors)

load("./ProcessedData/H1models.RData")
load("./ProcessedData/detect_data.RData")
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")
timeAtDepth <- read.csv("./Data/MM_dive_time_expand.csv")


### m1.1 figure ----------------------------------------------------------------

m1.1_predictions <- expand_grid(depth = 0:500, BestTaxon = as.factor(unique(detect_data$BestTaxon)))

m1.1preds <- predict(m1.1, m1.1_predictions, type = "response", se.fit = TRUE)
m1.1_sePreds <- data.frame(m1.1_predictions,
                           mu   = exp(m1.1preds$fit),
                           low  = exp(m1.1preds$fit - 1.96 * m1.1preds$se.fit),
                           high = exp(m1.1preds$fit + 1.96 * m1.1preds$se.fit)) %>% 
  left_join(mmEcoEvo, by = c("BestTaxon" = "Species"))

m1.1POD <- ggplot(m1.1_sePreds, aes(x = depth, color = Broad_taxa, fill = Broad_taxa)) +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity") +
  ylab("POD") +
  scale_fill_manual(values = c(pnw_palette("Cascades",2, type = "continuous"),
                               pnw_palette("Sunset",2, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Cascades",2, type = "continuous"),
                                pnw_palette("Sunset",2, type = "continuous"))) +
  facet_wrap(~abbrev, scales = "free_y") +
  geom_rug(data = detect_data, aes(x=depth), color = "grey")+
  geom_rug(data = filter(detect_data, Detected == 1), aes(x=depth))+
  theme_minimal() +
  theme(legend.position = "bottom")

png(file = "./Figures/m1.1POD.png")
m1.1POD
dev.off()

### m1.2 figure ----------------------------------------------------------------

m1.2_predictions <- expand_grid(depth = 0:500, BestTaxon = as.factor(unique(detect_data$BestTaxon)))

m1.2preds <- predict(m1.2, m1.2_predictions, type = "response", se.fit = TRUE)
m1.2_sePreds <- data.frame(m1.2_predictions,
                           mu   = exp(m1.2preds$fit),
                           low  = exp(m1.2preds$fit - 1.96 * m1.2preds$se.fit),
                           high = exp(m1.2preds$fit + 1.96 * m1.2preds$se.fit)) %>% 
  left_join(mmEcoEvo, by = c("BestTaxon" = "Species")) %>% 
  mutate(common_name = case_when(common_name == "killer whale"~"mammal eating killer whale",
                                 TRUE~common_name)) %>% 
  mutate(depth = case_when(depth == 0~1,
                           TRUE~depth)) %>%
  left_join(timeAtDepth, by = c("common_name" = "Species", "depth" = "depth")) %>% 
  left_join(maxDepth_species, by = c("common_name" = "Species")) %>% 
  mutate(time_10m = case_when(depth > maxDepth~0,
                              TRUE~time_10m))

# transform time_10m to plot with POD by depth
m1.2_scaled <- m1.2_sePreds %>%
  group_by(abbrev) %>%
  mutate(time_scaled = time_10m/545) %>% 
  mutate(
    a = diff(range(mu, na.rm = TRUE)) / diff(range(time_10m, na.rm = TRUE)),
    b = min(mu, na.rm = TRUE) - a * min(time_10m, na.rm = TRUE),
    time_scaled = a * time_10m + b
  )

m1.2POD <- ggplot(m1.2_scaled, aes(x = depth, color = Broad_taxa, fill = Broad_taxa)) +
  geom_line(aes(y = mu)) +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity") +
  geom_smooth(aes(y = time_scaled), linetype = 2) +
  scale_fill_manual(values = c(pnw_palette("Cascades",2, type = "continuous"),
                               pnw_palette("Sunset",2, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Cascades",2, type = "continuous"),
                                pnw_palette("Sunset",2, type = "continuous"))) +
  facet_wrap(~abbrev, scales = "free_y") +
  geom_rug(data = detect_data, aes(x=depth), color = "grey")+
  geom_rug(data = filter(detect_data, Detected == 1), aes(x=depth))+
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ylab("POD") +
  xlab("Depth")

png(file = "./Figures/m1.2POD.png")
m1.2POD
dev.off()

### m1.2b figure ---------------------------------------------------------------

m1.2b_predictions <- expand_grid(depth = 0:500, Prey.family = as.factor(unique(detect_data$Prey.family)))

m1.2bpreds <- predict(m1.2b, m1.2b_predictions, type = "response", se.fit = TRUE)
m1.2b_sePreds <- data.frame(m1.2b_predictions,
                            mu   = exp(m1.2bpreds$fit),
                            low  = exp(m1.2bpreds$fit - 1.96 * m1.2bpreds$se.fit),
                            high = exp(m1.2bpreds$fit + 1.96 * m1.2bpreds$se.fit))

m1.2bPOD <- ggplot(m1.2b_sePreds, aes(x = depth, color = Prey.family, fill = Prey.family)) +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity") +
  scale_fill_manual(values = c(pnw_palette("Moth",5, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Moth",5, type = "continuous"))) +
  #guides(color = "none", fill = "none") +
  facet_wrap(~Prey.family, scales = "free_y", ncol = 1) +
  geom_rug(data = detect_data, aes(x=depth), color = "grey")+
  geom_rug(data = filter(detect_data, Detected == 1), aes(x=depth))+
  theme_minimal() +
  theme_minimal() +
  xlab("Depth")+
  ylab("POD") +
  theme(legend.position = "none", legend.title = element_blank())

png(file = "./Figures/m1.2bPOD.png")
m1.2bPOD
dev.off()

### m1.2e figure ---------------------------------------------------------------

m1.2e_predictions <- expand_grid(time_10m = seq(min(detect_species_divetime$time_10m),
                                                max(detect_species_divetime$time_10m), 
                                                by = 5),
                                 BestTaxon = as.factor(unique(detect_species_divetime$BestTaxon)))

m1.2epreds <- predict(m1.2e, m1.2e_predictions, type = "response", se.fit = TRUE)

m1.2e_sePreds <- data.frame(m1.2e_predictions,
                            mu   = exp(m1.2epreds$fit),
                            low  = exp(m1.2epreds$fit - 1.96 * m1.2epreds$se.fit),
                            high = exp(m1.2epreds$fit + 1.96 * m1.2epreds$se.fit)) %>% 
  left_join(mmEcoEvo, by = c("BestTaxon" = "Species")) 

m1.2ePOD <- ggplot(m1.2e_sePreds, aes(x = time_10m, color = Broad_taxa, fill = Broad_taxa)) +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", 
              alpha = 0.2) +
  ylab("POD") +
  xlab("Time at depth") +
  facet_wrap(~abbrev, scales = "free_y") +
  scale_fill_manual(values = c(pnw_palette("Cascades",2, type = "continuous"),
                               pnw_palette("Sunset",2, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Cascades",2, type = "continuous"),
                                pnw_palette("Sunset",2, type = "continuous"))) +
  geom_rug(data = detect_species_divetime, aes(x=time_10m), color = "grey20")+
  geom_rug(data = filter(detect_species_divetime, Detected == 1), aes(x=time_10m))+
  theme_minimal() +
  theme(legend.position = "bottom")

png(file = "./Figures/m1.2ePOD.png")
m1.2ePOD
dev.off()

save(m1.1POD, m1.2POD, m1.2bPOD, m1.2ePOD, file = "./Figures/H1PODplots.Rdata")
