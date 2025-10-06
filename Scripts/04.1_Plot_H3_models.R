#### 3D Distribution
#### Figures for Q3 models
#### October 2025
#### AVC

library(tidyverse)
library(PNWColors)
library(mgcv)
library(patchwork)

load("./ProcessedData/m3.0c.RData")
load("./ProcessedData/detect_data.RData")
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")
metadata <- read.csv("./Data/Hake_2019_metadata.csv")

### m3.0c figure ---------------------------------------------------------------

m3.0c_predictions <- expand_grid(depth = 0:500, 
                                 lat = c(34.39501,41.47946,48.5639),
                                 lon = (min(metadata$lon) + max(metadata$lon))/2,
                                 BestTaxon = as.factor(c("Lagenorhynchus obliquidens",
                                                         "Megaptera novaeangliae",
                                                         "Berardius bairdii")))

m3.0cpreds <- predict(m3.0c, m3.0c_predictions, type = "response", se.fit = TRUE)
m3.0c_sePreds <- data.frame(m3.0c_predictions,
                           mu   = exp(m3.0cpreds$fit),
                           low  = exp(m3.0cpreds$fit - 1.96 * m3.0cpreds$se.fit),
                           high = exp(m3.0cpreds$fit + 1.96 * m3.0cpreds$se.fit)) %>% 
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
