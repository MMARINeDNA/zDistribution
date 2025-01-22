
# Requires running 01_Read_format_data.R
# to generate object detect_by_station

ggplot(detect_by_station, aes(y = BestTaxon, x = depth, 
                              fill = BestTaxon, color = BestTaxon)) +
  geom_count() +
  theme_minimal() + 
  coord_flip() +
  scale_x_reverse() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = c(pnw_palette("Cascades",12, type = "continuous"),
                               pnw_palette("Sunset",12, type = "continuous")[1:11])) +
  scale_color_manual(values = c(pnw_palette("Cascades",12, type = "continuous"),
                                pnw_palette("Sunset",12, type = "continuous")[1:11])) +
  theme(legend.position = "none")