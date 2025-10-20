# explore effect of freeze-thaw on probability of detection

load("./ProcessedData/detect_data.Rdata")

summarize_thaw <- detect_data %>%
  filter(!is.na(Thaw)) %>%
  group_by(primer, Thaw) %>%
  summarize(POD = sum(Detected)/length(Detected))

ggplot(summarize_thaw) +
  geom_col(aes(x= Thaw, y = POD)) +
  facet_wrap(~primer, nrow = 3) +
  theme_bw()
