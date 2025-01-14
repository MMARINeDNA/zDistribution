#### 3D distribution
#### Aggregate and format data
#### Jan 2025

#### Set up environment --------------------------------------------------------

library(tidyverse)
library(ggridges)
library(PNWColors)

#### Get data ------------------------------------------------------------------

detect_list <- list.files(path = "./Data", pattern = "taxon_table.csv", 
                          recursive = TRUE, full.names = TRUE)

detect_data <- readr::read_csv(detect_list, id = "file_name") %>% 
  filter(Class == "Mammalia") %>% 
  filter(!(BestTaxon %in% c("Homo sapiens", "Bos taurus", 
                            "Sus scrofa", "Bos",
                            "Ovis", "Felis catus",
                            "Artiodactyla","Canis lupus",
                            "Macropodidae", "Diprotodontia",
                            "Macropus fuliginosus",
                            "Macropus giganteus",
                            "Osphranter","Osphranter robustus",
                            "Osphranter rufus", "Dama dama",
                            "Canis"))) %>% 
  mutate(BestTaxon = case_when(BestTaxon == "Lagenorhynchus" ~ "Lagenorhynchus obliquidens",
                               TRUE ~ BestTaxon))

metadata <- read.csv("./Data/Hake_2019_metadata.csv")

rm(detect_list)

#### Merge and format detection data -------------------------------------------

detect_data_meta <- detect_data %>% 
  separate(Sample_name, 
           into = c("primer","NWFSCpopID","NWFSCsampNum","dilution","techRep"),
           sep = "-",
           remove = FALSE) %>% 
  unite(NWFSCsampleID, NWFSCpopID:NWFSCsampNum, sep = "-") %>% 
  separate(techRep, into = c("techRep",NA), sep = "_") %>% 
  separate(file_name, into = c(NA,NA,"data",NA,NA), sep = "/") %>% 
  separate(data, into = c(NA,"Plate",NA,NA), sep = "_") %>% #I'm not sure what LCA3 means in the filesnames
  left_join(metadata, by = c("NWFSCsampleID" = "sampleID"))

#### Make a fun little ridgeplot -----------------------------------------------

detect_by_station <- detect_data_meta %>% 
  group_by(station, depth, techRep) %>% 
  mutate(totReads = sum(nReads)) %>% 
  mutate(detect = case_when(totReads > 1 ~ 1,
                            TRUE ~ 0)) %>% 
  slice_head() %>% 
  ungroup() %>% 
  filter(detect == 1)

lowdetect_subset <- detect_by_station %>% 
  group_by(BestTaxon) %>% 
  mutate(nDetect = n()) %>% 
  filter(nDetect < 3)

ggplot(detect_by_station, aes(y = BestTaxon, x = depth, 
                              fill = BestTaxon, color = BestTaxon)) +
  geom_density_ridges(scale = 1.5, 
                      bandwidth = 50, 
                      jittered_points = TRUE,
                      point_alpha = 1,
                      point_shape = 21,
                      alpha = 0.6) +
  geom_point(data = lowdetect_subset, aes()) +
  theme_minimal() + 
  coord_flip() +
  scale_x_reverse() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = c(pnw_palette("Cascades",12, type = "continuous"),
                               pnw_palette("Sunset",12, type = "continuous")[1:11])) +
  scale_color_manual(values = c(pnw_palette("Cascades",12, type = "continuous"),
                               pnw_palette("Sunset",12, type = "continuous")[1:11])) +
  theme(legend.position = "none")

save(detect_data_meta, file = "./ProcessedData/detect_data_meta.RData")

                    