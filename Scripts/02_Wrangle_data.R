#### 3D Distribution
#### wrangling data for presence/absence modelling
#### January 2025
#### EKJ

library(tidyverse)

# First, load data and keep records from all samples, even if no cetaceans detected

metadata <- read.csv("./Data/Hake_2019_metadata.csv")

detect_list <- list.files(path = "./Data", pattern = "taxon_table.csv", 
                          recursive = TRUE, full.names = TRUE)

data_samples <- readr::read_csv(detect_list, id = "file_name") %>%
  select(-BestTaxon, -Class, -nReads) %>%
  unique()

# these are all of the samples that were run, with SampleUID as the unique identifier
samples_info <- data_samples %>% 
  separate(Sample_name, 
           into = c("primer","NWFSCpopID","NWFSCsampNum","dilution","techRep"),
           sep = "-",
           remove = FALSE) %>% 
  filter(NWFSCpopID == 52193) %>% 
  unite(NWFSCsampleID, NWFSCpopID:NWFSCsampNum, sep = "-") %>% 
  separate(techRep, into = c("techRep",NA), sep = "_") %>% 
  separate(file_name, into = c(NA,NA,"data",NA,NA), sep = "/") %>% 
  separate(data, into = c(NA,"Plate",NA,NA), sep = "_")  %>%
  mutate(SampleUID = paste0(Sample_name, "_", Plate))

# okay, now want a record for each possible species x sample

samples_info_species <- expand_grid(SampleUID = data_samples_info$SampleUID, 
                                    BestTaxon = unique(detect_data$BestTaxon)) %>%
  left_join(samples_info, by = "SampleUID")

# now want to fill in whether each record resulted in a detection or not

# first, modify AVC code to not include metadata yet, just detections
detect_data_nometa <- detect_data %>% 
  separate(Sample_name, 
           into = c("primer","NWFSCpopID","NWFSCsampNum","dilution","techRep"),
           sep = "-",
           remove = FALSE) %>% 
  filter(NWFSCpopID == 52193) %>% 
  unite(NWFSCsampleID, NWFSCpopID:NWFSCsampNum, sep = "-") %>% 
  separate(techRep, into = c("techRep",NA), sep = "_") %>% 
  separate(file_name, into = c(NA,NA,"data",NA,NA), sep = "/") %>% 
  separate(data, into = c(NA,"Plate",NA,NA), sep = "_") %>%
  mutate(SampleUID = paste0(Sample_name, "_", Plate))

# binary data for all possible sample x species combos

binary_detect_species <- left_join(samples_info_species, detect_data_nometa,
                                   by = c("Plate", "Sample_name", "primer", 
                                          "NWFSCsampleID", "dilution", "techRep", 
                                          "BestTaxon", "SampleUID")) %>%
  select(-Class) %>%
  replace_na(list(nReads = 0)) %>%
  mutate(Detected = ifelse(nReads>0, 1, 0)) %>%
  filter(techRep != "control" & techRep != "d1") %>%
  mutate(techRep = as.numeric(techRep))

# now want to collapse across techReps so that nTechReps is max(techRep)

detect_species_techreps <- binary_detect_species %>%
  group_by(BestTaxon, Plate, primer, NWFSCsampleID, dilution) %>%
  summarize(nTechReps = max(techRep), nReads = max(nReads), Detected = max(Detected))

# okay now we're ready to attach metadata

detect_species_meta <- left_join(detect_species_techreps, metadata, by = c("NWFSCsampleID" = "sampleID"))

# save the file and we're done!

save(detect_species_meta, file = "./ProcessedData/detect_species_meta.RData")

         
         