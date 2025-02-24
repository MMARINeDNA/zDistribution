#### 3D Distribution
#### wrangling data for presence/absence modelling
#### January 2025
#### EKJ&AVC

library(tidyverse)

# Load metadata
metadata <- read.csv("./Data/Hake_2019_metadata.csv")
timeAtDepth <- read.csv("./Data/MM_dive_time_expand.csv")
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")

# get all detections
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
                            "Canis", "Felis", "Petrogale xanthopus", 
                            "Homo", "Tursiops truncatus"))) %>% 
  mutate(BestTaxon = case_when(BestTaxon == "Lagenorhynchus" ~ "Lagenorhynchus obliquidens",
                               TRUE ~ BestTaxon)) %>% 
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

# All samples, even if no cetaceans detected
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
  filter(dilution != "pos", dilution != "poscontrol") %>%
  filter(techRep != "control") %>%
  unite(NWFSCsampleID, NWFSCpopID:NWFSCsampNum, sep = "-") %>% 
  separate(techRep, into = c("techRep",NA), sep = "_") %>% 
  separate(file_name, into = c(NA,NA,"data",NA,NA), sep = "/") %>% 
  separate(data, into = c(NA,"Plate",NA,NA), sep = "_")  %>%
  mutate(SampleUID = paste0(Sample_name, "_", Plate))

# Number of observations per dilution
check <- samples_info %>% 
  group_by(NWFSCsampleID, primer, techRep) %>% 
  filter(n() > 1) %>% 
  left_join(detect_data, by = c("NWFSCsampleID", "primer", "techRep", "dilution")) %>% 
  filter(!is.na(BestTaxon)) %>% 
  group_by(dilution) %>% 
  summarize(n.obs = n())

# Number of detects per species/primer/sample for multi-dilution samples
check2 <- samples_info %>% 
  group_by(NWFSCsampleID, primer, techRep) %>% 
  filter(n() > 1) %>% 
  left_join(detect_data, by = c("NWFSCsampleID", "primer", "techRep", "dilution")) %>% 
  group_by(NWFSCsampleID, primer, BestTaxon, techRep) %>% 
  filter(!is.na(BestTaxon)) %>% 
  summarize(n.obs.sample = n())

# Number of detects per species/primer/sample removing 1 multi-dilution replicate
check3 <- samples_info %>% 
  group_by(NWFSCsampleID, primer, techRep) %>% 
  filter(n() > 1) %>% 
  left_join(detect_data, by = c("NWFSCsampleID", "primer", "techRep", "dilution")) %>% 
  group_by(NWFSCsampleID, primer, dilution, techRep) %>% 
  mutate(onTarget_reads = sum(nReads)) %>% 
  group_by(NWFSCsampleID, primer, techRep) %>% 
  slice_max(onTarget_reads, na_rm = TRUE) %>% 
  group_by(NWFSCsampleID, primer, BestTaxon, techRep) %>% 
  summarize(n.obs.sample = n())

# okay, now want a record for each possible species x sample

samples_info_species <- expand_grid(SampleUID = samples_info$SampleUID, 
                                    BestTaxon = unique(detect_data$BestTaxon)) %>%
  left_join(samples_info, by = "SampleUID")

# now want to fill in whether each record resulted in a detection or not
# binary data for all possible sample x species combos

diluteProp <- function(d){
  dn <- as.numeric(str_extract_all(d, "\\d+"))
  prop <- 1/dn
  return(prop)
}

binary_detect_species <- left_join(samples_info_species, detect_data,
                                   by = c("Plate", "Sample_name", "primer", 
                                          "NWFSCsampleID", "dilution", "techRep", 
                                          "BestTaxon", "SampleUID")) %>%
  select(-Class) %>%
  replace_na(list(nReads = 0)) %>%
  mutate(Detected = ifelse(nReads>0, 1, 0)) %>%
#  mutate(techRep1 = techRep, .after = techRep) %>% 
#  mutate(techRep = ifelse(Plate == 313 | Sample_name == "MV1-52193-460-5-d10-1_S39", dilution, techRep)) %>% 
#  mutate(dilution = ifelse(Plate == 313 | Sample_name == "MV1-52193-460-5-d10-1_S39", techRep1, dilution)) %>% 
  filter(!grepl("control", Sample_name)) %>% 
#  select(-techRep1) %>% 
  mutate(techRep = as.numeric(techRep))

# now want to collapse across techReps so that nTechReps is max(techRep)
# also create a variable for dilution proportion

detect_species_techreps <- binary_detect_species %>%
  group_by(BestTaxon, Plate, primer, NWFSCsampleID, dilution) %>%
  summarize(nTechReps = max(techRep), nReads = max(nReads), Detected = max(Detected)) %>%
  mutate(DilutionP = diluteProp(dilution))

# okay now we're ready to attach metadata
detect_species_meta <- left_join(detect_species_techreps, metadata, 
                                 by = c("NWFSCsampleID" = "sampleID")) %>% 
  left_join(mmEcoEvo, by = c("BestTaxon" = "Species"))

# max dive depth by species, per the Navy reports
maxDepth_species <- timeAtDepth %>% 
  group_by(Species) %>% 
  summarize(maxDepth = max(depth))

detect_species_divetime <- detect_species_meta %>% 
  filter(!(BestTaxon %in% c("Balaenoptera", "Delphinidae"))) %>% 
  mutate(common_name = case_when(common_name == "killer whale"~"mammal eating killer whale",
                               TRUE~common_name)) %>% 
  left_join(timeAtDepth, by = c("common_name" = "Species", "depth" = "depth")) %>% 
  left_join(maxDepth_species, by = c("common_name" = "Species")) %>% 
  mutate(time_per_m = case_when(depth > maxDepth~0,
                                TRUE~time_per_m))
  
# save the file and we're done!
save(detect_data, detect_species_meta, detect_species_divetime,
     file = "./ProcessedData/detect_species_meta.RData")

         