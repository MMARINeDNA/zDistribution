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
                            "Homo", "Tursiops truncatus",
                            "Phoca vitulina", "Callorhinus ursinus",
                            "Zalophus californianus", "Eumetopias jubatus",
                            "Mirounga angustirostris",
                            "Tursiops truncatus", "Balaenoptera"))) %>% 
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
  mutate(SampleUID = paste0(Sample_name, "_", Plate)) %>% 
  filter(!(Plate == "314" & primer %in% c("DL", "DLL1"))) %>%
  filter(primer != "DLL1") %>% # take out all DLL1s %>%
  mutate(Plate = as.numeric(Plate))

# All samples MFU + MV1 positive samples, even if no cetaceans detected
data_samples <- readr::read_csv("./Data/MURIsampleInfo_2025-03-05.csv") 

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
  mutate(SampleUID = paste0(Sample_name, "_", Plate))

# Number of observations per dilution
check <- samples_info %>% 
  group_by(NWFSCsampleID, primer, techRep) %>% 
  filter(n() > 1) %>% 
  left_join(detect_data, by = c("NWFSCsampleID", "primer", "techRep", "dilution")) %>% 
  filter(!is.na(BestTaxon)) %>% 
  group_by(dilution) %>% 
  summarize(n.obs = n())

# Which samples have more than one dilution
check.double.dil <- samples_info %>% 
  group_by(NWFSCsampleID, primer, techRep) %>% 
  filter(n() > 1) %>% 
  slice_head()

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

# For samples with multiple dilutions, remove the dilution with the lower on-target read count
samples_info_onedil <- samples_info %>% 
  left_join(detect_data, by = c("NWFSCsampleID", "primer", "techRep", "dilution", "Plate", "Sample_name", "SampleUID")) %>%
  group_by(NWFSCsampleID, primer, dilution, techRep, Plate) %>% 
  mutate(onTarget_reads = sum(nReads)) %>% 
  replace_na(list(onTarget_reads = 0)) %>% 
  group_by(NWFSCsampleID, primer, techRep) %>% 
  slice_max(onTarget_reads, na_rm = TRUE) %>% 
  slice_head() %>%
  select(-c(BestTaxon, Class, nReads, onTarget_reads)) 

#check again
check.final <- samples_info_onedil %>% 
  group_by(NWFSCsampleID, primer, techRep) %>%
  filter(n() > 1) #should be an empty dataframe

# okay, now want a record for each possible species x sample

samples_info_species <- expand_grid(samples_info, BestTaxon = unique(detect_data$BestTaxon))

# now want to fill in whether each record resulted in a detection or not
# binary data for all possible sample x species combos



binary_detect_species <- left_join(samples_info_species, detect_data,
                                   by = c("Plate", "NWFSCsampleID", "primer", 
                                           "dilution", "techRep", 
                                          "BestTaxon")) %>%
  select(-Class) %>%
  replace_na(list(nReads = 0)) %>%
  mutate(Detected = ifelse(nReads>0, 1, 0)) %>%
#  mutate(techRep1 = techRep, .after = techRep) %>% 
#  mutate(techRep = ifelse(Plate == 313 | Sample_name == "MV1-52193-460-5-d10-1_S39", dilution, techRep)) %>% 
#  mutate(dilution = ifelse(Plate == 313 | Sample_name == "MV1-52193-460-5-d10-1_S39", techRep1, dilution)) %>% 
#  filter(!grepl("control", Sample_name)) %>% 
#  select(-techRep1) %>% 
  mutate(techRep = as.numeric(techRep))

# now want to collapse across techReps so that nTechReps is max(techRep)
# also create a variable for dilution proportion

diluteProp <- function(d){
  dn <- as.numeric(str_extract_all(d, "\\d+"))
  prop <- 1/dn
  return(prop)
}

detect_species_techreps <- binary_detect_species %>%
  group_by(BestTaxon, Plate, primer, NWFSCsampleID, dilution) %>%
  summarize(nTechReps = max(techRep), nReads = max(nReads), Detected = max(Detected)) %>%
  mutate(DilutionP = diluteProp(dilution))

### TODO:
### Thinking about this a bit more, I'm worried about combining a logistic relationship
### that we have a pretty good sense of (between tech reps and detection prob) with 
### a nonlinear relationship that we have very poor understanding of (with dilution).
### Would modeling each separately first tell us something? -AVC

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
save(detect_data, detect_species_meta, detect_species_divetime, samples_info_onedil,
     file = "./ProcessedData/detect_species_meta.RData")

write.csv(samples_info_onedil, file = paste("MURIsampleInfo_", as.character(Sys.Date()), ".csv"))       
