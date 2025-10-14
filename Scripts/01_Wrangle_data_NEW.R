#### 3D Distribution
#### wrangling data for presence/absence modelling
#### Summer 2025
#### EKJ&AVC

library(tidyverse)

## Get data --------------------------------------------------------------------

metadata <- read.csv("./Data/Hake_2019_metadata.csv")
timeAtDepth <- read.csv("./Data/MM_dive_time_expand.csv")
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")
freezethaw <- read.csv("./Data/HAKE2019_miseq_runs_thaw.csv")

detect_data_raw <- read.csv("./Data/M3_compiled_taxon_table_wide.csv") %>% 
  pivot_longer(-c(BestTaxon, Class), names_to = "SampleUID", values_to = "nReads") %>% 
  group_by(SampleUID) %>% 
  mutate(totalReads = sum(nReads)) %>% 
  separate(SampleUID, into = c("Sample_name", NA), remove = FALSE, sep = "_") %>% 
  separate(Sample_name, into = c("plate", "primer", "pop", "sample", "dilution", "techRep", "seqRep"), remove = FALSE, sep = "\\.") %>% 
  unite(pop:sample, col = "NWFSCsampleID", sep = "-") %>% 
  mutate(techRep = as.numeric(techRep)) %>% 
  mutate(Detected = ifelse(nReads>0, 1, 0)) %>% 
  filter(!(primer %in% c("MFU", "MV1") & totalReads == 0)) %>% 
  filter(Class == "Mammalia") %>% 
  filter(!BestTaxon %in% c("Moschus", "Equus caballus"))

## Add freeze/thaw info INCOMPLETE

freezethaw_mod <- freezethaw %>%
  select(-Plate) %>% # take this out to avoid confusion
  mutate("plate" = paste0("MURI", RunNo)) %>%
  rename("primer" = Markers) %>%
  mutate(plate = as.character(plate))

detect_data_thaw <- detect_data_raw %>%
  left_join(freezethaw_mod, by = c("plate", "primer"))

## Filter out DLL1, C16 primer, plate 309, and DL/DLL1 from plate 314 ----------

detect_data_filt <- detect_data_raw %>% 
  filter(plate != "MURI309") %>% 
  filter(!(primer %in% c("DLL1N", "C16"))) %>% 
  filter(!(primer == "DL" & plate == "MURI314"))



## Reduce sequencing reps ------------------------------------------------------

detect_data_1seq <- detect_data_filt %>% 
  group_by(plate, primer, NWFSCsampleID, dilution, techRep, seqRep) %>% 
  mutate(totReads = sum(nReads)) %>% 
  ungroup() %>% 
  group_by(primer, NWFSCsampleID, dilution, techRep) %>% 
  filter(totReads == max(totReads)) %>% 
  ungroup() %>% 
  mutate(seqRep = replace_na(seqRep, "sr1")) %>% 
  filter(!(totReads == 0 & seqRep %in% c("sr2", "sr3"))) %>% 
  select(-totReads)
  

## Reduce dilutions ------------------------------------------------------------

detect_data_1dil <- detect_data_1seq %>% 
  group_by(plate, primer, NWFSCsampleID, dilution, techRep) %>% 
  mutate(totReads = sum(nReads)) %>% 
  ungroup() %>% 
  group_by(primer, NWFSCsampleID, techRep) %>% 
  filter(totReads == max(totReads)) %>% 
  mutate(dilution = substr(dilution, 2, nchar(dilution))) %>% 
  mutate(diluti0n = as.numeric(dilution)) %>% 
  arrange(dilution, .by_group = TRUE) %>% 
  slice_head(n = 22) %>% 
  select(-totReads) 
          

## Reduce tech reps ------------------------------------------------------------
## EKJ Note I think we don't want to do this anymore, since
## each tech rep will have a different number of freeze-thaw cycles
## associated with it.
# detect_data_1rep <- detect_data_1dil %>% 
#   group_by(primer, NWFSCsampleID, BestTaxon) %>% 
#   mutate(nReps = n()) %>% 
#   mutate(Detected = max(Detected)) %>% 
#   slice_head()
  
## Count number of samples -----------------------------------------------------

nSamps_primer <- detect_data_1rep %>% 
  group_by(primer, NWFSCsampleID) %>% 
  n_groups()

length(unique(detect_data_1rep$NWFSCsampleID))

## Add metadata ----------------------------------------------------------------

detect_data_meta <- detect_data_1rep %>% 
  left_join(metadata, by = c("NWFSCsampleID" = "sampleID")) %>%
  left_join(mmEcoEvo, by = c("BestTaxon" = "Species"))

## Check species are all marine mammals
unique(detect_data_meta$BestTaxon)

## Remove Delphinidae family and all pinniped species, add common names
detect_data <- detect_data_meta %>% 
  filter(!(BestTaxon %in% c('Delphinidae', "Callorhinus ursinus",
                            "Eumetopias jubatus", "Phoca vitulina",
                            "Zalophus californianus", "Mirounga angustirostris")))  

  
## Check species are all cetaceans
unique(detect_data$BestTaxon)
unique(detect_data$common_name)


## Remove delphinid and baleen detections <100m from bottom (likely whalefall) -

detect_data_nowf <- detect_data %>%
  mutate(dist_to_bottom = bottom.depth.consensus - depth) %>%
  mutate(Detected = case_when(dist_to_bottom < 100 &
                  Detected == 1 &
                  bottom.depth.consensus > 200 &
                  Broad_taxa %in% c("Baleen whale", "Dolphin/Porpoise")~0,
                  TRUE~Detected))

## count number of detections by species ---------------------------------------

detect_per_species <- detect_data %>% 
  group_by(BestTaxon) %>% 
  summarize(nDetect = sum(Detected))

detect_per_species_nowf <- detect_data_nowf %>% 
  group_by(BestTaxon) %>% 
  summarize(nDetect = sum(Detected))

detect_per_family <- detect_data %>% 
  group_by(Family) %>% 
  summarize(nDetect = sum(Detected))

## add time at depth per species -----------------------------------------------

# max dive depth by species, per the Navy reports
maxDepth_species <- read.csv("./Data/MM_dive_time_expand.csv") %>% 
  group_by(Species) %>% 
  summarize(maxDepth = max(depth))

detect_species_divetime <- detect_data %>% 
  mutate(common_name = case_when(common_name == "killer whale"~"mammal eating killer whale",
                                 TRUE~common_name)) %>% 
  mutate(depth = case_when(depth == 0~1,
                           TRUE~depth)) %>% 
  left_join(timeAtDepth, by = c("common_name" = "Species", "depth" = "depth")) %>% 
  left_join(maxDepth_species, by = c("common_name" = "Species")) %>% 
  mutate(time_per_m = case_when(depth > maxDepth~0,
                                TRUE~time_per_m)) %>% 
  mutate(time_10m = case_when(depth > maxDepth~0,
                              TRUE~time_10m)) %>% 
  mutate(time_equalBin = case_when(depth > maxDepth~0,
                                   TRUE~time_equalBin))

save(detect_data, detect_species_divetime,
     detect_per_species, detect_per_family, 
     maxDepth_species, file = "./ProcessedData/detect_data.Rdata")
