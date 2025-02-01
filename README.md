# 3DDistribution
Modelling 3D distributions of marine mammal detections from eDNA samples, with the goal of answering the following questions:

1. Does the probability of a detecting cetaceans in eDNA samples vary with sample depth?
   - H0: Probability of detection does not vary with depth.
   - H1: Probability of detection varies across depth agnostic to species or functional group.
   - H2: Probability of detection varies across depth according to species or functional group.
2. Does the probability of a detecting cetaceans in eDNA samples change with the number of technical replicates?
   * NOTE should we wrap dilution into this question as well?
   - H0: Detection does not vary with # technical replicates.
   - H1: Detection varies with # technical replicated agnostic to species/fuctional group or depth.
   - H2: Detection varies with # technical replicates according to species/functional group, depth, or a combination of the two.
3. Does depth distribution of detections vary across xy spatial distribution?
   - H0: Depth distribution of detections does not vary across xy spatial distribution.
   - H1: Depth distribution of detection does vary across xy spatial distribution agnostic to oceanography (e.g. upwelling).
   - H2: Depth distribution of detection does vary across xy spatial distribution according to oceanography (e.g. upwelling).

# The Plan

Repo will contain code, data must be downloaded from [Google Drive](https://drive.google.com/drive/folders/1EZEfbxgRszwmN4RmaoQe7wh6S6zPgm5A?usp=drive_link) and placed in a folder named "Data".

Marine mammal metadata (taxonomy, species-specific time-at-depth and prey preference): contact AVC. Dive data come from [here](https://apps.dtic.mil/sti/tr/pdf/ADA560975.pdf) and [here](https://www.nepa.navy.mil/Portals/20/Documents/aftteis4/Dive%20Profile%20and%20Group%20Size_TR_2017_05_22.pdf).

EKJ TODO: clean up Analysis 3

# Analysis 1

- collapse all data across X and Y, ignore xy distribution
- model detection probability by depth ignoring species
- compare with models incorporating species, family, and prey category

# Analysis 2

- Using best model(s) from Analysis 1, retest varying number of replicates
- Potential alternative: run [Brice's replication model](https://github.com/BriceSemmens/eDNA_patch) without assuming species' presence and adding depth as a covariate

# Analysis 3

- develop 2D models of distribution for each species
- given that distribution, does sampling depth matter?
- does it vary by species or species type?
- develop explicit 3D models of distribution for each species
- assume that detections at depth reflect species z distribution
- incorporate spread and decay from Mod 1?
- does eDNA depth distribution interact with xy distribution due to differences in, e.g. oceanographic upwelling or downwelling?
- incorproate prey (Note for future consideration: this is one of the broad Mod 3 goals but not sure if it fits in this paper)
