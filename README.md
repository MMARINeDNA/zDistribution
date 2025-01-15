# 3DDistribution
Modelling 3D distributions of marine mammal detections from eDNA samples, with the goal of nswering the following questions:

1. Does detection probability vary across depth?
   - H0: Detection does not vary with depth.
   - H1: Detection varies across depth agnostic to species or functional group.
   - H2: Detection varies across depth according to species or functional group.
2. Does detection probability vary with technical replicates?
   - H0: Detection does not vary with # technical replicates.
   - H1: Detection varies with # technical replicated agnostic to species/fuctional group or depth.
   - H2: Detection varies with # technical replicates according to species/functional group, depth, or a combination of the two.
3. Does depth distribution of detections vary across xy spatial distribution?
   - H0: Depth distribution of detections does not vary across xy spatial distribution.
   - H1: Depth distribution of detection does vary across xy spatial distribution agnostic to oceanography (e.g. upwelling).
   - H2: Depth distribution of detection does vary across xy spatial distribution according to oceanography (e.g. upwelling).

# The Plan

Repo will contain code, data must be downloaded from [Google Drive](https://drive.google.com/drive/folders/1EZEfbxgRszwmN4RmaoQe7wh6S6zPgm5A?usp=drive_link) and placed in a folder named "Data".

EKJ TODO: clean this up and add connections with Brice's model https://github.com/BriceSemmens/eDNA_patch

# Analysis 1

- collapse all data across X and Y, ignore xy distribution
- model detection probability
- add species/species functional type
- add # of technical replicates (volume)
- add xy distrbution explicitly?

# Analysis 2

- develop 2D models of distribution for each species
- given that distribution, does sampling depth matter?
- does it vary by species or species type?
- by # of technical replicates (aka volume)?
- does eDNA depth distribution interact with xy distribution due to differences in, e.g. oceanographic upwelling or downwelling?

# Analysis 3

- develop explicit 3D models of distribution for each species
- assume that detections at depth reflect species z distribution
- incorporate spread and decay from Mod 1?
- incorproate prey (Note for future consideration: this is one of the broad Mod 3 goals but not sure if it fits in this paper)
