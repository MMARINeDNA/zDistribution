# Archive: species-specific depth splines in nimble

Unfinished work on fitting species-specific depth splines in nimble. None of these scripts is used in the manuscript, and some may not run as they are. They're kept here to return to later.

| File | What it is |
|---|---|
| `05_Model_w_Nimble_m1.2.R` | GAM `m1.2` (depth smooth by species) fitted in nimble, from its jagam export |
| `05_Model_w_Nimble_m1.2-slurmin.R` | Cluster version of the above (paths relative to the script's folder) |
| `05_Model_w_Nimble_m1.2+Occ.R` | Currently identical to `05_Model_w_Nimble_m1.2.R` (no occupancy layer yet, despite the name) |
| `05_Model_w_Nimble_m1.2+2LevelOcc-slurmin.R` | `m1.2` depth smooth inside the 2-level occupancy model, cluster version |
| `05_Model_w_Nimble_m1.0+2LevelOccModel_Bbai.R`, `_Lobl.R`, `_Mnov.R` | The 2-level occupancy model fitted separately to one species at a time (northern giant bottlenose whale, Pacific white-sided dolphin, humpback whale) |
| `nimble.def` | Apptainer/Singularity definition for a container with R, nimble and the tidyverse, for running these on a cluster |

Inputs, all produced by scripts in `Scripts/`:

- `ProcessedData/jagam_m1.2.RData`: the jagam export of `m1.2`, from `Scripts/02.5_Model_with_JAGS_Q1.R`
- `ProcessedData/jagam_m1.0.RData`: the depth-spline basis, from the same script
- `ProcessedData/detect_data.Rdata` and the per-species `detect_data_bbai/lobl/mnov.RData`, from `Scripts/01_Wrangle_data.R`

Most scripts read paths relative to the repo root (`./ProcessedData/...`); the `-slurmin` scripts read them relative to the script's own folder (`../ProcessedData/...`).
