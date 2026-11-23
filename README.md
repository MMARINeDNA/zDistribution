# Z-variability in POD of rare taxa

This repo contains the code and analyses used to examine depth variability in the probability of detecting cetaceans in eDNA samples. A living version of the manuscript can be found [here](https://mmarinedna.github.io/zDistribution/). We analyse marine mammal detections from water samples collected throughout the California Current at depths of 0–500 m to answer three questions:

1. How does observation error affect the probability of detection (POD) of cetacean DNA?
2. Does species-specific detection rate vary with depth?
3. Does depth variability in detection rate vary spatially?

Raw data are in `Data/` and intermediate data products in `ProcessedData/`; some data files are kept out of git. Dive data come from [here](https://apps.dtic.mil/sti/tr/pdf/ADA560975.pdf) and [here](https://www.nepa.navy.mil/Portals/20/Documents/aftteis4/Dive%20Profile%20and%20Group%20Size_TR_2017_05_22.pdf). Q3 needs the Becker et al. (2020) density predictions in `Data/CCE_model_run/`.

The manuscript is `index.qmd`, which includes `zDist_methods.qmd` and `zDist_results.qmd`. Render it from the repo root with:

```bash
LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 quarto render index.qmd
```

(The UTF-8 locale keeps the Greek letters in Table 2 from rendering as `<U+03B2>`.)

## Old question numbering

The order of the questions changed during the project, and some script, model and object names still use the old numbering:

- The manuscript's **Q1** (observation error, occupancy model) was the old Q2. Its scripts start with `05_`, and its figure objects start with `Q2_`.
- The manuscript's **Q2** (detection rate by depth, GAMs) was the old Q1. Its scripts start with `02_`, its models are `m1.x`, and its objects start with `H1`.
- **Q3** kept its number: scripts start with `04_`, models are `m3.x`, and objects start with `H3` or `m3`.

## Scripts by manuscript question

Run the scripts in the order listed within each block.

| Manuscript | Script | Produces |
|---|---|---|
| **Data** | `0.5_MM_dive_data_expand.R` | `Data/MM_dive_time_expand.csv` (time at depth, used in Fig. 3) |
| | `01_Wrangle_data.R` | `ProcessedData/detect_data*.RData` (514 samples) |
| | `01.1_Plot_detections.R` | Fig. 1 (left) |
| | `01.2_Sampling_map.R` | Fig. 1 (right) |
| | `01.2_Examine_deep_detections.R` | Deep and near-seafloor detections; Fig. S4. It saves to `detect_species_meta.RData`, but the manuscript loads `ProcessedData/potential_whalefall.RData`. |
| **Q1** | `02.5_Model_with_JAGS_Q1.R` | Depth-spline basis for the occupancy model (`ProcessedData/jagam_m1.0.RData`) |
| | `05_Model_w_Nimble_m1.0+2LevelOccModel_wThaw.R` | Occupancy model; Fig. 2 (left, right) |
| | `05_Calc_POD_by_sample.R` | Per-sample POD |
| | `05_Plot_TechRepsXPrimer.R` | Fig. 2 (middle) |
| **Q2** | `02_Model_with_GAMs_Q1.R` | Depth GAMs; Table 2 |
| | `02.1_Plot_H1_models.R` | Figs. 3, 4, S1, S3 |
| | `02.2_m1.2_performance.R` | Cross-validated AUC and TSS; Fig. S2 (top) |
| **Q3** | `04_Model_with_GAMs_Q3.R` | xyz and density GAMs |
| | `04.1_Plot_H3_models.R` | Figs. 5, 6 |
| | `04.2_m3.0c_performance.R` | Cross-validated AUC and TSS; Fig. S2 (bottom) |

Other scripts in `Scripts/` are exploratory and aren't used in the manuscript.

## Models and objects

| Manuscript | Model in code | Equation | Saved in | Shown in |
|---|---|---|---|---|
| **Q1** occupancy model | 2-level occupancy model with the `m1.0` depth spline (`05_Model_w_Nimble_m1.0+2LevelOccModel_wThaw.R`) | Eqs. 1–3 | `Results/nimbleOut_m1.0_2LevelOcc_wThaw.RData`, `Results/post_pdetect_m1.0+2LevelOcc_wThaw.RData` | Fig. 2 (`Figures/Q2_detect_probability_plots.Rdata`, `Figures/Q2_TechRepsXPrimer.RData`) |
| **Q2** depth smooth | `m1.0` | Eq. 4, line 1: β₀ + f(z) | `ProcessedData/H1models.RData` | Table 2; Fig. S1 |
| **Q2** depth smooth, species intercept | `m1.1` | Eq. 4, line 2: β_s + f(z) | `ProcessedData/H1models.RData` | Table 2 |
| **Q2** depth smooth by species (best model) | `m1.2` | Eq. 4, line 3: β_s + f_s(z) | `ProcessedData/H1models.RData` | Table 2; Figs. 3, S2 (top), S3 |
| **Q2** depth smooth by family | `m1.2a` | Eq. 4, line 4: β_F + f_F(z) | `ProcessedData/H1models.RData` | Table 2; Fig. 4 (left) |
| **Q2** depth smooth by prey category | `m1.2b` | Eq. 4, line 5: β_e + f_e(z) | `ProcessedData/H1models.RData` | Table 2; Fig. 4 (right) |
| **Q3** approach 1 (xyz) | `m3.0c` | Eq. 5 | `ProcessedData/m3.0models_preds_0.05degree.Rdata` | Fig. 5; Fig. S2 (bottom) |
| **Q3** approach 2 (density) | `m3.1` | Eq. 6 | `ProcessedData/m3.0models_preds_0.05degree.Rdata` | Fig. 6 |

Notes:

- In the occupancy model, `pdetect_intercept` and `thaw_par` are the primer intercepts (γ_m) and freeze/thaw slopes (γ_t,m) in Eq. 3, and `b_depth` holds the occupancy intercept and depth-spline coefficients in Eq. 1. Primers are indexed 1 = DL, 2 = MFU, 3 = MV1.
- Script 04 also fits `m3.0a` (depth + xy) and `m3.0b` (a tensor product of xyz); neither is reported in the manuscript.
- Equation numbers are as rendered; the labels in the `.qmd` files differ (`eq-4-8` renders as Eq. 4, `eq-9` as Eq. 5 and `eq-10` as Eq. 6).
- Cross-validation results are in `ProcessedData/m1.2clean_performance.Rdata` (Q2) and `ProcessedData/m3.0cperformance.Rdata` (Q3); the figure objects for Q2 and Q3 are in `Figures/H1PODplots.Rdata` and `Figures/H3.0c_map_clean.Rdata`.

## Archive

`Archive/` holds unfinished work on species-specific depth splines in nimble. None of it is used in the manuscript, and some of it may not run as is. See `Archive/README.md`.
