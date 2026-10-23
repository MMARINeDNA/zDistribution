# Critical review of methods, results and code

*Reviewer: Claude Code, 22 October 2026. Branch reviewed: `claude/manuscript-code-review-d6debd` (at `06d49d0`).*

## Scope and approach

I read the manuscript sources (`index.qmd`, `zDist_methods.qmd`, `zDist_results.qmd`) and every script that produces an object the manuscript loads. I then loaded the saved data and model objects in R to check each number in the text, rather than relying on code comments.

Two caveats affect what could be checked:

- Several inputs are not in the repository (see [§1.3](#13-the-pipeline-cannot-be-re-run-from-the-repository)). I could not re-run the pipeline end to end, so most checks are against the saved `.RData` outputs.
- The manuscript appears to have been rendered from files that exist only in the main checkout and are not committed: a modified, 250 MB `ProcessedData/m3.0models_preds_0.05degree.Rdata` and an untracked `ProcessedData/H3.0c_pred_flt_clean.Rdata`. I opened these read-only and noted where they differ from the committed versions.

Severity key: **High** = a reported result or conclusion is wrong, or can't be reproduced. **Medium** = the method departs from what the paper says, or a modelling choice materially affects inference. **Low** = wording, a stale comment, or fragile code.

---

## Summary

| # | Issue | Severity |
|---|---|---|
| 1.1 | The Q3 statistics, maps and cross-validation come from three different "m3.0c" models, none of which matches the methods | High |
| 1.2 | Q3 depth-of-maximum maps include depths below the seafloor; the area percentages in the text match no saved output | High |
| 1.3 | Pipeline can't be re-run from the repo (missing inputs, uncommitted outputs, scripts that error) | High |
| 2.1 | Occupancy-model priors are informative (nimble `dnorm(0, 10)` is precision, so SD = 0.32), not "weakly informative" | High |
| 2.2 | The depth spline's basis and priors don't match (`diagonalize` is ignored by `gam()`) | Medium |
| 2.3 | DL was never observed at 1 freeze/thaw cycle, so the headline "POD at Thaw = 1" comparison for DL is an extrapolation | Medium |
| 2.4 | Eq. 3 has one thaw slope; the code fits one per primer | Medium |
| 2.5 | All MCMC chains start from identical initial values | Low |
| 3.1 | Sequencing-rep and dilution selection maximises **mammal** reads, not total reads (selection on the outcome) | Medium–High |
| 3.2 | 129 failed (zero-read) MFU/MV1 replicates are counted as non-detections, contrary to the methods | Medium |
| 3.3 | Whale-fall filtering is computed but never used | Low–Medium |
| 4.1 | GAMs treat replicates, species rows and depths within a station as independent (pseudoreplication) | Medium–High |
| 4.2 | The text says the family and prey models had *lower* AIC than the species model; they are higher | High (text) |
| 4.3 | Cross-validation: station leakage, no seed, the wrong data or model in places, and a figure built from the wrong object | Medium |
| 4.4 | The time-at-depth overlay in Fig. 3 is min–max rescaled, not "on the same scale" | Low–Medium |
| 5 | Methods text vs. code mismatches for Q3 (coordinates, Eq. 10, density units/year, species subset) | Medium |
| 6 | Other numeric and text inconsistencies | Low |

---

## 1. Reproducibility and Q3 results (High)

### 1.1 Three different "m3.0c" models feed the Q3 results

The methods ([zDist_methods.qmd:99-118](../zDist_methods.qmd#L99-L118)) describe a model fitted to a **three-species subset**, with a thin-plate spline on **projected** coordinates. That is what [04_Model_with_GAMs_Q3.R:143-173](../Scripts/04_Model_with_GAMs_Q3.R#L143-L173) fits (UTM, `k = c(20,16)` / `c(5,16)`). However:

| Where it is used | Model actually used | Data | Coordinates | Species |
|---|---|---|---|---|
| Text: AIC = 1410, dev. expl. = 26.2% ([zDist_results.qmd:115](../zDist_results.qmd#L115)) | 3-species UTM model from script 04 | clean | UTM | 3 |
| Fig. 5 maps (`maxPOD_depth_clipped_clean`) | `m3.0c_clean` (only in the uncommitted local `.Rdata`) | clean, 41,280 rows | **lon/lat degrees** | **16** |
| Text: AUC = 0.77, TSS = 0.43 ([04.2_m3.0c_performance.R:41-71](../Scripts/04.2_m3.0c_performance.R#L41-L71)) | lon/lat model with `k = c(10,16)` on the space×taxon and depth×taxon terms | **non-clean** `detect_data` | lon/lat | 11 (those with >10 detections) |

The `m3.0c` object in both the committed and the local `m3.0models_preds_0.05degree.Rdata` is a 16-species lon/lat model on the **non-clean** data (48,656 rows; AIC 4985). Neither file contains `m3.1`, `detect_data_dens` or `D_grid`. Yet `index.qmd` calls `AIC(m3.1)`, and `Dmarg_plot` depends on `D_grid`. The values 1410 and 1571 in the committed render therefore came from an R session whose objects aren't saved anywhere.

The CV metrics should also use the clean data (`AUC_df_m3.0c_clean` = 0.78, `TSS_df_m3.0c_clean` = 0.45). Even then they evaluate a different model from the one described.

**Fix:** pick one Q3 model (the methods describe the 3-species, projected-coordinate version). Refit it, and generate the predictions, maps and CV from that same object and data. Save `m3.0c`, `m3.1`, `detect_data_dens` and `D_grid` together, from a committed script.

### 1.2 The depth-of-maximum maps include sub-seafloor depths; the area percentages don't match

- `m3.0c_clean_sePreds` (the source of Fig. 5) has **all 51 depths (0–500 m) at every grid cell**. No bathymetry mask was applied. On the shelf, the "depth of maximum detection" can be a depth below the seabed.
- As a result, in `maxPOD_depth_clean` the northern giant bottlenose whale's depth of maximum detection is **460 or 470 m in every cell** (5,996 and 8,563 cells). That contradicts the text: "greater than 150 m in 52% of the study area … deeper to the south and east" ([zDist_results.qmd:119](../zDist_results.qmd#L119)).
- The committed script's mask has a sign bug: `depth <= pmax(abs(seafloor), 0)` ([04_Model_with_GAMs_Q3.R:282-283](../Scripts/04_Model_with_GAMs_Q3.R#L282-L283); same at L391-392 and in `dave_futz.R`). In marmap, land has positive elevation, so `abs()` keeps land cells for depths up to the land's elevation. Use `depth <= pmax(-seafloor, 0)`, and treat `NA` deliberately rather than keeping it.
- The percentages in [zDist_results.qmd:117-121](../zDist_results.qmd#L117-L121) are hard-coded, and no script computes them. Recomputing them from the saved outputs gives:

| Statement in text | Text | `maxPOD_depth_clipped_clean` (local) | `maxPOD_depth` (committed) |
|---|---|---|---|
| Pacific white-sided dolphin: max depth < 100 m | 68% | 61.7% | 79.6% |
| Pacific white-sided dolphin: area masked (< 0.005) | 3% | 0.3% | 1.7% |
| Northern giant bottlenose whale: max depth > 150 m | 52% | 100% | 72.1% |
| Northern giant bottlenose whale: area masked | 31% | 0% | 18.0% |
| Humpback whale: max depth < 150 m | 40% | 39.3% | 37.8% |
| Humpback whale: max depth > 450 m | 33% | 57.6% | 54.5% |
| Humpback whale: area masked | 1.5% | 0% | 5.8% |

In the committed `H3.0c_pred_flt.Rdata`, humpback `max_mu` reaches **2.00**, which is impossible for a probability. That file must have been built with a different link transform (for example `exp()`). It is stale, but it shouldn't remain in the repo looking like a valid output.

**Fix:** mask predictions by bathymetry before taking the maximum over depth. Compute the area shares in code, as inline R in the `.qmd`, and weight by cell area if the grid is in degrees.

### 1.3 The pipeline cannot be re-run from the repository

**Missing inputs:**
- `duplicated_station_depth.csv` ([01_Wrangle_data.R:248](../Scripts/01_Wrangle_data.R#L248)). This defines which biological replicates are dropped, so `detect_data_clean` can't be regenerated without it.
- `Data/CCE_model_run/` (the Becker et al. density grids, [04_Model_with_GAMs_Q3.R:31](../Scripts/04_Model_with_GAMs_Q3.R#L31)).
- `Data/MFU_database.csv` and `Data/MV1_database.csv` ([01.2_Examine_deep_detections.R:12-13](../Scripts/01.2_Examine_deep_detections.R#L12-L13)).
- `Data/ddf/` and `Data/MASTER_CCE_1991_2018_GS.csv` (`extract_eiren.R`).

**Uncommitted outputs that the render depends on:**
- `ProcessedData/H3.0c_pred_flt_clean.Rdata` is loaded at [index.qmd:47](../index.qmd#L47) but is untracked.
- The local `m3.0models_preds_0.05degree.Rdata` is 250 MB, which exceeds GitHub's 100 MB file limit, so it can't be committed without LFS.
- No committed script creates `m3.0c_clean`, `m3.0c_clean_sePreds`, `maxPOD_depth_clean` or `H3.0c_map_clean.Rdata`. Only a partial trial exists, in `dave_futz.R` on rows 1:100.

**Scripts that error or write to the wrong place:**
- [02_Model_with_GAMs_Q1.R:301](../Scripts/02_Model_with_GAMs_Q1.R#L301) saves `m1.2clean`, which the script never creates. (In the saved file, `m1.2clean` is identical to `m1.2`.)
- [01.2_Examine_deep_detections.R:164-165](../Scripts/01.2_Examine_deep_detections.R#L164-L165) saves `potential_whalefall` and `whalefall_map` to `detect_species_meta.RData`, but the manuscript loads `potential_whalefall.RData`. The script also overwrites the object that `03_Model_with_GAMs_Q2.R` expects.
- [04.1_Plot_H3_models.R:114](../Scripts/04.1_Plot_H3_models.R#L114) has its save commented out.
- `index.qmd` loads `Figures/Q2_DetectabilityByPrimer.RData`, which is produced by now-commented code. It isn't used, but it's a stale dependency.

**Fix:** add the missing inputs (or document where to get them), commit a single `run_all` script or `Makefile` in execution order, and have the `.qmd` load only objects that a committed script creates.

---

## 2. Q1 occupancy model ([05_Model_w_Nimble_m1.0+2LevelOccModel_wThaw.R](../Scripts/05_Model_w_Nimble_m1.0+2LevelOccModel_wThaw.R))

### 2.1 The priors on POD are informative, not "uninformative or weakly informative" (High)

[L116-122](../Scripts/05_Model_w_Nimble_m1.0+2LevelOccModel_wThaw.R#L116-L122): `pdetect_intercept[k] ~ dnorm(0, 10)` and `thaw_par[k] ~ dnorm(0, 10)`. In nimble (BUGS parameterisation), the second argument is **precision**. I confirmed that `getParam(..., "sd")` returns 0.316.

On the logit scale, that prior puts about 95% of its mass on POD between 0.35 and 0.65 at Thaw = 0, and on thaw slopes between −0.62 and +0.62 per cycle. The posteriors show the prior is doing much of the work:

| Parameter | Posterior SD | Prior SD |
|---|---|---|
| `pdetect_intercept[DL]` | 0.265 | 0.316 |
| `pdetect_intercept[MFU]` | 0.202 | 0.316 |
| `thaw_par[MFU]` | 0.173 | 0.316 |

This contradicts [zDist_methods.qmd:71](../zDist_methods.qmd#L71). It likely shrinks every primer's POD towards 0.5 and the thaw effects towards 0, and so it affects every Q1 number: 0.44/0.35/0.31, the replicate counts needed for 90%, and the sample-level POD of 0.87.

**Fix:** use, for example, `dnorm(0, sd = 1.5)` (or `tau = 0.44`), refit, and report a prior-sensitivity check.

### 2.2 The depth smooth uses a basis that doesn't match its priors (Medium)

The latent occupancy smooth uses `X = predict(m1.0, type = "lpmatrix")` from a `gam(..., bs = "bs", select = TRUE, diagonalize = TRUE)` fit ([L25-28](../Scripts/05_Model_w_Nimble_m1.0+2LevelOccModel_wThaw.R#L25-L28), [L89](../Scripts/05_Model_w_Nimble_m1.0+2LevelOccModel_wThaw.R#L89)). `diagonalize` is a `jagam()` argument; `gam()` silently ignores it. The priors `b[2:4] ~ N(0, 1/λ1)` and `b[5] ~ N(0, 1/λ2)` were copied from the jagam template, which assumes the **diagonalised** basis.

I checked this directly:
- The `gam` model matrix differs from `q1Model_m1.0$jags.data$X`: the maximum absolute difference is 0.93.
- The `gam` penalty matrices are not diagonal.

So the prior is not the smoothing penalty implied by Eq. 1. It acts as an arbitrary ridge penalty on B-spline coefficients. The occupancy-by-depth curve is not shown in the paper, but the methods rely on it for the occupancy vs. GAM comparison ([zDist_methods.qmd:75](../zDist_methods.qmd#L75)).

**Fix:** use `q1Model_m1.0$jags.data$X`, and for prediction use `predict(q1Model_m1.0$pregam, newdata, type = "lpmatrix")` or the jagam basis. Alternatively, build the jagam object with the same `select`/`diagonalize` settings and take X from it.

### 2.3 The "POD at one freeze/thaw cycle" for DL is an extrapolation (Medium)

The observed thaw counts per technical replicate in `detect_data_clean`:

| Primer | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| DL | **0** | 62 | 3 | 192 | 88 | 171 | 0 |
| MFU | 502 | 9 | 0 | 0 | 4 | 1 | 0 |
| MV1 | 136 | 87 | 364 | 27 | 711 | 222 | 1 |

- DL was never run at Thaw = 1. The DL value of 0.35 at one cycle ([zDist_results.qmd:38](../zDist_results.qmd#L38)) is a linear-logit extrapolation from cycles 2–6. That value underpins two claims: the "MV1 > DL > MFU" ranking, and "~5–6 DL replicates for 90%". At the thaw levels DL actually experienced, its POD is 0.09–0.28, the lowest of the three primers (and its mean across samples, 0.14, is the lowest too).
- The MFU thaw slope rests on 14 replicates with Thaw > 1. With the tight prior from 2.1 it is largely prior-driven. Fig. 2 (right) extends every primer to Thaw = 7, where MFU's POD "rises" to 0.45 (the posterior slope is +0.09). That is extrapolation, not a finding.
- Thaw count is a plate- or run-level variable, so it is confounded with sequencing run, run date (2023–2025) and plate. The thaw effect can't be separated from run effects, and the paper should say so.

**Fix:** report each primer's POD at its observed thaw range, or at a common observed level. Restrict the plotted thaw range to the observed range for each primer. Discuss the confounding.

### 2.4 Eq. 3 doesn't match the code (Medium)

Eq. 3 ([zDist_methods.qmd:65-69](../zDist_methods.qmd#L65-L69)) has a single slope, γ_t. The code fits a slope for each primer, `thaw_par[Y_primer_index[i]]` ([L141](../Scripts/05_Model_w_Nimble_m1.0+2LevelOccModel_wThaw.R#L141)), and the results interpret primer-specific slopes. Eq. 3 should read `logit(Ω) = γ_m + γ_{t,m} T_r`.

### 2.5 Minor points on the occupancy model (Low)

- `inits = inits_fn_sitedepth()` is evaluated once ([L200](../Scripts/05_Model_w_Nimble_m1.0+2LevelOccModel_wThaw.R#L200)), so all four chains start from identical values. Gelman–Rubin assumes over-dispersed starts. The current R-hat values (≤ 1.03) and ESS (> 2,100) look fine, but pass the function itself (`inits = inits_fn_sitedepth`) and add jitter.
- `thin = 10` isn't mentioned in the methods.
- The "raw" points in Fig. 2 (right) (`raw_thaw_primer`) are unconditional detection proportions, roughly ψ × p. They sit on the same axis as POD *conditional on presence*, so they aren't comparable. Either drop them or state this in the caption.
- `N = q1Model_m1.0$jags.data$n` ([L155](../Scripts/05_Model_w_Nimble_m1.0+2LevelOccModel_wThaw.R#L155)) currently equals `nrow(mm.data)` (2,580). If `detect_data_allcet` is regenerated without re-running 02.5, the two will silently diverge. Use `nrow(mm.data)`.

---

## 3. Data wrangling ([01_Wrangle_data.R](../Scripts/01_Wrangle_data.R))

### 3.1 Replicate and dilution selection is on mammal reads, which selects on the outcome (Medium–High)

The methods say "keeping only the iteration with the greatest total read count" ([zDist_methods.qmd:43](../zDist_methods.qmd#L43)). In the code, the data are filtered to `Class == "Mammalia"` at L35 **before** `totReads <- sum(nReads)` is computed in the seqRep step ([L96-105](../Scripts/01_Wrangle_data.R#L96-L105)) and the dilution step ([L110-121](../Scripts/01_Wrangle_data.R#L110-L121)). So the code keeps whichever re-run or dilution had the most **mammal** reads. That biases the data toward detections.

I re-ran the selection on the raw table:
- Of 296 primer × sample × techRep groups with more than one sequencing rep or dilution, the choice differs in **86** (MFU 36, MV1 44, DL 6).
- Selecting by mammal reads gives **89** mammal detections in those groups; selecting by total reads gives 78.

This count includes pinnipeds and Delphinidae, and my re-run approximates the two-step logic, so treat the size as indicative. The direction is clear: detection rates are inflated in exactly the replicates that were re-run because they performed poorly.

**Fix:** use the existing `totalReads` column (computed over all taxa at L25, before the filter), or re-run selection before filtering to mammals.

### 3.2 Failed replicates are counted as non-detections (Medium)

The methods define non-detection as "0 reads for a given cetacean species in a station/depth location that was positive for non-target species … in either MFU or MV1" ([zDist_methods.qmd:43](../zDist_methods.qmd#L43)). The filter that would do this is commented out ([L34](../Scripts/01_Wrangle_data.R#L34)).

In `detect_data_clean`, **33 MFU and 96 MV1 replicates have zero reads of any taxon**, across 91 samples. In 2 samples, every MFU and MV1 replicate is empty. These are failed reactions, but they enter:
- the occupancy model as Bernoulli failures, which biases POD downward and conflates failure with non-detection; and
- the GAMs as 16 simultaneous zeros each.

**Fix:** implement the stated rule, or change the methods. Then re-fit.

### 3.3 Other wrangling points

- **Read threshold** (Low): the methods say "> 1 sequencing read", but the code uses `nReads > 0` ([L33](../Scripts/01_Wrangle_data.R#L33)). The raw table has 2 mammal cells with exactly 1 read; none survive into `detect_data_clean`, so there is no effect today. Align the code with the text anyway.
- **Whale-fall filter is unused** (Low–Medium): `detect_data_nowf` ([L193-199](../Scripts/01_Wrangle_data.R#L193-L199)) would zero **18** baleen and dolphin detections within 100 m of the seafloor, but it is never used. The discussion leans on whale fall to explain deep baleen detections. Decide whether to run the models with and without these detections as a sensitivity analysis, and document the choice.
- **Hard-coded counts** (Low): `slice_head(n = 22)` ([L119](../Scripts/01_Wrangle_data.R#L119)) relies on there being exactly 22 mammal taxa per replicate, and `n()/16` ([L258](../Scripts/01_Wrangle_data.R#L258)) on exactly 16 cetacean taxa. Both silently break if a taxon is added or removed. Use `n_distinct(BestTaxon)` or `distinct()`. `tooFewMV1` also counts MV1 techReps > 3 before they are dropped.
- **Typo** (Low): `mutate(diluti0n = ...)` ([L117](../Scripts/01_Wrangle_data.R#L117)) creates a stray column, so `arrange(dilution)` sorts the character column.
- **Stale comments** (Low): "515 samples, 41,200 rows, 2575 technical replicates" (L280). The current data have 516, 41,280 and 2,580.

---

## 4. Q2 GAMs ([02_Model_with_GAMs_Q1.R](../Scripts/02_Model_with_GAMs_Q1.R), [02.2](../Scripts/02.2_m1.2_performance.R), [02.1](../Scripts/02.1_Plot_H1_models.R))

### 4.1 Pseudoreplication (Medium–High)

Each row is a species × technical replicate (41,280 rows = 16 species × 2,580 replicates). The 5 replicates of a water sample, and the up-to-5 depths at a station, are treated as independent Bernoulli trials. There is no sample or station random effect; `ti(NWFSCsampleID, bs = "re")` is commented out at [L56](../Scripts/02_Model_with_GAMs_Q1.R#L56). There is also no spatial term in Q2.

This makes standard errors too small, the p-values and confidence bands in Figs. 3–4 too optimistic, and AIC too willing to pick complex models. The Pearson dispersion of 0.757 also suggests the model is misspecified for this correlated structure.

**Fix:** at minimum, add `s(NWFSCsampleID, bs = "re")`, or a station random effect, or aggregate to a binomial (successes out of 5) per sample × species with a sample random effect. Then re-check the model ranking.

### 4.2 Text vs. output (High for the text, easy to fix)

- [zDist_results.qmd:93](../zDist_results.qmd#L93): "Models stratifying detections by taxonomic family and preferred prey type both had **lower** AIC values than the model stratified by species." They are **higher**: m1.2 = 4547, m1.2a = 4854, m1.2b = 4960. The sentence and its "however" logic need rewriting.
- [zDist_results.qmd:56](../zDist_results.qmd#L56): "explained 14.3% of variance". The saved `m1.2` gives deviance explained = **13.13%**. Use inline R (`summary(m1.2)$dev.expl`) rather than a hard-coded value. Deviance explained is also not "variance explained".
- The AIC values in the code comments (5549, 5749, 5901, 6051) are stale. The actual values are m1.1 = 4689, m1.2a = 4854, m1.2b = 4960 and m1.0 = 5096.
- The `Stratification` labels in Table 2 ([L66](../zDist_results.qmd#L66)) are hard-coded by AIC rank. They happen to be correct now, but they will silently mislabel the table if the ranking changes. Join them by model name instead.
- m1.1 uses `s(depth)` (k = 10, thin-plate), while m1.0 uses `k = 5`. Minor, but the methods present these as a nested sequence.

### 4.3 Cross-validation (Medium)

- **Station leakage:** folds are assigned by `NWFSCsampleID` ([02.2:12-16](../Scripts/02.2_m1.2_performance.R#L12-L16), [04.2:13-18](../Scripts/04.2_m3.0c_performance.R#L13-L18)), so the other depths from the same station stay in the training set. For the xyz model (Q3), a flexible spatial smooth can effectively memorise stations, which inflates AUC and TSS. Use station-blocked, or spatially blocked, folds.
- **No `set.seed()`**, so the fold assignment is not reproducible.
- **Species subset:** the "clean" CV runs use `detect_per_species` from the **non-clean** data to choose species with >10 detections ([02.2:231](../Scripts/02.2_m1.2_performance.R#L231), [04.2:132](../Scripts/04.2_m3.0c_performance.R#L132)). The CV model is also fitted only to those 11 species, not to the 16-species model being reported.
- `pROC::roc()` defaults to `direction = "auto"`, which picks whichever direction gives AUC > 0.5 in each test fold. Set `direction = "<"`.
- **Supplementary Fig. S2** ([index.qmd:194](../index.qmd#L194)) plots `TSS_plot` and `TSS_plot_m3.0c`, both from the **non-clean** data. The text quotes the clean values for Q2. Use `TSS_plot_clean` and `TSS_plot_m3.0c_clean`.
- **Table S3** is referenced twice in the results, but its chunk is `eval: FALSE` and fully commented out ([index.qmd:161-180](../index.qmd#L161-L180)), so it doesn't exist in the rendered paper.
- **Per-species TSS bug:** `bin_label = breaks[names(breaks) == bins]` ([02.2:177](../Scripts/02.2_m1.2_performance.R#L177)) compares two vectors of different lengths, so R recycles them and assigns the wrong labels. Use `breaks[bins]`. Also, `ylab()` is called twice (L190-191), so the x-axis is never labelled. Neither is in the paper yet.

### 4.4 Time-at-depth overlay in Fig. 3 (Low–Medium)

The Fig. 3 caption says "Y-axis scale is the same for detection rate … and proportional time at depth." In the code, `time_scaled` ([02.1:119-120](../Scripts/02.1_Plot_H1_models.R#L119-L120)) min–max rescales time-at-depth to the range of `mu`, pooled across all species, while the panels use `scales = "free_y"`. The secondary-axis labels are also hidden. The dashed line therefore shows shape only, not magnitude, and its height is not comparable across panels. Fix the caption, or show the secondary axis for each panel.

---

## 5. Q3 methods text vs. code (Medium)

- **Coordinates:** the methods say the spline is on "projected coordinates (since they would obey the isotropy property)" ([zDist_methods.qmd:118](../zDist_methods.qmd#L118)). The models behind Fig. 5 and the AUC use **lon/lat in degrees**. At 38–48°N, 1° of longitude is about 0.67–0.79 of 1° of latitude, so an isotropic thin-plate spline on degrees is not isotropic in distance.
- **Eq. 10:** the equation writes `f_s(z)` (a depth smooth per species), but the text below it says "a smooth of depth, common across all species". The code (`s(depth, bs = "ts")`, [04:327](../Scripts/04_Model_with_GAMs_Q3.R#L327)) is common. Change the equation to `f(z)`. The equation is also missing a closing parenthesis: `logit(μ(i,z,s)`.
- **Density units:** "animals per kilometre". Becker et al. (2020) densities are animals per **km²**.
- **Density year:** the code filters the density predictions to `year == 2014` ([04:36](../Scripts/04_Model_with_GAMs_Q3.R#L36)), while the survey was in 2019. If that was intentional (for example, the latest year available), say so. Otherwise it is a bug.
- **Species subset:** the methods say the input data were subset to three species. The model used for the maps was fitted to all 16 (see 1.1).
- **CI band width:** the "50% CI depth range" ([04.1:30-41](../Scripts/04.1_Plot_H3_models.R#L30-L41)) collects the depths whose *point estimate* falls inside the 50% CI of the maximum, then takes max − min. For humpback whales, which have two separate peaks, that range spans the gap between them and overstates the uncertainty. Describe this construction in the methods, or compute the width differently (for example, from posterior simulation of the argmax depth).

---

## 6. Smaller numeric and text inconsistencies (Low)

| Location | Text | Data / output |
|---|---|---|
| [zDist_results.qmd:24](../zDist_results.qmd#L24), Discussion | 3 saber-toothed whale detections (1 at 300 m, 2 at 500 m) | 3 before cleaning (one was MV1 techRep 4, now dropped); **2 in the analysed data** (1 at 300 m, 1 at 500 m) |
| Discussion, whale-fall paragraph | "baleen whales were detected below 250 m at 28 locations"; 8 = 28.5% | 28 species × station-depth detections at **23** unique locations; 8/28 = **28.6%** |
| [zDist_results.qmd:6](../zDist_results.qmd#L6) | "1548 successfully amplified and sequenced samples" | 1,548 = sample × primer combinations; there are 2,580 sequenced replicates |
| [zDist_results.qmd:81](../zDist_results.qmd#L81) | "observed detections in Supplemental Figure S5" | That figure is **S3**; S5 is `eval: FALSE` |
| Table 2 labels vs. Eq. 4–8 | `β_s + f_F(z)`, `β_s + f_e(z)` | Methods use `β_F`, `β_e` |
| Fig. S1, S2 and S3 captions | "Q1 …", "Q@ model" | Q1 is now the occupancy model; these refer to Q2 |
| [zDist_results.qmd:115](../zDist_results.qmd#L115) | "interaction between species and $xys$" | Typo for $xyz$ |
| Discussion, Q1 | Occupancy POD is "higher … than GAM detection rate" because the occupancy model conditions on presence | The two also measure different things: P(any cetacean DNA in a bottle) vs. the per-species, per-replicate detection rate. The comparison figure isn't in the paper |
| [01.2_Examine_deep_detections.R:56](../Scripts/01.2_Examine_deep_detections.R#L56) | `crs = 4236` | Almost certainly meant 4326. EPSG:4236 is a different datum (Hu Tzu Shan 1950) |
| Sampling map ([01.2_Sampling_map.R:10](../Scripts/01.2_Sampling_map.R#L10)) | Fig. 1 map | Built from `detect_data` (non-clean), so it may show stations that were later dropped |

---

## 7. Checked and consistent

So the authors know what does *not* need attention, these checks passed against the saved outputs:

- **Sample and detection counts:** 516 samples, 175 stations, 2,580 replicates (516 DL, 516 MFU, 1,548 MV1); the depth breakdown 148/147/74/73/74; 468 detections of 16 species (6 baleen, 5 dolphin, 2 porpoise, 3 beaked); 11 species with >10 detections and 5 with fewer than 10; no station-depth with more than one sample after cleaning.
- **Q1 numbers** (conditional on the model as fitted):
  - POD at Thaw = 1: MV1 0.44, DL 0.35, MFU 0.31.
  - Replicates needed for 90%: MV1 4, DL 6 (0.885 at 5), MFU 7.
  - Mean per-primer POD: 0.14, 0.31, 0.39; lowest DL POD 0.087.
  - Sample-level POD: 0.87 (0.85–0.92).
  - "5 of 516 MFU samples with >2 thaws".
- **MCMC:** R-hat ≤ 1.03; ESS > 2,100 for every parameter.
- **Q2:** the AIC values in Table 2; the m1.0 deviance explained (0.7%) and p ≈ 1e-6; the clean-data AUC of 0.73 and TSS of 0.38.
- **Discussion counts:** "19 deep-water detections (68%) unexplained"; Delphinidae losses of 2 MFU and 17 MV1 (pre-cleaning; 16 MV1 after cleaning).
- **Fig. 5 caption:** maximum detection rates of about 0.30, 0.35 and 0.54 match `maxPOD_depth_clipped_clean`.

---

## 8. Suggested order of work

1. **Data first (3.1, 3.2):** fix the replicate selection and the handling of failed replicates, then regenerate `detect_data_clean`. Commit `duplicated_station_depth.csv`.
2. **Occupancy model (2.1–2.4):** fix the priors and the spline basis, update Eq. 3, and report POD within the observed thaw range.
3. **GAMs (4.1):** add sample and station random effects, and re-check the model ranking.
4. **Q3 (1.1, 1.2, 5):** settle on one model specification, fit on projected coordinates, apply the bathymetry mask, and compute the area statistics in code.
5. **Cross-validation (4.3):** use station-blocked folds, `set.seed()`, the same model and data as reported, and the correct S2 figure.
6. **Text (4.2, 6):** replace every hard-coded number in the `.qmd` files with inline R that reads saved objects, and fix the AIC sentence and the cross-references.
7. **Reproducibility (1.3):** add a `run_all.R` (or `targets`) pipeline, and keep large outputs out of git (LFS, or regenerate them on demand).
