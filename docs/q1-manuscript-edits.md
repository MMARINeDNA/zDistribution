# Q1 occupancy model: manuscript edits to make

*Notes from the Q1 rerun on branch `claude/q1-occupancy-fixes`, 26 October 2026. Scripts: `Scripts/01_Wrangle_data.R` and `Scripts/05_Model_w_Nimble_m1.0+2LevelOccModel_wThaw.R`.*

## Dataset changes (affect all three questions)

`01_Wrangle_data.R` now (agreed in [#23](https://github.com/MMARINeDNA/zDistribution/issues/23)):

1. **Drops failed MFU/MV1 runs** (zero reads of any taxon) before choosing among re-runs and dilutions, so a successful re-run is used where one exists. DL runs with zero reads are kept as genuine non-detections, since DL amplifies only cetaceans.
2. **Keeps each sample's first successful replicates:** 1 DL, 1 MFU and 3 MV1, by lowest replicate number. Spare replicates (MFU #2, MV1 #4–6) fill in for failed ones; 32 replicates come in this way.
3. **Keeps only samples with all five.** Previously the completeness check counted replicates before trimming to DL/MFU #1 and MV1 #1–3, and counted failed runs as successful.

| | Before | Now |
|---|---|---|
| Samples | 516 | **470** |
| Stations | 175 | **167** |
| Technical replicates | 2,580 | **2,350** (470 DL, 470 MFU, 1,410 MV1) |
| Detections | 468 | **447** |
| Samples by depth (0 / 50 / 150 / 300 / 500 m) | 148 / 147 / 74 / 73 / 74 | **129 / 132 / 70 / 69 / 70** |
| Species with > 10 detections | 11 | **10** (gray whale now has exactly 10) |
| Saber-toothed whale detections | 2 | **3** (matches the text again) |

These counts appear in the abstract, the methods (sample collection) and the first paragraph of the results, and most are inline R. **The Q2 and Q3 models and their outputs also need rerunning on the new data.**

## Q1 model changes

- **Depth spline:** uses the jagam (diagonalised) basis that matches its priors. It previously used a `gam()` basis that didn't match.
- **Thaw centred on 1 cycle:** each intercept is a primer's POD at one freeze/thaw cycle (previously at zero cycles, which no sample had).
- **Priors:**
  - Intercepts: weakly informative, `dnorm(0, sd = 1.5)`, so POD at one thaw is roughly 0.05–0.95 (95%). Previously `dnorm(0, tau = 10)`, i.e. SD 0.32, which pulled every primer's POD toward 0.5.
  - DL and MV1 thaw slopes: weakly informative, `dnorm(0, sd = 1)`, so the odds of detection change by ×0.14 to ×7.1 per cycle (95%). Previously SD 0.32.
  - **MFU thaw slope fixed at 0** (an assumption, not an estimate). MFU was run at more than one cycle in only 26 of 470 replicates, 20 of them dilution-10 re-runs on one sequencing run (MURI318), so its thaw effect can't be separated from run and dilution.
- **MCMC:** each chain starts from different, jittered values; a seed is set; monitors are listed explicitly.

## Q1 results (current run)

Converged: R-hat upper bound ≤ 1.03; minimum effective sample size 2,312.

| Primer | POD at 1 thaw (95% CI) | Thaw slope, logit per cycle (95% CI) | Replicates for 90% at 1 thaw |
|---|---|---|---|
| MV1 | 0.44 (0.36–0.52) | −0.06 (−0.15 to 0.04) | 4 (POD 0.90) |
| DL | 0.33 (0.17–0.52) | −0.32 (−0.58 to −0.08) | 6 (POD 0.91) |
| MFU | 0.27 (0.21–0.33) | 0 (fixed) | 8 (POD 0.92) |

- DL's value at one thaw is an extrapolation: DL was run only at 2–6 cycles. DL's POD at the thaws it was actually run at: 0.26 (2), 0.20 (3), 0.16 (4), 0.12 (5), 0.09 (6).
- Mean POD per replicate across the data: DL 0.14, MFU 0.27, MV1 0.40.
- Per-sample POD across all five replicates: mean 0.86, range 0.83–0.91.

## Sensitivity: fixing vs estimating the MFU slope (for the supplement)

Same 470-sample data; the only difference is whether MFU's slope is estimated (SD 1 prior).

| | MFU slope fixed at 0 (current) | MFU slope estimated |
|---|---|---|
| DL slope | −0.32 (−0.58 to −0.08) | −0.33 (−0.58 to −0.08) |
| MFU slope | 0 (fixed) | −0.27 (−0.93 to 0.21) |
| MV1 slope | −0.06 (−0.15 to 0.04) | −0.06 (−0.15 to 0.04) |
| DL POD, 1 thaw | 0.33 | 0.34 |
| MFU POD, 1 thaw | 0.27 | 0.27 |
| MV1 POD, 1 thaw | 0.44 | 0.44 |

With the complete-sets data, MFU's estimated slope is uncertain and its 95% interval includes 0, so fixing it at 0 is consistent with the data. Nothing else changes.

On the earlier 516-sample data (failed runs counted as non-detections and later dropped), MFU's estimated slope was positive (0.48 to 1.54 depending on the prior). That signal came from 4 positive re-run replicates at five thaws; three of those samples don't have complete successful sets and drop out of the new data. The swing shows how little information there is on MFU's thaw effect.

## Methods edits

- [ ] **Sample counts** (zDist_methods.qmd:11): 470 samples at 167 stations; 129 / 132 / 70 / 69 / 70 samples at 0 / 50 / 150 / 300 / 500 m. Check whether the latitude/longitude range changes.
- [ ] **Detection data QAQC** (zDist_methods.qmd:43): describe the failed-run rule and the complete-set rule. For example: "MFU and MV1 replicates with no reads of any taxon were treated as failed reactions and excluded. For each sample we retained the first successful DL, MFU and three MV1 replicates; samples without a complete set of five successful replicates were excluded."
- [ ] **Eq. 3:** change to primer-specific slopes: logit(Ω) = γ_m + γ_{t,m} · (T_r − 1), with γ_{t,MFU} = 0. Define T_r − 1 as thaw cycles beyond the first.
- [ ] **Priors** (zDist_methods.qmd:71): replace "uninformative or weakly informative" with an accurate description: weakly informative normal priors on the logit scale, SD 1.5 on the primer intercepts and SD 1 on the DL and MV1 thaw slopes.
- [ ] **MFU slope assumption:** state that MFU was run at more than one freeze/thaw cycle in only 26 of 470 replicates, mostly re-runs from a single sequencing run, so its thaw effect was fixed at zero; estimating it gave −0.27 (−0.93 to 0.21) and left the other estimates unchanged (supplementary table).
- [ ] **Thinning:** mention thin = 10 (2,500 retained draws per chain).
- [ ] **Starting values:** optionally mention that each chain started from dispersed values.

## Results edits

- [ ] **Opening paragraph and abstract:** update the counts in the table above (most are inline R and will update on re-render; check the abstract's hard-coded "516 samples", "175 unique geographic locations", "468 positive cetacean detections" and "2,580 sequencing replicates").
- [ ] **Primer ranking** (zDist_results.qmd:38): MV1 0.44, DL 0.33, MFU 0.27, the same order as before. The values are inline R.
- [ ] **DL at one thaw is an extrapolation:** DL was only run at 2–6 cycles. Say so wherever DL's one-thaw value or replicate count is quoted.
- [ ] **Replicates for 90%:** MV1 4, DL 6, MFU 8 (text says four, ~5–6, seven). Make these inline R rather than hard-coded.
- [ ] **Thaw sentences** (zDist_results.qmd:40): DL −0.32 (−0.58 to −0.08) per cycle; MV1 −0.06 (−0.15 to 0.04), indistinguishable from zero, so "weaker relationship for MV1" could become "no clear relationship for MV1". For MFU, replace "no relationship detected … only 5 of 516 total MFU samples sequenced in this study had > 2 thaw cycles" with the fixed-slope explanation above.
- [ ] **Hard-coded numbers** to replace with inline R: "0.14, 0.31, 0.39" (now DL 0.14, MFU 0.27, MV1 0.40), "lowest DL POD = 0.08" (now 0.09 at six thaws), "four / ~5–6 / seven technical replicates" (now 4 / 6 / 8). The per-sample POD (0.86, range 0.83–0.91) is already inline.
- [ ] **Fig. 2 (right):** restrict each primer's line to the thaw counts it was actually run at (DL 2–6; MFU 1, 2, 5, 6, 7; MV1 1–7). MFU's line is flat by construction, so the caption should say its slope is fixed at zero. In the caption, either drop the raw-proportion points or say they aren't conditional on presence.

## Discussion edits

- [ ] The argument that shorter fragments (MFU, MV1) are robust to freeze/thaw should rest on MV1, whose slope is near zero and well estimated. MFU's effect was assumed, not estimated.
- [ ] The best-practice recommendation (thaw extracts as few times as possible) still follows from DL.
- [ ] Re-check the deep-detection and whale-fall counts (28 detections below 250 m, 8 near the seafloor, 19 unexplained) on the new data.

## Downstream scripts (rerun on the 470-sample data)

- [x] `Scripts/05_Calc_POD_by_sample.R` (per-sample POD)
- [x] `Scripts/05_Plot_TechRepsXPrimer.R` (replicates figure; now loads tidyr, which it was missing)
- [ ] Q2 GAMs (`02_Model_with_GAMs_Q1.R`, `02.1`, `02.2`) and Q3 models (`04_*`) still use outputs built from the old data.
