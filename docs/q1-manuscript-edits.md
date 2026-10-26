# Q1 occupancy model: manuscript edits to make

*Notes from the Q1 rerun on branch `claude/q1-occupancy-fixes`, 26 October 2026. Script: `Scripts/05_Model_w_Nimble_m1.0+2LevelOccModel_wThaw.R`.*

## What changed in the model

- **Depth spline:** now uses the jagam (diagonalised) basis that matches its priors. It previously used a `gam()` basis that didn't match.
- **Failed replicates:** 129 MFU/MV1 replicates with zero reads of any taxon are dropped (2,451 replicates remain). DL replicates with zero reads are kept as non-detections. Pending Amy's confirmation in [#23](https://github.com/MMARINeDNA/zDistribution/issues/23).
- **Thaw centred on 1 cycle:** each intercept is now a primer's POD at one freeze/thaw cycle (previously at zero cycles, which no sample had).
- **Priors:**
  - Intercepts: weakly informative, `dnorm(0, sd = 1.5)`, so POD at one thaw is roughly 0.05–0.95 (95%). Previously `dnorm(0, tau = 10)`, i.e. SD 0.32, which pulled every primer's POD toward 0.5.
  - DL and MV1 thaw slopes: weakly informative, `dnorm(0, sd = 1)`, so the odds of detection change by ×0.14 to ×7.1 per cycle (95%). Previously SD 0.32.
  - **MFU thaw slope: fixed at 0** (an assumption, not an estimate). MFU was run at more than one cycle in only 14 replicates, all dilution-10 re-runs from essentially one sequencing run (MURI318), so its thaw effect can't be separated from run, plate and dilution.
- **MCMC:** each chain starts from different, jittered values; a seed is set; parameters to monitor are listed explicitly.

## Results (current run)

Converged: R-hat upper bound ≤ 1.004 for every detection parameter. The depth-spline coefficients are slightly higher (`b_depth[5]` 1.05, effective sample size 1,294), still within the usual 1.1 threshold. Minimum effective sample size 1,294.

| Primer | POD at 1 thaw (95% CI) | Thaw slope, logit per cycle (95% CI) | Replicates for 90% at 1 thaw |
|---|---|---|---|
| MV1 | 0.44 (0.36–0.52) | −0.04 (−0.13 to 0.05) | 4 (POD 0.901) |
| DL | 0.28 (0.14–0.45) | −0.26 (−0.50 to −0.02) | 8 (seven give 0.897) |
| MFU | 0.29 (0.23–0.35) | 0 (fixed) | 7 |

DL's value at one thaw is an extrapolation (DL was run only at 2–6 cycles). DL's POD at the thaws it was actually run at: 0.23 (2), 0.18 (3), 0.15 (4), 0.12 (5), 0.09 (6).

## Prior sensitivity (for a supplementary table)

| | Original model | **Current: MFU fixed, slope SD 1** | Slope SD 0.32 | Slope SD 1 | Slope SD 1.5 |
|---|---|---|---|---|---|
| DL slope | −0.35 (−0.50 to −0.21) | −0.26 (−0.50 to −0.02) | −0.23 (−0.46 to 0.00) | −0.27 (−0.51 to −0.03) | −0.27 (−0.51 to −0.02) |
| MFU slope | 0.09 (−0.26 to 0.43) | 0 (fixed) | 0.48 (0.08 to 0.93) | 1.18 (0.33 to 2.46) | 1.54 (0.41 to 3.48) |
| MV1 slope | −0.06 (−0.14 to 0.02) | −0.04 (−0.13 to 0.05) | −0.04 (−0.13 to 0.05) | −0.04 (−0.14 to 0.05) | −0.04 (−0.14 to 0.05) |
| DL POD, 1 thaw | 0.35 | 0.28 | 0.26 | 0.28 | 0.28 |
| MFU POD, 1 thaw | 0.31 | 0.29 | 0.28 | 0.28 | 0.28 |
| MV1 POD, 1 thaw | 0.44 | 0.44 | 0.44 | 0.44 | 0.44 |

The last three columns estimate MFU's slope and differ only in the slope prior (intercept SD 1.5 throughout). The "original model" column also differs in its spline basis, intercept prior, thaw reference and replicate set, so it isn't a pure prior comparison.

**Takeaways:**
- Fixing MFU's slope at 0 doesn't change the DL or MV1 results, or any primer's POD at one thaw.
- MV1's slope and every primer's POD at one thaw are insensitive to the slope prior.
- DL's slope is negative under every prior; under SD 0.32 its 95% interval just reaches 0.
- When MFU's slope is estimated, it is positive under every corrected prior and moves with the prior (0.48 to 1.54). That supports treating it as not estimable rather than as a finding. In the original model it was most likely held near zero by the tight *intercept* prior at zero thaws rather than by the slope prior (not tested directly).

## Methods edits

- [ ] **Eq. 3:** change to primer-specific slopes: logit(Ω) = γ_m + γ_{t,m} · (T_r − 1), with γ_{t,MFU} = 0. Define T_r − 1 as thaw cycles beyond the first.
- [ ] **Priors** (zDist_methods.qmd:71): replace "uninformative or weakly informative" with an accurate description: weakly informative normal priors on the logit scale, SD 1.5 on the primer intercepts and SD 1 on the DL and MV1 thaw slopes.
- [ ] **MFU slope assumption:** state that MFU was run at more than one freeze/thaw cycle in only 14 replicates, all dilution-10 re-runs from a single sequencing run, so its thaw effect could not be estimated and was fixed at zero. Point to the sensitivity table.
- [ ] **Thinning:** mention thin = 10 (2,500 retained draws per chain).
- [ ] **Failed replicates:** say that MFU/MV1 replicates with zero reads of any taxon were treated as failed reactions and excluded from the occupancy model (once Amy confirms, #23).
- [ ] **Starting values:** optionally mention that each chain started from dispersed values.

## Results edits

- [ ] **Primer ranking** (zDist_results.qmd:38): MV1 is still highest (0.44); MFU (0.29) and DL (0.28) are now essentially tied, so "followed by DL and then MFU" no longer holds.
- [ ] **DL at one thaw is an extrapolation:** DL was only run at 2–6 cycles. Say so wherever DL's one-thaw value or replicate count is quoted.
- [ ] **Replicates for 90%:** MV1 4, MFU 7, DL 8 (was 4, 7, ~5–6). Make these inline R rather than hard-coded.
- [ ] **MFU thaw sentence** (zDist_results.qmd:40): replace "no relationship detected … noting limited statistical power" with: MFU's thaw effect could not be estimated, because only 14 MFU replicates were run at more than one cycle, all re-runs from a single sequencing run, so it was assumed to be zero.
- [ ] **DL and MV1 thaw sentences:** DL's slope is −0.26 (−0.50 to −0.02) per cycle; MV1's is −0.04 (−0.13 to 0.05), indistinguishable from zero. "Weaker relationship for MV1" could become "no clear relationship for MV1".
- [ ] **Hard-coded numbers** to replace with inline R: "0.14, 0.31, 0.39" (mean POD by primer), "lowest DL POD = 0.08", "four / ~5–6 / seven technical replicates", and the per-sample POD mean and range once `05_Calc_POD_by_sample.R` is rerun.
- [ ] **Fig. 2 (right):** restrict each primer's line to the thaw counts it was actually run at (DL 2–6; MFU 1, 2, 5, 6; MV1 1–7). MFU's line is now flat by construction, so the caption should say its slope is fixed at zero. In the caption, either drop the raw-proportion points or say they aren't conditional on presence.

## Discussion edits

- [ ] The argument that shorter fragments (MFU, MV1) are robust to freeze/thaw should rest on MV1, whose slope is near zero and well estimated. MFU's effect was assumed, not estimated.
- [ ] The best-practice recommendation (thaw extracts as few times as possible) still follows from DL.

## Downstream scripts (rerun)

- [x] `Scripts/05_Calc_POD_by_sample.R` (per-sample POD). Now drops the same failed MFU/MV1 replicates as the model.
  - Mean POD per replicate by primer: DL 0.13, MFU 0.29, MV1 0.41 (text currently says 0.14, 0.31, 0.39).
  - Per-sample POD across all five replicates: mean 0.84, range 0.09–0.90 (text currently says 0.87, range 0.85–0.92).
  - The low end of the range comes from the 2 samples whose MFU and MV1 replicates all failed, leaving only their DL replicate. 64 samples have 4 usable replicates, 18 have 3, 7 have 2, and 2 have 1. The text could give the range for samples with all five replicates, or say why the minimum is low.
- [x] `Scripts/05_Plot_TechRepsXPrimer.R` (replicates figure). Fixed a missing `library(tidyr)`, which had stopped the script from running on its own. POD for k replicates at one thaw: MV1 0.90 at 4; MFU 0.91 at 7; DL 0.90 at 7 (just under the threshold) and 0.93 at 8.
