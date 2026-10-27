# Q1 occupancy model: manuscript edits to make

*Revised 27 October 2026 for the minor revision. Scripts: `Scripts/01_Wrangle_data.R`, `Scripts/05_Model_w_Nimble_m1.0+2LevelOccModel_wThaw.R` and `Scripts/05_Plot_TechRepsXPrimer.R`. The text edits are tracked in [#26](https://github.com/MMARINeDNA/zDistribution/issues/26).*

## Dataset

`01_Wrangle_data.R` applies the non-detection rule as the methods state it (confirmed by Amy in [#23](https://github.com/MMARINeDNA/zDistribution/issues/23)). It works on the station/depth sample:

- If **every** MFU and MV1 run for a sample returned zero reads of any taxon, the sample is treated as failed and dropped. DL is assumed to have failed too; in both such samples it also had zero reads.
- Otherwise, a zero for a cetacean species is a non-detection, including in a failed MFU or MV1 run.

The replicate selection is the original one: DL #1, MFU #1 and MV1 #1–3.

| | Paper | Now |
|---|---|---|
| Samples | 516 | **514** |
| Stations | 175 | 175 |
| Technical replicates | 2,580 | **2,570** (514 DL, 514 MFU, 1,542 MV1) |
| Detections | 468 | 468 |
| Samples by depth (0 / 50 / 150 / 300 / 500 m) | 148 / 147 / 74 / 73 / 74 | 148 / **146** / 74 / 73 / **73** |

The two dropped samples (52193-1353 at 50 m and 52193-1682 at 500 m) had no detections. Apart from them, the data are identical to the data behind the submitted paper. An earlier change (#25) dropped every failed run individually and required complete sets of successful replicates (470 samples); that was stricter than the methods and has been reverted.

## Q1 model

The priors and thaw coding are the paper's: `dnorm(0, 10)` on every primer intercept and thaw slope, and thaw cycles counted from 0. Three bug fixes from #25 are kept; none of them changes the detection estimates:

- **Depth spline:** uses the jagam basis that matches its priors (`gam()` ignored `diagonalize = TRUE`).
- **Starting values:** each chain starts from different, jittered values, and seeds are set.
- **`N`:** taken from the data rather than from a stored object.

## Q1 results

They match the submitted paper. Converged: R-hat ≤ 1.03; minimum effective sample size 1,394.

| Primer | POD at 1 thaw (95% CI) | Thaw slope, logit per cycle (95% CI) | Replicates for 90% at 1 thaw |
|---|---|---|---|
| MV1 | 0.44 (0.38–0.51) | −0.06 (−0.14 to 0.02) | 4 (POD 0.90) |
| DL | 0.35 (0.26–0.45) | −0.35 (−0.50 to −0.21) | 6 (POD 0.93; 0.89 at 5) |
| MFU | 0.31 (0.26–0.37) | +0.09 (−0.26 to 0.43) | 7 (POD 0.93) |

- Mean POD per replicate across the data: DL 0.14, MFU 0.31, MV1 0.39. DL's lowest POD (6 thaws) is 0.087.
- Per-sample POD across all five replicates: mean 0.87, range 0.85–0.92.
- MFU was run at more than one thaw in 14 of 514 replicates (13 of them dilution-10 re-runs on one sequencing run, MURI318), and at more than two thaws in 5.

## How much the priors matter

The priors are informative: in nimble, `dnorm(0, 10)` has precision 10, i.e. SD 0.32 on the logit scale, so about 95% of the prior mass for POD at zero thaws lies between 0.35 and 0.65. Refitting the same data with weakly informative priors (SD 1.5 on the intercepts, SD 1 on the slopes) gives:

| | Informative priors (paper) | Weakly informative priors |
|---|---|---|
| MFU thaw slope | +0.09 (−0.26 to 0.43) | +1.21 (0.39 to 2.47) |
| DL thaw slope | −0.35 (−0.50 to −0.21) | −0.26 (−0.51 to −0.02) |
| MV1 thaw slope | −0.06 (−0.14 to 0.02) | −0.05 (−0.14 to 0.04) |
| POD at 1 thaw, DL / MFU / MV1 | 0.35 / 0.31 / 0.44 | 0.28 / 0.26 / 0.42 |

The MV1 results and the direction of the DL thaw effect don't depend on the prior. Two results do:

- **MFU's thaw slope.** The data can't separate a thaw effect from sequencing run and dilution, so the near-zero estimate reflects the prior. The text shouldn't present it as evidence of no relationship.
- **DL's POD at one thaw.** DL was only run at 2–6 thaws, so its value at one thaw is an extrapolation, and it leans on the intercept prior (posterior SD 0.27 vs. prior SD 0.32).

The weakly informative runs were done in scratch and aren't saved in the repo.

## Fig. 2

Dashed lines mark extrapolation in every panel:

- **Left** (POD at one thaw) and **middle** (POD vs. number of replicates, at one thaw): DL is dashed, because it was never run at one thaw. The 90% reference line in the middle panel is dotted.
- **Right** (POD vs. freeze/thaw cycles): each primer is solid over the thaw range it was run at (DL 2–6, MFU 1–6, MV1 1–7) and dashed, with a lighter ribbon, outside it. Raw detection proportions are not shown, since they are not conditional on presence.
