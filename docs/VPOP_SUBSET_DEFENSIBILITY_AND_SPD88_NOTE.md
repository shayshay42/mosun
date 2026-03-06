# VPop Subset Defensibility and New SPD Dataset Note

Date: 2026-03-05

## Scope

This note answers:
1. Whether the currently varied VPop parameter subset is defensible for the current calibration targets.
2. How to record and use the new SPD dataset target for 21-day-cycle step-up dosing.

## What is varied today

Current VPop generation varies 17 patient-level columns total:
- Core tumor-design variables (6):
  - `kBtumorprolif`, `KBptumor`, `KTrptumor`, `BT_ratio_tumor_init`, `Btumor_perml_init`, `tumor_burden_factor`
- Extra kinetic variables (11) from `sample_and_prune_vpop.py` default list:
  - `kBapop`, `kBkill`, `fBkill`, `kBprolif`, `kTaexit`, `kTact`, `fTadeact`, `fTap`, `kTaapop`, `fTaprolif`, `fTrapop`

Implementation references:
- Default extra list: `scripts/sample_and_prune_vpop.py` (`default_extra_names`)
- Priors source: `generated/vpop_targets/parameter_prior_bounds_from_params_sheet.csv`
- Run-level record of selected names: `pruning_summary.json -> inputs.sampled_extra_params`

## Current target set actually constraining VPop

`generated/vpop_targets/vpop_target_evaluation_plan.json` currently binds to:
- 2 regimens: `step_1_2_13_5mg`, `step_1_2_60_30mg`
- Horizon: `84` days
- Features: day-42 tumor proxy and day-0..21 IL6/CRS proxy
- Target modes: parameter bounds, population fractions, pairwise gap constraints

Not directly evaluable (already marked in plan):
- PK dose proportionality window
- Baseline SPD distribution mapping
- Long-term outcomes (PFS/OS/DOR)
- Baseline residual rituximab detectability

## Is the current subset defensible?

Short answer: **yes, but only as a first-pass proxy-calibrated subset**.

Why it is defensible for the current objective:
- All 11 extra variables are bounded and available in the prior table (Cyno LB/UB extraction).
- They are mechanistically central to B/T turnover and activation dynamics that drive the two current observables (tumor day-42 proxy and IL6 day0-21 proxy).
- Dimensionality is moderate, which keeps random-subset pruning computationally manageable.

Why it is not fully sufficient clinically:
- The active target set does not yet constrain PK behavior or long-term SPD-based disease trajectories.
- Several bounded estimated parameters likely relevant to exposure-response are currently fixed (for example: `fTact`, `KmBT_act`, `ndrugactT`, `KdrugactT`, `KmTB_kill`, `fKmTB_kill`, `nkill`, `KdrugB`, `kTrexit`, `fa0`, `fAICD`, `fTa0deact`, `fTa0apop`).
- Therefore, identifiability and realism can be limited when calibrating to richer endpoints.

Conclusion:
- Keep this subset for current 84-day proxy calibration workflows.
- Expand subset in phases once SPD long-horizon targets are added.

## New dataset note (user-provided)

User-provided clinical scenario to add as calibration target:
- 21-day cycles
- Cycle 1 step-up: C1D1 `1 mg`, C1D8 `2 mg`, C1D15 `60 mg`
- Cycle 2: C2D1 `60 mg`
- Cycle 3+: `30 mg` Q3W
- Endpoint: best SPD % change over 8 cycles

Likely event-day map (days from treatment start):
- C1D1: day 0, 1 mg
- C1D8: day 7, 2 mg
- C1D15: day 14, 60 mg
- C2D1: day 21, 60 mg
- C3D1: day 42, 30 mg
- C4D1: day 63, 30 mg
- C5D1: day 84, 30 mg
- C6D1: day 105, 30 mg
- C7D1: day 126, 30 mg
- C8D1: day 147, 30 mg
- Suggested simulation horizon for 8 cycles: day 168

## Important data-format caveat

Current file `assets/musun_2022_88patients_SPD.csv` in repo has 41 rows with columns `x` and ` y` (not a clear 88-row patient table).

Before using as a hard target, confirm whether it is:
- full patient-level SPD (% best change) values, or
- digitized plot coordinates from a published waterfall/curve.

## Recommended integration path for this new target

1. Add regimen to evaluation plan:
- `go29781_stepup_c1_1_2_60_c2_60_c3plus_30_q3w`
- with dose events above and `horizon_day = 168`.

2. Add feature binding for SPD proxy:
- `spd_pct(t) = 100 * (Btumor(t)/Btumor(0) - 1)`
- `best_spd_pct_8cycles = min_t spd_pct(t), t in [0,168]`

3. Add target bindings:
- distribution quantile anchors (p10, p25, p50, p75, p90) from patient-level SPD if available
- response-threshold fractions (for example <= -50%, <= -30%, >= +20%)

4. After SPD target activation, expand varied parameter subset (phased):
- Phase A (current 11 + core 6): keep as baseline
- Phase B: add ER potency/slope terms (`fTact`, `KmBT_act`, `KdrugactT`, `KmTB_kill`, `fKmTB_kill`, `nkill`, `KdrugB`, `ndrugactT`)
- Phase C: add T-cell homeostasis/adaptation terms (`kTrexit`, `fa0`, `fAICD`, `fTa0deact`, `fTa0apop`) only if needed by calibration diagnostics

## Practical recommendation

Treat the current subset as a defensible **MVP subset** for the existing proxy targets, not as the final clinical VPop definition. The new SPD-over-8-cycles dataset should become a first-class calibration target, and then subset expansion should be data-driven from calibration residuals and sensitivity analysis.
