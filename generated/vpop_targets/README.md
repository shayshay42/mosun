# VPop Target Bundle

This folder contains machine-readable targets for virtual-patient sampling/pruning.

## Files

- `vpop_calibration_targets.json`
  - Curated target spec with:
    - hard priors,
    - cohort-composition targets,
    - dose-schedule constraints,
    - clinical endpoint targets,
    - PK/baseline distribution anchors,
    - known data gaps.
- `parameter_prior_bounds_from_params_sheet.csv`
  - Extracted bounds from `assets/params_41540_2020_145_MOESM2_ESM.xlsx` (rows with explicit `Cyno_LB/Cyno_UB`).
- `vpop_calibration_targets_flat.csv`
  - Flat scalar target table derived from `vpop_calibration_targets.json` for quick ingestion in pruning code.
- `vpop_target_evaluation_plan.json`
  - Executable mapping of each target to:
    - dosing scenario,
    - simulation horizon/config,
    - model observable or derived proxy function.

## Rebuild

Run:

```bash
python3 scripts/build_vpop_target_bundle.py
```

## Sample-and-prune workflow

Run:

```bash
python3 scripts/sample_and_prune_vpop.py
```

Outputs go to `generated/vpop_pruning/<run_tag>/`.

If Julia is not on `PATH`, pass:

```bash
python3 scripts/sample_and_prune_vpop.py --julia-bin /path/to/julia --julia-project /path/to/project
```

## Important note

The follow-up paper extraction has a tumor-burden spread inconsistency:
- `SPD +/-10%` in extracted Table S1 text.
- `+/-50%` in extracted supplementary methods prose.

Both are retained in JSON (`hard` and `sensitivity` entries).
