# VPop Generation + Sensitivity (Windows Parallel Runbook)

Date: 2026-03-05

This runbook covers:
- VPop target bundle usage,
- candidate generation and pruning,
- Windows parallel subset-search execution,
- eFAST sensitivity for scalar output `best %SPD`.

## 1) What is included

Targets bundle:
- `generated/vpop_targets/vpop_calibration_targets.json`
- `generated/vpop_targets/vpop_target_evaluation_plan.json`
- `generated/vpop_targets/parameter_prior_bounds_from_params_sheet.csv`
- `generated/vpop_targets/vpop_calibration_targets_flat.csv`

Core scripts:
- `scripts/sample_and_prune_vpop.py`
- `scripts/generate_phase1_design.py`
- `scripts/generate_phase1_standard_design.py`
- `scripts/windows/run_vpop.ps1`
- `scripts/windows/run_vpop_parallel.ps1`

Sensitivity scripts:
- `julia/run_efast_best_spd.jl`
- `scripts/run_efast_best_spd.sh` (Unix runner)
- `scripts/plot_efast_best_spd.py`

Docs:
- `docs/VPOP_SUBSET_DEFENSIBILITY_AND_SPD88_NOTE.md`
- `docs/VPOP_GENERATION_WINDOWS_PARALLEL_RUNBOOK.md` (this file)

Additional dataset:
- `assets/musun_2022_88patients_SPD.csv`

## 2) One-time setup on Windows

From repository root in PowerShell:

```powershell
pwsh ./scripts/windows/setup_windows.ps1 -RepoRoot . -PythonExe python -JuliaExe julia
```

## 3) Single-run VPop generation (baseline)

```powershell
pwsh ./scripts/windows/run_vpop.ps1 `
  -RepoRoot . `
  -JuliaExe julia `
  -RunTag vpop_windows_single `
  -NCandidates 384 `
  -NSelect 140 `
  -NRandomSubsets 100000
```

## 4) Parallel VPop pruning on Windows

This mode:
1. Runs one base simulation pass to create a metrics table.
2. Launches multiple worker jobs in parallel, each performing random subset search on the same metrics table with different `subset_seed`.
3. Picks best worker by `objective.best_score`.
4. Optionally publishes best selected cohort to `assets/generated_vpop/`.

```powershell
pwsh ./scripts/windows/run_vpop_parallel.ps1 `
  -RepoRoot . `
  -JuliaExe julia `
  -RunTagPrefix vpop_parallel_win `
  -NCandidates 384 `
  -NSelect 140 `
  -NRandomSubsets 100000 `
  -NumWorkers 8 `
  -Seed 20260228 `
  -SpreadMode primary `
  -PublishBestToAssets
```

Outputs:
- Base run:
  - `generated/vpop_pruning/<RunTagPrefix>_base/`
- Worker runs:
  - `generated/vpop_pruning/<RunTagPrefix>_w*/`
- Combined summary:
  - `generated/vpop_pruning/<RunTagPrefix>_parallel_summary.json`
- Optional published best cohort:
  - `assets/generated_vpop/selected_patients.csv`
  - `assets/generated_vpop/pruning_summary.json`

## 5) Important seed behavior for parallel pruning

`sample_and_prune_vpop.py` now supports:
- `--seed`: candidate sampling seed,
- `--subset-seed`: subset-search seed.

For valid `--skip-sim` reuse, workers keep the same `--seed` (same candidate design) and vary `--subset-seed` (different random subset draws).

## 6) eFAST sensitivity for `best %SPD`

Output definition in `julia/run_efast_best_spd.jl`:
- `best_spd_pct = min_t 100*(Btumor(t)/Btumor(0)-1)`, `t ∈ [0, horizon]`.

Regimen used (8-cycle step-up):
- Day 0: 1 mg
- Day 7: 2 mg
- Day 14: 60 mg
- Day 21: 60 mg
- Days 42,63,84,105,126,147: 30 mg
- Horizon: 168 days

Run from repo root (PowerShell equivalent of env + Julia call):

```powershell
$env:EFAST_OUT_DIR = "generated/figures/sensitivity/efast_best_spd"
$env:EFAST_N = "65"
$env:EFAST_M = "4"
$env:EFAST_SAVE_DT = "1.0"
$env:TCE_SOLVER = "qndf"
julia --project=./julia julia/run_efast_best_spd.jl
python scripts/plot_efast_best_spd.py --in-dir generated/figures/sensitivity/efast_best_spd
```

Outputs:
- `generated/figures/sensitivity/efast_best_spd/efast_best_spd_indices.csv`
- `generated/figures/sensitivity/efast_best_spd/efast_best_spd_samples.csv`
- `generated/figures/sensitivity/efast_best_spd/efast_best_spd_meta.json`
- `generated/figures/sensitivity/efast_best_spd/efast_best_spd_indices.png`

## 7) Notes

- The eFAST implementation follows SALib-compatible FAST sampling and analysis equations.
- If you switch to `TCE_SOLVER=cvode_bdf`, runtime may increase but robustness can improve for some parameter corners.
- `assets/musun_2022_88patients_SPD.csv` should be format-checked before using as a hard clinical calibration target.
