# Technical Summary: Genentech TCE QSP Translation and Workflows

This document summarizes the technical architecture, design decisions, and execution workflows in this repository.

## 1) Model Lineage and Scope

## Source model

- Original model context: SimBiology project from Genentech supplementary code (`assets/Supp Matlab Code/TDBr26_6_paper.sbproj`) with supporting MATLAB scripts and datasets.
- Source publication context:
  - `assets/s41540-020-00145-7.pdf`
  - `assets/iraj_tce_supp.pdf`
- Follow-up context:
  - `assets/followup_paper/Systems-based_Digital_Twins_to_Help_Characterize_C.pdf`
  - supplementary tables/docs in `assets/followup_paper/`.

## Translation objective

- Build a Julia implementation that is semantically faithful to the SimBiology model.
- Preserve equations, repeated-assignment rules, variant semantics, and dose event behavior.
- Support numerical simulation parity, optimization, and VPop workflows in open scripting tools.

## Export/reconstruction artifacts

- Flattened model tables: `generated/model_tables/`
  - species, parameters, rules, reactions, variants, doses, stoichiometry.
- Full reconstruction bundle: `generated/model_reconstruction_bundle/`
  - includes SBML and equation manifests.
- MATLAB generate-code bundle: `generated/model_generatecode_bundle_2025b/`.

## 2) Julia Model Architecture

Primary implementation: `julia/src/TCellEngagerQSP.jl`

## Core pieces

1. Model loader/parsers
- Reads exported TSV/CSV artifacts from `generated/model_tables`.
- Constructs state, parameter, rule, reaction, and dose structures.

2. Two RHS modes
- `legacy` mode: evaluates rule/reaction functions similar to exported dynamic semantics.
- `canonical` mode: compiles repeated assignments + reaction rates into a canonical Julia RHS closure.

3. Simulation context
- `SimContext` for legacy path.
- `CanonicalSimContext` for canonical path.

4. Solver abstraction
- Supported: `cvode_bdf`, `qndf`, `rodas4p`, `kencarp4`, `vcabm`, `tsit5`.
- Optional MTK Jacobian path for stiff solvers.

5. Dosing semantics
- Predominantly event-based dosing via preset-time callback in sweep scripts.
- Also includes segmented solve approach for canonical callback-vs-segmented comparisons.

## Important design note

The ODE core is canonical, but some observable evaluation still uses repeated-assignment rules in evaluation helpers to ensure exact SimBiology observable semantics for outputs such as `IL6combo`.

## 3) Dosing Scenarios and Simulation Entry Points

## Phase-1 scenario simulation

Script: `julia/run_phase1_callback_sweep.jl`

Inputs:
- `generated/phase1_design*/patients.csv`
- `generated/phase1_design*/regimen_events.csv`
- model tables under `generated/model_tables/`

Outputs (`phase1_metrics_julia.csv`):
- `il6_peak_0_2`
- `il6_peak_0_21`
- `il6_auc_0_2`
- `il6_auc_0_21`
- `btumor_t0`, `btumor_day42`, `tumor_resid_day42`, `tumor_cfbl_day42`
- `loss_raw`, plus metadata (`regimen`, `patient_id`, `status`)

## Standard phase-1 design generation

- `scripts/generate_phase1_standard_design.py`
- `scripts/generate_phase1_design.py`

These define fixed-dose and step-up regimen event schedules used by simulation workflows.

## 4) Proxies, Endpoints, and Distribution Plots

## Endpoint proxies used in workflows

- Toxicity proxy:
  - Early IL6 peak windows (`il6_peak_0_2` and now `il6_peak_0_21`) and AUC windows.
- Tumor burden proxy:
  - `tumor_resid_day42 = Btumor(day42) / Btumor(day0)`.

## Plotting scripts

- `scripts/plot_phase1_proxy_kde_fixed_stepup.py`
  - 4-panel fixed-vs-step-up overlays for scale-adjusted toxicity/tumor terms.
- `scripts/plot_phase1_standard_trace_overlay.py`
  - MATLAB/Julia overlays.
- `scripts/plot_phase1_tumor_trace_overlay_threeway.py`
  - MATLAB + Julia legacy + Julia canonical overlays.
- `scripts/plot_standard_endpoint_distributions.py`
  - endpoint distribution overlays.

## 5) Fidelity/Parity Testing History

Parity outputs are stored under `generated/`.

## Representative parity artifacts

1. Three-case MATLAB vs canonical Julia parity
- File: `generated/julia_reference/parity_summary_canonical.csv`
- Representative max relative errors are small (on the order of 0.25% to 1.6%).

2. Standard all-regimen endpoint parity
- File: `generated/phase1_standard_julia/parity_summary_standard_all_regimens.csv`
- Includes endpoint-level absolute/relative summaries.

3. Trace-level parity after variant semantic fix
- File: `generated/phase1_standard_julia_traces/parity_trace_summary_after_variant_semantics_fix.csv`

4. Canonical segmented-vs-callback internal consistency
- File: `generated/phase1_standard_julia_traces_canonical/canonical_segmented_vs_callback_summary.csv`
- Shows very small normalized discrepancies, confirming dosing-mode equivalence in canonical Julia path.

## Regression guard

- `julia/test_canonical_not_worse.jl`
- Intended to ensure canonical translation quality does not regress versus established legacy references.

## 6) AD Gradient and Optimization Workflows

## AD gradient checks

Script: `julia/gradient_check_dose3_ad.jl`

Compares dose gradients for a cycle-1 objective using:
- OtD forward
- OtD reverse
- DtO forward
- DtO reverse
- finite-difference references

Purpose:
- verify gradient consistency and conditioning before optimization.

## Optimization scripts

1. Single-patient cycle-1 optimization method comparison
- `julia/optimize_cycle1_doses_methods.jl`
- Supports methods including:
  - simulated annealing,
  - coordinate-descent/Nelder-Mead style search,
  - finite-diff + LBFGS,
  - forward AD + LBFGS (DtO and OtD variants).

2. Two-cycle population optimization
- `julia/optimize_two_cycle_population.jl`
- Produces cohort-level fixed-vs-optimized endpoint and dose distribution outputs.

Plot scripts:
- `scripts/plot_cycle1_optimization_results.py`
- `scripts/plot_two_cycle_population_optimization.py`
- `scripts/plot_two_cycle_optimized_dose_heatmap.py`

## 7) VPop Generation and Pruning Design

Main script: `scripts/sample_and_prune_vpop.py`

## Current (v2) features

1. Candidate VP sampling
- Samples core tumor parameters and additional kinetic parameters from bounds in:
  - `generated/vpop_targets/parameter_prior_bounds_from_params_sheet.csv`

2. Multi-regimen objective support
- Plan can be expanded to include standard phase-1 regimens.
- Supports target modes:
  - `population_fraction`
  - `population_mean`
  - `pairwise_difference`

3. Toxicity thresholds by feature
- CRS proxy thresholds are calibrated from reference regimen distributions and can use day0-21 IL6 features.

4. Shape penalties
- Optional monotonic/order penalties across regimen ladders (within regimen type) for response and toxicity trends.

## Target bundle files

- `generated/vpop_targets/vpop_calibration_targets.json`
- `generated/vpop_targets/vpop_target_evaluation_plan.json`
- `generated/vpop_targets/vpop_calibration_targets_flat.csv`

## Known current limitations

Some targets remain marked as not directly evaluable in current outputs, e.g.:
- explicit PK proportionality metrics,
- SPD-distribution mapping,
- long-term outcomes (PFS/OS/DOR),
- baseline residual rituximab detectability.

These are explicitly tracked in the evaluation plan summary fields.

## 8) Windows Execution Runbook

## One-time setup (PowerShell)

```powershell
pwsh ./scripts/windows/setup_windows.ps1 -RepoRoot . -PythonExe python -JuliaExe julia
```

This does:
- creates `.venv` and installs `requirements.txt`,
- instantiates/precompiles `julia/Project.toml`.

## Run VPop generation

```powershell
pwsh ./scripts/windows/run_vpop.ps1 `
  -RepoRoot . `
  -JuliaExe julia `
  -RunTag vpop_windows_run `
  -NCandidates 72 `
  -NSelect 50 `
  -NRandomSubsets 1500
```

Output summary:
- `generated/vpop_pruning/<run_tag>/pruning_summary.json`

## 9) Design Decisions (Why this shape)

1. Keep translation lossless first, then optimize structure
- Canonical RHS built from exported equations,
- repeated-assignment observables preserved to match SimBiology semantics.

2. Separate simulation from selection
- expensive ODE simulation writes reusable metrics,
- pruning can iterate quickly with `--skip-sim` on existing metrics.

3. Keep dosing explicit and auditable
- regimen events are CSV-defined, not hardcoded in solver internals.

4. Keep artifacts inspectable
- all parity summaries, metrics, and objective traces are written to files for auditability.

## 10) Suggested GitHub Layout for Team Use

Track:
- source scripts (`julia/`, `scripts/`),
- environment files (`julia/Project.toml`, `julia/Manifest.toml`, `requirements.txt`),
- target bundles (`generated/vpop_targets/`),
- model export tables required by Julia runtime (`generated/model_tables/`).

Avoid committing heavy runtime outputs repeatedly:
- `generated/vpop_pruning/`, `generated/figures/`, and temporary run folders (covered by `.gitignore`).

## 11) Quick Command Index

Run canonical phase-1 sweep:

```bash
julia --project=./julia julia/run_phase1_callback_sweep.jl
```

Run VPop pruning:

```bash
python scripts/sample_and_prune_vpop.py --run-tag vpop_run --julia-project ./julia --julia-bin "$(which julia)"
```

Plot 4-panel scaled proxy distributions:

```bash
python scripts/plot_phase1_proxy_kde_fixed_stepup.py
```

Run gradient checks:

```bash
julia --project=./julia julia/gradient_check_dose3_ad.jl
```

Run cycle-1 optimization comparison:

```bash
julia --project=./julia julia/optimize_cycle1_doses_methods.jl
```
