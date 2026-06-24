# Genentech TCE VPop Translation (MATLAB SimBiology -> Julia)

This repository contains a reproducible Julia translation of the Genentech T-cell engager (mosunetuzumab) QSP model, plus workflows for:

- phase-1 dosing scenario simulations,
- MATLAB vs Julia fidelity checks,
- endpoint distribution plotting,
- AD gradient checks (OtD/DtO; forward/reverse),
- dose optimization experiments,
- virtual population (VPop) sampling and pruning.

## Repository goals

1. Preserve model semantics from SimBiology while exposing a canonical Julia ODE RHS.
2. Verify parity against MATLAB/SimBiology reference outputs.
3. Support optimization and VPop generation in a scriptable Julia/Python pipeline.

## Environment Setup

### 1) Python

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -U pip
pip install -r requirements.txt
```

### 2) Julia (repo-local environment)

```bash
julia --project=./julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

The Julia environment is tracked in:

- `julia/Project.toml`
- `julia/Manifest.toml`

## Windows One-Shot Setup

PowerShell:

```powershell
pwsh ./scripts/windows/setup_windows.ps1 -RepoRoot . -PythonExe python -JuliaExe julia
```

## MOSUN Julia VPop Pipeline

The current production path for mosunetuzumab work is `MosunModelCore`, a checked-in Julia ODE implementation designed for validation, sensitivity analysis, and dose/VPop workflows. Generated campaign outputs should be written under ignored `generated/` subdirectories; only lightweight reference manifests and source inputs are tracked.

### Validate the Julia model

```bash
julia --project=./julia julia/test/runtests.jl
julia --project=./julia julia/test_canonical_not_worse.jl
```

For implementation notes and validation layers, see `docs/MOSUN_MODEL_CORE_VALIDATION.md`.

### Run a small eFAST smoke test

This exercises the Fig. 5 eFAST sensitivity runner with a tiny parameter/regimen subset and writes ignored outputs.

```bash
SUSILO_FIG5_EFAST_OUT_DIR=generated/tmp_publish_smoke/susilo_fig5_efast \
SUSILO_FIG5_EFAST_N=65 \
SUSILO_FIG5_EFAST_M=4 \
SUSILO_FIG5_EFAST_PARAM_LIMIT=3 \
SUSILO_FIG5_EFAST_REGIMEN_LIMIT=2 \
julia --project=./julia julia/run_susilo_fig5_il6_dummy_efast.jl
```

PowerShell equivalent:

```powershell
$env:SUSILO_FIG5_EFAST_OUT_DIR = "generated/tmp_publish_smoke/susilo_fig5_efast"
$env:SUSILO_FIG5_EFAST_N = "65"
$env:SUSILO_FIG5_EFAST_M = "4"
$env:SUSILO_FIG5_EFAST_PARAM_LIMIT = "3"
$env:SUSILO_FIG5_EFAST_REGIMEN_LIMIT = "2"
julia --project=./julia julia/run_susilo_fig5_il6_dummy_efast.jl
```

Full eFAST runs use the same entrypoints without the smoke limits. The main Fig. 5 sensitivity scripts are:

- `julia/run_susilo_fig5_il6_dummy_efast.jl`
- `julia/run_susilo_fig5_il6_day84_checkpointed_efast.jl`
- `scripts/plot_susilo_fig5_il6_dummy_efast.py`
- `scripts/plot_susilo_fig5_il6_day84_checkpointed_efast.py`

The default parameter-range input is `generated/figures/reference/dlbcl_stack_parameter_ranges_all/dlbcl_stack_parameter_ranges_all.csv`.

### Generate and prune the Hosseini Fig. 5 VPop

After eFAST candidate/resimulation outputs are available, build the candidate set and prune to a VPop using the digitized Hosseini Fig. 5 targets in `assets/digitization_hosseini_2020_vpop/`:

```bash
python scripts/build_susilo_efast_hosseini2020_vpop250_candidate_set.py
python scripts/prune_susilo_efast_hosseini2020_fig5_vpop250.py
```

The pruning script writes selected-parameter, trajectory, waterfall, diagnostics, and manifest files under `generated/figures/vpop_pruning/...` by default. Those outputs are intentionally ignored unless a specific lightweight manifest is promoted for publication.

## Core Workflows

### Run phase-1 sweep (canonical Julia)

```bash
julia --project=./julia julia/run_phase1_callback_sweep.jl
```

Outputs are written under `generated/phase1_julia/` (or `PHASE1_JULIA_OUT_DIR` if set).

### Plot 4-panel proxy distributions (fixed vs step-up)

```bash
python scripts/plot_phase1_proxy_kde_fixed_stepup.py
```

### Run VPop sampling + pruning

```bash
python scripts/sample_and_prune_vpop.py \
  --run-tag vpop_prune_run1 \
  --n-candidates 72 \
  --n-select 50 \
  --n-random-subsets 1500 \
  --julia-bin "$(which julia)" \
  --julia-project ./julia
```

Windows PowerShell wrapper:

```powershell
pwsh ./scripts/windows/run_vpop.ps1 -RepoRoot . -JuliaExe julia -RunTag vpop_windows_run
```

### Run parity checks / fidelity artifacts

```bash
julia --project=./julia julia/run_three_case_parity.jl
julia --project=./julia julia/test_canonical_not_worse.jl
```

### Run gradient checks

```bash
julia --project=./julia julia/gradient_check_dose3_ad.jl
```

### Run optimization experiments

```bash
julia --project=./julia julia/optimize_cycle1_doses_methods.jl
julia --project=./julia julia/optimize_two_cycle_population.jl
```

## Key Directories

- `assets/`: source papers, SimBiology supplementary MATLAB code, datasets, parameter sheet.
- `julia/src/TCellEngagerQSP.jl`: translated model and solver wiring.
- `julia/*.jl`: simulation, parity, gradient, and optimization entry points.
- `scripts/*.py`: design generation, VPop pipeline, plotting.
- `generated/model_tables/`: exported model tables used by the Julia build path.
- `generated/model_reconstruction_bundle/`: full reconstruction artifacts.

## Technical Overview

For model history, translation decisions, fidelity test context, dosing semantics, AD/optimization setup, and VPop design details, see:

- `docs/TECHNICAL_SUMMARY.md`
- `docs/TRACKING_LOSS_OPTIMIZATION_HANDOFF.md` (branch handoff + resume runbook for tracking-loss optimization)

## GitHub Push Checklist

```bash
git init
git add README.md docs requirements.txt julia scripts assets generated/model_tables generated/vpop_targets
git commit -m "Initial packaged translation + workflows"
git branch -M main
git remote add origin <your-github-repo-url>
git push -u origin main
```

If you already have a git repo, skip `git init` and only add/commit the new files.
