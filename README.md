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
