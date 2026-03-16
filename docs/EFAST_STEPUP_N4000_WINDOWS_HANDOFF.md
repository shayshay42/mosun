# eFAST Step-Up `N=4000` Windows Handoff

This note is for handing the step-up DLBCL eFAST runs to a collaborator with a larger Windows workstation.

## Purpose

Run the same two scalar-output eFAST analyses we already ran locally, but at a larger sample size:

- `AUC(Btumor)`
- `AUC(IL6combo)`

using the existing step-up regimen and the same `56` DLBCL-range parameters.

## Why `N=4000`

The current eFAST runner is [run_efast_aucsum_dlbcl_ranges.jl](/Users/shayanhajhashemi/genentech_tce_vpop_translation/julia/run_efast_aucsum_dlbcl_ranges.jl).

Its FAST sampler uses:

- `D = 56` sampled parameters
- `M = 4`
- primary frequency `omega0 = floor((N - 1) / (2M))`
- secondary frequency budget `m = floor(omega0 / (2M))`

The formal minimum validity condition in the code is:

- `N > 4*M^2`

which is only `N > 64` for `M=4`.

That is not a convergence threshold. For this implementation, a more meaningful structural threshold is the point where the secondary frequency budget can cover the other `D-1` parameters without reusing frequencies:

- require `m >= D - 1`
- approximately `N >= 4*M^2*(D-1) + 1`
- here: `N >= 3521`

So `N=4000` is the first clean target above that threshold.

At `N=4000`:

- total evaluations per endpoint = `56 * 4000 = 224000`
- two endpoints = `448000` solves total

## Regimen

This uses the same 8-cycle step-up regimen as the earlier local runs:

- dose days: `0, 7, 14, 21, 42, 63, 84, 105, 126, 147`
- doses mg: `1, 2, 60, 60, 30, 30, 30, 30, 30, 30`
- horizon: `168` days
- solver: `QNDF`

## Required Files On The Handoff Branch

Do not push only the new Windows wrapper. The collaborator needs the runnable Julia entrypoint and its support file too.

Minimum set that must be present on the branch:

- [run_efast_aucsum_dlbcl_ranges.jl](/Users/shayanhajhashemi/genentech_tce_vpop_translation/julia/run_efast_aucsum_dlbcl_ranges.jl)
- [MosunModelCoreSupport.jl](/Users/shayanhajhashemi/genentech_tce_vpop_translation/julia/src/MosunModelCoreSupport.jl)
- [run_efast_stepup_n4000.ps1](/Users/shayanhajhashemi/genentech_tce_vpop_translation/scripts/windows/run_efast_stepup_n4000.ps1)
- [plot_efast_best_spd.py](/Users/shayanhajhashemi/genentech_tce_vpop_translation/scripts/plot_efast_best_spd.py)
- `generated/figures/reference/dlbcl_stack_parameter_ranges_all/dlbcl_stack_parameter_ranges_all.csv`
- tracked model tables in `generated/model_reconstruction_bundle/`

Without the bounds CSV or the untracked Julia runner/support file, the Windows wrapper is not enough.

## Recommended Git Workflow

Yes. Push this as a new branch.

Reason:

- this is a compute handoff, not a stable `main` workflow yet
- the eFAST runner and support code are still local-only in this worktree
- the collaborator should get an isolated, reproducible branch containing the exact runner, bounds file, and launch script used for the handoff

Suggested branch name:

- `efast_stepup_n4000_handoff`

## Windows Setup

From the repo root in PowerShell:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/windows/setup_windows.ps1 -RepoRoot .
```

That creates the Python venv and instantiates the Julia project.

## Run Command

Default run:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/windows/run_efast_stepup_n4000.ps1 -RepoRoot .
```

This launches both endpoints in parallel by default.

Useful overrides:

```powershell
# Example: 32 logical cores, leave headroom, use 10 threads per Julia job
powershell -ExecutionPolicy Bypass -File scripts/windows/run_efast_stepup_n4000.ps1 `
  -RepoRoot . `
  -JuliaExe julia `
  -ThreadsPerJob 10
```

```powershell
# Run only the tumor endpoint
powershell -ExecutionPolicy Bypass -File scripts/windows/run_efast_stepup_n4000.ps1 `
  -RepoRoot . `
  -Only auc_btumor `
  -ThreadsPerJob 12
```

```powershell
# Run sequentially if the machine or scheduler dislikes two simultaneous Julia processes
powershell -ExecutionPolicy Bypass -File scripts/windows/run_efast_stepup_n4000.ps1 `
  -RepoRoot . `
  -Sequential `
  -ThreadsPerJob 20
```

```powershell
# Skip plotting if Python is unavailable
powershell -ExecutionPolicy Bypass -File scripts/windows/run_efast_stepup_n4000.ps1 `
  -RepoRoot . `
  -SkipPlots
```

## Output Locations

The script writes to:

- `generated/figures/sensitivity/auc_btumor_dlbcl_ranges_n4000/`
- `generated/figures/sensitivity/auc_il6combo_dlbcl_ranges_n4000/`

Each directory should contain:

- `*_indices.csv`
- `*_samples.csv`
- `*_meta.json`
- `*_indices.png` if plotting is enabled
- `run.log`

## Threading Guidance

The launcher defaults to:

- `ThreadsPerJob = floor((logical_cpu_count - 2) / 2)`

so that two Julia processes can run in parallel while leaving some headroom.

Examples:

- `16` logical cores: use `ThreadsPerJob 6` or `7`
- `24` logical cores: use `ThreadsPerJob 10` or `11`
- `32` logical cores: use `ThreadsPerJob 10` to `14`

If memory is tight or Julia starts swapping, reduce `ThreadsPerJob` before changing anything else.

## Current Stability Context

Earlier local runs showed:

- `N=129 -> N=1000` changed rankings materially
- `N=1000 -> N=1500` still changed rankings materially

So `N=4000` is not “guaranteed converged”. It is the first target that gets out of the obvious secondary-frequency reuse regime in the current eFAST implementation.

The goal of the Windows handoff is:

- get a much better-quality ranking estimate
- then compare `N=1500` vs `N=4000` before interpreting the result as stable

## What To Report Back

At minimum, ask the collaborator to send back:

- both `*_indices.csv`
- both `*_meta.json`
- both `*_indices.png`

If possible, also keep:

- `run.log` from each directory

Those files are enough to compare `N=1500` against `N=4000` locally.
