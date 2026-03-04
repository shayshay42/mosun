# Tracking-Loss Optimization Handoff

This note is the continuation runbook for branch `tracking_loss_optimization`.

## Scope of this branch

This branch packages the tracking-loss optimization workflow for multi-cycle dose control in the Julia-translated TCE model.

Main additions used for this workflow:

- `julia/optimize_patient_multicycle_tracking.jl`
- `julia/optimize_tracking_loss_vpop_nmcd.jl`
- `julia/export_vpop_tracking_nmcd_trajectories.jl`
- `julia/eval_tracking_loss_vpop_scenarios.jl`
- `scripts/run_vpop_tracking_nmcd.sh`
- `scripts/plot_vpop_tracking_nmcd_spaghetti.py`
- `scripts/run_vpop_tracking_loss_scenarios.sh`
- `scripts/plot_vpop_tracking_loss_scenarios.py`
- `scripts/run_multicycle_tracking_methods_patient1_2cycles.sh`
- `scripts/plot_multicycle_tracking_trajectories.py`
- `assets/generated_vpop/selected_patients.csv`
- `assets/generated_vpop/pruning_summary.json`

## Important implementation notes

1. `scripts/run_*` launchers are now machine-portable by default:
   - `JULIA_BIN` defaults to `julia`
   - `JULIA_PROJECT` defaults to `<repo>/julia`
   - `PYTHON_BIN` defaults to `python3`

2. Resume behavior for VPop NM+CD was fixed:
   - `julia/optimize_tracking_loss_vpop_nmcd.jl` now correctly detects `patient_id` in existing summary CSVs and skips completed patients when `OPTVPOP_RESUME=1`.
   - Existing duplicate rows are compacted to the latest row per patient at resume load.

3. Julia entry-point files in this workflow now use guarded `main()` calls (`if abspath(PROGRAM_FILE) == @__FILE__`).

## Quick setup on a new machine

From repo root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -U pip
python3 -m pip install -r requirements.txt

julia --project=./julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

Optional override (if you want a different Julia environment):

```bash
export JULIA_PROJECT=/path/to/your/julia_env
```

## Main run command (full VPop, NM+CD)

This is the high-budget command we discussed (400 NM iterations):

```bash
OPTVPOP_NM_ITERS=400 \
OPTVPOP_NM_F_CALLS=20000 \
OPTVPOP_CD_SWEEPS=10 \
OPTVPOP_RESUME=1 \
./scripts/run_vpop_tracking_nmcd.sh
```

What it writes:

- `vpop_tracking_nmcd_summary.csv`
- `vpop_tracking_nmcd_meta.json`
- `vpop_tracking_nmcd_trajectories.csv`
- `vpop_tracking_nmcd_trajectories_meta.json`
- `vpop_tracking_nmcd_spaghetti.png`

Default output directory:

- `generated/figures/optimization/vpop_tracking_nmcd/`

## Resume / checkpoint controls

You can tune checkpointing and progress logs without editing code:

- `OPTVPOP_RESUME=1` (default): resume from existing summary
- `OPTVPOP_PROGRESS_EVERY=10`: progress print frequency
- `OPTVPOP_CHECKPOINT_EVERY=5`: write summary every N newly completed patients
- `OPTVPOP_MAX_PATIENTS=<n>`: subset runs for calibration/debug

Example:

```bash
OPTVPOP_RESUME=1 \
OPTVPOP_PROGRESS_EVERY=1 \
OPTVPOP_CHECKPOINT_EVERY=1 \
OPTVPOP_MAX_PATIENTS=20 \
./scripts/run_vpop_tracking_nmcd.sh
```

## Where we left off

- The resume fix is in place and verified to skip completed patients.
- A calibration run with `NM_ITERS=400` measured about **270 seconds (4.5 min) per patient** for patient 1.
- No completed 140-patient, 400-iter NM+CD output is committed yet.

## Push this branch

```bash
git checkout tracking_loss_optimization
git add -A
git commit -m "Add tracking-loss optimization pipeline, resume-safe VPop NM+CD, and handoff runbook"
git push -u origin tracking_loss_optimization
```

If your remote is not set:

```bash
git remote add origin <your-repo-url>
git push -u origin tracking_loss_optimization
```
