#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -n "${JULIA_BIN:-}" ]]; then
  JULIA_BIN="$JULIA_BIN"
elif command -v julia >/dev/null 2>&1; then
  JULIA_BIN="$(command -v julia)"
elif [[ -x "$HOME/.juliaup/bin/julia" ]]; then
  JULIA_BIN="$HOME/.juliaup/bin/julia"
else
  JULIA_BIN="julia"
fi
JULIA_PROJECT="${JULIA_PROJECT:-$REPO_ROOT/julia}"
if [[ -n "${PYTHON_BIN:-}" ]]; then
  PYTHON_BIN="$PYTHON_BIN"
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN="$(command -v python3)"
else
  PYTHON_BIN="python3"
fi

OUT_DIR="${OUT_DIR:-$REPO_ROOT/generated/figures/optimization/vpop_tracking_nmcd}"
mkdir -p "$OUT_DIR"

echo "[run] VPop NM+CD optimization on tracking loss"
echo "[run] output dir: $OUT_DIR"

env \
  OPTMC_PATIENTS_CSV="${OPTMC_PATIENTS_CSV:-assets/generated_vpop/selected_patients.csv}" \
  OPTMC_OVERRIDES_CSV="${OPTMC_OVERRIDES_CSV:-generated/phase1_design/dlbcl_param_overrides.csv}" \
  OPTMC_N_CYCLES="${OPTMC_N_CYCLES:-2}" \
  OPTMC_SCHEDULE_MODE="${OPTMC_SCHEDULE_MODE:-three_per_cycle}" \
  OPTMC_SOLVER="${OPTMC_SOLVER:-tsit5}" \
  OPTMC_ABSTOL="${OPTMC_ABSTOL:-1e-6}" \
  OPTMC_RELTOL="${OPTMC_RELTOL:-1e-4}" \
  OPTMC_MAXITERS_SOLVE="${OPTMC_MAXITERS_SOLVE:-300000}" \
  OPTMC_TRACK_TUMOR_END_FRAC="${OPTMC_TRACK_TUMOR_END_FRAC:-0.05}" \
  OPTMC_TRACK_IL6_BASE="${OPTMC_TRACK_IL6_BASE:-0.0}" \
  OPTMC_TRACK_IL6_PULSE="${OPTMC_TRACK_IL6_PULSE:-120.0}" \
  OPTMC_TRACK_IL6_SIGMA_DAYS="${OPTMC_TRACK_IL6_SIGMA_DAYS:-0.40}" \
  OPTMC_TRACK_IL6_SCALE="${OPTMC_TRACK_IL6_SCALE:-300.0}" \
  OPTMC_TRACK_W_TUMOR="${OPTMC_TRACK_W_TUMOR:-0.70}" \
  OPTMC_TRACK_W_TOX="${OPTMC_TRACK_W_TOX:-0.25}" \
  OPTMC_TRACK_W_DOSE="${OPTMC_TRACK_W_DOSE:-0.05}" \
  OPTVPOP_MAX_PATIENTS="${OPTVPOP_MAX_PATIENTS:-999999}" \
  OPTVPOP_NM_ITERS="${OPTVPOP_NM_ITERS:-20}" \
  OPTVPOP_NM_F_CALLS="${OPTVPOP_NM_F_CALLS:-250}" \
  OPTVPOP_CD_SWEEPS="${OPTVPOP_CD_SWEEPS:-10}" \
  OPTVPOP_CD_STEP0_MG="${OPTVPOP_CD_STEP0_MG:-2.0}" \
  OPTVPOP_CD_TOL_MG="${OPTVPOP_CD_TOL_MG:-0.05}" \
  OPTVPOP_OUT_DIR="$OUT_DIR" \
  "$JULIA_BIN" --project="$JULIA_PROJECT" "$REPO_ROOT/julia/optimize_tracking_loss_vpop_nmcd.jl"

echo "[post] export trajectories (fixed + optimized) for spaghetti plot"
env \
  OPTMC_PATIENTS_CSV="${OPTMC_PATIENTS_CSV:-assets/generated_vpop/selected_patients.csv}" \
  OPTMC_OVERRIDES_CSV="${OPTMC_OVERRIDES_CSV:-generated/phase1_design/dlbcl_param_overrides.csv}" \
  OPTMC_N_CYCLES="${OPTMC_N_CYCLES:-2}" \
  OPTMC_SCHEDULE_MODE="${OPTMC_SCHEDULE_MODE:-three_per_cycle}" \
  OPTMC_SOLVER="${OPTMC_SOLVER:-tsit5}" \
  OPTMC_ABSTOL="${OPTMC_ABSTOL:-1e-6}" \
  OPTMC_RELTOL="${OPTMC_RELTOL:-1e-4}" \
  OPTMC_MAXITERS_SOLVE="${OPTMC_MAXITERS_SOLVE:-300000}" \
  OPTMC_TRACK_TUMOR_END_FRAC="${OPTMC_TRACK_TUMOR_END_FRAC:-0.05}" \
  OPTMC_TRACK_IL6_BASE="${OPTMC_TRACK_IL6_BASE:-0.0}" \
  OPTMC_TRACK_IL6_PULSE="${OPTMC_TRACK_IL6_PULSE:-120.0}" \
  OPTMC_TRACK_IL6_SIGMA_DAYS="${OPTMC_TRACK_IL6_SIGMA_DAYS:-0.40}" \
  OPTMC_TRACK_IL6_SCALE="${OPTMC_TRACK_IL6_SCALE:-300.0}" \
  OPTMC_TRACK_W_TUMOR="${OPTMC_TRACK_W_TUMOR:-0.70}" \
  OPTMC_TRACK_W_TOX="${OPTMC_TRACK_W_TOX:-0.25}" \
  OPTMC_TRACK_W_DOSE="${OPTMC_TRACK_W_DOSE:-0.05}" \
  OPTVPOP_SUMMARY_CSV="$OUT_DIR/vpop_tracking_nmcd_summary.csv" \
  OPTVPOP_OUT_DIR="$OUT_DIR" \
  "$JULIA_BIN" --project="$JULIA_PROJECT" "$REPO_ROOT/julia/export_vpop_tracking_nmcd_trajectories.jl"

echo "[plot] spaghetti trajectories"
OPTVPOP_RESULTS_DIR="$OUT_DIR" \
OPTVPOP_FIG_OUT_DIR="$OUT_DIR" \
  "$PYTHON_BIN" "$REPO_ROOT/scripts/plot_vpop_tracking_nmcd_spaghetti.py"

echo "[done] wrote:"
ls -1 "$OUT_DIR"
