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

OUT_DIR="${OUT_DIR:-$REPO_ROOT/generated/figures/optimization/multicycle_tracking_patient1_2cycles}"
mkdir -p "$OUT_DIR"

echo "[run] quadratic trajectory-tracking objective benchmark (patient=1, 2 cycles, 6 controls)"
echo "[run] output dir: $OUT_DIR"

env \
  OPTMC_PATIENTS_CSV="assets/generated_vpop/selected_patients.csv" \
  OPTMC_OVERRIDES_CSV="generated/phase1_design/dlbcl_param_overrides.csv" \
  OPTMC_PATIENT_ID="${OPTMC_PATIENT_ID:-1}" \
  OPTMC_N_CYCLES="${OPTMC_N_CYCLES:-2}" \
  OPTMC_SCHEDULE_MODE="${OPTMC_SCHEDULE_MODE:-three_per_cycle}" \
  OPTMC_SOLVER="${OPTMC_SOLVER:-tsit5}" \
  OPTMC_ABSTOL="${OPTMC_ABSTOL:-1e-6}" \
  OPTMC_RELTOL="${OPTMC_RELTOL:-1e-4}" \
  OPTMC_MAXITERS_SOLVE="${OPTMC_MAXITERS_SOLVE:-300000}" \
  OPTMC_OUT_DIR="$OUT_DIR" \
  OPTMC_SCALE_SAMPLES="${OPTMC_SCALE_SAMPLES:-4}" \
  OPTMC_SCALE_SEED="${OPTMC_SCALE_SEED:-20260303}" \
  OPTMC_SA_ITERS="${OPTMC_SA_ITERS:-120}" \
  OPTMC_NM_ITERS="${OPTMC_NM_ITERS:-40}" \
  OPTMC_CD_ITERS="${OPTMC_CD_ITERS:-20}" \
  OPTMC_LBFGS_ITERS="${OPTMC_LBFGS_ITERS:-30}" \
  OPTMC_SA_F_CALLS="${OPTMC_SA_F_CALLS:-500}" \
  OPTMC_NM_F_CALLS="${OPTMC_NM_F_CALLS:-400}" \
  OPTMC_LBFGS_F_CALLS="${OPTMC_LBFGS_F_CALLS:-600}" \
  OPTMC_RANDOM_SEED="${OPTMC_RANDOM_SEED:-20260303}" \
  OPTMC_TRACK_TUMOR_END_FRAC="${OPTMC_TRACK_TUMOR_END_FRAC:-0.05}" \
  OPTMC_TRACK_IL6_BASE="${OPTMC_TRACK_IL6_BASE:-0.0}" \
  OPTMC_TRACK_IL6_PULSE="${OPTMC_TRACK_IL6_PULSE:-120.0}" \
  OPTMC_TRACK_IL6_SIGMA_DAYS="${OPTMC_TRACK_IL6_SIGMA_DAYS:-0.40}" \
  OPTMC_TRACK_IL6_SCALE="${OPTMC_TRACK_IL6_SCALE:-300.0}" \
  OPTMC_TRACK_W_TUMOR="${OPTMC_TRACK_W_TUMOR:-0.70}" \
  OPTMC_TRACK_W_TOX="${OPTMC_TRACK_W_TOX:-0.25}" \
  OPTMC_TRACK_W_DOSE="${OPTMC_TRACK_W_DOSE:-0.05}" \
  "$JULIA_BIN" --project="$JULIA_PROJECT" "$REPO_ROOT/julia/optimize_patient_multicycle_tracking.jl"

echo "[plot] method convergence"
OPTMC_RESULTS_DIR="$OUT_DIR" \
OPTMC_FIG_OUT_DIR="$OUT_DIR" \
OPTMC_AD_XJITTER_FRAC="${OPTMC_AD_XJITTER_FRAC:-0.0015}" \
OPTMC_FIG_BASENAME="multicycle_tracking_method_convergence.png" \
  "$PYTHON_BIN" "$REPO_ROOT/scripts/plot_multicycle_method_comparison.py"

echo "[plot] target vs optimized trajectories"
OPTMC_RESULTS_DIR="$OUT_DIR" \
OPTMC_FIG_OUT_DIR="$OUT_DIR" \
OPTMC_TRACK_FIG_BASENAME="multicycle_tracking_target_vs_best.png" \
  "$PYTHON_BIN" "$REPO_ROOT/scripts/plot_multicycle_tracking_trajectories.py"

echo "[done] wrote:"
ls -1 "$OUT_DIR"
