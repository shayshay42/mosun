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

OUT_DIR="${OUT_DIR:-$REPO_ROOT/generated/figures/optimization/multicycle_methods_patient1_4cycles_12ctrl_budget}"
mkdir -p "$OUT_DIR"

echo "[run] multicycle all-method benchmark (patient=1, 4 cycles, 12 controls)"
echo "[run] output dir: $OUT_DIR"

env \
  OPTMC_PATIENTS_CSV="assets/generated_vpop/selected_patients.csv" \
  OPTMC_OVERRIDES_CSV="generated/phase1_design/dlbcl_param_overrides.csv" \
  OPTMC_PATIENT_ID="${OPTMC_PATIENT_ID:-1}" \
  OPTMC_N_CYCLES="${OPTMC_N_CYCLES:-4}" \
  OPTMC_SCHEDULE_MODE="${OPTMC_SCHEDULE_MODE:-three_per_cycle}" \
  OPTMC_OPT_METHOD="${OPTMC_OPT_METHOD:-all_methods}" \
  OPTMC_SOLVER="${OPTMC_SOLVER:-tsit5}" \
  OPTMC_ABSTOL="${OPTMC_ABSTOL:-1e-6}" \
  OPTMC_RELTOL="${OPTMC_RELTOL:-1e-4}" \
  OPTMC_MAXITERS_SOLVE="${OPTMC_MAXITERS_SOLVE:-300000}" \
  OPTMC_OUT_DIR="$OUT_DIR" \
  OPTMC_SCALE_SAMPLES="${OPTMC_SCALE_SAMPLES:-4}" \
  OPTMC_SA_ITERS="${OPTMC_SA_ITERS:-10}" \
  OPTMC_NM_ITERS="${OPTMC_NM_ITERS:-10}" \
  OPTMC_CD_ITERS="${OPTMC_CD_ITERS:-5}" \
  OPTMC_LBFGS_ITERS="${OPTMC_LBFGS_ITERS:-10}" \
  OPTMC_SA_F_CALLS="${OPTMC_SA_F_CALLS:-40}" \
  OPTMC_NM_F_CALLS="${OPTMC_NM_F_CALLS:-60}" \
  OPTMC_LBFGS_F_CALLS="${OPTMC_LBFGS_F_CALLS:-80}" \
  OPTMC_RANDOM_SEED="${OPTMC_RANDOM_SEED:-20260303}" \
  "$JULIA_BIN" --project="$JULIA_PROJECT" "$REPO_ROOT/julia/optimize_patient_multicycle_balanced.jl"

echo "[plot] multicycle method comparison"
OPTMC_RESULTS_DIR="$OUT_DIR" \
OPTMC_FIG_OUT_DIR="$OUT_DIR" \
  "$PYTHON_BIN" "$REPO_ROOT/scripts/plot_multicycle_method_comparison.py"

echo "[done] wrote:"
ls -1 "$OUT_DIR"
