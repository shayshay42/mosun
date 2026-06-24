#!/usr/bin/env bash
set -euo pipefail

DEFAULT_PROJECTS_ROOT="${PROJECTS_ROOT:-/projects/${USER:-${LOGNAME:-mosun}}}"
CAMPAIGN_ROOT="${CAMPAIGN_ROOT:-${DEFAULT_PROJECTS_ROOT}/mosun_susilo_fig5_efast}"
JULIA_VERSION="${JULIA_VERSION:-1.12.5}"
JULIA_BIN="${CAMPAIGN_ROOT}/tools/julia-${JULIA_VERSION}/bin/julia"
REPO_DIR="${CAMPAIGN_ROOT}/repo"
JULIA_PROJECT="${REPO_DIR}/julia"
LOG_ROOT="${LOG_ROOT:-${CAMPAIGN_ROOT}/logs}"
SMOKE_OUT="${SMOKE_OUT:-${REPO_DIR}/generated/figures/sensitivity/susilo_fig5_il6_dummy_efast_smoke_20260502}"
FULL_OUT="${FULL_OUT:-${REPO_DIR}/generated/figures/sensitivity/susilo_fig5_il6_dummy_efast_20260502}"
RUN_SMOKE="${RUN_SMOKE:-1}"
RUN_FULL="${RUN_FULL:-1}"

mkdir -p "${LOG_ROOT}"

export JULIA_DEPOT_PATH="${CAMPAIGN_ROOT}/julia_depot"
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-100}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_PKG_PRECOMPILE_AUTO=0

run_stage() {
  local stage_name="$1"
  local out_dir="$2"
  local n="$3"
  local seeds="$4"
  local log_file="${LOG_ROOT}/susilo_fig5_efast_${stage_name}_$(date +%Y%m%d_%H%M%S).log"

  echo "Running Susilo/Fig.5 eFAST ${stage_name}: N=${n}, seeds=${seeds}, out=${out_dir}"
  SUSILO_FIG5_EFAST_OUT_DIR="${out_dir}" \
  SUSILO_FIG5_EFAST_N="${n}" \
  SUSILO_FIG5_EFAST_M="${SUSILO_FIG5_EFAST_M:-4}" \
  SUSILO_FIG5_EFAST_SEEDS="${seeds}" \
  SUSILO_FIG5_EFAST_SOLVER="${SUSILO_FIG5_EFAST_SOLVER:-rodas4p}" \
  SUSILO_FIG5_EFAST_HORIZON_DAYS="${SUSILO_FIG5_EFAST_HORIZON_DAYS:-84}" \
  SUSILO_FIG5_EFAST_SAVEAT_DT="${SUSILO_FIG5_EFAST_SAVEAT_DT:-0.125}" \
  SUSILO_FIG5_EFAST_POST_EVENT_DT="${SUSILO_FIG5_EFAST_POST_EVENT_DT:-0.01}" \
  SUSILO_FIG5_EFAST_BW_KG="${SUSILO_FIG5_EFAST_BW_KG:-70}" \
  "${JULIA_BIN}" --project="${JULIA_PROJECT}" "${REPO_DIR}/julia/run_susilo_fig5_il6_dummy_efast.jl" \
    > "${log_file}" 2>&1

  echo "Plotting Susilo/Fig.5 eFAST ${stage_name}"
  if python3 -c "import matplotlib, pandas, numpy" >/dev/null 2>&1; then
    python3 "${REPO_DIR}/scripts/plot_susilo_fig5_il6_dummy_efast.py" --in-dir "${out_dir}" \
      >> "${log_file}" 2>&1
  else
    echo "Skipping server-side plotting for ${stage_name}: matplotlib/pandas/numpy not available in python3." \
      >> "${log_file}" 2>&1
  fi
  echo "Completed ${stage_name}; log=${log_file}"
}

if [[ "${RUN_SMOKE}" == "1" ]]; then
  run_stage "smoke" "${SMOKE_OUT}" "${SUSILO_FIG5_EFAST_SMOKE_N:-256}" "${SUSILO_FIG5_EFAST_SMOKE_SEEDS:-20260502}"
fi

if [[ "${RUN_FULL}" == "1" ]]; then
  run_stage "full" "${FULL_OUT}" "${SUSILO_FIG5_EFAST_FULL_N:-4096}" "${SUSILO_FIG5_EFAST_FULL_SEEDS:-20260502,20260503,20260504}"
fi
