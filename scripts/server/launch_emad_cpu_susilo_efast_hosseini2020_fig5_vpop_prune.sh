#!/usr/bin/env bash
set -euo pipefail

DEFAULT_PROJECTS_ROOT="${PROJECTS_ROOT:-/projects/${USER:-${LOGNAME:-mosun}}}"
CAMPAIGN_ROOT="${CAMPAIGN_ROOT:-${DEFAULT_PROJECTS_ROOT}/mosun_susilo_fig5_efast}"
REPO_DIR="${CAMPAIGN_ROOT}/repo"
LOG_ROOT="${LOG_ROOT:-${CAMPAIGN_ROOT}/logs}"
OUT_DIR="${SUSILO_EFAST_HOSSEINI_VPOP_OUT_DIR:-${REPO_DIR}/generated/figures/vpop_pruning/susilo_efast_hosseini2020_fig5_vpop250_20260505}"
NICE_LEVEL="${NICE_LEVEL:-5}"
PYTHON_BIN="${PYTHON_BIN:-${CAMPAIGN_ROOT}/pyenvs/vpop-prune/bin/python}"

mkdir -p "${LOG_ROOT}"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export SUSILO_EFAST_HOSSEINI_VPOP_OUT_DIR="${OUT_DIR}"
export SUSILO_EFAST_VPOP_PRUNE_WORKERS="${SUSILO_EFAST_VPOP_PRUNE_WORKERS:-24}"
export SUSILO_EFAST_VPOP_PRUNE_INITIAL_CANDIDATES="${SUSILO_EFAST_VPOP_PRUNE_INITIAL_CANDIDATES:-240}"
export SUSILO_EFAST_VPOP_PRUNE_STARTS="${SUSILO_EFAST_VPOP_PRUNE_STARTS:-6}"
export SUSILO_EFAST_VPOP_PRUNE_SWAP_ITERATIONS="${SUSILO_EFAST_VPOP_PRUNE_SWAP_ITERATIONS:-900}"
export SUSILO_EFAST_VPOP_RANDOM_BASELINES="${SUSILO_EFAST_VPOP_RANDOM_BASELINES:-32}"

log_file="${LOG_ROOT}/susilo_efast_hosseini2020_fig5_vpop_prune_$(date +%Y%m%d_%H%M%S).log"

echo "Running eFAST-derived Hosseini Fig.5 VPop250 pruning: out=${OUT_DIR}, workers=${SUSILO_EFAST_VPOP_PRUNE_WORKERS}, starts=${SUSILO_EFAST_VPOP_PRUNE_STARTS}, swaps=${SUSILO_EFAST_VPOP_PRUNE_SWAP_ITERATIONS}, log=${log_file}"

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "Python environment not found at ${PYTHON_BIN}; cannot run cluster-side pruning." | tee "${log_file}"
  exit 2
fi

if ! "${PYTHON_BIN}" -c "import matplotlib, numpy, pandas" >/dev/null 2>&1; then
  echo "${PYTHON_BIN} is missing matplotlib/numpy/pandas; cannot run cluster-side pruning." | tee "${log_file}"
  exit 2
fi

cd "${REPO_DIR}"
nice -n "${NICE_LEVEL}" "${PYTHON_BIN}" "${REPO_DIR}/scripts/prune_susilo_efast_hosseini2020_fig5_vpop250.py" \
  > "${log_file}" 2>&1

echo "Completed eFAST-derived Hosseini Fig.5 VPop250 pruning; log=${log_file}"
