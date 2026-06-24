#!/usr/bin/env bash
set -euo pipefail

DEFAULT_PROJECTS_ROOT="${PROJECTS_ROOT:-/projects/${USER:-${LOGNAME:-mosun}}}"
CAMPAIGN_ROOT="${CAMPAIGN_ROOT:-${DEFAULT_PROJECTS_ROOT}/mosun_campaign}"
JULIA_VERSION="${JULIA_VERSION:-1.12.5}"
JULIA_SERIES="$(printf '%s' "$JULIA_VERSION" | awk -F. '{print $1 "." $2}')"
JULIA_DIR="${CAMPAIGN_ROOT}/tools/julia-${JULIA_VERSION}"
JULIA_TARBALL="julia-${JULIA_VERSION}-linux-x86_64.tar.gz"
JULIA_URL="${JULIA_URL:-https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_SERIES}/${JULIA_TARBALL}}"
REPO_DIR="${CAMPAIGN_ROOT}/repo"
JULIA_PROJECT="${REPO_DIR}/julia"
DEPOT_DIR="${CAMPAIGN_ROOT}/julia_depot"

mkdir -p "${CAMPAIGN_ROOT}"/{repo,tools,julia_depot,runs,logs}

if [[ ! -d "${JULIA_DIR}" ]]; then
  tmp_dir="$(mktemp -d)"
  trap 'rm -rf "${tmp_dir}"' EXIT
  archive="${tmp_dir}/${JULIA_TARBALL}"
  echo "Downloading Julia ${JULIA_VERSION} from ${JULIA_URL}"
  if command -v curl >/dev/null 2>&1; then
    curl -fsSL "${JULIA_URL}" -o "${archive}"
  elif command -v wget >/dev/null 2>&1; then
    wget -qO "${archive}" "${JULIA_URL}"
  else
    echo "Neither curl nor wget is available" >&2
    exit 1
  fi
  tar -xzf "${archive}" -C "${tmp_dir}"
  mv "${tmp_dir}/julia-${JULIA_VERSION}" "${JULIA_DIR}"
fi

if [[ ! -f "${JULIA_PROJECT}/Project.toml" ]]; then
  echo "Expected synced repo at ${JULIA_PROJECT}, but Project.toml was not found" >&2
  exit 1
fi

export JULIA_DEPOT_PATH="${DEPOT_DIR}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

"${JULIA_DIR}/bin/julia" --version
"${JULIA_DIR}/bin/julia" --project="${JULIA_PROJECT}" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

echo "Bootstrap complete"
echo "CAMPAIGN_ROOT=${CAMPAIGN_ROOT}"
echo "JULIA_DIR=${JULIA_DIR}"
echo "JULIA_DEPOT_PATH=${JULIA_DEPOT_PATH}"
