#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
JULIA_BIN="${JULIA_BIN:-julia}"
JULIA_PROJECT="${JULIA_PROJECT:-$REPO_ROOT/julia}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
OUT_DIR="${EFAST_OUT_DIR:-$REPO_ROOT/generated/figures/sensitivity/efast_best_spd}"

mkdir -p "$OUT_DIR"

echo "[run] eFAST best%SPD"
echo "[run] out dir: $OUT_DIR"
echo "[run] julia project: $JULIA_PROJECT"

"$JULIA_BIN" --project="$JULIA_PROJECT" "$REPO_ROOT/julia/run_efast_best_spd.jl" 2>&1 | tee "$OUT_DIR/run.log"

echo "[plot] eFAST sensitivity bars"
"$PYTHON_BIN" "$REPO_ROOT/scripts/plot_efast_best_spd.py" --in-dir "$OUT_DIR"

echo "[done] wrote:"
ls -1 "$OUT_DIR"
