#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


def resolve_repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    repo = resolve_repo_root()
    parser = argparse.ArgumentParser(
        description="Plot distributions of varied VPop parameters with vertical reference from params xlsx.",
    )
    parser.add_argument(
        "--patients-csv",
        type=Path,
        default=repo / "assets" / "generated_vpop" / "selected_patients.csv",
    )
    parser.add_argument(
        "--params-xlsx",
        type=Path,
        default=repo / "assets" / "params_41540_2020_145_MOESM2_ESM.xlsx",
    )
    parser.add_argument("--params-sheet", type=str, default="Sheet1")
    parser.add_argument(
        "--reference-column",
        type=str,
        default="DLBCL",
        help="Column in params workbook for vertical reference values (e.g., DLBCL, Human).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=repo / "generated" / "figures" / "optimization" / "vpop_parameter_distributions",
    )
    parser.add_argument(
        "--fig-name",
        type=str,
        default="vpop_varied_parameter_distributions_with_xlsx_reference.png",
    )
    return parser.parse_args()


def candidate_reference_name(col_name: str) -> str:
    if col_name.endswith("_init"):
        return col_name[: -len("_init")]
    return col_name


def should_use_log_x(values: np.ndarray) -> bool:
    vals = values[np.isfinite(values)]
    if len(vals) < 2:
        return False
    if np.any(vals <= 0):
        return False
    vmin = float(np.min(vals))
    vmax = float(np.max(vals))
    if vmin <= 0:
        return False
    return vmax / vmin >= 100.0


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    patients = pd.read_csv(args.patients_csv)
    params = pd.read_excel(args.params_xlsx, sheet_name=args.params_sheet)
    if "NAME" not in params.columns:
        raise ValueError(f"Expected NAME column in {args.params_xlsx}")
    if args.reference_column not in params.columns:
        raise ValueError(
            f"Reference column {args.reference_column} not found in {args.params_xlsx}. "
            f"Available: {list(params.columns)}"
        )

    params["NAME"] = params["NAME"].astype(str)
    ref_map = (
        params[params[args.reference_column].notna()][["NAME", args.reference_column]]
        .drop_duplicates(subset=["NAME"], keep="first")
        .set_index("NAME")[args.reference_column]
        .to_dict()
    )

    varied_cols = [c for c in patients.columns if c != "patient_id" and patients[c].nunique(dropna=False) > 1]
    plotted: list[tuple[str, str, float]] = []
    for c in varied_cols:
        ref_name = candidate_reference_name(c)
        if ref_name in ref_map:
            plotted.append((c, ref_name, float(ref_map[ref_name])))

    if not plotted:
        raise ValueError("No varied parameter columns matched NAME entries in params workbook.")

    n = len(plotted)
    ncols = 4
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.4 * ncols, 3.2 * nrows), constrained_layout=True)
    axes_arr = np.array(axes).reshape(-1)

    for ax in axes_arr[n:]:
        ax.axis("off")

    for i, (col, ref_name, ref_val) in enumerate(plotted):
        ax = axes_arr[i]
        x = pd.to_numeric(patients[col], errors="coerce").to_numpy(dtype=float)
        x = x[np.isfinite(x)]
        if len(x) == 0:
            ax.set_title(f"{col}\n(no finite values)")
            continue

        bins = max(12, int(np.sqrt(len(x))))
        ax.hist(x, bins=bins, density=True, color="#4C78A8", alpha=0.45, edgecolor="none")

        if len(np.unique(x)) > 1:
            kde = gaussian_kde(x)
            grid = np.linspace(float(np.min(x)), float(np.max(x)), 300)
            ax.plot(grid, kde(grid), color="#1F4E79", lw=1.8)

        ax.axvline(ref_val, color="#D62728", linestyle="--", lw=1.6)

        if should_use_log_x(x):
            ax.set_xscale("log")

        title = col if col == ref_name else f"{col} ({ref_name})"
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("value")
        ax.set_ylabel("density")

    fig.suptitle(
        f"Generated VPop Varied Parameter Distributions (n={len(patients)})\n"
        f"Dashed line = {args.reference_column} value from {args.params_xlsx.name}",
        fontsize=13,
    )

    out_path = args.out_dir / args.fig_name
    fig.savefig(out_path, dpi=170)
    plt.close(fig)
    print(out_path)
    print(f"n_varied_mapped={len(plotted)}")
    print("mapped_params=" + ",".join([f"{c}->{r}" if c != r else c for c, r, _ in plotted]))


if __name__ == "__main__":
    main()
