#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_endpoint_grid(
    reg_order: list[str],
    matlab_map: dict[str, Path],
    julia_map: dict[str, Path],
    ycol: str,
    title: str,
    ylabel: str,
    out_path: Path,
) -> None:
    n = len(reg_order)
    ncols = 4
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 3.0 * nrows), squeeze=False)
    axes_flat = axes.flatten()

    for i, reg in enumerate(reg_order):
        ax = axes_flat[i]
        mdf = pd.read_csv(matlab_map[reg])
        jdf = pd.read_csv(julia_map[reg])
        merged = mdf.merge(jdf, on=["time", "regimen", "patient_id"], suffixes=("_matlab", "_julia"))

        ym = merged[f"{ycol}_matlab"]
        yj = merged[f"{ycol}_julia"]
        abs_err = (yj - ym).abs().max()
        amp = max(ym.abs().max(), 1e-12)
        norm_err = abs_err / amp

        ax.plot(merged["time"], ym, color="black", linewidth=1.8, label="MATLAB")
        ax.plot(merged["time"], yj, color="#D55E00", linestyle="--", linewidth=1.8, label="Julia")
        ax.set_title(f"{reg}\nmax|Δ|={abs_err:.3g}, max|Δ|/max|MATLAB|={norm_err:.3g}", fontsize=9)
        ax.set_xlabel("Time (days)")
        ax.grid(alpha=0.25)

    for i in range(n, len(axes_flat)):
        axes_flat[i].axis("off")

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False)
    fig.suptitle(title, y=0.995, fontsize=12)
    fig.text(0.01, 0.5, ylabel, va="center", rotation="vertical")
    fig.tight_layout(rect=[0.02, 0.0, 1.0, 0.95])
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    design_path = repo / "generated" / "phase1_design_standard" / "regimen_events.csv"
    matlab_manifest_path = repo / "generated" / "phase1_standard_matlab_traces" / "manifest.csv"
    julia_manifest_path = repo / "generated" / "phase1_standard_julia_traces" / "manifest.csv"
    out_dir = repo / "generated" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    reg_df = pd.read_csv(design_path)
    reg_order = reg_df["regimen"].drop_duplicates().tolist()

    mm = pd.read_csv(matlab_manifest_path)
    jm = pd.read_csv(julia_manifest_path)
    mm = mm[mm["status"] == "ok"].copy()
    jm = jm[jm["status"] == "ok"].copy()

    matlab_map = {r.regimen: Path(r.trace_csv) for _, r in mm.iterrows()}
    julia_map = {r.regimen: Path(r.trace_csv) for _, r in jm.iterrows()}

    missing = [r for r in reg_order if r not in matlab_map or r not in julia_map]
    if missing:
        raise RuntimeError(f"Missing trace CSV for regimens: {missing}")

    tox_out = out_dir / "phase1_standard_toxicity_trace_overlay_all_regimens.png"
    tum_out = out_dir / "phase1_standard_tumor_trace_overlay_all_regimens.png"

    plot_endpoint_grid(
        reg_order=reg_order,
        matlab_map=matlab_map,
        julia_map=julia_map,
        ycol="IL6combo",
        title="MATLAB vs Julia Trace Parity: Toxicity Proxy (IL6combo)",
        ylabel="IL6combo",
        out_path=tox_out,
    )
    plot_endpoint_grid(
        reg_order=reg_order,
        matlab_map=matlab_map,
        julia_map=julia_map,
        ycol="Btumor",
        title="MATLAB vs Julia Trace Parity: Tumor Burden Proxy (Btumor)",
        ylabel="Btumor",
        out_path=tum_out,
    )

    print(tox_out)
    print(tum_out)


if __name__ == "__main__":
    main()
