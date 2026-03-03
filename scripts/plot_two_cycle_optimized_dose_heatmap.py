#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    in_dir = Path(os.environ.get("OPT2C_RESULTS_DIR", repo / "generated" / "figures" / "optimization" / "two_cycle_population"))
    out_dir = Path(os.environ.get("OPT2C_FIG_OUT_DIR", in_dir))
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_path = in_dir / "two_cycle_population_optimization_summary.csv"
    meta_path = in_dir / "two_cycle_population_optimization_meta.json"
    if not summary_path.exists() or not meta_path.exists():
        missing = [str(p) for p in [summary_path, meta_path] if not p.exists()]
        raise FileNotFoundError(f"Missing required inputs: {missing}")

    df = pd.read_csv(summary_path)
    meta = json.loads(meta_path.read_text())

    dose_cols = [c for c in df.columns if c.startswith("dose_day") and c.endswith("_mg")]
    dose_cols = sorted(dose_cols, key=lambda c: int(c.split("dose_day", 1)[1].split("_mg", 1)[0]))
    if not dose_cols:
        raise RuntimeError("No optimized dose columns found in summary CSV.")

    sort_col = "loss_opt" if "loss_opt" in df.columns else "loss_fixed"
    df = df.sort_values(sort_col, ascending=True).reset_index(drop=True)

    dose_mat = df[dose_cols].to_numpy(dtype=float)
    patient_ids = df["patient_id"].to_numpy(dtype=int)
    loss_vals = df[sort_col].to_numpy(dtype=float)

    upper = float(np.nanmax(dose_mat))
    if not np.isfinite(upper) or upper <= 0.0:
        upper = 1.0

    dose_days = [int(c.split("dose_day", 1)[1].split("_mg", 1)[0]) for c in dose_cols]
    xticks = np.arange(len(dose_cols))

    plt.style.use("default")
    fig = plt.figure(figsize=(13.0, 14.0), constrained_layout=True)
    gs = GridSpec(1, 2, width_ratios=[4.8, 2.2], wspace=0.04, figure=fig)
    ax = fig.add_subplot(gs[0, 0])
    ax_txt = fig.add_subplot(gs[0, 1], sharey=ax)
    im = ax.imshow(
        dose_mat,
        aspect="auto",
        interpolation="nearest",
        cmap="Greys",
        vmin=0.0,
        vmax=upper,
    )

    ax.set_xlabel("Dose day")
    ax.set_ylabel(f"Patients (sorted by {sort_col}, ascending)")
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"d{d}" for d in dose_days])

    n = len(patient_ids)
    step = max(1, n // 12)
    yt = np.arange(0, n, step)
    ax.set_yticks(yt)
    ax.set_yticklabels([str(patient_ids[i]) for i in yt])

    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label("Optimized dose (mg)")

    # Right-side per-row text labels with exact optimized loss values.
    ax_txt.set_xlim(0, 1)
    ax_txt.set_ylim(n - 0.5, -0.5)
    ax_txt.axis("off")
    ax_txt.text(0.0, -1.2, "Row Labels: patient_id | loss", fontsize=9, fontweight="bold")
    for i, (pid, lossv) in enumerate(zip(patient_ids, loss_vals)):
        ax_txt.text(0.0, i, f"{pid:>3d} | {lossv:.6g}", fontsize=6.0, family="monospace", va="center")

    best = float(loss_vals[0]) if len(loss_vals) else float("nan")
    worst = float(loss_vals[-1]) if len(loss_vals) else float("nan")
    fig.suptitle(
        f"Optimized Dose Heatmap (2 cycles, 3 doses/cycle)\nRows sorted by {sort_col}: best={best:.3g}, worst={worst:.3g}",
        fontsize=12,
    )
    out_path = out_dir / "two_cycle_optimized_dose_heatmap_sorted_by_loss.png"
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(out_path)


if __name__ == "__main__":
    main()
