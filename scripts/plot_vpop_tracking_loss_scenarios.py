#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


def plot_kde(ax: plt.Axes, vals: np.ndarray, label: str, color: str, ls: str = "-") -> None:
    vals = vals[np.isfinite(vals)]
    if len(vals) < 3:
        return
    grid = np.linspace(float(np.min(vals)), float(np.max(vals)), 300)
    if np.allclose(grid[0], grid[-1]):
        grid = np.linspace(float(np.min(vals)) - 1e-6, float(np.max(vals)) + 1e-6, 300)
    kde = gaussian_kde(vals)
    ax.plot(grid, kde(grid), color=color, linestyle=ls, linewidth=2.2, label=label)


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    in_dir = Path(
        os.environ.get(
            "OPTVPOP_RESULTS_DIR",
            repo / "generated" / "figures" / "optimization" / "vpop_tracking_loss_scenarios",
        )
    )
    out_dir = Path(os.environ.get("OPTVPOP_FIG_OUT_DIR", in_dir))
    out_dir.mkdir(parents=True, exist_ok=True)
    out_name = os.environ.get("OPTVPOP_FIG_BASENAME", "vpop_tracking_loss_kde_overlay.png")
    log_space = os.environ.get("OPTVPOP_LOG_SPACE", "0").strip().lower() in {"1", "true", "yes", "y"}

    data_path = in_dir / "vpop_tracking_loss_by_scenario.csv"
    meta_path = in_dir / "vpop_tracking_loss_meta.json"
    df = pd.read_csv(data_path)
    meta = json.loads(meta_path.read_text())

    order = ["fixed_min", "random_shared", "fixed_max"]
    color = {"fixed_min": "#1f77b4", "random_shared": "#2ca02c", "fixed_max": "#d62728"}
    label = {"fixed_min": "Fixed Min Doses", "random_shared": "Random Shared Doses", "fixed_max": "Fixed Max Doses"}

    plt.style.use("default")
    fig, ax = plt.subplots(figsize=(11, 6), facecolor="white")

    for sc in order:
        vals_raw = df.loc[df["scenario"] == sc, "loss"].to_numpy(float)
        vals = np.log10(np.maximum(vals_raw, 1e-12)) if log_space else vals_raw
        if len(vals) == 0:
            continue
        med_raw = float(np.median(vals_raw))
        plot_kde(ax, vals, f"{label[sc]} (median={med_raw:.4g})", color[sc])

    ax.set_title("VPop Tracking Loss Distribution by Dosing Scenario" + (" (log10 space)" if log_space else ""))
    ax.set_xlabel("log10(Tracking loss)" if log_space else "Tracking loss")
    ax.set_ylabel("density")
    ax.grid(False)
    ax.legend(frameon=False)
    fig.suptitle(
        f"n={meta['n_patients']} patients, cycles={meta['n_cycles']}, random_seed={meta['random_seed']}",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out_path = out_dir / out_name
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(out_path)


if __name__ == "__main__":
    main()
