#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    in_dir = Path(
        os.environ.get(
            "OPTVPOP_RESULTS_DIR",
            repo / "generated" / "figures" / "optimization" / "vpop_tracking_nmcd",
        )
    )
    out_dir = Path(os.environ.get("OPTVPOP_FIG_OUT_DIR", in_dir))
    out_dir.mkdir(parents=True, exist_ok=True)
    out_name = os.environ.get("OPTVPOP_SPAGHETTI_FIG_BASENAME", "vpop_tracking_nmcd_spaghetti.png")

    traj_path = in_dir / "vpop_tracking_nmcd_trajectories.csv"
    if not traj_path.exists():
        raise FileNotFoundError(f"Missing trajectory file: {traj_path}")
    df = pd.read_csv(traj_path)

    plt.style.use("default")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), facecolor="white")
    color = {"fixed": "#1f77b4", "optimized": "#d62728"}
    alpha = 0.09

    for sc in ["fixed", "optimized"]:
        sub = df[df["scenario"] == sc]
        for pid, g in sub.groupby("patient_id", sort=False):
            g = g.sort_values("time_day")
            axes[0].plot(g["time_day"], g["Btumor"], color=color[sc], alpha=alpha, linewidth=1.0)
            axes[1].plot(g["time_day"], g["IL6combo"], color=color[sc], alpha=alpha, linewidth=1.0)

    # median envelopes for readability
    for sc in ["fixed", "optimized"]:
        sub = df[df["scenario"] == sc].copy()
        med = sub.groupby("time_day", as_index=False)[["Btumor", "IL6combo"]].median()
        axes[0].plot(med["time_day"], med["Btumor"], color=color[sc], linewidth=2.4, label=f"{sc} median")
        axes[1].plot(med["time_day"], med["IL6combo"], color=color[sc], linewidth=2.4, label=f"{sc} median")

    axes[0].set_title("Tumor Spaghetti (all patients)")
    axes[0].set_xlabel("Time (days)")
    axes[0].set_ylabel("Btumor")
    axes[0].grid(False)

    axes[1].set_title("Toxicity Spaghetti (all patients)")
    axes[1].set_xlabel("Time (days)")
    axes[1].set_ylabel("IL6combo")
    axes[1].grid(False)

    handles0, labels0 = axes[0].get_legend_handles_labels()
    fig.legend(handles0, labels0, loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=2, frameon=False)
    fig.suptitle("VPop Tracking Loss Optimization: Fixed vs Optimized Trajectories", fontsize=13)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    out_path = out_dir / out_name
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(out_path)


if __name__ == "__main__":
    main()
