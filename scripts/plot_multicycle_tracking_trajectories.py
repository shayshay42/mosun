#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    in_dir = Path(
        os.environ.get(
            "OPTMC_RESULTS_DIR",
            repo / "generated" / "figures" / "optimization" / "multicycle_tracking_patient_1_2cycles",
        )
    )
    out_dir = Path(os.environ.get("OPTMC_FIG_OUT_DIR", in_dir))
    out_dir.mkdir(parents=True, exist_ok=True)
    out_name = os.environ.get("OPTMC_TRACK_FIG_BASENAME", "multicycle_tracking_target_vs_best.png")

    fixed_path = in_dir / "trajectory_fixed.csv"
    best_path = in_dir / "trajectory_best.csv"
    target_path = in_dir / "trajectory_target.csv"
    meta_path = in_dir / "optimization_meta.json"

    fixed = pd.read_csv(fixed_path)
    best = pd.read_csv(best_path)
    target = pd.read_csv(target_path)
    meta = json.loads(meta_path.read_text())

    plt.style.use("default")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), facecolor="white")

    axes[0].plot(target["time_day"], target["Btumor_target"], "k--", linewidth=2.0, label="Target")
    axes[0].plot(fixed["time_day"], fixed["Btumor"], color="#1f77b4", linewidth=2.0, label="Fixed regimen")
    axes[0].plot(best["time_day"], best["Btumor"], color="#d62728", linewidth=2.0, label=f"Best ({meta['best_method']})")
    axes[0].set_title("Tumor Trajectory")
    axes[0].set_xlabel("Time (days)")
    axes[0].set_ylabel("Btumor")
    axes[0].grid(False)

    axes[1].plot(target["time_day"], target["IL6_target"], "k--", linewidth=2.0, label="Target")
    axes[1].plot(fixed["time_day"], fixed["IL6combo"], color="#1f77b4", linewidth=2.0, label="Fixed regimen")
    axes[1].plot(best["time_day"], best["IL6combo"], color="#d62728", linewidth=2.0, label=f"Best ({meta['best_method']})")
    axes[1].set_title("Toxicity Proxy Trajectory (IL6)")
    axes[1].set_xlabel("Time (days)")
    axes[1].set_ylabel("IL6combo")
    axes[1].grid(False)

    for td in meta["dose_times_days"]:
        axes[0].axvline(td, color="#cccccc", linestyle=":", linewidth=1.0, zorder=0)
        axes[1].axvline(td, color="#cccccc", linestyle=":", linewidth=1.0, zorder=0)

    fig.suptitle(
        f"Quadratic Tracking Objective: patient {meta['patient_id']}, {meta['n_cycles']} cycles, best={meta['best_method']}",
        fontsize=13,
    )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.02), frameon=False, ncol=3)
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    out_path = out_dir / out_name
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(out_path)


if __name__ == "__main__":
    main()
