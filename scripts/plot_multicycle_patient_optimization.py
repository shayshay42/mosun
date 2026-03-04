#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    in_dir = Path(
        os.environ.get(
            "OPTMC_RESULTS_DIR",
            repo / "generated" / "figures" / "optimization" / "multicycle_patient_1_8cycles",
        )
    )
    out_dir = Path(os.environ.get("OPTMC_FIG_OUT_DIR", in_dir))
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_path = in_dir / "optimization_summary.json"
    dose_path = in_dir / "dose_schedule.csv"
    fixed_path = in_dir / "trajectory_fixed.csv"
    opt_path = in_dir / "trajectory_optimized.csv"
    trace_path = in_dir / "loss_trace.csv"

    required = [summary_path, dose_path, fixed_path, opt_path]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required files: {missing}")

    summary = json.loads(summary_path.read_text())
    dose_df = pd.read_csv(dose_path)
    fixed_df = pd.read_csv(fixed_path)
    opt_df = pd.read_csv(opt_path)
    trace_df = pd.read_csv(trace_path) if trace_path.exists() else pd.DataFrame()

    fixed_loss = float(summary["objective_fixed"])
    opt_loss = float(summary["objective_optimized"])
    pid = int(summary["patient_id"])
    cycles = int(summary["n_cycles"])

    plt.style.use("default")
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    for ax in axes.flat:
        ax.grid(False)

    ax = axes[0, 0]
    ax.plot(
        dose_df["time_day"],
        dose_df["fixed_dose_mg"],
        color="#1f77b4",
        marker="o",
        linewidth=1.8,
        label=f"Fixed (loss={fixed_loss:.4g})",
    )
    ax.plot(
        dose_df["time_day"],
        dose_df["optimized_dose_mg"],
        color="#d62728",
        marker="s",
        linewidth=1.8,
        linestyle="--",
        label=f"Optimized (loss={opt_loss:.4g})",
    )
    ax.set_title("Dose Schedule")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Dose (mg)")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ax.plot(
        fixed_df["t_day"],
        fixed_df["btumor"],
        color="#1f77b4",
        linewidth=2.0,
        label=f"Fixed (loss={fixed_loss:.4g})",
    )
    ax.plot(
        opt_df["t_day"],
        opt_df["btumor"],
        color="#d62728",
        linewidth=2.0,
        linestyle="--",
        label=f"Optimized (loss={opt_loss:.4g})",
    )
    ax.set_title("Tumor Burden Trajectory")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Btumor")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(
        fixed_df["t_day"],
        fixed_df["il6combo"],
        color="#1f77b4",
        linewidth=2.0,
        label=f"Fixed (loss={fixed_loss:.4g})",
    )
    ax.plot(
        opt_df["t_day"],
        opt_df["il6combo"],
        color="#d62728",
        linewidth=2.0,
        linestyle="--",
        label=f"Optimized (loss={opt_loss:.4g})",
    )
    ax.set_title("Toxicity Proxy Trajectory (IL6combo)")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("IL6combo")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    term_names = [
        "tox_peak_mean_s",
        "tox_peak_max_s",
        "tox_auc_s",
        "tumor_terminal_s",
        "tumor_auc_s",
    ]
    labels = ["tox_mean", "tox_max", "tox_auc", "tum_end", "tum_auc"]
    f = [float(summary["fixed_scaled_terms"][k]) for k in term_names]
    o = [float(summary["optimized_scaled_terms"][k]) for k in term_names]
    x = np.arange(len(labels))
    w = 0.36
    ax.bar(x - w / 2, f, width=w, color="#1f77b4", alpha=0.85, label=f"Fixed (loss={fixed_loss:.4g})")
    ax.bar(x + w / 2, o, width=w, color="#d62728", alpha=0.85, label=f"Optimized (loss={opt_loss:.4g})")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20)
    ax.set_title("Scaled Loss Components")
    ax.set_ylabel("scaled term value")
    ax.legend(frameon=False, fontsize=8)

    fig.suptitle(f"Patient {pid}: Multi-cycle Balanced Optimization ({cycles} cycles)", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig_path = out_dir / "multicycle_patient_balanced_optimization.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)

    if not trace_df.empty:
        fig2, axes2 = plt.subplots(1, 2, figsize=(13, 4.2))
        for ax in axes2:
            ax.grid(False)
        axes2[0].plot(trace_df["elapsed_s"], trace_df["loss"], color="#2ca02c", marker="o", markersize=3, linewidth=1.2)
        axes2[0].set_title("Optimization Trace: Wallclock vs Loss")
        axes2[0].set_xlabel("Elapsed time (s)")
        axes2[0].set_ylabel("Loss")
        axes2[1].plot(trace_df["eval"], trace_df["loss"], color="#9467bd", marker="o", markersize=3, linewidth=1.2)
        axes2[1].set_title("Optimization Trace: Evaluations vs Loss")
        axes2[1].set_xlabel("Objective evaluation")
        axes2[1].set_ylabel("Loss")
        fig2.tight_layout()
        trace_fig = out_dir / "multicycle_patient_balanced_loss_trace.png"
        fig2.savefig(trace_fig, dpi=180)
        plt.close(fig2)
        print(trace_fig)

    print(fig_path)


if __name__ == "__main__":
    main()
