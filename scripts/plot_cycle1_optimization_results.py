#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def prettify(name: str) -> str:
    mapping = {
        "simulated_annealing": "Simulated Annealing",
        "nelder_mead_coordinate_descent": "Nelder-Mead Coordinate Descent",
        "finite_diff_lbfgs": "Finite Diff + L-BFGS",
        "forward_ad_dto_lbfgs": "Forward-AD DTO + L-BFGS",
        "forward_ad_otd_lbfgs": "Forward-AD OTD + L-BFGS",
    }
    return mapping.get(name, name)


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    in_dir_default = repo / "generated" / "optimization_cycle1"
    opt_results_dir = os.environ.get("OPT_RESULTS_DIR", "").strip()
    in_dir = Path(opt_results_dir) if opt_results_dir else in_dir_default
    trace_path = in_dir / "optimization_traces.csv"
    summary_path = in_dir / "optimization_summary.csv"
    traj_path = in_dir / "best_dose_trajectory.csv"
    meta_path = in_dir / "optimization_meta.json"
    opt_fig_out = os.environ.get("OPT_FIG_OUT_DIR", "").strip()
    out_dir = Path(opt_fig_out) if opt_fig_out else in_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if not all(p.exists() for p in [trace_path, summary_path, traj_path, meta_path]):
        missing = [str(p) for p in [trace_path, summary_path, traj_path, meta_path] if not p.exists()]
        raise FileNotFoundError(f"Missing optimization outputs: {missing}")

    traces = pd.read_csv(trace_path)
    summary = pd.read_csv(summary_path)
    traj = pd.read_csv(traj_path)
    with open(meta_path) as f:
        meta = json.load(f)

    method_order = [
        "simulated_annealing",
        "nelder_mead_coordinate_descent",
        "finite_diff_lbfgs",
        "forward_ad_dto_lbfgs",
        "forward_ad_otd_lbfgs",
    ]
    method_order = [m for m in method_order if m in set(traces["method"].unique())]
    method_to_best = dict(zip(summary["method"], summary["best_loss"]))

    plt.style.use("default")
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="white")
    colors = plt.cm.tab20(np.linspace(0, 1, max(len(method_order), 1)))

    for i, m in enumerate(method_order):
        sub = traces[traces["method"] == m].sort_values("eval").copy()
        if sub.empty:
            continue
        sub["best_loss_so_far"] = sub["loss"].cummin()
        best_loss = float(method_to_best.get(m, sub["best_loss_so_far"].iloc[-1]))
        label = f"{prettify(m)} (L*={best_loss:.4f})"
        mark_every = max(len(sub) // 14, 1)
        axes[0].plot(
            sub["elapsed_s"],
            sub["best_loss_so_far"],
            label=label,
            color=colors[i],
            linewidth=2.0,
            marker="o",
            markersize=4,
            markevery=mark_every,
            alpha=0.95,
        )
        axes[1].plot(
            sub["eval"],
            sub["best_loss_so_far"],
            label=label,
            color=colors[i],
            linewidth=2.0,
            marker="o",
            markersize=4,
            markevery=mark_every,
            alpha=0.95,
        )

    axes[0].set_title("A) Wallclock Time vs Loss")
    axes[0].set_xlabel("Wallclock time (s)")
    axes[0].set_ylabel("Best loss so far")
    axes[1].set_title("B) Iterations vs Loss")
    axes[1].set_xlabel("Iteration (objective evaluations)")
    axes[1].set_ylabel("Best loss so far")
    for ax in axes:
        ax.grid(False)
        ax.set_facecolor("white")
        for spine in ax.spines.values():
            spine.set_color("#222222")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.02), fontsize=9, ncol=1, frameon=False)
    fig.suptitle("Cycle-1 Dose Optimization Comparison (3 doses @ day 0,7,14)", fontsize=14)
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])

    fig1_path = out_dir / "cycle1_dose_optimization_convergence_overlay.png"
    fig.savefig(fig1_path, dpi=180, bbox_inches="tight")

    best_method = meta["best_method"]
    best_loss = float(meta["best_loss"])
    best_doses = meta["best_doses_mg"]
    traj_label = (
        f"Best: {prettify(best_method)} doses={np.round(best_doses, 4).tolist()} "
        f"(L*={best_loss:.6f})"
    )

    fig2, ax2 = plt.subplots(1, 2, figsize=(16, 6), sharex=True, facecolor="white")
    ax2[0].plot(
        traj["time_day"],
        traj["Btumor"],
        color="#1f77b4",
        linewidth=2.2,
        marker="o",
        markersize=3.5,
        markevery=max(len(traj) // 20, 1),
        label=traj_label,
    )
    ax2[0].set_title("Tumor Burden Trajectory")
    ax2[0].set_xlabel("Time (days)")
    ax2[0].set_ylabel("Btumor")
    for td in [0.0, 7.0, 14.0]:
        ax2[0].axvline(td, color="#555555", linestyle="--", alpha=0.4, linewidth=1.5)
    ax2[0].legend(fontsize=9, frameon=False, loc="best")

    ax2[1].plot(
        traj["time_day"],
        traj["IL6combo"],
        color="#d62728",
        linewidth=2.2,
        marker="o",
        markersize=3.5,
        markevery=max(len(traj) // 20, 1),
        label=traj_label,
    )
    ax2[1].set_title("Toxicity Proxy Trajectory")
    ax2[1].set_xlabel("Time (days)")
    ax2[1].set_ylabel("IL6combo")
    for td in [0.0, 7.0, 14.0]:
        ax2[1].axvline(td, color="#555555", linestyle="--", alpha=0.4, linewidth=1.5)
    ax2[1].legend(fontsize=9, frameon=False, loc="best")
    for ax in ax2:
        ax.grid(False)
        ax.set_facecolor("white")
        for spine in ax.spines.values():
            spine.set_color("#222222")

    fig2.suptitle("Forward Simulation at Best Optimized Doses", fontsize=14)
    fig2.tight_layout(rect=[0, 0, 1, 0.95])
    fig2_path = out_dir / "cycle1_best_dose_forward_trajectories.png"
    fig2.savefig(fig2_path, dpi=180)

    print(fig1_path)
    print(fig2_path)


if __name__ == "__main__":
    main()
