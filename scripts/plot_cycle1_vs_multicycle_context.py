#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def prettify(name: str) -> str:
    mapping = {
        "simulated_annealing": "Simulated Annealing",
        "nelder_mead_coordinate_descent": "Nelder-Mead + CD",
        "finite_diff_lbfgs": "Finite Diff + L-BFGS",
        "forward_ad_dto_lbfgs": "Forward DTO + L-BFGS",
        "forward_ad_otd_lbfgs": "Forward OTD + L-BFGS",
        "coordinate_descent": "Coordinate Descent",
    }
    return mapping.get(name, name)


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    cycle1_dir = repo / "generated" / "figures" / "optimization" / "patient_1_cycle1_allmethods"
    multicycle_dir = repo / "generated" / "figures" / "optimization" / "multicycle_patient_1_8cycles_vpop_cd"
    out_path = repo / "generated" / "figures" / "optimization" / "patient_1_cycle1_vs_multicycle_context.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    tr1 = pd.read_csv(cycle1_dir / "optimization_traces.csv")
    sm1 = pd.read_csv(cycle1_dir / "optimization_summary.csv")
    trm = pd.read_csv(multicycle_dir / "loss_trace.csv")

    method_order = [
        "simulated_annealing",
        "nelder_mead_coordinate_descent",
        "finite_diff_lbfgs",
        "forward_ad_dto_lbfgs",
        "forward_ad_otd_lbfgs",
    ]
    best_map = dict(zip(sm1["method"], sm1["best_loss"]))

    plt.style.use("default")
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    colors = plt.cm.tab10(np.linspace(0, 1, len(method_order)))

    for i, m in enumerate(method_order):
        sub = tr1[tr1["method"] == m].sort_values("eval").copy()
        if sub.empty:
            continue
        sub["best"] = sub["loss"].cummin()
        lbl = f"{prettify(m)} (L*={best_map.get(m, sub['best'].iloc[-1]):.4g})"
        axes[0].plot(
            sub["eval"],
            sub["best"],
            color=colors[i],
            linewidth=2.0,
            marker="o",
            markersize=3.5,
            markevery=max(len(sub) // 14, 1),
            label=lbl,
        )

    subm = trm.sort_values("eval").copy()
    subm["best"] = subm["loss"].cummin()
    axes[1].plot(
        subm["eval"],
        subm["best"],
        color="#2ca02c",
        linewidth=2.2,
        marker="o",
        markersize=3.5,
        markevery=max(len(subm) // 20, 1),
        label=f"{prettify('coordinate_descent')} (L*={subm['best'].iloc[-1]:.4g})",
    )

    axes[0].set_title("Cycle-1 (3-dose): all benchmarked methods")
    axes[0].set_xlabel("Objective evaluations")
    axes[0].set_ylabel("Best loss so far")
    axes[1].set_title("Multi-cycle robust objective: current run")
    axes[1].set_xlabel("Objective evaluations")
    axes[1].set_ylabel("Best loss so far")
    for ax in axes:
        ax.grid(False)
        ax.legend(frameon=False, fontsize=8)

    fig.suptitle("Method Context Comparison", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(out_path)


if __name__ == "__main__":
    main()
