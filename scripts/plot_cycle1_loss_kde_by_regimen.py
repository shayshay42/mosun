#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


def kde_or_hist(ax, values: np.ndarray, label: str, color):
    vals = values[np.isfinite(values)]
    if len(vals) == 0:
        return
    if len(vals) >= 3 and np.std(vals) > 1e-12:
        xs = np.linspace(vals.min() - 0.1 * np.std(vals), vals.max() + 0.1 * np.std(vals), 300)
        ys = gaussian_kde(vals)(xs)
        ax.plot(xs, ys, label=label, color=color, linewidth=1.8)
    else:
        ax.hist(vals, bins=10, density=True, alpha=0.3, label=label, color=color)


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    metrics_path = repo / "generated" / "cycle1_loss" / "cycle1_loss_metrics_julia.csv"
    design_path = repo / "generated" / "phase1_design" / "regimen_events.csv"
    out_dir = repo / "generated" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not metrics_path.exists():
        raise FileNotFoundError(f"Missing metrics file: {metrics_path}")
    if not design_path.exists():
        raise FileNotFoundError(f"Missing design file: {design_path}")

    df = pd.read_csv(metrics_path)
    design = pd.read_csv(design_path)
    reg_order = (
        design[["regimen", "regimen_type"]]
        .drop_duplicates()
        .sort_values(["regimen_type", "regimen"])
    )
    fixed_regs = reg_order.loc[reg_order["regimen_type"] == "fixed_q3w", "regimen"].tolist()
    step_regs = reg_order.loc[reg_order["regimen_type"] == "stepup_q3w", "regimen"].tolist()

    ok = df[df["status"] == "ok"].copy()
    ok = ok[np.isfinite(ok["loss_cycle1"])]

    plt.style.use("ggplot")
    fig, axes = plt.subplots(1, 2, figsize=(17, 6), sharey=True)

    fixed_colors = plt.cm.tab10(np.linspace(0, 1, max(len(fixed_regs), 1)))
    for i, reg in enumerate(fixed_regs):
        vals = ok.loc[ok["regimen"] == reg, "loss_cycle1"].to_numpy(dtype=float)
        kde_or_hist(axes[0], vals, reg, fixed_colors[i % len(fixed_colors)])
    axes[0].set_title("Cycle-1 Loss Distribution - Fixed q3w")
    axes[0].set_xlabel("loss_cycle1")
    axes[0].set_ylabel("density")
    axes[0].legend(fontsize=8, ncol=2)

    step_colors = plt.cm.Set2(np.linspace(0, 1, max(len(step_regs), 1)))
    for i, reg in enumerate(step_regs):
        vals = ok.loc[ok["regimen"] == reg, "loss_cycle1"].to_numpy(dtype=float)
        kde_or_hist(axes[1], vals, reg, step_colors[i % len(step_colors)])
    axes[1].set_title("Cycle-1 Loss Distribution - Step-up q3w")
    axes[1].set_xlabel("loss_cycle1")
    axes[1].legend(fontsize=8, ncol=2)

    fig.suptitle(
        "Cycle-1 Balanced Loss KDE by Dosing Scenario\n"
        "(loss=0.5*tox/tox_scale + 0.5*tumor/tumor_scale; all patients)",
        fontsize=14,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = out_dir / "cycle1_loss_kde_all_regimens_fixed_vs_stepup.png"
    fig.savefig(out_path, dpi=180)
    print(out_path)


if __name__ == "__main__":
    main()
