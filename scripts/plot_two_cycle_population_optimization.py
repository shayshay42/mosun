#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


def kde_or_hist(ax: plt.Axes, x: np.ndarray, label: str, color: str, linestyle: str = "-") -> None:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 2:
        ax.hist(x, bins=max(1, min(20, x.size + 1)), density=True, alpha=0.35, color=color, label=label)
        return
    xmin, xmax = float(np.min(x)), float(np.max(x))
    if np.isclose(xmin, xmax):
        ax.axvline(xmin, color=color, linestyle=linestyle, linewidth=2.0, label=label)
        return
    pad = 0.06 * (xmax - xmin)
    grid = np.linspace(xmin - pad, xmax + pad, 300)
    kde = gaussian_kde(x)
    ax.plot(grid, kde(grid), color=color, linestyle=linestyle, linewidth=2.2, label=label)


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
    fixed_doses = np.asarray(meta["fixed_doses_mg"], dtype=float)
    dose_times = np.asarray(meta["dose_times_days"], dtype=float)

    dose_cols = [c for c in df.columns if c.startswith("dose_day") and c.endswith("_mg")]
    dose_cols = sorted(dose_cols, key=lambda c: int(c.split("dose_day", 1)[1].split("_mg", 1)[0]))

    plt.style.use("default")

    # Figure 1: fixed-label regimen distributions over 2 cycles.
    fig1, axes1 = plt.subplots(1, 3, figsize=(15, 4.8))
    for ax in axes1:
        ax.grid(False)

    kde_or_hist(axes1[0], df["tox_fixed"].to_numpy(), "Fixed label", "#1f77b4")
    axes1[0].set_title("Toxicity Proxy (IL6 smoothmax day 0-2)")
    axes1[0].set_xlabel("toxicity proxy")
    axes1[0].set_ylabel("density")
    axes1[0].legend(frameon=False)

    kde_or_hist(axes1[1], df["tumor_fixed"].to_numpy(), "Fixed label", "#2ca02c")
    axes1[1].set_title("Tumor Proxy (Btumor day42/day0)")
    axes1[1].set_xlabel("tumor proxy")
    axes1[1].set_ylabel("density")
    axes1[1].legend(frameon=False)

    kde_or_hist(axes1[2], df["loss_fixed"].to_numpy(), "Fixed label", "#d62728")
    axes1[2].set_title("Loss")
    axes1[2].set_xlabel("loss")
    axes1[2].set_ylabel("density")
    axes1[2].legend(frameon=False)

    label_reg = meta.get("label_regimen_name", "step_0.8_2_6mg")
    fig1.suptitle(f"Fixed Regimen Distribution Over 2 Cycles ({label_reg})", fontsize=14)
    fig1.tight_layout(rect=(0, 0, 1, 0.93))
    fig1_path = out_dir / "two_cycle_fixed_label_distribution_endpoints.png"
    fig1.savefig(fig1_path, dpi=180)
    plt.close(fig1)

    # Figure 2: optimized-per-patient vs fixed distributions (all virtual patients).
    fig2, axes2 = plt.subplots(1, 3, figsize=(15, 4.8))
    for ax in axes2:
        ax.grid(False)

    kde_or_hist(axes2[0], df["tox_fixed"].to_numpy(), "Fixed", "#1f77b4", "-")
    kde_or_hist(axes2[0], df["tox_opt"].to_numpy(), "Optimized", "#ff7f0e", "--")
    axes2[0].set_title("Toxicity Proxy")
    axes2[0].set_xlabel("IL6 smoothmax day 0-2")
    axes2[0].set_ylabel("density")
    axes2[0].legend(frameon=False)

    kde_or_hist(axes2[1], df["tumor_fixed"].to_numpy(), "Fixed", "#1f77b4", "-")
    kde_or_hist(axes2[1], df["tumor_opt"].to_numpy(), "Optimized", "#ff7f0e", "--")
    axes2[1].set_title("Tumor Proxy")
    axes2[1].set_xlabel("Btumor day42/day0")
    axes2[1].set_ylabel("density")
    axes2[1].legend(frameon=False)

    kde_or_hist(axes2[2], df["loss_fixed"].to_numpy(), "Fixed", "#1f77b4", "-")
    kde_or_hist(axes2[2], df["loss_opt"].to_numpy(), "Optimized", "#ff7f0e", "--")
    axes2[2].set_title("Loss")
    axes2[2].set_xlabel("loss")
    axes2[2].set_ylabel("density")
    axes2[2].legend(frameon=False)

    fig2.suptitle("All Virtual Patients: Fixed vs Individualized Optimal Outcomes", fontsize=14)
    fig2.tight_layout(rect=(0, 0, 1, 0.93))
    fig2_path = out_dir / "two_cycle_fixed_vs_optimized_endpoint_distributions.png"
    fig2.savefig(fig2_path, dpi=180)
    plt.close(fig2)

    # Figure 3: optimized dose distributions by day with fixed-regimen overlay.
    fig3, ax3 = plt.subplots(figsize=(12.5, 5.4))
    ax3.grid(False)
    pos = np.arange(1, len(dose_cols) + 1)
    dose_data = [df[c].to_numpy() for c in dose_cols]

    bp = ax3.boxplot(
        dose_data,
        positions=pos,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.3},
    )
    for patch in bp["boxes"]:
        patch.set(facecolor="#9ecae1", edgecolor="#4a4a4a", alpha=0.75, linewidth=1.0)
    for whisker in bp["whiskers"]:
        whisker.set(color="#4a4a4a", linewidth=1.0)
    for cap in bp["caps"]:
        cap.set(color="#4a4a4a", linewidth=1.0)

    ax3.plot(
        pos,
        fixed_doses,
        color="#d62728",
        marker="o",
        markersize=6,
        linewidth=1.8,
        linestyle="--",
        label="Fixed regimen doses",
    )

    rng = np.random.default_rng(20260302)
    for i, c in enumerate(dose_cols):
        vals = df[c].to_numpy()
        xj = pos[i] + rng.uniform(-0.08, 0.08, size=vals.size)
        ax3.scatter(xj, vals, s=13, alpha=0.18, color="#1f77b4", edgecolors="none")

    xticklabels = [f"d{int(d)}" for d in dose_times]
    ax3.set_xticks(pos)
    ax3.set_xticklabels(xticklabels)
    ax3.set_xlabel("Dose day")
    ax3.set_ylabel("Dose (mg)")
    ax3.set_title("Per-Patient Optimized Dose Distribution (3 days/cycle over 2 cycles)")
    ax3.legend(frameon=False)
    fig3.tight_layout()
    fig3_path = out_dir / "two_cycle_optimized_dose_distribution_vs_fixed.png"
    fig3.savefig(fig3_path, dpi=180)
    plt.close(fig3)

    print(fig1_path)
    print(fig2_path)
    print(fig3_path)


if __name__ == "__main__":
    main()
