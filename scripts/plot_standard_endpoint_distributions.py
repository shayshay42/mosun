#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


def plot_overlay(
    matlab_vals: np.ndarray,
    julia_vals: np.ndarray,
    xlabel: str,
    title: str,
    out_path: Path,
) -> None:
    all_vals = np.concatenate([matlab_vals, julia_vals])
    use_log = np.all(all_vals > 0.0) and (all_vals.max() / all_vals.min() > 100.0)
    if use_log:
        matlab_vals = np.log10(matlab_vals)
        julia_vals = np.log10(julia_vals)
        xlabel = f"log10({xlabel})"

    lo = min(matlab_vals.min(), julia_vals.min())
    hi = max(matlab_vals.max(), julia_vals.max())
    if hi <= lo:
        hi = lo + 1.0
    bins = np.linspace(lo, hi, 10)
    grid = np.linspace(lo, hi, 300)

    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    ax.hist(matlab_vals, bins=bins, density=True, alpha=0.35, color="#1f77b4", label="MATLAB")
    ax.hist(julia_vals, bins=bins, density=True, alpha=0.35, color="#ff7f0e", label="Julia")

    if len(matlab_vals) >= 3:
        ax.plot(grid, gaussian_kde(matlab_vals)(grid), color="#1f77b4", linewidth=2.0)
    if len(julia_vals) >= 3:
        ax.plot(grid, gaussian_kde(julia_vals)(grid), color="#ff7f0e", linewidth=2.0)

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    matlab_path = repo / "generated" / "phase1_standard_matlab" / "phase1_metrics_matlab.csv"
    julia_path = repo / "generated" / "phase1_standard_julia" / "phase1_metrics_julia.csv"
    out_dir = repo / "generated" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    m = pd.read_csv(matlab_path)
    j = pd.read_csv(julia_path)
    m = m[m["status"] == "ok"].copy()
    j = j[j["status"] == "ok"].copy()

    merged = m.merge(j, on=["regimen", "patient_id"], suffixes=("_m", "_j"))
    merged = merged.sort_values("regimen")

    plot_overlay(
        matlab_vals=merged["il6_peak_0_2_m"].to_numpy(float),
        julia_vals=merged["il6_peak_0_2_j"].to_numpy(float),
        xlabel="IL6 peak day 0-2",
        title="Toxicity Proxy Distribution Across Dose Scenarios",
        out_path=out_dir / "phase1_standard_toxicity_distribution_overlay.png",
    )

    plot_overlay(
        matlab_vals=merged["tumor_resid_day42_m"].to_numpy(float),
        julia_vals=merged["tumor_resid_day42_j"].to_numpy(float),
        xlabel="Tumor residual at day 42 (Btumor_day42 / Btumor_day0)",
        title="Tumor Burden Proxy Distribution Across Dose Scenarios",
        out_path=out_dir / "phase1_standard_tumor_distribution_overlay.png",
    )

    print(out_dir / "phase1_standard_toxicity_distribution_overlay.png")
    print(out_dir / "phase1_standard_tumor_distribution_overlay.png")


if __name__ == "__main__":
    main()

