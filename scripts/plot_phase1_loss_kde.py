#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


def add_norm_loss(df: pd.DataFrame, il6_min: float, il6_max: float, tum_min: float, tum_max: float) -> pd.DataFrame:
    out = df.copy()
    out["il6_norm"] = (out["il6_peak_0_2"] - il6_min) / max(il6_max - il6_min, 1e-12)
    out["tumor_norm"] = (out["tumor_resid_day42"] - tum_min) / max(tum_max - tum_min, 1e-12)
    out["loss_norm"] = 0.5 * out["il6_norm"] + 0.5 * out["tumor_norm"]
    return out


def plot_kde(ax: plt.Axes, x: np.ndarray, label: str, color: str, ls: str = "-") -> None:
    x = x[np.isfinite(x)]
    if len(x) < 3:
        return
    grid = np.linspace(max(0.0, np.min(x) - 0.05), np.max(x) + 0.05, 300)
    kde = gaussian_kde(x)
    ax.plot(grid, kde(grid), label=label, color=color, linestyle=ls, linewidth=2.0)


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    m_path = repo / "generated" / "phase1_matlab" / "phase1_metrics_matlab.csv"
    j_path = repo / "generated" / "phase1_julia" / "phase1_metrics_julia.csv"
    out_dir = repo / "generated" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    m = pd.read_csv(m_path)
    j = pd.read_csv(j_path)
    m = m[m["status"] == "ok"].copy()
    j = j[j["status"] == "ok"].copy()
    m["engine"] = "MATLAB"
    j["engine"] = "Julia"

    both = pd.concat([m, j], ignore_index=True)
    il6_min, il6_max = both["il6_peak_0_2"].min(), both["il6_peak_0_2"].max()
    tum_min, tum_max = both["tumor_resid_day42"].min(), both["tumor_resid_day42"].max()

    m = add_norm_loss(m, il6_min, il6_max, tum_min, tum_max)
    j = add_norm_loss(j, il6_min, il6_max, tum_min, tum_max)
    both = pd.concat([m, j], ignore_index=True)
    both.to_csv(repo / "generated" / "phase1_loss_combined.csv", index=False)

    # MATLAB-only toxicity KDE across all regimens.
    fig, ax = plt.subplots(figsize=(11, 5.5))
    regs = sorted(m["regimen"].unique())
    cmap = plt.colormaps.get_cmap("tab20")
    for i, reg in enumerate(regs):
        vals = m.loc[m["regimen"] == reg, "il6_peak_0_2"].to_numpy()
        plot_kde(ax, vals, reg, color=cmap(i / max(len(regs) - 1, 1)))
    ax.set_title("MATLAB phase-1 sweep: toxicity proxy KDE by regimen")
    ax.set_xlabel("IL6 peak day 0-2")
    ax.set_ylabel("density")
    ax.grid(alpha=0.25)
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "phase1_toxicity_kde_matlab_all_regimens.png", dpi=180)
    plt.close(fig)

    # MATLAB-only tumor burden KDE across all regimens.
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for i, reg in enumerate(regs):
        vals = m.loc[m["regimen"] == reg, "tumor_resid_day42"].to_numpy()
        plot_kde(ax, vals, reg, color=cmap(i / max(len(regs) - 1, 1)))
    ax.set_title("MATLAB phase-1 sweep: tumor burden proxy KDE by regimen")
    ax.set_xlabel("Tumor residual day 42 (Btumor_day42 / Btumor_day0)")
    ax.set_ylabel("density")
    ax.grid(alpha=0.25)
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "phase1_tumor_burden_kde_matlab_all_regimens.png", dpi=180)
    plt.close(fig)

    # MATLAB-only overview across all regimens.
    fig, ax = plt.subplots(figsize=(11, 5.5))
    cmap = plt.colormaps.get_cmap("tab20")
    for i, reg in enumerate(regs):
        vals = m.loc[m["regimen"] == reg, "loss_norm"].to_numpy()
        plot_kde(ax, vals, reg, color=cmap(i / max(len(regs) - 1, 1)))
    ax.set_title("MATLAB phase-1 sweep: normalized loss KDE by regimen")
    ax.set_xlabel("loss_norm (0.5*IL6_norm + 0.5*tumor_residual_norm)")
    ax.set_ylabel("density")
    ax.grid(alpha=0.25)
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "phase1_loss_kde_matlab_all_regimens.png", dpi=180)
    plt.close(fig)

    # MATLAB vs Julia for representative regimens.
    selected = [
        "fixed_0.4mg",
        "fixed_2.8mg",
        "step_0.4_1_2.8mg",
        "step_1_2_13.5mg",
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    axes = axes.flatten()
    colors = {"MATLAB": "#1b9e77", "Julia": "#d95f02"}
    for ax, reg in zip(axes, selected):
        for eng in ["MATLAB", "Julia"]:
            vals = both.loc[(both["regimen"] == reg) & (both["engine"] == eng), "loss_norm"].to_numpy()
            ls = "-" if eng == "MATLAB" else "--"
            plot_kde(ax, vals, eng, color=colors[eng], ls=ls)
        ax.set_title(reg)
        ax.set_xlabel("loss_norm")
        ax.grid(alpha=0.25)
        ax.legend()
    fig.suptitle("Phase-1 sweep loss KDE: MATLAB vs Julia")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_dir / "phase1_loss_kde_matlab_vs_julia_selected.png", dpi=180)
    plt.close(fig)

    # Quick numeric summary for scanability.
    summary = (
        both.groupby(["engine", "regimen"], as_index=False)["loss_norm"]
        .agg(["median", "mean", "std"])
        .reset_index()
    )
    summary.to_csv(repo / "generated" / "phase1_loss_summary_by_regimen.csv", index=False)

    print(out_dir / "phase1_toxicity_kde_matlab_all_regimens.png")
    print(out_dir / "phase1_tumor_burden_kde_matlab_all_regimens.png")
    print(out_dir / "phase1_loss_kde_matlab_all_regimens.png")
    print(out_dir / "phase1_loss_kde_matlab_vs_julia_selected.png")
    print(repo / "generated" / "phase1_loss_summary_by_regimen.csv")


if __name__ == "__main__":
    main()
