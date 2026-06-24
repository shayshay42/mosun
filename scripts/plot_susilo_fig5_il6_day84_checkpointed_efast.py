#!/usr/bin/env python3
"""Plot checkpointed IL6/day84 eFAST outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REGIMEN_LABELS = {
    "c1d1_20_only_then_20_q3w": "20 mg",
    "c1d1_d8_d15_6p7_then_20_q3w": "6.7/6.7/6.7",
    "c1d1_d8_d15_1p6_10_10_then_20_q3w": "1.6/10/10",
    "c1d1_d8_1p6_20_then_20_q3w": "1.6/20",
}

FAMILY_COLORS = {
    "IL6": "#d62728",
    "tumor_burden_growth_infiltration": "#2ca02c",
    "t_cell_activation_infiltration": "#1f77b4",
    "b_cell_killing": "#9467bd",
    "PK_exposure": "#ff7f0e",
    "dummy_null": "#7f7f7f",
    "other_model": "#8c564b",
}


def summarize_indices(indices: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        indices.groupby(
            ["endpoint", "endpoint_family", "regimen", "metric", "parameter", "parameter_family"],
            as_index=False,
        )
        .agg(median_ST=("ST", "median"), mean_ST=("ST", "mean"), median_S1=("S1", "median"))
        .sort_values(["endpoint", "median_ST"], ascending=[True, False])
    )
    grouped["regimen_label_short"] = grouped["regimen"].map(REGIMEN_LABELS).fillna(grouped["regimen"])
    return grouped


def plot_top_parameters(root: Path, summary: pd.DataFrame, thresholds: pd.DataFrame) -> list[str]:
    endpoints = list(summary["endpoint"].drop_duplicates())
    n = len(endpoints)
    ncols = 2
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(13, max(4, 4.4 * nrows)), squeeze=False)
    outputs: list[str] = []

    threshold_map = dict(zip(thresholds["endpoint"], thresholds["dummy_threshold_ST"]))
    for ax, endpoint in zip(axes.flat, endpoints):
        sub = summary[(summary["endpoint"] == endpoint) & (summary["parameter_family"] != "dummy_null")].head(15)
        sub = sub.iloc[::-1]
        colors = [FAMILY_COLORS.get(f, "#333333") for f in sub["parameter_family"]]
        ax.barh(sub["parameter"], sub["median_ST"], color=colors, alpha=0.88)
        thr = threshold_map.get(endpoint)
        if pd.notna(thr):
            ax.axvline(thr, color="#111111", linestyle="--", linewidth=1.1, label="max dummy median ST")
        label = sub["regimen_label_short"].iloc[0] if len(sub) else endpoint
        metric = sub["metric"].iloc[0] if len(sub) else ""
        ax.set_title(f"{label}: {metric}", fontsize=10)
        ax.set_xlabel("median ST")
        ax.grid(axis="x", alpha=0.25)
        ax.legend(loc="lower right", fontsize=7, frameon=False)
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.suptitle("Checkpointed eFAST: top parameters for dose-problem scalar", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    for suffix in ("png", "pdf"):
        out = root / f"susilo_fig5_il6_day84_checkpointed_efast_top_parameters.{suffix}"
        fig.savefig(out, dpi=220 if suffix == "png" else None)
        outputs.append(str(out))
    plt.close(fig)
    return outputs


def plot_dummy_threshold(root: Path, summary: pd.DataFrame, thresholds: pd.DataFrame) -> list[str]:
    joined = summary.merge(thresholds[["endpoint", "dummy_threshold_ST"]], on="endpoint", how="left")
    fig, ax = plt.subplots(figsize=(11, 6))
    for family, sub in joined.groupby("parameter_family"):
        ax.scatter(
            sub["dummy_threshold_ST"],
            sub["median_ST"],
            s=18 if family != "dummy_null" else 34,
            alpha=0.65 if family != "dummy_null" else 0.95,
            label=family,
            color=FAMILY_COLORS.get(family, "#333333"),
        )
    lim = max(float(np.nanmax(joined[["dummy_threshold_ST", "median_ST"]].to_numpy())), 0.01) * 1.05
    ax.plot([0, lim], [0, lim], color="#111111", linestyle="--", linewidth=1.0)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("endpoint dummy threshold ST")
    ax.set_ylabel("parameter median ST")
    ax.set_title("Dummy-null threshold comparison")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, ncols=2, frameon=False)
    fig.tight_layout()
    outputs = []
    for suffix in ("png", "pdf"):
        out = root / f"susilo_fig5_il6_day84_checkpointed_efast_dummy_thresholds.{suffix}"
        fig.savefig(out, dpi=220 if suffix == "png" else None)
        outputs.append(str(out))
    plt.close(fig)
    return outputs


def plot_family_summary(root: Path, summary: pd.DataFrame) -> list[str]:
    fam = (
        summary.groupby(["regimen_label_short", "parameter_family"], as_index=False)
        .agg(max_median_ST=("median_ST", "max"), median_median_ST=("median_ST", "median"))
        .sort_values(["regimen_label_short", "max_median_ST"], ascending=[True, False])
    )
    pivot = fam.pivot(index="parameter_family", columns="regimen_label_short", values="max_median_ST").fillna(0.0)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="viridis")
    ax.set_xticks(range(pivot.shape[1]), labels=pivot.columns, rotation=25, ha="right")
    ax.set_yticks(range(pivot.shape[0]), labels=pivot.index)
    ax.set_title("Maximum median ST by parameter family and Fig. 5 regimen")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("max median ST")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            ax.text(j, i, f"{pivot.iloc[i, j]:.2f}", ha="center", va="center", color="white" if pivot.iloc[i, j] > pivot.to_numpy().max() * 0.45 else "black", fontsize=8)
    fig.tight_layout()
    outputs = []
    for suffix in ("png", "pdf"):
        out = root / f"susilo_fig5_il6_day84_checkpointed_efast_family_summary.{suffix}"
        fig.savefig(out, dpi=220 if suffix == "png" else None)
        outputs.append(str(out))
    plt.close(fig)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--in-dir", type=Path, default=Path("generated/figures/sensitivity/susilo_fig5_il6_day84_checkpointed_efast_20260504"))
    args = parser.parse_args()
    root = args.in_dir
    indices = pd.read_csv(root / "efast_indices_by_endpoint.csv")
    thresholds = pd.read_csv(root / "dummy_null_thresholds.csv")
    summary = summarize_indices(indices)
    summary.to_csv(root / "efast_indices_endpoint_parameter_summary.csv", index=False)

    outputs: list[str] = []
    outputs += plot_top_parameters(root, summary, thresholds)
    outputs += plot_dummy_threshold(root, summary, thresholds)
    outputs += plot_family_summary(root, summary)
    meta = {
        "created_by": "scripts/plot_susilo_fig5_il6_day84_checkpointed_efast.py",
        "input_dir": str(root),
        "output_files": outputs,
    }
    (root / "susilo_fig5_il6_day84_checkpointed_efast_plot_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
