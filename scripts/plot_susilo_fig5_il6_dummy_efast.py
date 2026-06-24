#!/usr/bin/env python3
"""Plot the Susilo/Fig. 5 IL6 dummy-null eFAST sensitivity outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FAMILY_ORDER = [
    "tumor_burden_growth_infiltration",
    "T_cell_activation",
    "t_cell_activation_infiltration",
    "b_cell_killing",
    "PK_exposure",
    "IL6",
    "other_model",
    "dummy_null",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--in-dir",
        type=Path,
        default=Path("generated/figures/sensitivity/susilo_fig5_il6_dummy_efast_20260502"),
        help="Directory containing the eFAST CSV outputs.",
    )
    return parser.parse_args()


def short_endpoint(endpoint: str) -> str:
    parts = endpoint.split("__")
    if len(parts) >= 3 and parts[0] == "delta":
        return f"delta {parts[1]} {parts[-1]}"
    if len(parts) >= 2:
        return f"{parts[0]} {parts[-1]}"
    return endpoint


def save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def prepare_screen(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["median_ST"] = pd.to_numeric(df["median_ST"], errors="coerce")
    df["dummy_threshold_ST"] = pd.to_numeric(df["dummy_threshold_ST"], errors="coerce")
    df["endpoint_short"] = df["endpoint"].map(short_endpoint)
    return df


def plot_heatmap(screen: pd.DataFrame, out_dir: Path) -> Path:
    real = screen[screen["parameter_family"] != "dummy_null"].copy()
    top_params = (
        real.groupby("parameter")["median_ST"].max().sort_values(ascending=False).head(35).index.tolist()
    )
    top_endpoints = (
        real.groupby("endpoint")["median_ST"].max().sort_values(ascending=False).head(45).index.tolist()
    )
    heat = (
        real[real["parameter"].isin(top_params) & real["endpoint"].isin(top_endpoints)]
        .pivot_table(index="endpoint", columns="parameter", values="median_ST", aggfunc="max")
        .reindex(index=top_endpoints, columns=top_params)
        .fillna(0.0)
    )

    fig_w = max(12.0, 0.34 * len(top_params))
    fig_h = max(10.0, 0.23 * len(top_endpoints))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(heat.values, aspect="auto", cmap="viridis", vmin=0.0)
    ax.set_xticks(np.arange(len(top_params)))
    ax.set_xticklabels(top_params, rotation=70, ha="right", fontsize=7)
    ax.set_yticks(np.arange(len(top_endpoints)))
    ax.set_yticklabels([short_endpoint(x) for x in top_endpoints], fontsize=6)
    ax.set_title("Susilo/Fig. 5 eFAST total-order sensitivity (median ST)")
    ax.set_xlabel("Parameter")
    ax.set_ylabel("Endpoint")
    cbar = fig.colorbar(im, ax=ax, shrink=0.85)
    cbar.set_label("median ST across seeds")
    out = out_dir / "susilo_fig5_il6_dummy_efast_st_heatmap.png"
    save(fig, out)
    return out


def plot_dummy_threshold(screen: pd.DataFrame, out_dir: Path) -> Path:
    real = screen[screen["parameter_family"] != "dummy_null"].copy()
    top = (
        real.groupby("endpoint", as_index=False)
        .agg(top_real_ST=("median_ST", "max"), dummy_threshold_ST=("dummy_threshold_ST", "max"))
        .sort_values("top_real_ST", ascending=False)
        .head(50)
    )
    x = np.arange(len(top))
    fig, ax = plt.subplots(figsize=(15, 6))
    ax.bar(x - 0.18, top["top_real_ST"], width=0.36, color="#2563eb", label="Top real parameter")
    ax.bar(x + 0.18, top["dummy_threshold_ST"], width=0.36, color="#9ca3af", label="Max dummy-null threshold")
    ax.set_xticks(x)
    ax.set_xticklabels([short_endpoint(e) for e in top["endpoint"]], rotation=75, ha="right", fontsize=7)
    ax.set_ylabel("median ST")
    ax.set_title("Dummy-null threshold comparison by endpoint")
    ax.grid(axis="y", color="#e5e7eb", linewidth=0.8)
    ax.legend(frameon=False)
    out = out_dir / "susilo_fig5_il6_dummy_threshold_comparison.png"
    save(fig, out)
    return out


def plot_top_parameter_barplots(screen: pd.DataFrame, out_dir: Path) -> Path:
    real = screen[screen["parameter_family"] != "dummy_null"].copy()
    families = ["IL6", "T_cell_activation", "tumor_response", "paired_delta"]
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    for ax, family in zip(axes.flat, families):
        sub = real[real["endpoint_family"] == family]
        if sub.empty:
            ax.axis("off")
            continue
        top = (
            sub.groupby(["parameter", "parameter_family"], as_index=False)["median_ST"]
            .max()
            .sort_values("median_ST", ascending=True)
            .tail(15)
        )
        colors = ["#dc2626" if pf == "IL6" else "#2563eb" for pf in top["parameter_family"]]
        ax.barh(top["parameter"], top["median_ST"], color=colors)
        ax.set_title(f"Top parameters for {family}")
        ax.set_xlabel("max median ST")
        ax.grid(axis="x", color="#e5e7eb", linewidth=0.8)
    fig.suptitle("Top eFAST parameters by endpoint family", y=0.995)
    out = out_dir / "susilo_fig5_il6_top_parameter_barplots_by_family.png"
    save(fig, out)
    return out


def plot_susilo_family_summary(screen: pd.DataFrame, out_dir: Path) -> Path:
    real = screen[screen["parameter_family"] != "dummy_null"].copy()
    summary = (
        real.groupby(["endpoint_family", "parameter_family"], as_index=False)["median_ST"]
        .max()
        .rename(columns={"median_ST": "max_median_ST"})
    )
    endpoint_families = summary.groupby("endpoint_family")["max_median_ST"].max().sort_values(ascending=False).index.tolist()
    parameter_families = [f for f in FAMILY_ORDER if f in set(summary["parameter_family"])]
    pivot = (
        summary.pivot_table(index="endpoint_family", columns="parameter_family", values="max_median_ST", aggfunc="max")
        .reindex(index=endpoint_families, columns=parameter_families)
        .fillna(0.0)
    )
    fig, ax = plt.subplots(figsize=(max(10, 1.3 * len(parameter_families)), max(5, 0.6 * len(endpoint_families))))
    im = ax.imshow(pivot.values, aspect="auto", cmap="magma", vmin=0.0)
    ax.set_xticks(np.arange(len(parameter_families)))
    ax.set_xticklabels(parameter_families, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(endpoint_families)))
    ax.set_yticklabels(endpoint_families)
    ax.set_title("Susilo-family sensitivity summary")
    ax.set_xlabel("Parameter family")
    ax.set_ylabel("Endpoint family")
    cbar = fig.colorbar(im, ax=ax, shrink=0.85)
    cbar.set_label("max median ST")
    out = out_dir / "susilo_fig5_il6_susilo_family_summary.png"
    save(fig, out)
    return out


def main() -> None:
    args = parse_args()
    root = args.in_dir.resolve()
    screen_path = root / "significant_parameter_screen.csv"
    if not screen_path.exists():
        raise FileNotFoundError(screen_path)

    screen = prepare_screen(screen_path)
    outputs = [
        plot_heatmap(screen, root),
        plot_dummy_threshold(screen, root),
        plot_top_parameter_barplots(screen, root),
        plot_susilo_family_summary(screen, root),
    ]
    meta = {
        "source_root": str(root),
        "input": str(screen_path),
        "outputs": [str(path) for path in outputs] + [str(path.with_suffix(".pdf")) for path in outputs],
    }
    with (root / "susilo_fig5_il6_dummy_efast_plot_meta.json").open("w", encoding="utf-8") as handle:
        json.dump(meta, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
