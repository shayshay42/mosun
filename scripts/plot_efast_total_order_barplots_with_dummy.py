#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


REGIMEN_LABELS = {
    "no_dose": "No dose",
    "c1d1_20_only_then_20_q3w": "20 mg",
    "c1d1_d8_d15_6p7_then_20_q3w": "6.7/6.7/6.7",
    "c1d1_d8_d15_1p6_10_10_then_20_q3w": "1.6/10/10",
    "c1d1_d8_1p6_20_then_20_q3w": "1.6/20",
}

FAMILY_COLORS = {
    "IL6": "#CC79A7",
    "tumor_burden_growth_infiltration": "#4C78A8",
    "t_cell_activation_infiltration": "#7F3C8D",
    "T_cell_activation": "#7F3C8D",
    "b_cell_killing": "#D55E00",
    "PK_exposure": "#009E73",
    "other_model": "#6E6E6E",
    "dummy_null": "#C7C7C7",
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_path(raw: str | Path) -> Path:
    path = Path(raw)
    return path if path.is_absolute() else repo_root() / path


def fs_path(path: Path) -> str:
    resolved = str(path.resolve())
    if os.name == "nt" and not resolved.startswith("\\\\?\\"):
        return "\\\\?\\" + resolved
    return resolved


def default_full_root() -> Path:
    return repo_root() / "generated/server_results/susilo_fig5_il6_dummy_efast_20260502/susilo_fig5_il6_dummy_efast_20260502"


def default_checkpointed_root() -> Path:
    return (
        repo_root()
        / "generated/server_results/susilo_fig5_il6_day84_checkpointed_efast_20260504/susilo_fig5_il6_day84_checkpointed_efast_20260504"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot eFAST total-order ST barplots with dummy-null controls visible.")
    parser.add_argument("--full-root", default=os.environ.get("SUSILO_FULL_EFAST_ROOT", str(default_full_root())))
    parser.add_argument("--checkpointed-root", default=os.environ.get("SUSILO_CHECKPOINTED_EFAST_ROOT", str(default_checkpointed_root())))
    parser.add_argument("--top-full-endpoints", type=int, default=12)
    parser.add_argument("--top-real-parameters", type=int, default=18)
    return parser.parse_args()


def short_endpoint(row: pd.Series) -> str:
    regimen = REGIMEN_LABELS.get(str(row.get("regimen", "")), str(row.get("regimen", "")))
    metric = str(row.get("metric", ""))
    metric = metric.replace("first_peak_il6_0_7d", "IL6 peak 0-7d")
    metric = metric.replace("peak_il6_0_21d", "IL6 peak 0-21d")
    metric = metric.replace("peak_il6_0_42d", "IL6 peak 0-42d")
    metric = metric.replace("il6combo_auc_0_7d", "IL6 AUC 0-7d")
    metric = metric.replace("il6combo_auc_0_21d", "IL6 AUC 0-21d")
    metric = metric.replace("il6combo_auc_0_42d", "IL6 AUC 0-42d")
    metric = metric.replace("tafraction_peak_0_42d", "T activation peak")
    metric = metric.replace("day84_tumor_size_change_pct", "Day-84 tumor change")
    metric = metric.replace("best_spd_pct", "Best SPD")
    metric = metric.replace("dose_problem_scalar", "Dose-tradeoff scalar")
    return f"{regimen}: {metric}"


def summarize_indices(indices: pd.DataFrame) -> pd.DataFrame:
    summary = (
        indices.groupby(
            ["endpoint", "endpoint_family", "regimen", "metric", "parameter", "parameter_family"],
            as_index=False,
        )
        .agg(
            median_ST=("ST", "median"),
            mean_ST=("ST", "mean"),
            median_S1=("S1", "median"),
            n_seeds=("seed", "nunique"),
        )
        .sort_values(["endpoint", "median_ST"], ascending=[True, False])
    )
    endpoint_meta = summary.drop_duplicates("endpoint").set_index("endpoint")
    summary["endpoint_label"] = summary["endpoint"].map(lambda e: short_endpoint(endpoint_meta.loc[e]))
    summary["is_dummy"] = summary["parameter_family"].eq("dummy_null") | summary["parameter"].astype(str).str.startswith("dummy_null")
    return summary


def choose_endpoints(summary: pd.DataFrame, top_n: int | None) -> list[str]:
    real = summary[~summary["is_dummy"]].copy()
    scores = real.groupby("endpoint")["median_ST"].max().sort_values(ascending=False)
    if top_n is None or len(scores) <= top_n:
        return scores.index.tolist()
    return scores.head(top_n).index.tolist()


def make_plot_table(summary: pd.DataFrame, thresholds: pd.DataFrame, endpoints: list[str], top_real: int) -> pd.DataFrame:
    threshold_map = thresholds.set_index("endpoint").to_dict(orient="index")
    rows = []
    for endpoint in endpoints:
        sub = summary[summary["endpoint"].eq(endpoint)].copy()
        real = sub[~sub["is_dummy"]].nlargest(top_real, "median_ST")
        dummy = sub[sub["is_dummy"]].copy()
        plot_sub = pd.concat([real, dummy], ignore_index=True).sort_values("median_ST", ascending=False)
        threshold = threshold_map.get(endpoint, {})
        for _, row in plot_sub.iterrows():
            rows.append(
                {
                    **row.to_dict(),
                    "dummy_threshold_ST": threshold.get("dummy_threshold_ST", np.nan),
                    "dummy_parameter_argmax": threshold.get("dummy_parameter_argmax", ""),
                    "n_dummy_params": threshold.get("n_dummy_params", np.nan),
                    "shown_as_top_real_or_dummy": bool(row["is_dummy"] or row["parameter"] in set(real["parameter"])),
                }
            )
    return pd.DataFrame(rows)


def plot_bundle(root: Path, prefix: str, title: str, top_endpoints: int | None, top_real: int) -> list[Path]:
    indices_path = root / "efast_indices_by_endpoint.csv"
    thresholds_path = root / "dummy_null_thresholds.csv"
    if not indices_path.exists():
        raise FileNotFoundError(indices_path)
    if not thresholds_path.exists():
        raise FileNotFoundError(thresholds_path)

    indices = pd.read_csv(indices_path)
    thresholds = pd.read_csv(thresholds_path)
    summary = summarize_indices(indices)
    endpoints = choose_endpoints(summary, top_endpoints)
    plot_table = make_plot_table(summary, thresholds, endpoints, top_real)

    n = len(endpoints)
    ncols = 2 if n <= 4 else 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.0 * ncols, max(4.6 * nrows, 5.5)), squeeze=False)
    threshold_map = thresholds.set_index("endpoint").to_dict(orient="index")

    for ax, endpoint in zip(axes.ravel(), endpoints):
        sub = plot_table[plot_table["endpoint"].eq(endpoint)].copy()
        sub = sub.sort_values("median_ST", ascending=True)
        y = np.arange(len(sub))
        labels = [
            f"{p} (dummy)" if bool(is_dummy) else str(p)
            for p, is_dummy in zip(sub["parameter"], sub["is_dummy"])
        ]
        colors = [FAMILY_COLORS.get(str(fam), "#999999") for fam in sub["parameter_family"]]
        bars = ax.barh(y, sub["median_ST"], color=colors, alpha=0.88, edgecolor="#333333", linewidth=0.35)
        for bar, is_dummy in zip(bars, sub["is_dummy"]):
            if bool(is_dummy):
                bar.set_hatch("///")
                bar.set_alpha(0.72)

        threshold = threshold_map.get(endpoint, {})
        thr = threshold.get("dummy_threshold_ST", np.nan)
        if pd.notna(thr):
            ax.axvline(float(thr), color="#111111", linestyle=(0, (4, 3)), linewidth=1.25)
            ymax = max(len(sub) - 0.5, 1)
            ax.text(
                float(thr),
                ymax,
                f" max dummy\n {threshold.get('dummy_parameter_argmax', '')}",
                ha="center",
                va="top",
                fontsize=7,
                color="#111111",
                bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.86},
            )
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel("Median total-order index, ST")
        label = sub["endpoint_label"].iloc[0] if len(sub) else endpoint
        n_seeds = int(sub["n_seeds"].max()) if len(sub) else 0
        ax.set_title(f"{label}\n{n_seeds} eFAST seed(s)", loc="left", fontsize=10)
        xmax = max(float(sub["median_ST"].max()) if len(sub) else 0.0, float(thr) if pd.notna(thr) else 0.0)
        ax.set_xlim(0.0, max(xmax * 1.12, 0.05))
        ax.grid(True, axis="x", alpha=0.22)

    for ax in axes.ravel()[n:]:
        ax.axis("off")

    families = [f for f in FAMILY_COLORS if f in set(plot_table["parameter_family"].astype(str))]
    family_handles = [
        Patch(facecolor=FAMILY_COLORS[f], edgecolor="#333333", label=f.replace("_", " "))
        for f in families
        if f != "dummy_null"
    ]
    dummy_handle = Patch(facecolor=FAMILY_COLORS["dummy_null"], edgecolor="#333333", hatch="///", label="dummy null parameter")
    threshold_handle = Line2D([0], [0], color="#111111", linestyle=(0, (4, 3)), label="max dummy threshold")
    fig.legend(
        handles=[*family_handles, dummy_handle, threshold_handle],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=4,
        frameon=False,
        fontsize=8,
    )
    fig.suptitle(title, y=0.992, fontsize=14)
    fig.text(
        0.01,
        0.012,
        f"Each panel shows the top {top_real} real parameters by median ST plus all dummy_null controls. Parameters above the dashed line exceed the endpoint-specific max dummy-null ST.",
        fontsize=8,
        color="#333333",
    )
    fig.tight_layout(rect=(0.0, 0.05, 1.0, 0.965))

    outputs = []
    for ext in [".png", ".pdf"]:
        out = root / f"{prefix}{ext}"
        fig.savefig(fs_path(out), dpi=220, bbox_inches="tight")
        outputs.append(out)
    plt.close(fig)

    table_path = root / f"{prefix}_plotted_parameters.csv"
    plot_table.to_csv(table_path, index=False)
    outputs.append(table_path)

    meta_path = root / f"{prefix}_meta.json"
    meta = {
        "created_by": "scripts/plot_efast_total_order_barplots_with_dummy.py",
        "input_indices": str(indices_path),
        "input_thresholds": str(thresholds_path),
        "title": title,
        "n_endpoints_plotted": int(len(endpoints)),
        "endpoints": endpoints,
        "top_real_parameters_per_endpoint": int(top_real),
        "dummy_parameters": sorted(summary.loc[summary["is_dummy"], "parameter"].astype(str).unique().tolist()),
        "outputs": [str(p) for p in outputs],
    }
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    outputs.append(meta_path)
    return outputs


def main() -> None:
    args = parse_args()
    full_root = resolve_path(args.full_root)
    checkpointed_root = resolve_path(args.checkpointed_root)
    outputs = []
    outputs += plot_bundle(
        full_root,
        "susilo_fig5_full_efast_total_order_st_barplots_with_dummy_top_endpoints",
        "Full Susilo Fig. 5 eFAST: total-order sensitivity with dummy-null controls",
        args.top_full_endpoints,
        args.top_real_parameters,
    )
    outputs += plot_bundle(
        checkpointed_root,
        "susilo_fig5_checkpointed_efast_total_order_st_barplots_with_dummy",
        "Checkpointed Susilo Fig. 5 eFAST: dose-tradeoff total-order sensitivity with dummy-null controls",
        None,
        args.top_real_parameters,
    )
    print("\n".join(str(p) for p in outputs if p.suffix in {".png", ".pdf"}))


if __name__ == "__main__":
    main()
