#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
VPOP250_ROOT = Path(
    os.environ.get(
        "SUSILO_EFAST_HOSSEINI_VPOP250_ROOT",
        REPO_ROOT / "generated/figures/vpop_pruning/susilo_efast_hosseini2020_fig5_vpop250_20260505",
    )
)
VPOP1000_ROOT = Path(
    os.environ.get(
        "SUSILO_EFAST_HOSSEINI_VPOP1000_ROOT",
        REPO_ROOT / "generated/figures/vpop_pruning/susilo_efast_hosseini2020_fig5_vpop1000_20260516",
    )
)
FULL1000_TRAJ = Path(
    os.environ.get(
        "TOPPARAM10_FULL1000_FIG5_TRAJECTORIES",
        REPO_ROOT / "generated/figures/optimization/topparam10_mosun_fig5_20260327/topparam1000_mosun_fig5_trajectories.csv",
    )
)
FULL1000_WATERFALL = Path(
    os.environ.get(
        "TOPPARAM10_FULL1000_FIG5_WATERFALL",
        REPO_ROOT / "generated/figures/optimization/topparam10_mosun_fig5_20260327/topparam1000_mosun_fig5_waterfall.csv",
    )
)
HYBRID_TRAJ = Path(
    os.environ.get(
        "TOPPARAM10_HYBRID_OLD1000_SINKHORN_MISSING_FIG5_TRAJECTORIES",
        REPO_ROOT
        / "generated/figures/vpop_pruning/topparam10_old1000_plus_sinkhorn_missing_hosseini2020_fig5_vpop1000_20260517/candidate_trajectory_summary.csv",
    )
)
HYBRID_WATERFALL = Path(
    os.environ.get(
        "TOPPARAM10_HYBRID_OLD1000_SINKHORN_MISSING_FIG5_WATERFALL",
        REPO_ROOT
        / "generated/figures/vpop_pruning/topparam10_old1000_plus_sinkhorn_missing_hosseini2020_fig5_vpop1000_20260517/candidate_waterfall.csv",
    )
)
OUT_DIR = Path(os.environ.get("SUSILO_EFAST_HOSSEINI_VPOP_VIOLIN_OUT_DIR", VPOP1000_ROOT))
DIGITIZED_ROOT = Path(
    os.environ.get(
        "HOSSEINI2020_DIGITIZED_ROOT",
        REPO_ROOT / "assets/digitization_hosseini_2020_vpop",
    )
)

REGIMEN_ORDER = [
    "c1d1_20_only_then_20_q3w",
    "c1d1_d8_d15_6p7_then_20_q3w",
    "c1d1_d8_d15_1p6_10_10_then_20_q3w",
    "c1d1_d8_1p6_20_then_20_q3w",
]

DATASET_ORDER = ["VPop250", "Full1000", "Hybrid", "Sinkhorn VPop1000"]
DATASET_COLORS = {
    "VPop250": "#4C78A8",
    "Full1000": "#6A3D9A",
    "Hybrid": "#009E73",
    "Sinkhorn VPop1000": "#D55E00",
}
DATASET_POSITIONS = {
    "VPop250": 0.68,
    "Full1000": 0.93,
    "Hybrid": 1.18,
    "Sinkhorn VPop1000": 1.43,
}
REFERENCE_POSITION = 1.68
REGIMEN_LABELS = {
    "c1d1_20_only_then_20_q3w": "20 mg C1D1",
    "c1d1_d8_d15_6p7_then_20_q3w": "6.7/6.7/6.7 mg",
    "c1d1_d8_d15_1p6_10_10_then_20_q3w": "1.6/10/10 mg",
    "c1d1_d8_1p6_20_then_20_q3w": "1.6/20 mg",
}
DIGITIZED_REGIMEN_TO_CODE = {
    "20 mg C1D1; 20 mg Day 1 subsequent cycles": "c1d1_20_only_then_20_q3w",
    "6.7 mg C1D1/C1D8/C1D15; 20 mg Day 1 subsequent cycles": "c1d1_d8_d15_6p7_then_20_q3w",
    "1.6/10/10 mg C1D1/C1D8/C1D15; 20 mg Day 1 subsequent cycles": "c1d1_d8_d15_1p6_10_10_then_20_q3w",
    "1.6/20 mg C1D1/C1D8; 20 mg Day 1 subsequent cycles": "c1d1_d8_1p6_20_then_20_q3w",
}
REGIMEN_TO_DIGITIZED = {v: k for k, v in DIGITIZED_REGIMEN_TO_CODE.items()}

METRICS = [
    {
        "metric": "profile_peak_il6_pg_ml",
        "title": "IL6 peak",
        "ylabel": "IL6 peak (pg/mL)",
        "yscale": "log",
        "ylim": (1e-3, 3e4),
        "biomarker": "IL6",
    },
    {
        "metric": "profile_peak_cd69_cd8_pct",
        "title": "CD69+ CD8+ T-cell peak",
        "ylabel": "Peak CD69+ CD8+ T cells (%)",
        "yscale": "linear",
        "ylim": (0, 105),
        "biomarker": "Tcell",
    },
    {
        "metric": "day84_tumor_change_pct",
        "title": "Day-84 tumor change",
        "ylabel": "Day-84 tumor change from baseline (%)",
        "yscale": "linear",
        "ylim": (-120, 320),
        "biomarker": "tumor",
    },
]


def fs_path(path: Path) -> str:
    resolved = str(path.resolve())
    if os.name == "nt" and not resolved.startswith("\\\\?\\"):
        return "\\\\?\\" + resolved
    return resolved


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def load_vpop_metrics(root: Path, stem: str, label: str, id_col: str = "candidate_id") -> pd.DataFrame:
    traj_path = root / f"{stem}_trajectories.csv"
    wf_path = root / f"{stem}_waterfall.csv"
    require(traj_path)
    require(wf_path)

    traj = pd.read_csv(traj_path, usecols=[id_col, "regimen", "time_day", "il6combo", "tafraction_pb_pct"])
    traj = traj.rename(columns={id_col: "candidate_id"})
    traj["candidate_id"] = traj["candidate_id"].astype(str)
    traj = traj[traj["regimen"].isin(REGIMEN_ORDER)].copy()
    peaks = (
        traj.groupby(["candidate_id", "regimen"], sort=False)
        .agg(
            profile_peak_il6_pg_ml=("il6combo", "max"),
            profile_peak_cd69_cd8_pct=("tafraction_pb_pct", "max"),
        )
        .reset_index()
    )
    wf = pd.read_csv(wf_path, usecols=[id_col, "regimen", "tumor_size_change_pct"])
    wf = wf.rename(columns={id_col: "candidate_id"})
    wf["candidate_id"] = wf["candidate_id"].astype(str)
    wf = wf[wf["regimen"].isin(REGIMEN_ORDER)].rename(columns={"tumor_size_change_pct": "day84_tumor_change_pct"})
    merged = peaks.merge(wf, on=["candidate_id", "regimen"], how="inner", validate="one_to_one")
    merged["vpop"] = label
    return merged


def load_direct_metrics(traj_path: Path, wf_path: Path, label: str, id_col: str) -> pd.DataFrame:
    require(traj_path)
    require(wf_path)
    traj = pd.read_csv(traj_path, usecols=[id_col, "regimen", "time_day", "il6combo", "tafraction_pb_pct"])
    traj = traj.rename(columns={id_col: "candidate_id"})
    traj["candidate_id"] = traj["candidate_id"].astype(str)
    traj = traj[traj["regimen"].isin(REGIMEN_ORDER)].copy()
    peaks = (
        traj.groupby(["candidate_id", "regimen"], sort=False)
        .agg(
            profile_peak_il6_pg_ml=("il6combo", "max"),
            profile_peak_cd69_cd8_pct=("tafraction_pb_pct", "max"),
        )
        .reset_index()
    )
    wf = pd.read_csv(wf_path, usecols=[id_col, "regimen", "tumor_size_change_pct"])
    wf = wf.rename(columns={id_col: "candidate_id", "tumor_size_change_pct": "day84_tumor_change_pct"})
    wf["candidate_id"] = wf["candidate_id"].astype(str)
    wf = wf[wf["regimen"].isin(REGIMEN_ORDER)].copy()
    merged = peaks.merge(wf, on=["candidate_id", "regimen"], how="inner", validate="one_to_one")
    merged["vpop"] = label
    return merged


def load_precomputed_metrics(source_path: Path, source_label: str, output_label: str) -> pd.DataFrame:
    require(source_path)
    long = pd.read_csv(source_path)
    long = long[(long["vpop"] == source_label) & (long["regimen"].isin(REGIMEN_ORDER))].copy()
    if long.empty:
        raise ValueError(f"No rows for {source_label} in {source_path}")
    wide = (
        long.pivot_table(
            index=["candidate_id", "regimen"],
            columns="metric",
            values="value",
            aggfunc="first",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    required = [m["metric"] for m in METRICS]
    missing = [col for col in required if col not in wide.columns]
    if missing:
        raise ValueError(f"Missing precomputed metric columns for {source_label}: {missing}")
    wide["candidate_id"] = wide["candidate_id"].astype(str)
    wide["vpop"] = output_label
    return wide[["candidate_id", "regimen", *required, "vpop"]]


def load_sinkhorn_vpop1000_metrics() -> pd.DataFrame:
    traj_path = VPOP1000_ROOT / "selected_vpop1000_trajectories.csv"
    wf_path = VPOP1000_ROOT / "selected_vpop1000_waterfall.csv"
    if traj_path.exists() and wf_path.exists():
        return load_vpop_metrics(VPOP1000_ROOT, "selected_vpop1000", "Sinkhorn VPop1000")
    source_path = OUT_DIR / "vpop250_vpop1000_hosseini_fig5_endpoint_violin_source.csv"
    return load_precomputed_metrics(source_path, "VPop1000", "Sinkhorn VPop1000")


def load_digitized_reference() -> pd.DataFrame:
    traj = pd.read_csv(DIGITIZED_ROOT / "hosseini2020_fig5_digitized_IL6_Tcell_green_gray_trajectories_anchor_scaled.csv")
    wf = pd.read_csv(DIGITIZED_ROOT / "hosseini2020_fig5_digitized_tumor_waterfall_shapes.csv")
    rows: list[dict[str, object]] = []
    for digitized_regimen, regimen in DIGITIZED_REGIMEN_TO_CODE.items():
        for biomarker, metric in [("IL6", "profile_peak_il6_pg_ml"), ("Tcell", "profile_peak_cd69_cd8_pct")]:
            sub = traj[(traj["regimen"] == digitized_regimen) & (traj["biomarker"] == biomarker)]
            curve_peaks = sub.groupby("curve")["value"].max()
            if {"p5_green", "median_green", "p95_green"}.issubset(curve_peaks.index):
                rows.append(
                    {
                        "regimen": regimen,
                        "metric": metric,
                        "lower": float(curve_peaks["p5_green"]),
                        "median": float(curve_peaks["median_green"]),
                        "upper": float(curve_peaks["p95_green"]),
                        "reference_kind": "digitized_curve_peak_p5_median_p95",
                    }
                )
        tumor = wf[wf["regimen"] == digitized_regimen]["tumor_change_percent_day84"].to_numpy(float)
        if tumor.size:
            rows.append(
                {
                    "regimen": regimen,
                    "metric": "day84_tumor_change_pct",
                    "lower": float(np.nanmin(tumor)),
                    "median": float(np.nanmedian(tumor)),
                    "upper": float(np.nanmax(tumor)),
                    "reference_kind": "digitized_waterfall_min_median_max",
                }
            )
    return pd.DataFrame(rows)


def write_metric_summaries(metrics: pd.DataFrame, reference: pd.DataFrame) -> list[Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    metric_long = metrics.melt(
        id_vars=["vpop", "candidate_id", "regimen"],
        value_vars=[m["metric"] for m in METRICS],
        var_name="metric",
        value_name="value",
    )
    metric_long["regimen_label"] = metric_long["regimen"].map(REGIMEN_LABELS)
    metric_path = OUT_DIR / "vpop250_full1000_hybrid_sinkhorn1000_hosseini_fig5_endpoint_violin_source.csv"
    metric_long.to_csv(metric_path, index=False)

    summary = (
        metric_long.groupby(["vpop", "regimen", "regimen_label", "metric"], sort=False)["value"]
        .agg(n="size", median="median", p05=lambda x: np.nanpercentile(x, 5), p95=lambda x: np.nanpercentile(x, 95))
        .reset_index()
    )
    summary_path = OUT_DIR / "vpop250_full1000_hybrid_sinkhorn1000_hosseini_fig5_endpoint_violin_summary.csv"
    summary.to_csv(summary_path, index=False)

    ref_path = OUT_DIR / "hosseini2020_fig5_digitized_endpoint_reference_for_violins.csv"
    reference.to_csv(ref_path, index=False)
    return [metric_path, summary_path, ref_path]


def style_violin(parts: dict[str, object], color: str) -> None:
    for body in parts["bodies"]:
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.32)
        body.set_linewidth(0.8)
    for key in ["cbars", "cmins", "cmaxes", "cmedians"]:
        if key in parts:
            parts[key].set_color(color)
            parts[key].set_linewidth(1.0)


def plot_violins(metrics: pd.DataFrame, reference: pd.DataFrame) -> list[Path]:
    outputs: list[Path] = []
    colors = DATASET_COLORS
    width = 0.17

    fig, axes = plt.subplots(len(METRICS), len(REGIMEN_ORDER), figsize=(17.4, 8.8), sharex="col")
    for row, metric_info in enumerate(METRICS):
        metric = str(metric_info["metric"])
        for col, regimen in enumerate(REGIMEN_ORDER):
            ax = axes[row, col]
            for vpop in DATASET_ORDER:
                values = metrics[(metrics["vpop"] == vpop) & (metrics["regimen"] == regimen)][metric].to_numpy(float)
                values = values[np.isfinite(values)]
                if metric_info["yscale"] == "log":
                    values = np.maximum(values, 1e-12)
                parts = ax.violinplot(
                    [values],
                    positions=[DATASET_POSITIONS[vpop]],
                    widths=width,
                    showmeans=False,
                    showmedians=True,
                    showextrema=False,
                )
                style_violin(parts, colors[vpop])
                q05, q50, q95 = np.nanpercentile(values, [5, 50, 95])
                ax.scatter([DATASET_POSITIONS[vpop]], [q50], s=18, color=colors[vpop], zorder=4)
                ax.vlines(DATASET_POSITIONS[vpop], q05, q95, color=colors[vpop], linewidth=1.2, alpha=0.95, zorder=3)

            ref = reference[(reference["regimen"] == regimen) & (reference["metric"] == metric)]
            if not ref.empty:
                lower = float(ref.iloc[0]["lower"])
                median = float(ref.iloc[0]["median"])
                upper = float(ref.iloc[0]["upper"])
                if metric_info["yscale"] == "log":
                    lower = max(lower, 1e-12)
                    median = max(median, 1e-12)
                    upper = max(upper, 1e-12)
                ax.vlines(REFERENCE_POSITION, lower, upper, color="#111111", linewidth=1.25, zorder=5)
                ax.hlines(median, REFERENCE_POSITION - 0.09, REFERENCE_POSITION + 0.09, color="#111111", linewidth=2.0, zorder=6)
                ax.scatter([REFERENCE_POSITION], [median], marker="_", s=90, color="#111111", zorder=7)

            ax.set_yscale(str(metric_info["yscale"]))
            ax.set_ylim(metric_info["ylim"])
            ax.set_xlim(0.47, 1.85)
            ax.set_xticks(
                [
                    DATASET_POSITIONS["VPop250"],
                    DATASET_POSITIONS["Full1000"],
                    DATASET_POSITIONS["Hybrid"],
                    DATASET_POSITIONS["Sinkhorn VPop1000"],
                    REFERENCE_POSITION,
                ],
                ["250", "old\n1000", "hybrid", "Sinkhorn\n1000", "Hosseini\nref"],
                fontsize=8,
            )
            if row == 0:
                ax.set_title(REGIMEN_LABELS[regimen], fontsize=10, fontweight="bold")
            if col == 0:
                ax.set_ylabel(str(metric_info["ylabel"]), fontsize=9)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(axis="both", labelsize=8)
            ax.grid(False)

    handles = [
        plt.Line2D([0], [0], color=colors["VPop250"], marker="s", linestyle="", markersize=8, alpha=0.65, label="VPop250"),
        plt.Line2D([0], [0], color=colors["Full1000"], marker="s", linestyle="", markersize=8, alpha=0.65, label="Old full1000"),
        plt.Line2D([0], [0], color=colors["Hybrid"], marker="s", linestyle="", markersize=8, alpha=0.65, label="Hybrid old1000 + eFAST missing"),
        plt.Line2D([0], [0], color=colors["Sinkhorn VPop1000"], marker="s", linestyle="", markersize=8, alpha=0.65, label="Sinkhorn VPop1000"),
        plt.Line2D([0], [0], color="#111111", marker="_", linestyle="-", markersize=10, label="Digitized Hosseini reference"),
    ]
    fig.legend(handles=handles, frameon=False, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.015))
    fig.suptitle("Side-by-side Fig. 5 endpoint distributions", y=1.055, fontsize=15, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    for ext in [".png", ".pdf", ".svg"]:
        path = OUT_DIR / f"vpop250_full1000_hybrid_sinkhorn1000_hosseini_fig5_endpoint_side_by_side_violins{ext}"
        fig.savefig(fs_path(path), dpi=260, bbox_inches="tight")
        outputs.append(path)
    plt.close(fig)
    return outputs


def main() -> None:
    for path in [
        FULL1000_TRAJ,
        FULL1000_WATERFALL,
        HYBRID_TRAJ,
        HYBRID_WATERFALL,
        DIGITIZED_ROOT / "hosseini2020_fig5_digitized_IL6_Tcell_green_gray_trajectories_anchor_scaled.csv",
        DIGITIZED_ROOT / "hosseini2020_fig5_digitized_tumor_waterfall_shapes.csv",
    ]:
        require(path)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    metrics = pd.concat(
        [
            load_precomputed_metrics(OUT_DIR / "vpop250_vpop1000_hosseini_fig5_endpoint_violin_source.csv", "VPop250", "VPop250")
            if not (VPOP250_ROOT / "selected_vpop250_trajectories.csv").exists()
            else load_vpop_metrics(VPOP250_ROOT, "selected_vpop250", "VPop250"),
            load_direct_metrics(FULL1000_TRAJ, FULL1000_WATERFALL, "Full1000", "sample_id"),
            load_direct_metrics(HYBRID_TRAJ, HYBRID_WATERFALL, "Hybrid", "candidate_id"),
            load_sinkhorn_vpop1000_metrics(),
        ],
        ignore_index=True,
    )
    reference = load_digitized_reference()
    outputs = []
    outputs.extend(write_metric_summaries(metrics, reference))
    outputs.extend(plot_violins(metrics, reference))

    meta = {
        "vpop250_root": str(VPOP250_ROOT),
        "full1000_trajectories": str(FULL1000_TRAJ),
        "full1000_waterfall": str(FULL1000_WATERFALL),
        "hybrid_trajectories": str(HYBRID_TRAJ),
        "hybrid_waterfall": str(HYBRID_WATERFALL),
        "vpop1000_root": str(VPOP1000_ROOT),
        "digitized_root": str(DIGITIZED_ROOT),
        "metrics": [m["metric"] for m in METRICS],
        "digitized_reference_caveat": "IL6 and T-cell reference markers are p5/median/p95 curve peaks from digitized published trajectories, not raw Hosseini patient-level simulations. Tumor reference uses digitized waterfall min/median/max.",
        "outputs": [str(p) for p in outputs],
    }
    meta_path = OUT_DIR / "vpop250_full1000_hybrid_sinkhorn1000_hosseini_fig5_endpoint_side_by_side_violins_meta.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"Wrote side-by-side violin outputs to {OUT_DIR}")


if __name__ == "__main__":
    main()
