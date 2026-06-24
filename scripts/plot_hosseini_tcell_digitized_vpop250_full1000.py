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
OUT_ROOT = (
    REPO_ROOT
    / "generated/server_results/hosseini_efast_vpop250_rp2d_stepup_clinical_to_individualized_improvement_20260512"
)
VPOP250_TRAJ = (
    REPO_ROOT
    / "generated/figures/vpop_pruning/susilo_efast_hosseini2020_fig5_vpop250_20260505/selected_vpop250_trajectories.csv"
)
VPOP250_SUMMARY = (
    REPO_ROOT
    / "generated/figures/vpop_pruning/susilo_efast_hosseini2020_fig5_vpop250_20260505/selected_vpop250_summary.csv"
)
FULL1000_TRAJ = (
    REPO_ROOT
    / "generated/figures/optimization/topparam10_mosun_fig5_20260327/topparam1000_mosun_fig5_trajectories.csv"
)
FULL1000_SUMMARY = (
    REPO_ROOT
    / "generated/figures/optimization/topparam10_mosun_fig5_20260327/topparam1000_mosun_fig5_summary.csv"
)
DIGITIZED_TCELL = (
    REPO_ROOT
    / "assets/digitization_hosseini_2020_vpop/hosseini2020_fig5_digitized_IL6_Tcell_green_gray_trajectories_anchor_scaled.csv"
)

REGIMEN_ORDER = [
    "hosseini_20_only_then_20q3w",
    "hosseini_6p7_triplet_then_20q3w",
    "paper_best",
    "hosseini_1p6_20_then_20q3w",
]
REGIMEN_LABELS = {
    "hosseini_20_only_then_20q3w": "20 q3w",
    "hosseini_6p7_triplet_then_20q3w": "6.7 x3 +20",
    "paper_best": "Hosseini step-up",
    "hosseini_1p6_20_then_20q3w": "1.6/20 +20",
}
SIM_REGIMEN_MAP = {
    "c1d1_20_only_then_20_q3w": "hosseini_20_only_then_20q3w",
    "c1d1_d8_d15_6p7_then_20_q3w": "hosseini_6p7_triplet_then_20q3w",
    "c1d1_d8_d15_1p6_10_10_then_20_q3w": "paper_best",
    "c1d1_d8_1p6_20_then_20_q3w": "hosseini_1p6_20_then_20q3w",
}
DIGITIZED_REGIMEN_MAP = {
    "20 mg C1D1; 20 mg Day 1 subsequent cycles": "hosseini_20_only_then_20q3w",
    "6.7 mg C1D1/C1D8/C1D15; 20 mg Day 1 subsequent cycles": "hosseini_6p7_triplet_then_20q3w",
    "1.6/10/10 mg C1D1/C1D8/C1D15; 20 mg Day 1 subsequent cycles": "paper_best",
    "1.6/20 mg C1D1/C1D8; 20 mg Day 1 subsequent cycles": "hosseini_1p6_20_then_20q3w",
}

BLUE = "#4C78A8"
ORANGE = "#D55E00"
BLACK = "#111111"


def fs_path(path: Path) -> str:
    resolved = str(path.resolve())
    if os.name == "nt" and not resolved.startswith("\\\\?\\"):
        return "\\\\?\\" + resolved
    return resolved


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def save_figure(fig: plt.Figure, out_base: Path, outputs: list[Path]) -> None:
    out_base.parent.mkdir(parents=True, exist_ok=True)
    for ext in [".png", ".pdf"]:
        path = out_base.with_suffix(ext)
        fig.savefig(fs_path(path), dpi=220, bbox_inches="tight")
        outputs.append(path)
    plt.close(fig)


def write_manifest(outputs: list[Path], out_dir: Path) -> None:
    rows = []
    for path in outputs:
        rows.append(
            {
                "path": str(path),
                "name": path.name,
                "kind": path.suffix.lstrip("."),
                "bytes": int(path.stat().st_size) if path.exists() else 0,
            }
        )
    pd.DataFrame(rows).to_csv(out_dir / "tcell_digitized_vpop250_full1000_manifest.csv", index=False)


def mapped_summary(path: Path, dataset: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["dashboard_regimen"] = df["regimen"].map(SIM_REGIMEN_MAP)
    df = df[df["dashboard_regimen"].isin(REGIMEN_ORDER)].copy()
    df["dataset"] = dataset
    keep = [
        "dataset",
        "dashboard_regimen",
        "regimen",
        "regimen_label",
        "time_day",
        "median_tafraction_pb_pct",
        "p05_tafraction_pb_pct",
        "p95_tafraction_pb_pct",
    ]
    return df[keep].sort_values(["dataset", "dashboard_regimen", "time_day"])


def mapped_peak_table(path: Path, dataset: str, id_col: str) -> pd.DataFrame:
    df = pd.read_csv(path, usecols=[id_col, "regimen", "time_day", "tafraction_pb_pct"])
    df["dashboard_regimen"] = df["regimen"].map(SIM_REGIMEN_MAP)
    df = df[df["dashboard_regimen"].isin(REGIMEN_ORDER)].copy()
    df = df[pd.to_numeric(df["time_day"], errors="coerce").between(0.0, 42.0)]
    peaks = (
        df.groupby([id_col, "regimen", "dashboard_regimen"], as_index=False)["tafraction_pb_pct"]
        .max()
        .rename(columns={id_col: "sample_key", "tafraction_pb_pct": "peak_tafraction_pb_pct_0_42d"})
    )
    peaks["dataset"] = dataset
    return peaks[["dataset", "sample_key", "regimen", "dashboard_regimen", "peak_tafraction_pb_pct_0_42d"]]


def digitized_tcell_curves(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df[df["biomarker"].astype(str) == "Tcell"].copy()
    df["dashboard_regimen"] = df["regimen"].map(DIGITIZED_REGIMEN_MAP)
    df = df[df["dashboard_regimen"].isin(REGIMEN_ORDER)].copy()
    return df.sort_values(["dashboard_regimen", "curve", "time_day"])


def digitized_tcell_reference(curves: pd.DataFrame) -> pd.DataFrame:
    rows = []
    required = {"p5_green", "median_green", "p95_green"}
    for digitized_regimen, dashboard_regimen in DIGITIZED_REGIMEN_MAP.items():
        sub = curves[(curves["regimen"] == digitized_regimen) & curves["time_day"].between(0.0, 42.0)]
        peaks = pd.to_numeric(sub["value"], errors="coerce").groupby(sub["curve"].astype(str)).max()
        missing = sorted(required - set(peaks.index))
        if missing:
            raise ValueError(f"Missing digitized T-cell curves for {digitized_regimen}: {missing}")
        rows.append(
            {
                "metric": "peak_tafraction_pb_pct_0_42d",
                "dashboard_regimen": dashboard_regimen,
                "digitized_regimen": digitized_regimen,
                "available": True,
                "lower_value": float(peaks["p5_green"]),
                "median_value": float(peaks["median_green"]),
                "upper_value": float(peaks["p95_green"]),
                "max_gray_envelope_value": float(peaks.get("max_gray_envelope", np.nan)),
                "lower_kind": "p5_green_peak_0_42d",
                "median_kind": "median_green_peak_0_42d",
                "upper_kind": "p95_green_peak_0_42d",
                "units": "CD69+ CD8+ T cells (%)",
                "source_file": str(DIGITIZED_TCELL),
                "digitization_basis": "Hosseini 2020 Fig. 5 T-cell panels E-H, days 0-42; green p5/median/p95 curves and visible gray upper envelope.",
            }
        )
    return pd.DataFrame(rows)


def plot_trajectory_panel(summary: pd.DataFrame, digitized: pd.DataFrame, out_dir: Path, outputs: list[Path]) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(17.0, 4.3), sharey=True)
    for ax, regimen in zip(axes, REGIMEN_ORDER):
        for dataset, color, alpha in [("eFAST VPop250", BLUE, 0.22), ("top10 full1000", ORANGE, 0.18)]:
            sub = summary[(summary["dataset"] == dataset) & (summary["dashboard_regimen"] == regimen)].sort_values("time_day")
            ax.fill_between(
                sub["time_day"],
                sub["p05_tafraction_pb_pct"],
                sub["p95_tafraction_pb_pct"],
                color=color,
                alpha=alpha,
                linewidth=0,
            )
            ax.plot(sub["time_day"], sub["median_tafraction_pb_pct"], color=color, linewidth=2.0, label=dataset)
        target = digitized[digitized["dashboard_regimen"] == regimen]
        for curve, linestyle, linewidth in [
            ("median_green", "-", 2.0),
            ("p5_green", (0, (4, 3)), 1.3),
            ("p95_green", (0, (4, 3)), 1.3),
        ]:
            t = target[target["curve"] == curve].sort_values("time_day")
            ax.plot(
                t["time_day"],
                t["value"],
                color=BLACK,
                linestyle=linestyle,
                linewidth=linewidth,
                label="Hosseini digitized" if curve == "median_green" else None,
            )
        ax.set_title(REGIMEN_LABELS[regimen], fontsize=11)
        ax.set_xlim(0, 42)
        ax.set_ylim(0, 100)
        ax.set_xlabel("Day")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("CD69+ CD8+ T cells (%)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", frameon=False, bbox_to_anchor=(0.985, 0.99))
    fig.suptitle("Activated T-cell trajectories: VPop250, full1000, and digitized Hosseini 2020 Fig. 5", y=1.03)
    fig.tight_layout()
    save_figure(fig, out_dir / "tcell_trajectory_panel_vpop250_vs_full1000_hosseini4500_digitized_overlay", outputs)


def plot_peak_violin(peaks: pd.DataFrame, reference: pd.DataFrame, out_dir: Path, outputs: list[Path]) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    positions = np.arange(1, len(REGIMEN_ORDER) + 1, dtype=float)
    offsets = {"eFAST VPop250": -0.17, "top10 full1000": 0.17}
    colors = {"eFAST VPop250": BLUE, "top10 full1000": ORANGE}
    for dataset in ["eFAST VPop250", "top10 full1000"]:
        data = [
            peaks.loc[
                (peaks["dataset"] == dataset) & (peaks["dashboard_regimen"] == regimen),
                "peak_tafraction_pb_pct_0_42d",
            ].to_numpy(float)
            for regimen in REGIMEN_ORDER
        ]
        parts = ax.violinplot(data, positions=positions + offsets[dataset], widths=0.30, showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor(colors[dataset])
            body.set_edgecolor(colors[dataset])
            body.set_alpha(0.30)
            body.set_linewidth(0.9)
        for x_pos, vals in zip(positions + offsets[dataset], data):
            ax.scatter([x_pos], [np.median(vals)], s=36, color="white", edgecolor=colors[dataset], linewidth=1.1, zorder=5)
    for _, row in reference.iterrows():
        if row["dashboard_regimen"] not in REGIMEN_ORDER:
            continue
        x = float(REGIMEN_ORDER.index(row["dashboard_regimen"]) + 1)
        low = float(row["lower_value"])
        med = float(row["median_value"])
        high = float(row["upper_value"])
        ax.vlines(x, low, high, color=BLACK, linewidth=1.45, zorder=8)
        ax.hlines([low, high], x - 0.09, x + 0.09, color=BLACK, linewidth=1.0, zorder=9)
        ax.hlines(med, x - 0.19, x + 0.19, color=BLACK, linewidth=2.0, zorder=10)
        max_gray = float(row["max_gray_envelope_value"])
        if np.isfinite(max_gray):
            ax.scatter([x], [max_gray], marker="^", s=34, facecolors="none", edgecolors=BLACK, alpha=0.45, zorder=11)
    ax.set_xticks(positions, [REGIMEN_LABELS[r] for r in REGIMEN_ORDER])
    ax.set_ylabel("Peak CD69+ CD8+ T cells over 0-42 days (%)")
    ax.set_ylim(0, 100)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    handles = [
        plt.Line2D([0], [0], color=BLUE, lw=8, alpha=0.45, label="eFAST VPop250"),
        plt.Line2D([0], [0], color=ORANGE, lw=8, alpha=0.38, label="top10 full1000"),
        plt.Line2D([0], [0], color=BLACK, lw=1.8, marker="_", markersize=12, label="Hosseini digitized p5/median/p95"),
    ]
    ax.legend(handles=handles, loc="upper right", frameon=False)
    ax.set_title("Activated T-cell peak distributions with digitized Hosseini 4500-VPop reference")
    fig.tight_layout()
    save_figure(fig, out_dir / "tcell_peak_violin_vpop250_vs_full1000_hosseini4500_whisker_overlay", outputs)


def main() -> None:
    for path in [VPOP250_TRAJ, VPOP250_SUMMARY, FULL1000_TRAJ, FULL1000_SUMMARY, DIGITIZED_TCELL]:
        require(path)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []

    summary = pd.concat(
        [
            mapped_summary(VPOP250_SUMMARY, "eFAST VPop250"),
            mapped_summary(FULL1000_SUMMARY, "top10 full1000"),
        ],
        ignore_index=True,
    )
    digitized = digitized_tcell_curves(DIGITIZED_TCELL)
    reference = digitized_tcell_reference(digitized)
    peaks = pd.concat(
        [
            mapped_peak_table(VPOP250_TRAJ, "eFAST VPop250", "candidate_id"),
            mapped_peak_table(FULL1000_TRAJ, "top10 full1000", "sample_id"),
        ],
        ignore_index=True,
    )

    counts = peaks.groupby(["dataset", "dashboard_regimen"]).size().unstack(fill_value=0)
    if not (counts.loc["eFAST VPop250", REGIMEN_ORDER] == 250).all():
        raise ValueError(f"Unexpected VPop250 T-cell peak counts:\n{counts}")
    if not (counts.loc["top10 full1000", REGIMEN_ORDER] == 1000).all():
        raise ValueError(f"Unexpected full1000 T-cell peak counts:\n{counts}")
    if reference["dashboard_regimen"].nunique() != 4:
        raise ValueError("Digitized T-cell reference did not map all four Hosseini regimens.")

    summary_path = OUT_ROOT / "tcell_trajectory_summary_vpop250_full1000.csv"
    digitized_path = OUT_ROOT / "hosseini4500_digitized_tcell_curves.csv"
    ref_path = OUT_ROOT / "hosseini4500_digitized_tcell_reference_metrics.csv"
    peaks_path = OUT_ROOT / "tcell_peak_metric_vpop250_full1000.csv"
    summary.to_csv(summary_path, index=False)
    digitized.to_csv(digitized_path, index=False)
    reference.to_csv(ref_path, index=False)
    peaks.to_csv(peaks_path, index=False)
    outputs.extend([summary_path, digitized_path, ref_path, peaks_path])

    plot_trajectory_panel(summary, digitized, OUT_ROOT, outputs)
    plot_peak_violin(peaks, reference, OUT_ROOT, outputs)

    meta = {
        "created_by": "scripts/plot_hosseini_tcell_digitized_vpop250_full1000.py",
        "output_root": str(OUT_ROOT),
        "inputs": {
            "vpop250_trajectories": str(VPOP250_TRAJ),
            "vpop250_summary": str(VPOP250_SUMMARY),
            "full1000_trajectories": str(FULL1000_TRAJ),
            "full1000_summary": str(FULL1000_SUMMARY),
            "digitized_tcell": str(DIGITIZED_TCELL),
        },
        "regimen_order": REGIMEN_ORDER,
        "counts": counts.to_dict(),
        "caveat": "Digitized Hosseini T-cell curves are raster-derived Fig. 5 references, not raw 4500-patient simulation outputs.",
        "outputs": [str(p) for p in outputs],
    }
    meta_path = OUT_ROOT / "tcell_digitized_vpop250_full1000_meta.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    outputs.append(meta_path)
    write_manifest(outputs, OUT_ROOT)
    print(OUT_ROOT)


if __name__ == "__main__":
    main()
