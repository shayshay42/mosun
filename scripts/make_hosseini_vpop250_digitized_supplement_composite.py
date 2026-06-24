from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


REPO_ROOT = Path(__file__).resolve().parents[1]
VPOP_ROOT = (
    REPO_ROOT
    / "generated"
    / "figures"
    / "vpop_pruning"
    / "susilo_efast_hosseini2020_fig5_vpop250_20260505"
)
DIGITIZED_ROOT = REPO_ROOT / "assets" / "digitization_hosseini_2020_vpop"
OUT_ROOT = VPOP_ROOT / "supplement_vpop250_hosseini_digitized_validation_20260520"

TRAJ_PATH = VPOP_ROOT / "selected_vpop250_trajectories.csv"
SUMMARY_PATH = VPOP_ROOT / "selected_vpop250_summary.csv"
WATERFALL_PATH = VPOP_ROOT / "selected_vpop250_waterfall.csv"
DIGITIZED_TRAJ_PATH = DIGITIZED_ROOT / "hosseini2020_fig5_digitized_IL6_Tcell_green_gray_trajectories_anchor_scaled.csv"
DIGITIZED_WATERFALL_PATH = DIGITIZED_ROOT / "hosseini2020_fig5_digitized_tumor_waterfall_shapes.csv"

REGIMEN_ORDER = [
    "c1d1_20_only_then_20_q3w",
    "c1d1_d8_d15_6p7_then_20_q3w",
    "c1d1_d8_d15_1p6_10_10_then_20_q3w",
    "c1d1_d8_1p6_20_then_20_q3w",
]
REGIMEN_LABELS_SHORT = {
    "c1d1_20_only_then_20_q3w": "20 mg\nC1D1",
    "c1d1_d8_d15_6p7_then_20_q3w": "6.7/6.7/6.7 mg\nC1D1/D8/D15",
    "c1d1_d8_d15_1p6_10_10_then_20_q3w": "1.6/10/10 mg\nC1D1/D8/D15",
    "c1d1_d8_1p6_20_then_20_q3w": "1.6/20 mg\nC1D1/D8",
}
REGIMEN_ROW_LABELS = {
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

VPOP_FILL = "#BFC7BC"
VPOP_LINE = "#5E8E54"
VPOP_DOT = "#7E877E"
DIGITIZED = "#D97706"
DIGITIZED_DARK = "#7C2D12"
GRAY_TRACE = "#C8CDD3"
BLACK = "#222222"


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 450,
            "font.family": "DejaVu Sans",
            "font.size": 7.8,
            "axes.titlesize": 9.0,
            "axes.titleweight": "bold",
            "axes.labelsize": 8.0,
            "axes.labelweight": "bold",
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "legend.fontsize": 7.0,
            "axes.linewidth": 0.75,
            "xtick.major.width": 0.75,
            "ytick.major.width": 0.75,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    for path in [TRAJ_PATH, SUMMARY_PATH, WATERFALL_PATH, DIGITIZED_TRAJ_PATH, DIGITIZED_WATERFALL_PATH]:
        require(path)
    traj = pd.read_csv(TRAJ_PATH, usecols=["candidate_id", "regimen", "time_day", "il6combo", "tafraction_pb_pct"])
    summary = pd.read_csv(SUMMARY_PATH)
    wf = pd.read_csv(WATERFALL_PATH)
    dtraj = pd.read_csv(DIGITIZED_TRAJ_PATH)
    dwf = pd.read_csv(DIGITIZED_WATERFALL_PATH)
    traj = traj[traj["regimen"].isin(REGIMEN_ORDER)].copy()
    summary = summary[summary["regimen"].isin(REGIMEN_ORDER)].copy()
    wf = wf[wf["regimen"].isin(REGIMEN_ORDER)].copy()
    dtraj["regimen_code"] = dtraj["regimen"].map(DIGITIZED_REGIMEN_TO_CODE)
    dwf["regimen_code"] = dwf["regimen"].map(DIGITIZED_REGIMEN_TO_CODE)
    dtraj = dtraj[dtraj["regimen_code"].isin(REGIMEN_ORDER)].copy()
    dwf = dwf[dwf["regimen_code"].isin(REGIMEN_ORDER)].copy()
    return traj, summary, wf, dtraj, dwf


def build_endpoint_tables(
    traj: pd.DataFrame, wf: pd.DataFrame, dtraj: pd.DataFrame, dwf: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    peaks = (
        traj.groupby(["candidate_id", "regimen"], as_index=False)
        .agg(
            profile_peak_il6_pg_ml=("il6combo", "max"),
            profile_peak_cd69_cd8_pct=("tafraction_pb_pct", "max"),
        )
    )
    tumor = wf[["candidate_id", "regimen", "tumor_size_change_pct"]].rename(
        columns={"tumor_size_change_pct": "day84_tumor_change_pct"}
    )
    endpoints = peaks.merge(tumor, on=["candidate_id", "regimen"], how="inner", validate="one_to_one")
    endpoint_long = endpoints.melt(
        id_vars=["candidate_id", "regimen"],
        value_vars=["profile_peak_il6_pg_ml", "profile_peak_cd69_cd8_pct", "day84_tumor_change_pct"],
        var_name="metric",
        value_name="value",
    )
    endpoint_long["regimen_label"] = endpoint_long["regimen"].map(REGIMEN_LABELS_SHORT)
    endpoint_long["source"] = "VPop250"

    ref_rows: list[dict[str, object]] = []
    for regimen in REGIMEN_ORDER:
        for biomarker, metric in [
            ("IL6", "profile_peak_il6_pg_ml"),
            ("Tcell", "profile_peak_cd69_cd8_pct"),
        ]:
            sub = dtraj[(dtraj["regimen_code"] == regimen) & (dtraj["biomarker"] == biomarker)]
            peaks_by_curve = sub.groupby("curve")["value"].max()
            if {"p5_green", "median_green", "p95_green"}.issubset(set(peaks_by_curve.index)):
                ref_rows.append(
                    {
                        "regimen": regimen,
                        "metric": metric,
                        "lower": float(peaks_by_curve["p5_green"]),
                        "median": float(peaks_by_curve["median_green"]),
                        "upper": float(peaks_by_curve["p95_green"]),
                        "reference_kind": "digitized p5/median/p95 curve peak",
                    }
                )
        t = dwf[dwf["regimen_code"] == regimen]["tumor_change_percent_day84"].to_numpy(float)
        if t.size:
            ref_rows.append(
                {
                    "regimen": regimen,
                    "metric": "day84_tumor_change_pct",
                    "lower": float(np.nanpercentile(t, 5)),
                    "median": float(np.nanmedian(t)),
                    "upper": float(np.nanpercentile(t, 95)),
                    "reference_kind": "digitized waterfall p5/median/p95",
                }
            )
    reference = pd.DataFrame(ref_rows)
    reference["regimen_label"] = reference["regimen"].map(REGIMEN_LABELS_SHORT)
    return endpoint_long, reference


def panel_label(fig: plt.Figure, ax: plt.Axes, label: str, dx: float = -0.025, dy: float = 0.01) -> None:
    bbox = ax.get_position()
    fig.text(bbox.x0 + dx, bbox.y1 + dy, label, fontsize=11, fontweight="bold", ha="left", va="bottom")


def clean_axis(ax: plt.Axes) -> None:
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_endpoint_panel(
    ax: plt.Axes,
    endpoint_long: pd.DataFrame,
    reference: pd.DataFrame,
    metric: str,
    title: str,
    ylabel: str,
    yscale: str,
    ylim: tuple[float, float],
    rng: np.random.Generator,
) -> None:
    positions = np.arange(len(REGIMEN_ORDER), dtype=float)
    values = [
        endpoint_long[(endpoint_long["regimen"] == regimen) & (endpoint_long["metric"] == metric)]["value"].to_numpy(float)
        for regimen in REGIMEN_ORDER
    ]
    plot_values = values
    if yscale == "log":
        plot_values = [np.maximum(v, 1e-3) for v in values]

    violins = ax.violinplot(plot_values, positions=positions, widths=0.66, showmeans=False, showmedians=False, showextrema=False)
    for body in violins["bodies"]:
        body.set_facecolor(VPOP_FILL)
        body.set_edgecolor("#6D756B")
        body.set_linewidth(0.6)
        body.set_alpha(0.55)

    for x, vals in zip(positions, plot_values):
        jitter = rng.normal(0, 0.035, size=len(vals))
        ax.scatter(np.full(len(vals), x) + jitter, vals, s=4.5, color=VPOP_DOT, alpha=0.22, linewidths=0, zorder=2)
        q25, med, q75 = np.nanpercentile(vals, [25, 50, 75])
        ax.vlines(x, q25, q75, color=BLACK, lw=1.0, zorder=4)
        ax.scatter([x], [med], s=18, color=BLACK, zorder=5)

    ref = reference[reference["metric"] == metric].set_index("regimen")
    for x, regimen in zip(positions, REGIMEN_ORDER):
        if regimen not in ref.index:
            continue
        row = ref.loc[regimen]
        lower = max(float(row["lower"]), 1e-3) if yscale == "log" else float(row["lower"])
        median = max(float(row["median"]), 1e-3) if yscale == "log" else float(row["median"])
        upper = max(float(row["upper"]), 1e-3) if yscale == "log" else float(row["upper"])
        ax.vlines(x + 0.25, lower, upper, color=DIGITIZED, lw=1.2, zorder=6)
        ax.scatter([x + 0.25], [median], s=22, marker="D", color=DIGITIZED_DARK, edgecolor="white", linewidth=0.35, zorder=7)

    ax.set_title(title, loc="left", pad=6)
    ax.set_ylabel(ylabel)
    ax.set_xticks(positions)
    ax.set_xticklabels([REGIMEN_LABELS_SHORT[r] for r in REGIMEN_ORDER])
    ax.set_yscale(yscale)
    ax.set_ylim(*ylim)
    ax.set_xlim(-0.55, len(REGIMEN_ORDER) - 0.35)
    clean_axis(ax)


def draw_il6_overlay(ax: plt.Axes, traj: pd.DataFrame, summary: pd.DataFrame, dtraj: pd.DataFrame, regimen: str, show_ylabel: bool) -> None:
    sub = traj[traj["regimen"] == regimen]
    for _, patient in sub.groupby("candidate_id", sort=False):
        ax.plot(patient["time_day"], np.maximum(patient["il6combo"], 1.0), color=GRAY_TRACE, lw=0.22, alpha=0.18)
    summ = summary[summary["regimen"] == regimen].sort_values("time_day")
    x = summ["time_day"].to_numpy(float)
    ax.fill_between(x, np.maximum(summ["p05_il6combo"].to_numpy(float), 1.0), np.maximum(summ["p95_il6combo"].to_numpy(float), 1.0), color=VPOP_LINE, alpha=0.16, linewidth=0)
    ax.plot(x, np.maximum(summ["median_il6combo"].to_numpy(float), 1.0), color=VPOP_LINE, lw=1.4)
    target = dtraj[(dtraj["regimen_code"] == regimen) & (dtraj["biomarker"] == "IL6")]
    for curve, style, color, lw in [
        ("median_green", "-", DIGITIZED_DARK, 1.25),
        ("p5_green", "--", DIGITIZED, 0.95),
        ("p95_green", "--", DIGITIZED, 0.95),
    ]:
        t = target[target["curve"] == curve].sort_values("time_day")
        ax.plot(t["time_day"], np.maximum(t["value"], 1.0), linestyle=style, color=color, lw=lw)
    ax.set_yscale("log")
    ax.set_ylim(1, 1e4)
    ax.set_xlim(0, 42)
    if show_ylabel:
        ax.set_ylabel("IL6\n(pg/mL)")
    clean_axis(ax)


def draw_tcell_overlay(ax: plt.Axes, traj: pd.DataFrame, summary: pd.DataFrame, dtraj: pd.DataFrame, regimen: str, show_ylabel: bool) -> None:
    sub = traj[traj["regimen"] == regimen]
    for _, patient in sub.groupby("candidate_id", sort=False):
        ax.plot(patient["time_day"], patient["tafraction_pb_pct"], color=GRAY_TRACE, lw=0.22, alpha=0.18)
    summ = summary[summary["regimen"] == regimen].sort_values("time_day")
    x = summ["time_day"].to_numpy(float)
    ax.fill_between(x, summ["p05_tafraction_pb_pct"], summ["p95_tafraction_pb_pct"], color=VPOP_LINE, alpha=0.16, linewidth=0)
    ax.plot(x, summ["median_tafraction_pb_pct"], color=VPOP_LINE, lw=1.4)
    target = dtraj[(dtraj["regimen_code"] == regimen) & (dtraj["biomarker"] == "Tcell")]
    for curve, style, color, lw in [
        ("median_green", "-", DIGITIZED_DARK, 1.25),
        ("p5_green", "--", DIGITIZED, 0.95),
        ("p95_green", "--", DIGITIZED, 0.95),
    ]:
        t = target[target["curve"] == curve].sort_values("time_day")
        ax.plot(t["time_day"], t["value"], linestyle=style, color=color, lw=lw)
    ax.set_ylim(0, 100)
    ax.set_xlim(0, 42)
    if show_ylabel:
        ax.set_ylabel("CD69+ CD8+\nT cells (%)")
    clean_axis(ax)


def draw_waterfall_overlay(ax: plt.Axes, wf: pd.DataFrame, dwf: pd.DataFrame, regimen: str, show_ylabel: bool) -> None:
    sub = wf[wf["regimen"] == regimen].sort_values("tumor_size_change_pct", ascending=False)
    ax.plot(sub["rank_fraction"], sub["tumor_size_change_pct"], color=VPOP_LINE, lw=1.5, label="VPop250")
    ax.scatter(sub["rank_fraction"], sub["tumor_size_change_pct"], color=GRAY_TRACE, s=5, alpha=0.55, linewidths=0)
    target = dwf[dwf["regimen_code"] == regimen].sort_values("patient_rank_fraction")
    ax.plot(target["patient_rank_fraction"], target["tumor_change_percent_day84"], color=DIGITIZED_DARK, lw=1.35, label="Digitized Hosseini")
    ax.axhline(-50, color="#777777", lw=0.65, ls="--")
    ax.set_xlim(0, 1)
    ax.set_ylim(-110, 260)
    if show_ylabel:
        ax.set_ylabel("Day-84 tumor\nchange (%)")
    clean_axis(ax)


def make_figure() -> list[Path]:
    setup_style()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    traj, summary, wf, dtraj, dwf = load_data()
    endpoint_long, reference = build_endpoint_tables(traj, wf, dtraj, dwf)

    endpoint_long.to_csv(OUT_ROOT / "vpop250_endpoint_distribution_source.csv", index=False)
    reference.to_csv(OUT_ROOT / "digitized_endpoint_targets_for_overlay.csv", index=False)

    rng = np.random.default_rng(20260520)
    fig = plt.figure(figsize=(13.4, 12.3), constrained_layout=False)
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 2.2], left=0.07, right=0.985, top=0.965, bottom=0.065, hspace=0.28)
    top = outer[0].subgridspec(1, 3, wspace=0.32)
    bottom = outer[1].subgridspec(4, 3, hspace=0.22, wspace=0.25)

    ax_a = fig.add_subplot(top[0, 0])
    draw_endpoint_panel(
        ax_a,
        endpoint_long,
        reference,
        "profile_peak_il6_pg_ml",
        "Endpoint distribution: IL6 peak",
        "Peak IL6 (pg/mL)",
        "log",
        (1, 2.5e4),
        rng,
    )
    ax_b = fig.add_subplot(top[0, 1])
    draw_endpoint_panel(
        ax_b,
        endpoint_long,
        reference,
        "profile_peak_cd69_cd8_pct",
        "Endpoint distribution: activated T cells",
        "Peak CD69+ CD8+ T cells (%)",
        "linear",
        (0, 105),
        rng,
    )
    ax_c = fig.add_subplot(top[0, 2])
    draw_endpoint_panel(
        ax_c,
        endpoint_long,
        reference,
        "day84_tumor_change_pct",
        "Endpoint distribution: day-84 tumor change",
        "Tumor change from baseline (%)",
        "linear",
        (-120, 280),
        rng,
    )

    bottom_axes: list[plt.Axes] = []
    for row, regimen in enumerate(REGIMEN_ORDER):
        ax1 = fig.add_subplot(bottom[row, 0])
        ax2 = fig.add_subplot(bottom[row, 1])
        ax3 = fig.add_subplot(bottom[row, 2])
        bottom_axes.extend([ax1, ax2, ax3])
        draw_il6_overlay(ax1, traj, summary, dtraj, regimen, show_ylabel=True)
        draw_tcell_overlay(ax2, traj, summary, dtraj, regimen, show_ylabel=True)
        draw_waterfall_overlay(ax3, wf, dwf, regimen, show_ylabel=True)
        ax1.text(-0.34, 0.5, REGIMEN_ROW_LABELS[regimen], transform=ax1.transAxes, ha="right", va="center", fontsize=7.2, fontweight="bold")
        if row == 0:
            ax1.set_title("IL6 trajectories", loc="left", pad=5)
            ax2.set_title("Activated T-cell trajectories", loc="left", pad=5)
            ax3.set_title("Day-84 tumor waterfall", loc="left", pad=5)
        if row < len(REGIMEN_ORDER) - 1:
            for ax in [ax1, ax2, ax3]:
                ax.tick_params(axis="x", labelbottom=False)
        else:
            ax1.set_xlabel("Day")
            ax2.set_xlabel("Day")
            ax3.set_xlabel("Patient rank")

    handles = [
        Patch(facecolor=VPOP_FILL, edgecolor="#6D756B", alpha=0.55, label="VPop250 distribution"),
        Line2D([0], [0], color=VPOP_LINE, lw=1.5, label="VPop250 median / 5th-95th"),
        Line2D([0], [0], marker="D", color=DIGITIZED_DARK, markerfacecolor=DIGITIZED_DARK, lw=0, markersize=4.5, label="Digitized Hosseini median"),
        Line2D([0], [0], color=DIGITIZED, lw=1.0, ls="--", label="Digitized 5th/95th curves"),
    ]
    fig.legend(handles=handles, frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.55, 0.018), columnspacing=1.5, handlelength=2.0)

    panel_label(fig, ax_a, "A", dx=-0.035, dy=0.008)
    panel_label(fig, ax_b, "B", dx=-0.035, dy=0.008)
    panel_label(fig, ax_c, "C", dx=-0.035, dy=0.008)
    panel_label(fig, bottom_axes[0], "D", dx=-0.035, dy=0.01)

    outputs: list[Path] = []
    base = OUT_ROOT / "vpop250_hosseini_digitized_violin_spaghetti_supplement"
    for suffix in [".png", ".pdf", ".svg"]:
        out = base.with_suffix(suffix)
        fig.savefig(out, bbox_inches="tight", facecolor="white")
        outputs.append(out)
    plt.close(fig)

    manifest = pd.DataFrame({"artifact": [p.suffix.lstrip(".") for p in outputs], "path": [str(p) for p in outputs]})
    manifest.to_csv(OUT_ROOT / "manifest.csv", index=False)
    meta = {
        "created_by": "scripts/make_hosseini_vpop250_digitized_supplement_composite.py",
        "purpose": "Supplement multi-panel validation figure comparing eFAST-derived VPop250 to digitized Hosseini et al. Fig. 5 targets.",
        "inputs": {
            "selected_vpop250_trajectories": str(TRAJ_PATH),
            "selected_vpop250_summary": str(SUMMARY_PATH),
            "selected_vpop250_waterfall": str(WATERFALL_PATH),
            "digitized_trajectories": str(DIGITIZED_TRAJ_PATH),
            "digitized_waterfall": str(DIGITIZED_WATERFALL_PATH),
        },
        "outputs": [str(p) for p in outputs],
        "n_vpop": int(traj["candidate_id"].nunique()),
        "regimens": REGIMEN_ORDER,
        "notes": [
            "Trajectory overlays use VPop250 individual gray traces, VPop250 median/5-95% green summaries, and digitized Hosseini median/5-95% orange curves.",
            "Tumor panels use day-84 waterfall data because the digitized Hosseini source is a day-84 waterfall rather than a continuous tumor trajectory.",
            "Endpoint distribution panels show VPop250 violins with digitized target medians and 5-95% intervals.",
        ],
    }
    (OUT_ROOT / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    caption = (
        "Supplementary Fig. X. eFAST-derived VPop250 calibration to digitized Hosseini et al. Fig. 5 behavior. "
        "Panels A-C compare endpoint distributions for first-profile IL6 peak, activated CD69+ CD8+ T-cell peak, "
        "and day-84 tumor change across four Hosseini dosing scenarios. Gray/green violins show the selected VPop250; "
        "orange markers and intervals show digitized Hosseini targets. Panel D shows the corresponding trajectory "
        "and waterfall overlays: thin gray curves are individual VPop members, green curves/bands are VPop median and "
        "5th-95th percentiles, and orange curves are digitized Hosseini median and 5th-95th percentile targets."
    )
    (OUT_ROOT / "notes.md").write_text(caption + "\n", encoding="utf-8")
    return outputs


if __name__ == "__main__":
    make_figure()
