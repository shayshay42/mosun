from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


STRATEGY_ORDER = [
    "min_dose",
    "max_dose",
    "individual_random",
    "initial_random",
    "cohort_optimum",
    "individualized_optimum",
]

STRATEGY_LABELS = {
    "min_dose": "Min dose",
    "max_dose": "Max dose",
    "individual_random": "Per-individual random",
    "initial_random": "Random initial",
    "cohort_optimum": "Cohort optimum",
    "individualized_optimum": "Optimal individualized",
}

STRATEGY_COLORS = {
    "min_dose": "#7f8c8d",
    "max_dose": "#cb4335",
    "individual_random": "#f39c12",
    "initial_random": "#7d3c98",
    "cohort_optimum": "#2874a6",
    "individualized_optimum": "#138d75",
}


def kde_curve(values: np.ndarray, n_grid: int = 256):
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size < 2:
        return None, None
    vmin = vals.min()
    vmax = vals.max()
    if np.isclose(vmin, vmax):
        pad = max(abs(vmin) * 0.05, 1e-6)
        grid = np.linspace(vmin - pad, vmax + pad, n_grid)
        dens = np.zeros_like(grid)
        dens[len(grid) // 2] = 1.0
        return grid, dens
    pad = 0.08 * (vmax - vmin)
    grid = np.linspace(vmin - pad, vmax + pad, n_grid)
    dens = gaussian_kde(vals)(grid)
    return grid, dens


def kde_curve_log10(values: np.ndarray, n_grid: int = 256):
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if vals.size < 2:
        return None, None
    logv = np.log10(vals)
    return kde_curve(logv, n_grid=n_grid)


def should_draw_vline(values: np.ndarray) -> bool:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return False
    if vals.size == 1:
        return True
    if np.allclose(vals, vals[0], rtol=1e-8, atol=1e-10):
        return True
    q05, q95 = np.quantile(vals, [0.05, 0.95])
    spread = float(q95 - q05)
    loc = float(np.median(np.abs(vals)))
    scale = max(loc, 1.0)
    return spread / scale < 1e-4


def plot_strategy(ax, values: np.ndarray, mode: str, color: str, label: str) -> None:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return
    if mode == "log10":
        vals = vals[vals > 0]
        if vals.size == 0:
            return
        xvals = np.log10(vals)
    else:
        xvals = vals

    if should_draw_vline(xvals):
        xpos = float(np.median(xvals))
        ax.axvline(xpos, color=color, linewidth=2.0, linestyle="--", alpha=0.95, label=label)
        return

    x, y = kde_curve(xvals)
    if x is None:
        return
    ax.plot(x, y, color=color, linewidth=2.0, label=label)
    ax.fill_between(x, y, color=color, alpha=0.10)


def plot_strategy_with_hist(ax, values: np.ndarray, mode: str, color: str, label: str) -> None:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return
    if mode == "log10":
        vals = vals[vals > 0]
        if vals.size == 0:
            return
        xvals = np.log10(vals)
    else:
        xvals = vals

    if should_draw_vline(xvals):
        xpos = float(np.median(xvals))
        ax.axvline(xpos, color=color, linewidth=2.0, linestyle="--", alpha=0.95, label=label)
        return

    ax.hist(
        xvals,
        bins=18,
        density=True,
        color=color,
        alpha=0.10,
        edgecolor="none",
    )
    x, y = kde_curve(xvals)
    if x is None:
        return
    ax.plot(x, y, color=color, linewidth=2.0, label=label)


def solver_short(name: str) -> str:
    text = str(name)
    for candidate in ["QNDF", "CVODE_BDF", "Rodas4P", "Tsit5", "KenCarp4", "TRBDF2"]:
        if candidate in text:
            return candidate
    return text.split("{", 1)[0]


def scenario_title(summary: dict) -> str:
    return "DLBCL-like virtual cohort under an 8-cycle mosunetuzumab step-up regimen"


def schedule_label(schedule_df: pd.DataFrame, strategy: str) -> str:
    if strategy == "individualized_optimum":
        mean_rows = schedule_df.loc[schedule_df["strategy"] == "individualized_optimum_mean"].copy()
        parts = [f"{row['decision_label']}={row['decision_dose_mg']:.1f}" for _, row in mean_rows.iterrows()]
        return f"{STRATEGY_LABELS[strategy]} [seeded local search, no AD] (mean: " + ", ".join(parts) + ")"
    if strategy == "individual_random":
        mean_rows = schedule_df.loc[schedule_df["strategy"] == "individual_random_mean"].copy()
        parts = [f"{row['decision_label']}={row['decision_dose_mg']:.1f}" for _, row in mean_rows.iterrows()]
        return f"{STRATEGY_LABELS[strategy]} (mean: " + ", ".join(parts) + ")"
    if strategy == "cohort_optimum":
        rows = schedule_df.loc[schedule_df["strategy"] == strategy].copy()
        parts = [f"{row['decision_label']}={row['decision_dose_mg']:.1f}" for _, row in rows.iterrows()]
        return f"{STRATEGY_LABELS[strategy]} [LBFGS, finite-diff grad] (" + ", ".join(parts) + ")"
    rows = schedule_df.loc[schedule_df["strategy"] == strategy].copy()
    parts = [f"{row['decision_label']}={row['decision_dose_mg']:.1f}" for _, row in rows.iterrows()]
    return STRATEGY_LABELS[strategy] + " (" + ", ".join(parts) + ")"


def build_specs():
    return [
        ("best_spd_pct", "Best %SPD", "linear"),
        ("peak_il6", "Peak IL6 (log10 scale)", "log10"),
        ("day_of_global_peak_il6", "Day of Global IL6 Peak", "linear"),
        ("auc_tdbc", "AUC of Central Drug Concentration", "linear"),
        ("tumor_auc_abs", "Tumor Trajectory AUC (log10 scale)", "log10"),
        ("tox_peak_mean", "Mean Dose-Window IL6 Peak (log10 scale)", "log10"),
        ("tox_peak_max", "Max Dose-Window IL6 Peak (log10 scale)", "log10"),
        ("tox_auc", "IL6 AUC / Horizon (log10 scale)", "log10"),
        ("tumor_terminal", "Terminal Tumor Ratio (log10 scale)", "log10"),
    ]


def draw_figure(df: pd.DataFrame, schedule_df: pd.DataFrame, summary: dict, out_path: Path, with_hist: bool) -> None:
    specs = build_specs()
    fig, axes = plt.subplots(len(specs), 1, figsize=(14, 3.2 * len(specs)))
    axes = np.atleast_1d(axes)

    for ax, (col, title, mode) in zip(axes, specs):
        for strategy in STRATEGY_ORDER:
            vals = df.loc[df["strategy"] == strategy, col].to_numpy(dtype=float)
            color = STRATEGY_COLORS[strategy]
            label = schedule_label(schedule_df, strategy)
            if with_hist:
                plot_strategy_with_hist(ax, vals, mode, color, label)
            else:
                plot_strategy(ax, vals, mode, color, label)
        ax.set_title(title)
        ax.set_ylabel("Density")
        ax.grid(False)
        if mode == "log10":
            ax.set_xlabel("log10(value)")
        else:
            ax.set_xlabel("Value")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 0.995))
    variant = "Histogram + density overlays" if with_hist else "Density overlays"
    fig.suptitle(
        scenario_title(summary)
        + "\n"
        + f"{variant}; 168-day horizon; {solver_short(summary['solver'])} solver; cohort n={summary['cohort_size']}",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 0.78, 0.985])
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    out_dir_env = os.environ.get(
        "STRAT_EVAL_OUT_DIR",
        str(repo_root / "generated" / "figures" / "optimization" / "clinical_lhs_100_dose_strategy_overlay_20260310"),
    )
    out_dir = Path(out_dir_env)

    df = pd.read_csv(out_dir / "strategy_endpoint_distributions.csv")
    schedule_df = pd.read_csv(out_dir / "strategy_dose_schedules.csv")
    summary = json.loads((out_dir / "strategy_overlay_summary.json").read_text())

    plt.rcParams.update({"figure.dpi": 160})
    draw_figure(df, schedule_df, summary, out_dir / "dose_strategy_endpoint_overlay.png", with_hist=False)
    draw_figure(df, schedule_df, summary, out_dir / "dose_strategy_endpoint_overlay_hist.png", with_hist=True)

    print(out_dir / "dose_strategy_endpoint_overlay.png")
    print(out_dir / "dose_strategy_endpoint_overlay_hist.png")


if __name__ == "__main__":
    main()
