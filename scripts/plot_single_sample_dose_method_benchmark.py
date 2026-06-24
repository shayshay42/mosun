#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D


OBJECTIVES = ["simple", "auc_combo", "remission_il6", "clinical", "tracking"]
METHOD_ORDER = [
    "metropolis_anneal",
    "random_search",
    "nelder_mead",
    "cma_es",
    "finite_diff_bfgs",
    "finite_diff_lbfgs",
    "dto_forward_ad_lbfgs",
    "dto_forward_ad_bfgs",
    "dto_forward_ad_ipnewton",
    "dto_reverse_ad_lbfgs",
    "otd_reverse_ad_lbfgs",
]

METHOD_STYLES = {
    "metropolis_anneal": dict(color="#1f77b4", marker="o", linestyle="-"),
    "random_search": dict(color="#d55e00", marker="s", linestyle="-"),
    "nelder_mead": dict(color="#009e73", marker="^", linestyle="-"),
    "cma_es": dict(color="#3b5b92", marker="h", linestyle="-"),
    "finite_diff_bfgs": dict(color="#4d4d4d", marker="D", linestyle="--"),
    "finite_diff_lbfgs": dict(color="#4d4d4d", marker="D", linestyle="--"),
    "dto_forward_ad_lbfgs": dict(color="#7b61ff", marker="P", linestyle="-"),
    "dto_forward_ad_bfgs": dict(color="#6d28d9", marker="X", linestyle="-"),
    "dto_forward_ad_ipnewton": dict(color="#c2410c", marker="d", linestyle="-"),
    "dto_reverse_ad_lbfgs": dict(color="#8c564b", marker="X", linestyle=":"),
    "otd_reverse_ad_lbfgs": dict(color="#17becf", marker="*", linestyle=":"),
}

OBJECTIVE_SUMMARY_STYLES = {
    "simple": dict(marker="o", linestyle="-"),
    "auc_combo": dict(marker="s", linestyle="--"),
    "remission_il6": dict(marker="^", linestyle=":"),
    "clinical": dict(marker="D", linestyle="-."),
    "tracking": dict(marker="P", linestyle=(0, (3, 1, 1, 1))),
}

REFERENCE_STYLES = {
    "paper_best": dict(color="#111827", marker="*", linestyle="-", label="Paper 1.6/10/10/20/20/20"),
    "rp2d": dict(color="#7c3aed", marker="D", linestyle="--", label="RP2D 1/2/60/30/30/30"),
}


def active_methods() -> list[str]:
    include_env = os.environ.get("SINGLE_BENCH_PLOT_METHODS", "").strip()
    if include_env:
        requested = [m.strip() for m in include_env.split(",") if m.strip()]
        return [m for m in METHOD_ORDER if m in requested]
    exclude_reverse = os.environ.get("SINGLE_BENCH_EXCLUDE_REVERSE", "0").strip().lower() in {"1", "true", "yes"}
    if exclude_reverse:
        return [m for m in METHOD_ORDER if "reverse" not in m]
    return list(METHOD_ORDER)


def prettify(name: str) -> str:
    mapping = {
        "metropolis_anneal": "Metropolis annealing",
        "random_search": "Random search",
        "nelder_mead": "Nelder-Mead",
        "cma_es": "CMA-ES",
        "finite_diff_bfgs": "FiniteDiff + BFGS",
        "finite_diff_lbfgs": "FiniteDiff + L-BFGS",
        "dto_forward_ad_lbfgs": "DtO forward AD + L-BFGS",
        "dto_forward_ad_bfgs": "DtO forward AD + BFGS",
        "dto_forward_ad_ipnewton": "DtO forward AD + IPNewton",
        "dto_reverse_ad_lbfgs": "DtO reverse AD + L-BFGS",
        "otd_reverse_ad_lbfgs": "OtD reverse AD + L-BFGS",
    }
    return mapping.get(name, name)


def short_status(status: str) -> str:
    if status.startswith("ok"):
        return ""
    if status.startswith("stopped_partial"):
        return " [partial]"
    if status.startswith("capped_eval_"):
        return f" [{status}]"
    return " [error]"


def load_objective_dir(root: Path, objective: str):
    obj_dir = root / objective
    summary_frames = []
    summary_path = obj_dir / "optimization_summary.csv"
    if summary_path.exists() and summary_path.stat().st_size > 0:
        summary_frames.append(pd.read_csv(summary_path).assign(_source_priority=0))
    for shard in sorted(obj_dir.glob("optimization_summary__*.csv")):
        if shard.stat().st_size > 0:
            summary_frames.append(pd.read_csv(shard).assign(_source_priority=1))
    if not summary_frames:
        raise FileNotFoundError(f"No optimization summary files found in {obj_dir}")
    summary = pd.concat(summary_frames, ignore_index=True)
    summary = summary.sort_values(["method", "_source_priority"]).drop_duplicates(subset=["method"], keep="last").drop(columns=["_source_priority"])

    shard_methods = set()
    for shard in obj_dir.glob("optimization_traces__*.csv"):
        stem = shard.stem
        if "__" in stem:
            shard_methods.add(stem.split("__", 1)[1])
    trace_frames = []
    trace_path = obj_dir / "optimization_traces.csv"
    if trace_path.exists() and trace_path.stat().st_size > 0:
        try:
            combined = pd.read_csv(trace_path)
            if shard_methods and "method" in combined.columns:
                combined = combined[~combined["method"].astype(str).str.lower().str.replace(r"[^a-z0-9]+", "_", regex=True).isin(shard_methods)]
            trace_frames.append(combined)
        except EmptyDataError:
            pass
    for shard in sorted(obj_dir.glob("optimization_traces__*.csv")):
        if shard.stat().st_size > 0:
            try:
                trace_frames.append(pd.read_csv(shard))
            except EmptyDataError:
                pass
    if trace_frames:
        traces = pd.concat(trace_frames, ignore_index=True)
    else:
        traces = pd.DataFrame(columns=["objective", "method", "eval", "elapsed_s", "loss"])
    meta = json.loads((obj_dir / "optimization_meta.json").read_text())
    return obj_dir, traces, summary, meta


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    root = Path(
        os.environ.get(
            "SINGLE_BENCH_RESULTS_ROOT",
            repo / "generated" / "figures" / "optimization" / "single_sample_dose_method_benchmark_20260311",
        )
    )
    obj_env = os.environ.get("SINGLE_BENCH_PLOT_OBJECTIVES", "").strip()
    if obj_env:
        objectives = [obj.strip() for obj in obj_env.split(",") if obj.strip()]
    else:
        objectives = [obj for obj in OBJECTIVES if (root / obj).exists()]
    if not objectives:
        raise FileNotFoundError(f"No benchmark objective directories found under {root}")
    loaded = {obj: load_objective_dir(root, obj) for obj in objectives}
    first_meta = loaded[objectives[0]][3]
    reference_points = None
    reference_path = root / "reference_regimen_points.csv"
    if reference_path.exists() and reference_path.stat().st_size > 0:
        reference_points = pd.read_csv(reference_path)
    sample_id = first_meta["sample_id"]
    solver = first_meta["solver"]
    scenario = first_meta["scenario_label"]
    callback_mode = first_meta["callback_mode"]
    log_y = os.environ.get("SINGLE_BENCH_LOG_Y", "0").strip().lower() in {"1", "true", "yes"}
    add_optima_column = os.environ.get("SINGLE_BENCH_ADD_OPTIMA_COLUMN", "0").strip().lower() in {"1", "true", "yes"}
    max_eval_x_env = os.environ.get("SINGLE_BENCH_MAX_EVAL_X", "").strip()
    max_time_x_env = os.environ.get("SINGLE_BENCH_MAX_TIME_X", "").strip()
    max_eval_x = float(max_eval_x_env) if max_eval_x_env else None
    max_time_x = float(max_time_x_env) if max_time_x_env else None
    suffix = "_logy" if log_y else ""
    out_path = Path(
        os.environ.get(
            "SINGLE_BENCH_FIG_OUT",
            root / f"single_sample_dose_method_benchmark_grid_sample{sample_id}{suffix}.png",
        )
    )

    methods = active_methods()
    plt.style.use("default")
    ncols = len(objectives) + (1 if add_optima_column else 0)
    fig_width = 6.2 * len(objectives) + (5.8 if add_optima_column else 0)
    fig, axes = plt.subplots(2, ncols, figsize=(fig_width, 9.5), facecolor="white")
    if ncols == 1:
        axes = np.array(axes).reshape(2, 1)
    best_dose_notes: list[tuple[int, str]] = []

    for col, objective in enumerate(objectives):
        _, traces, summary, meta = loaded[objective]
        ax_eval = axes[0, col]
        ax_time = axes[1, col]
        eval_xmax = 0.0
        time_xmax = 0.0

        best_map = dict(zip(summary["method"], summary["best_loss"]))
        status_map = dict(zip(summary["method"], summary["status"]))
        label = meta["objective_label"]

        active_summary = summary[summary["method"].isin(methods)].copy()
        best_row = active_summary.sort_values("best_loss").iloc[0]
        dose_note = (
            f"Best overall: {prettify(best_row['method'])}\n"
            f"C1D1 {best_row['dose1_mg']:.2f} mg, C1D8 {best_row['dose2_mg']:.2f} mg\n"
            f"C1D15/C2D1 {best_row['dose3_mg']:.2f} mg, C3+ q3w {best_row['dose4_mg']:.2f} mg"
        )
        best_dose_notes.append((col, dose_note))

        for method in methods:
            sub = traces[traces["method"] == method].sort_values("eval").copy()
            if sub.empty:
                continue
            sub["best_loss_so_far"] = sub["loss"].cummin()
            best = float(best_map.get(method, sub["best_loss_so_far"].iloc[-1]))
            status = str(status_map.get(method, ""))
            suffix = short_status(status)
            plot_label = f"{prettify(method)}{suffix}"
            style = METHOD_STYLES[method]
            markevery = max(1, len(sub) // 10)
            common = dict(
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=2.6,
                marker=style["marker"],
                markersize=5.5,
                markerfacecolor="white",
                markeredgewidth=1.2,
                markevery=markevery,
                alpha=0.96,
                label=plot_label,
                path_effects=[pe.Stroke(linewidth=4.2, foreground="white"), pe.Normal()],
            )
            ax_eval.plot(sub["eval"], sub["best_loss_so_far"], **common)
            ax_time.plot(sub["elapsed_s"], sub["best_loss_so_far"], **common)
            ax_eval.scatter(
                sub["eval"].iloc[-1],
                sub["best_loss_so_far"].iloc[-1],
                s=42,
                color=style["color"],
                edgecolors="white",
                linewidths=1.0,
                zorder=5,
            )
            ax_time.scatter(
                sub["elapsed_s"].iloc[-1],
                sub["best_loss_so_far"].iloc[-1],
                s=42,
                color=style["color"],
                edgecolors="white",
                linewidths=1.0,
                zorder=5,
            )
            if method not in {"finite_diff_lbfgs", "finite_diff_bfgs"}:
                eval_xmax = max(eval_xmax, float(sub["eval"].max()))
                time_xmax = max(time_xmax, float(sub["elapsed_s"].max()))

        ax_eval.set_title(label)
        ax_eval.set_xlabel("Objective evaluations")
        ax_eval.set_ylabel("Best loss so far")
        ax_eval.grid(False)
        ax_time.set_xlabel("Wallclock time (s)")
        ax_time.set_ylabel("Best loss so far")
        ax_time.grid(False)
        for ax in (ax_eval, ax_time):
            ax.set_facecolor("white")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(colors="#374151", labelsize=9)
            ax.title.set_fontsize(11)
            ax.xaxis.label.set_size(10)
            ax.yaxis.label.set_size(10)
            if log_y:
                ax.set_yscale("log")
        if eval_xmax > 0:
            ax_eval.set_xlim(0, eval_xmax * 1.05)
        if time_xmax > 0:
            ax_time.set_xlim(0, time_xmax * 1.05)
        if max_eval_x is not None:
            ax_eval.set_xlim(0, max_eval_x)
        if max_time_x is not None:
            ax_time.set_xlim(0, max_time_x)

    if add_optima_column:
        ax_scatter = axes[0, -1]
        ax_sched = axes[1, -1]
        dose_days = np.array(first_meta["dose_times_days"], dtype=float)
        max_dose = 60.0
        for objective in objectives:
            _, _, summary, _ = loaded[objective]
            if "tumor_auc" not in summary.columns or "peak_il6" not in summary.columns:
                continue
            ostyle = OBJECTIVE_SUMMARY_STYLES.get(objective, dict(marker="o", linestyle="-"))
            active_summary = summary[summary["method"].isin(methods)].copy()
            for _, row in active_summary.iterrows():
                method = row["method"]
                style = METHOD_STYLES.get(method, dict(color="#374151"))
                x_peak = float(row["peak_il6"])
                y_tumor = float(row["tumor_auc"])
                if np.isfinite(x_peak) and x_peak > 0 and np.isfinite(y_tumor) and y_tumor > 0:
                    ax_scatter.scatter(
                        x_peak,
                        y_tumor,
                        color=style["color"],
                        marker=ostyle["marker"],
                        s=54,
                        alpha=0.95,
                        edgecolors="white",
                        linewidths=0.9,
                        zorder=4,
                    )
                doses = np.array(
                    [
                        float(row["dose1_mg"]),
                        float(row["dose2_mg"]),
                        float(row["dose3_mg"]),
                        float(row["dose4_mg"]),
                        float(row["dose4_mg"]),
                        float(row["dose4_mg"]),
                    ],
                    dtype=float,
                )
                max_dose = max(max_dose, float(np.nanmax(doses)))
                ax_sched.plot(
                    dose_days,
                    doses,
                    color=style["color"],
                    linestyle=ostyle["linestyle"],
                    marker=ostyle["marker"],
                    markersize=5.5,
                    linewidth=2.0,
                    alpha=0.85,
                    markerfacecolor="white",
                    markeredgewidth=1.0,
                    path_effects=[pe.Stroke(linewidth=3.2, foreground="white"), pe.Normal()],
                )

        ax_scatter.set_title("Optima in Endpoint Space")
        ax_scatter.set_xlabel("Peak IL6 magnitude")
        ax_scatter.set_ylabel("BTumor AUC")
        ax_scatter.set_xscale("log")
        ax_scatter.set_yscale("log")
        ax_scatter.grid(False)

        ax_sched.set_title("Optimal Dose Schedules")
        ax_sched.set_xlabel("Dose day")
        ax_sched.set_ylabel("Dose (mg)")
        ax_sched.set_ylim(0, max(60.0, max_dose * 1.05))
        ax_sched.set_xlim(float(dose_days.min()) - 1.0, float(dose_days.max()) + 1.0)
        ax_sched.set_xticks(dose_days)
        ax_sched.grid(False)

        for ax in (ax_scatter, ax_sched):
            ax.set_facecolor("white")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(colors="#374151", labelsize=9)
            ax.title.set_fontsize(11)
            ax.xaxis.label.set_size(10)
            ax.yaxis.label.set_size(10)

        if reference_points is not None and not reference_points.empty:
            for _, row in reference_points.iterrows():
                ref_name = str(row["reference"])
                style = REFERENCE_STYLES.get(ref_name)
                if style is None:
                    continue
                x_peak = float(row["peak_il6"])
                y_tumor = float(row["tumor_auc"])
                if np.isfinite(x_peak) and x_peak > 0 and np.isfinite(y_tumor) and y_tumor > 0:
                    ax_scatter.scatter(
                        x_peak,
                        y_tumor,
                        color=style["color"],
                        marker=style["marker"],
                        s=120,
                        alpha=0.95,
                        edgecolors="white",
                        linewidths=1.2,
                        zorder=6,
                    )
                doses = np.array(
                    [
                        float(row["dose1_mg"]),
                        float(row["dose2_mg"]),
                        float(row["dose3_mg"]),
                        float(row["dose4_mg"]),
                        float(row["dose4_mg"]),
                        float(row["dose4_mg"]),
                    ],
                    dtype=float,
                )
                ax_sched.plot(
                    dose_days,
                    doses,
                    color=style["color"],
                    linestyle=style["linestyle"],
                    marker=style["marker"],
                    markersize=7.0,
                    linewidth=2.8,
                    alpha=0.95,
                    markerfacecolor="white",
                    markeredgewidth=1.2,
                    path_effects=[pe.Stroke(linewidth=4.0, foreground="white"), pe.Normal()],
                    zorder=5,
                )

    sample_prefix = "Randomly selected sample" if first_meta.get("sample_selection_mode", "random") == "random" else "Sample"
    fig.suptitle(
        f"{scenario}\n{sample_prefix} {sample_id}; {solver}; callback dosing via {callback_mode}; reverse methods omitted" +
        ("; log-y scale" if log_y else ""),
        y=0.98,
        fontsize=13,
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.01, 0.62), frameon=False, fontsize=9, title="Methods")
    if add_optima_column:
        objective_handles = [
            Line2D(
                [0],
                [0],
                color="#4b5563",
                marker=OBJECTIVE_SUMMARY_STYLES[obj]["marker"],
                linestyle=OBJECTIVE_SUMMARY_STYLES[obj]["linestyle"],
                linewidth=2.0,
                markersize=6,
                markerfacecolor="white",
                markeredgewidth=1.0,
            )
            for obj in objectives
        ]
        objective_labels = [loaded[obj][3]["objective_label"] for obj in objectives]
        fig.legend(
            objective_handles,
            objective_labels,
            loc="center left",
            bbox_to_anchor=(1.01, 0.37),
            frameon=False,
            fontsize=9,
            title="Objectives",
        )
        if reference_points is not None and not reference_points.empty:
            reference_handles = [
                Line2D(
                    [0],
                    [0],
                    color=style["color"],
                    marker=style["marker"],
                    linestyle=style["linestyle"],
                    linewidth=2.4,
                    markersize=7,
                    markerfacecolor="white",
                    markeredgewidth=1.1,
                )
                for style in REFERENCE_STYLES.values()
            ]
            reference_labels = [style["label"] for style in REFERENCE_STYLES.values()]
            fig.legend(
                reference_handles,
                reference_labels,
                loc="center left",
                bbox_to_anchor=(1.01, 0.18),
                frameon=False,
                fontsize=9,
                title="References",
            )
    fig.tight_layout(rect=(0, 0.13, 0.83, 0.94))
    fig.subplots_adjust(wspace=0.25, hspace=0.28)
    for col, note in best_dose_notes:
        bbox = axes[1, col].get_position()
        xmid = (bbox.x0 + bbox.x1) / 2
        fig.text(
            xmid,
            0.045,
            note,
            ha="center",
            va="bottom",
            fontsize=9,
            color="#1f2937",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#f9fafb", edgecolor="#d1d5db"),
        )
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(out_path)


if __name__ == "__main__":
    main()
