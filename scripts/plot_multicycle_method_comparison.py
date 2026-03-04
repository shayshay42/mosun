#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def prettify(name: str) -> str:
    mapping = {
        "simulated_annealing": "Simulated Annealing",
        "nelder_mead_coordinate_descent": "Nelder-Mead + CD",
        "finite_diff_lbfgs": "Finite Diff + L-BFGS",
        "forward_ad_dto_lbfgs": "Forward DTO + L-BFGS",
        "forward_ad_otd_lbfgs": "Forward OTD + L-BFGS",
        "random_search": "Random Search (global+local)",
        "metropolis_anneal": "Metropolis Annealing",
        "mcmc_anneal": "MCMC Annealing",
    }
    return mapping.get(name, name)


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    in_dir = Path(
        os.environ.get(
            "OPTMC_RESULTS_DIR",
            repo / "generated" / "figures" / "optimization" / "multicycle_methods_patient1_8cycles",
        )
    )
    out_dir = Path(os.environ.get("OPTMC_FIG_OUT_DIR", in_dir))
    out_dir.mkdir(parents=True, exist_ok=True)
    # Display-only horizontal jitter to separate overlapping AD curves.
    # Backward-compatible: if OPTMC_AD_XJITTER_FRAC is unset, fall back to OPTMC_AD_JITTER_FRAC.
    ad_xjitter_frac = float(os.environ.get("OPTMC_AD_XJITTER_FRAC", os.environ.get("OPTMC_AD_JITTER_FRAC", "0.0")))
    fig_basename = os.environ.get("OPTMC_FIG_BASENAME", "multicycle_method_comparison_convergence.png")

    trace_path = in_dir / "optimization_traces.csv"
    summary_path = in_dir / "optimization_summary.csv"
    meta_path = in_dir / "optimization_meta.json"
    if not trace_path.exists() or not summary_path.exists() or not meta_path.exists():
        missing = [str(p) for p in [trace_path, summary_path, meta_path] if not p.exists()]
        raise FileNotFoundError(f"Missing required files: {missing}")

    traces = pd.read_csv(trace_path)
    summary = pd.read_csv(summary_path)
    meta = json.loads(meta_path.read_text())

    method_order = [
        "simulated_annealing",
        "random_search",
        "metropolis_anneal",
        "mcmc_anneal",
        "nelder_mead_coordinate_descent",
        "finite_diff_lbfgs",
        "forward_ad_dto_lbfgs",
        "forward_ad_otd_lbfgs",
    ]
    ad_jitter_sign = {
        "forward_ad_dto_lbfgs": -1.0,
        "forward_ad_otd_lbfgs": 1.0,
    }
    available = set(traces["method"].unique())
    method_order = [m for m in method_order if m in available]
    best_map = dict(zip(summary["method"], summary["best_loss"]))
    status_map = dict(zip(summary["method"], summary["status"]))
    t_span = max(float(traces["elapsed_s"].max() - traces["elapsed_s"].min()), 1e-12)
    e_span = max(float(traces["eval"].max() - traces["eval"].min()), 1e-12)

    plt.style.use("default")
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="white")
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(method_order), 1)))

    for i, m in enumerate(method_order):
        sub = traces[traces["method"] == m].sort_values("eval").copy()
        if sub.empty:
            continue
        sub["best_loss_so_far"] = sub["loss"].cummin()
        yvals = sub["best_loss_so_far"].to_numpy(copy=True)
        x_time = sub["elapsed_s"].to_numpy(copy=True)
        x_eval = sub["eval"].to_numpy(copy=True)
        if ad_xjitter_frac > 0.0 and m in ad_jitter_sign:
            x_time = np.maximum(0.0, x_time + ad_jitter_sign[m] * ad_xjitter_frac * t_span)
            x_eval = np.maximum(1.0, x_eval + ad_jitter_sign[m] * ad_xjitter_frac * e_span)
        best = float(best_map.get(m, sub["best_loss_so_far"].iloc[-1]))
        status = str(status_map.get(m, ""))
        suffix = "" if status.startswith("ok") else f", {status}"
        label = f"{prettify(m)} (L*={best:.4g}{suffix})"
        markevery = max(len(sub) // 18, 1)
        axes[0].plot(
            x_time,
            yvals,
            color=colors[i],
            linewidth=2.0,
            marker="o",
            markersize=3.5,
            markevery=markevery,
            label=label,
        )
        axes[1].plot(
            x_eval,
            yvals,
            color=colors[i],
            linewidth=2.0,
            marker="o",
            markersize=3.5,
            markevery=markevery,
            label=label,
        )

    axes[0].set_title("A) Wallclock Time vs Best Loss")
    axes[0].set_xlabel("Wallclock (s)")
    axes[0].set_ylabel("Best loss so far")
    axes[1].set_title("B) Objective Evaluations vs Best Loss")
    axes[1].set_xlabel("Objective evaluations")
    axes[1].set_ylabel("Best loss so far")
    for ax in axes:
        ax.grid(False)

    title = (
        f"Multi-cycle Dose Optimization Method Comparison: patient {meta['patient_id']}, "
        f"{meta['n_cycles']} cycles, {len(meta['dose_times_days'])} dose controls"
    )
    if ad_xjitter_frac > 0.0:
        title += f" (AD x-jitter={ad_xjitter_frac:.2e})"
    fig.suptitle(title, fontsize=14)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.02), frameon=False, ncol=1, fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    out_path = out_dir / fig_basename
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(out_path)


if __name__ == "__main__":
    main()
