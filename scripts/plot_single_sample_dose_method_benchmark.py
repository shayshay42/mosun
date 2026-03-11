#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError


OBJECTIVES = ["simple", "clinical", "tracking"]
METHOD_ORDER = [
    "metropolis_anneal",
    "random_search",
    "nelder_mead",
    "finite_diff_lbfgs",
    "dto_forward_ad_lbfgs",
    "dto_reverse_ad_lbfgs",
    "otd_forward_ad_lbfgs",
    "otd_reverse_ad_lbfgs",
]


def prettify(name: str) -> str:
    mapping = {
        "metropolis_anneal": "Metropolis annealing",
        "random_search": "Random search",
        "nelder_mead": "Nelder-Mead",
        "finite_diff_lbfgs": "FiniteDiff + L-BFGS",
        "dto_forward_ad_lbfgs": "DtO forward AD + L-BFGS",
        "dto_reverse_ad_lbfgs": "DtO reverse AD + L-BFGS",
        "otd_forward_ad_lbfgs": "OtD forward AD + L-BFGS",
        "otd_reverse_ad_lbfgs": "OtD reverse AD + L-BFGS",
    }
    return mapping.get(name, name)


def load_objective_dir(root: Path, objective: str):
    obj_dir = root / objective
    summary = pd.read_csv(obj_dir / "optimization_summary.csv")
    try:
        traces = pd.read_csv(obj_dir / "optimization_traces.csv")
    except EmptyDataError:
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
    out_path = Path(os.environ.get("SINGLE_BENCH_FIG_OUT", root / "single_sample_dose_method_benchmark_grid.png"))

    loaded = {obj: load_objective_dir(root, obj) for obj in OBJECTIVES}
    sample_id = loaded["simple"][3]["sample_id"]
    solver = loaded["simple"][3]["solver"]
    scenario = loaded["simple"][3]["scenario_label"]
    callback_mode = loaded["simple"][3]["callback_mode"]

    colors = {m: c for m, c in zip(METHOD_ORDER, plt.cm.tab10(np.linspace(0, 1, len(METHOD_ORDER))))}

    plt.style.use("default")
    fig, axes = plt.subplots(2, 3, figsize=(18, 8), facecolor="white")

    for col, objective in enumerate(OBJECTIVES):
        _, traces, summary, meta = loaded[objective]
        ax_eval = axes[0, col]
        ax_time = axes[1, col]

        best_map = dict(zip(summary["method"], summary["best_loss"]))
        status_map = dict(zip(summary["method"], summary["status"]))
        label = meta["objective_label"]

        for method in METHOD_ORDER:
            sub = traces[traces["method"] == method].sort_values("eval").copy()
            if sub.empty:
                continue
            sub["best_loss_so_far"] = sub["loss"].cummin()
            best = float(best_map.get(method, sub["best_loss_so_far"].iloc[-1]))
            status = str(status_map.get(method, ""))
            suffix = "" if status.startswith("ok") else f" [{status}]"
            plot_label = f"{prettify(method)} (L*={best:.4g}){suffix}"
            ax_eval.plot(sub["eval"], sub["best_loss_so_far"], color=colors[method], linewidth=2, marker="o", markersize=3, label=plot_label)
            ax_time.plot(sub["elapsed_s"], sub["best_loss_so_far"], color=colors[method], linewidth=2, marker="o", markersize=3, label=plot_label)

        ax_eval.set_title(label)
        ax_eval.set_xlabel("Objective evaluations")
        ax_eval.set_ylabel("Best loss so far")
        ax_eval.grid(False)
        ax_time.set_xlabel("Wallclock time (s)")
        ax_time.set_ylabel("Best loss so far")
        ax_time.grid(False)

    fig.suptitle(
        f"{scenario}\nRandomly selected sample {sample_id}; {solver}; callback dosing via {callback_mode}",
        y=0.98,
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0, 0.84, 0.94))
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(out_path)


if __name__ == "__main__":
    main()
