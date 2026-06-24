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
from matplotlib.lines import Line2D


FAMILY_COLORS = {
    "tumor_burden_growth_infiltration": "#4C78A8",
    "b_cell_killing": "#D55E00",
    "IL6": "#CC79A7",
    "PK_exposure": "#009E73",
    "t_cell_activation_infiltration": "#7F3C8D",
    "other_model": "#6E6E6E",
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_path(raw: str | Path) -> Path:
    path = Path(raw)
    return path if path.is_absolute() else repo_root() / path


def default_root() -> Path:
    return repo_root() / "generated/figures/vpop_pruning/susilo_efast_hosseini2020_fig5_vpop250_20260505"


def fs_path(path: Path) -> str:
    resolved = str(path.resolve())
    if os.name == "nt" and not resolved.startswith("\\\\?\\"):
        return "\\\\?\\" + resolved
    return resolved


def main() -> None:
    root = resolve_path(os.environ.get("HOSSEINI_EFAST_VPOP250_ROOT", str(default_root())))
    in_csv = root / "high_importance_parameters.csv"
    out_base = root / "hosseini2020_fig5_vpop250_efast_ranked_high_importance_parameters"
    if not in_csv.exists():
        raise FileNotFoundError(in_csv)

    df = pd.read_csv(in_csv)
    df = df[df["selected_high_importance"].astype(bool)].sort_values("importance_rank").copy()
    if len(df) != 33:
        raise RuntimeError(f"Expected 33 selected high-importance parameters, found {len(df)}")

    df["label"] = df["importance_rank"].astype(int).astype(str) + ". " + df["parameter"].astype(str)
    df_plot = df.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(df_plot))
    colors = [FAMILY_COLORS.get(str(fam), "#999999") for fam in df_plot["parameter_family"]]

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(13.2, max(8.0, 0.28 * len(df_plot) + 1.8)),
        gridspec_kw={"width_ratios": [1.15, 0.85]},
        sharey=True,
    )

    ax = axes[0]
    ax.barh(y, df_plot["full_n_endpoints"], color=colors, alpha=0.88, height=0.72)
    ax.scatter(
        df_plot["ckpt_n_endpoints"],
        y,
        s=34,
        color="#111111",
        marker="D",
        label="checkpointed endpoints",
        zorder=4,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(df_plot["label"], fontsize=8)
    ax.set_xlabel("Number of significant endpoints")
    ax.set_title("A. eFAST-selected parameters by endpoint coverage", loc="left", fontsize=11)
    ax.grid(True, axis="x", alpha=0.22)
    ax.set_xlim(0, max(float(df_plot["full_n_endpoints"].max()) * 1.08, 5.0))
    ax.legend(frameon=False, loc="lower right", fontsize=8)

    ax = axes[1]
    ax.scatter(df_plot["full_max_ST"], y, s=46, color="#4C78A8", label="full screen max ST")
    ckpt = pd.to_numeric(df_plot["ckpt_max_ST"], errors="coerce")
    mask = ckpt > 0
    ax.scatter(ckpt[mask], y[mask], s=36, color="#D55E00", marker="s", label="checkpoint max ST")
    for idx, row in df_plot.iterrows():
        if bool(row["forced_susilo_il6"]):
            ax.text(1.015, idx, "forced", va="center", ha="left", fontsize=7, color="#555555")
    ax.set_xlabel("Maximum total-order index, ST")
    ax.set_title("B. Total-order sensitivity strength", loc="left", fontsize=11)
    ax.set_xlim(0.78, 1.07)
    ax.grid(True, axis="x", alpha=0.22)
    ax.legend(frameon=False, loc="lower left", fontsize=8)

    family_handles = [
        Line2D([0], [0], marker="s", linestyle="", color=color, label=family.replace("_", " "))
        for family, color in FAMILY_COLORS.items()
        if family in set(df["parameter_family"].astype(str))
    ]
    fig.legend(handles=family_handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.01), fontsize=8)
    fig.suptitle(
        "Hosseini Fig. 5 VPop250 eFAST high-importance parameter set",
        y=0.995,
        fontsize=13,
    )
    fig.text(
        0.01,
        0.012,
        "Rank order follows the VPop250 prefilter rule: checkpointed significant endpoint count, full-screen significant endpoint count, then max total-order ST.",
        fontsize=8,
        color="#333333",
    )
    fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.965))

    outputs = []
    for ext in [".png", ".pdf"]:
        out = out_base.with_suffix(ext)
        fig.savefig(fs_path(out), dpi=220, bbox_inches="tight")
        outputs.append(out)
    plt.close(fig)

    ranked_csv = root / "hosseini2020_fig5_vpop250_efast_ranked_high_importance_parameters.csv"
    df.to_csv(ranked_csv, index=False)
    outputs.append(ranked_csv)

    meta = {
        "created_by": "scripts/plot_hosseini_efast_high_importance_parameters.py",
        "input_csv": str(in_csv),
        "output_base": str(out_base),
        "n_parameters": int(len(df)),
        "rank_rule": [
            "checkpointed significant endpoint count descending",
            "full-screen significant endpoint count descending",
            "checkpoint max total-order ST descending",
            "full-screen max total-order ST descending",
            "parameter name ascending",
        ],
        "outputs": [str(p) for p in outputs],
    }
    meta_path = root / "hosseini2020_fig5_vpop250_efast_ranked_high_importance_parameters_meta.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    outputs.append(meta_path)

    print(out_base.with_suffix(".png"))


if __name__ == "__main__":
    main()
