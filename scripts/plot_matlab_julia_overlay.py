#!/usr/bin/env python3
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    manifest_path = repo / "generated" / "matlab_reference" / "manifest_three_cases.csv"
    julia_dir = repo / "generated" / "julia_reference"
    out_dir = repo / "generated" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(manifest_path)
    written = []

    for _, row in manifest.iterrows():
        sim_key = str(row["sim_key"])
        case_no = int(row["case_no"])
        sid = int(row["id"])
        schedule = str(row["schedule"])
        doses = json.loads(row["selected_doses_json"])
        outvec = json.loads(row["outvec_json"])

        matlab_path = repo / str(row["sim_csv_relpath"])
        julia_path = julia_dir / f"{sim_key}.csv"

        mdf = pd.read_csv(matlab_path)
        jdf = pd.read_csv(julia_path)
        merged = mdf.merge(jdf, on="time", suffixes=("_matlab", "_julia"))

        n = len(outvec)
        ncols = 2
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(12, 3.8 * nrows), squeeze=False)
        axes_flat = axes.flatten()

        for i, name in enumerate(outvec):
            ax = axes_flat[i]
            ym = merged[f"{name}_matlab"]
            yj = merged[f"{name}_julia"]

            abs_err = (yj - ym).abs()
            denom = ym.abs().clip(lower=1e-12)
            max_abs = abs_err.max()
            max_rel = (abs_err / denom).max()

            ax.plot(merged["time"], ym, color="black", linewidth=2.0, label="MATLAB")
            ax.plot(
                merged["time"],
                yj,
                color="#D55E00",
                linestyle="--",
                linewidth=2.0,
                label="Julia",
            )
            ax.set_title(f"{name}\nmax|Δ|={max_abs:.3g}, max rel={max_rel:.3g}", fontsize=10)
            ax.set_xlabel("Time (days)")
            ax.grid(alpha=0.25)

        for i in range(n, len(axes_flat)):
            axes_flat[i].axis("off")

        handles, labels = axes_flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False)
        dose_str = ", ".join(doses)
        fig.suptitle(
            f"{sim_key}  | case={case_no} id={sid} schedule={schedule}\nDoses: {dose_str}",
            y=0.995,
            fontsize=12,
        )
        fig.tight_layout(rect=[0, 0, 1, 0.93])

        out_path = out_dir / f"{sim_key}_overlay.png"
        fig.savefig(out_path, dpi=180)
        plt.close(fig)
        written.append(out_path)

    # Summary bar chart of max errors for quick triage.
    summary = pd.read_csv(julia_dir / "parity_summary.csv")
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    ax.bar(summary["sim_key"], summary["max_abs_err"], color="#4C78A8")
    ax.set_ylabel("max absolute error")
    ax.set_xlabel("simulation key")
    ax.set_title("MATLAB vs Julia parity (max abs error)")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    summary_path = out_dir / "parity_max_abs_summary.png"
    fig.savefig(summary_path, dpi=180)
    plt.close(fig)
    written.append(summary_path)

    for p in written:
        print(p)


if __name__ == "__main__":
    main()
