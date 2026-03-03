#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def build_maps(manifest_path: Path) -> dict[str, Path]:
    df = pd.read_csv(manifest_path)
    df = df[df["status"] == "ok"].copy()
    return {r.regimen: Path(r.trace_csv) for _, r in df.iterrows()}


def compare_mode_traces(
    reg_order: list[str],
    seg_map: dict[str, Path],
    cb_map: dict[str, Path],
    ycol: str,
    title: str,
    ylabel: str,
    out_path: Path,
) -> pd.DataFrame:
    n = len(reg_order)
    ncols = 4
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 3.2 * nrows), squeeze=False)
    axes_flat = axes.flatten()

    rows: list[dict[str, float | str]] = []

    for i, reg in enumerate(reg_order):
        ax = axes_flat[i]
        sdf = pd.read_csv(seg_map[reg])
        cdf = pd.read_csv(cb_map[reg])
        merged = sdf.merge(cdf, on=["time", "regimen", "patient_id"], suffixes=("_seg", "_cb"))

        ys = merged[f"{ycol}_seg"]
        yc = merged[f"{ycol}_cb"]
        abs_err = (yc - ys).abs().max()
        amp = max(ys.abs().max(), 1e-12)
        norm_err = abs_err / amp
        rows.append({"regimen": reg, "observable": ycol, "max_abs": float(abs_err), "max_norm_abs": float(norm_err)})

        ax.plot(merged["time"], ys, color="#D55E00", linestyle="--", linewidth=1.7, label="Canonical segmented")
        ax.plot(merged["time"], yc, color="#0072B2", linestyle=":", linewidth=1.9, label="Canonical callback")
        ax.set_title(f"{reg}\nmax|Δ|/max|seg|={norm_err:.3g}", fontsize=9)
        ax.set_xlabel("Time (days)")
        ax.grid(alpha=0.25)

    for i in range(n, len(axes_flat)):
        axes_flat[i].axis("off")

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False)
    fig.suptitle(title, y=0.996, fontsize=12)
    fig.text(0.01, 0.5, ylabel, va="center", rotation="vertical")
    fig.tight_layout(rect=[0.02, 0.0, 1.0, 0.95])
    fig.savefig(out_path, dpi=180)
    plt.close(fig)

    return pd.DataFrame(rows)


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    out_dir = repo / "generated" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    reg_df = pd.read_csv(repo / "generated" / "phase1_design_standard" / "regimen_events.csv")
    reg_order = reg_df["regimen"].drop_duplicates().tolist()

    seg_map = build_maps(repo / "generated" / "phase1_standard_julia_traces_canonical_segmented" / "manifest.csv")
    cb_map = build_maps(repo / "generated" / "phase1_standard_julia_traces_canonical" / "manifest.csv")

    missing = [r for r in reg_order if r not in seg_map or r not in cb_map]
    if missing:
        raise RuntimeError(f"Missing traces for regimens: {missing}")

    tox_fig = out_dir / "phase1_standard_toxicity_trace_overlay_canonical_segmented_vs_callback.png"
    tum_fig = out_dir / "phase1_standard_tumor_trace_overlay_canonical_segmented_vs_callback.png"

    tox_df = compare_mode_traces(
        reg_order,
        seg_map,
        cb_map,
        ycol="IL6combo",
        title="Canonical Dosing Comparison: Segmented vs PresetTimeCallback (IL6combo)",
        ylabel="IL6combo",
        out_path=tox_fig,
    )
    tum_df = compare_mode_traces(
        reg_order,
        seg_map,
        cb_map,
        ycol="Btumor",
        title="Canonical Dosing Comparison: Segmented vs PresetTimeCallback (Btumor)",
        ylabel="Btumor",
        out_path=tum_fig,
    )

    summary = pd.concat([tox_df, tum_df], ignore_index=True)
    summary_path = repo / "generated" / "phase1_standard_julia_traces_canonical" / "canonical_segmented_vs_callback_summary.csv"
    summary.to_csv(summary_path, index=False)

    print(tox_fig)
    print(tum_fig)
    print(summary_path)


if __name__ == "__main__":
    main()

