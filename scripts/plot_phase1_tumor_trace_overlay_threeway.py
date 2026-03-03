#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    out_dir = repo / "generated" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    reg_df = pd.read_csv(repo / "generated" / "phase1_design_standard" / "regimen_events.csv")
    reg_order = reg_df["regimen"].drop_duplicates().tolist()

    m_manifest = pd.read_csv(repo / "generated" / "phase1_standard_matlab_traces" / "manifest.csv")
    l_manifest = pd.read_csv(repo / "generated" / "phase1_standard_julia_traces" / "manifest.csv")
    c_manifest = pd.read_csv(repo / "generated" / "phase1_standard_julia_traces_canonical" / "manifest.csv")

    m_manifest = m_manifest[m_manifest["status"] == "ok"].copy()
    l_manifest = l_manifest[l_manifest["status"] == "ok"].copy()
    c_manifest = c_manifest[c_manifest["status"] == "ok"].copy()

    m_map = {r.regimen: Path(r.trace_csv) for _, r in m_manifest.iterrows()}
    l_map = {r.regimen: Path(r.trace_csv) for _, r in l_manifest.iterrows()}
    c_map = {r.regimen: Path(r.trace_csv) for _, r in c_manifest.iterrows()}

    missing = [r for r in reg_order if r not in m_map or r not in l_map or r not in c_map]
    if missing:
        raise RuntimeError(f"Missing traces for regimens: {missing}")

    n = len(reg_order)
    ncols = 4
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 3.2 * nrows), squeeze=False)
    axes_flat = axes.flatten()

    for i, reg in enumerate(reg_order):
        ax = axes_flat[i]
        mdf = pd.read_csv(m_map[reg])[["time", "Btumor"]].rename(columns={"Btumor": "Btumor_matlab"})
        ldf = pd.read_csv(l_map[reg])[["time", "Btumor"]].rename(columns={"Btumor": "Btumor_legacy"})
        cdf = pd.read_csv(c_map[reg])[["time", "Btumor"]].rename(columns={"Btumor": "Btumor_canonical"})

        merged = mdf.merge(ldf, on="time").merge(cdf, on="time")

        ym = merged["Btumor_matlab"]
        yl = merged["Btumor_legacy"]
        yc = merged["Btumor_canonical"]

        amp = max(float(ym.abs().max()), 1e-12)
        err_l = float((yl - ym).abs().max() / amp)
        err_c = float((yc - ym).abs().max() / amp)

        ax.plot(merged["time"], ym, color="black", linewidth=1.7, label="MATLAB")
        ax.plot(merged["time"], yl, color="#D55E00", linestyle="--", linewidth=1.6, label="Julia legacy")
        ax.plot(merged["time"], yc, color="#0072B2", linestyle=":", linewidth=1.9, label="Julia canonical")
        ax.set_title(f"{reg}\nlegacy={err_l:.3g}, canonical={err_c:.3g}", fontsize=9)
        ax.set_xlabel("Time (days)")
        ax.grid(alpha=0.25)

    for i in range(n, len(axes_flat)):
        axes_flat[i].axis("off")

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)
    fig.suptitle("Phase-1 Tumor Trace Overlay: MATLAB vs Julia Legacy vs Julia Canonical", y=0.996, fontsize=12)
    fig.text(0.01, 0.5, "Btumor", va="center", rotation="vertical")
    fig.tight_layout(rect=[0.02, 0.0, 1.0, 0.95])

    out_path = out_dir / "phase1_standard_tumor_trace_overlay_matlab_legacy_canonical_all_regimens.png"
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(out_path)


if __name__ == "__main__":
    main()

