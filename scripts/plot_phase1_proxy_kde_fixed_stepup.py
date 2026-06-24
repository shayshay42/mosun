#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import os


def plot_kde(ax: plt.Axes, x: np.ndarray, label: str, color: str) -> None:
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return
    if len(x) < 3 or np.allclose(np.std(x), 0.0):
        ax.axvline(float(np.median(x)), color=color, linewidth=1.8, label=label)
        return
    try:
        lo = float(np.min(x))
        hi = float(np.max(x))
        span = max(hi - lo, 1e-12)
        grid = np.linspace(lo - 0.05 * span, hi + 0.05 * span, 300)
        kde = gaussian_kde(x)
        ax.plot(grid, kde(grid), color=color, linewidth=1.8, label=label)
    except Exception:
        ax.axvline(float(np.median(x)), color=color, linewidth=1.8, label=label)


def draw_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    regimens: list[str],
    value_col: str,
    legend_map: dict[str, str],
    title: str,
    xlabel: str,
    cmap_name: str,
) -> None:
    cmap = plt.colormaps.get_cmap(cmap_name)
    denom = max(len(regimens) - 1, 1)
    for i, reg in enumerate(regimens):
        vals = df.loc[df["regimen"] == reg, value_col].to_numpy()
        plot_kde(ax, vals, legend_map.get(reg, reg), color=cmap(i / denom))
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("density")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=6, frameon=False)


def _robust_scale(x: np.ndarray, *, prefer_positive: bool = False) -> float:
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return 1.0
    x = np.maximum(x, 0.0)
    med = float(np.nanmedian(x))
    if med > 1e-10 and not prefer_positive:
        return med
    pos = x[x > 1e-8]
    if len(pos):
        pos_med = float(np.nanmedian(pos))
        if pos_med > 1e-10:
            return pos_med
    q90 = float(np.nanquantile(x, 0.9))
    if q90 > 1e-10:
        return q90
    return 1.0


def resolve_scales(df: pd.DataFrame, reference_regimen: str) -> tuple[float, float, str]:
    ref = df[df["regimen"] == reference_regimen].copy()
    source = reference_regimen
    if len(ref) == 0:
        # Fallback to the labeled standard step-up regimen if present.
        fallback = "step_0.8_2_4.2mg"
        ref = df[df["regimen"] == fallback].copy()
        source = fallback if len(ref) else "all_regimens"
    if len(ref) == 0:
        ref = df.copy()
    tox_scale = _robust_scale(ref["il6_peak_0_2"].to_numpy(dtype=float), prefer_positive=False)
    tumor_scale = _robust_scale(ref["tumor_resid_day42"].to_numpy(dtype=float), prefer_positive=True)
    return max(tox_scale, 1e-12), max(tumor_scale, 1e-12), source


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    metrics_default = repo / "generated" / "phase1_julia" / "phase1_metrics_julia.csv"
    metrics_path = metrics_default
    metrics_env = os.getenv("PHASE1_METRICS_PATH", "").strip()
    if metrics_env:
        metrics_path = Path(metrics_env) if Path(metrics_env).is_absolute() else (repo / metrics_env)

    design_path = repo / "generated" / "phase1_design_standard" / "regimen_events.csv"
    out_path_env = os.getenv("PHASE1_PROXY_KDE_OUT_PATH", "").strip()
    if out_path_env:
        out_path = Path(out_path_env)
        out_path = out_path if out_path.is_absolute() else (repo / out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        out_dir = repo / "generated" / "figures"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "phase1_proxy_kde_scaled_fixed_stepup_4panel.png"

    if not metrics_path.exists():
        raise FileNotFoundError(f"Missing metrics file: {metrics_path}")

    df = pd.read_csv(metrics_path)
    if "status" in df.columns:
        df = df[df["status"] == "ok"].copy()
    df["il6_peak_0_2"] = np.maximum(df["il6_peak_0_2"].to_numpy(dtype=float), 0.0)
    df["tumor_resid_day42"] = np.maximum(df["tumor_resid_day42"].to_numpy(dtype=float), 0.0)

    req_cols = {"regimen", "regimen_type", "il6_peak_0_2", "tumor_resid_day42"}
    missing = sorted(req_cols - set(df.columns))
    if missing:
        raise RuntimeError(f"Missing required columns in {metrics_path}: {missing}")

    w_tox = float(os.getenv("PHASE1_LOSS_W_TOX", "0.5"))
    w_tum = float(os.getenv("PHASE1_LOSS_W_TUMOR", "0.5"))
    ref_regimen = os.getenv("PHASE1_SCALE_REGIMEN", "fixed_0.8mg").strip()
    tox_scale_env = os.getenv("PHASE1_TOX_SCALE", "").strip()
    tumor_scale_env = os.getenv("PHASE1_TUMOR_SCALE", "").strip()
    if tox_scale_env and tumor_scale_env:
        tox_scale = max(float(tox_scale_env), 1e-12)
        tumor_scale = max(float(tumor_scale_env), 1e-12)
        scale_source = "env_override"
    else:
        tox_scale, tumor_scale, scale_source = resolve_scales(df, ref_regimen)

    df["tox_scaled"] = df["il6_peak_0_2"] / tox_scale
    df["tumor_scaled"] = df["tumor_resid_day42"] / tumor_scale
    df["loss_scaled"] = w_tox * df["tox_scaled"] + w_tum * df["tumor_scaled"]

    reg_df = pd.read_csv(design_path)
    reg_order = reg_df["regimen"].drop_duplicates().tolist()
    fixed_regs = [r for r in reg_order if r.startswith("fixed_")]
    step_regs = [r for r in reg_order if r.startswith("step_")]

    dff = df[df["regimen_type"] == "fixed_q3w"].copy()
    dfs = df[df["regimen_type"] == "stepup_q3w"].copy()

    summary = (
        df.groupby("regimen", as_index=False)[["tox_scaled", "tumor_scaled", "loss_scaled"]]
        .mean(numeric_only=True)
        .sort_values("regimen")
    )
    summary_path_env = os.getenv("PHASE1_PROXY_KDE_SUMMARY_PATH", "").strip()
    if summary_path_env:
        summary_path = Path(summary_path_env)
        summary_path = summary_path if summary_path.is_absolute() else (repo / summary_path)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(summary_path, index=False)

    meta_path_env = os.getenv("PHASE1_PROXY_KDE_META_PATH", "").strip()
    if meta_path_env:
        import json

        meta_path = Path(meta_path_env)
        meta_path = meta_path if meta_path.is_absolute() else (repo / meta_path)
        meta_path.parent.mkdir(parents=True, exist_ok=True)
        meta_path.write_text(
            json.dumps(
                {
                    "created_by": "scripts/plot_phase1_proxy_kde_fixed_stepup.py",
                    "metrics_csv": str(metrics_path),
                    "output_png": str(out_path),
                    "summary_csv": str(summary_path) if summary_path_env else None,
                    "loss_weight_tox": w_tox,
                    "loss_weight_tumor": w_tum,
                    "tox_scale": tox_scale,
                    "tumor_scale": tumor_scale,
                    "scale_source": scale_source,
                    "n_rows": int(len(df)),
                    "n_regimens": int(df["regimen"].nunique()),
                },
                indent=2,
            ),
            encoding="utf-8",
        )
    legend_map = {
        row.regimen: (
            f"{row.regimen} "
            f"(tox={row.tox_scaled:.3f}, tum={row.tumor_scaled:.3f}, L={row.loss_scaled:.3f})"
        )
        for row in summary.itertuples(index=False)
    }

    fig, axes = plt.subplots(2, 2, figsize=(16, 10), squeeze=False)

    draw_panel(
        axes[0, 0],
        dff,
        fixed_regs,
        value_col="tox_scaled",
        legend_map=legend_map,
        title="Scaled Toxicity Term (IL6 day0-2 / tox_scale) - Fixed q3w",
        xlabel="toxicity term",
        cmap_name="tab10",
    )
    draw_panel(
        axes[0, 1],
        dfs,
        step_regs,
        value_col="tox_scaled",
        legend_map=legend_map,
        title="Scaled Toxicity Term (IL6 day0-2 / tox_scale) - Step-up q3w",
        xlabel="toxicity term",
        cmap_name="tab20",
    )
    draw_panel(
        axes[1, 0],
        dff,
        fixed_regs,
        value_col="tumor_scaled",
        legend_map=legend_map,
        title="Scaled Tumor Term (Btumor day42/day0 / tumor_scale) - Fixed q3w",
        xlabel="tumor term",
        cmap_name="tab10",
    )
    draw_panel(
        axes[1, 1],
        dfs,
        step_regs,
        value_col="tumor_scaled",
        legend_map=legend_map,
        title="Scaled Tumor Term (Btumor day42/day0 / tumor_scale) - Step-up q3w",
        xlabel="tumor term",
        cmap_name="tab20",
    )

    fig.suptitle(
        (
            "Phase-1 Scale-Adjusted VPop Proxy Distributions by Dosing Scenario\n"
            f"loss = {w_tox:.2f}*tox_term + {w_tum:.2f}*tumor_term; "
            f"tox_scale={tox_scale:.3g}, tumor_scale={tumor_scale:.3g}, source={scale_source}"
        ),
        y=0.99,
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(out_path)


if __name__ == "__main__":
    main()
