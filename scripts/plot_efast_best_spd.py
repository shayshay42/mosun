#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[1]
    default_dir = repo / "generated" / "figures" / "sensitivity" / "efast_best_spd"
    ap = argparse.ArgumentParser(description="Plot eFAST first-order and total sensitivity indices.")
    ap.add_argument("--in-dir", type=Path, default=default_dir)
    ap.add_argument("--indices-csv", type=Path, default=None)
    ap.add_argument("--meta-json", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    in_dir = args.in_dir
    idx_path = args.indices_csv or (in_dir / "efast_best_spd_indices.csv")
    if not idx_path.is_absolute():
        idx_path = (Path.cwd() / idx_path).resolve()
    if not idx_path.is_file():
        raise FileNotFoundError(f"Missing indices csv: {idx_path}")

    out_path = args.out or (idx_path.parent / "efast_best_spd_indices.png")
    if not out_path.is_absolute():
        out_path = (Path.cwd() / out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    meta_path = args.meta_json
    if meta_path is None:
        cand = idx_path.parent / idx_path.name.replace("_indices.csv", "_meta.json")
        meta_path = cand if cand.is_file() else None

    df = pd.read_csv(idx_path).copy()
    if "rank_ST" in df.columns:
        df = df.sort_values("rank_ST", ascending=True)
    else:
        df = df.sort_values("ST", ascending=False)
    df = df.reset_index(drop=True)

    y = range(len(df))
    h = 0.38

    fig_h = max(5.5, 0.42 * len(df) + 1.8)
    fig, ax = plt.subplots(figsize=(10.0, fig_h))
    ax.barh([yy + h / 2 for yy in y], df["ST"], height=h, color="#222222", label="ST (total)")
    ax.barh([yy - h / 2 for yy in y], df["S1"], height=h, color="#9a9a9a", label="S1 (first-order)")

    ax.set_yticks(list(y))
    ax.set_yticklabels(df["parameter"].tolist())
    ax.invert_yaxis()
    ax.set_xlabel("Sensitivity index")
    title = "eFAST Sensitivity"
    if meta_path is not None and meta_path.is_file():
        with meta_path.open("r") as f:
            meta = json.load(f)
        output = meta.get("output")
        horizon = meta.get("horizon_days")
        if output is not None:
            title = f"eFAST Sensitivity for {output}"
        if horizon is not None:
            title = f"{title} ({horizon:g}-day regimen)"
    ax.set_title(title)
    ax.legend(loc="lower right", frameon=False)
    ax.grid(False)

    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
