#!/usr/bin/env python3
"""Build an eFAST-derived candidate set for Hosseini 2020 Fig. 5 VPop pruning."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
FULL_ROOT = REPO_ROOT / "generated" / "server_results" / "susilo_fig5_il6_dummy_efast_20260502" / "susilo_fig5_il6_dummy_efast_20260502"
CKPT_ROOT = REPO_ROOT / "generated" / "server_results" / "susilo_fig5_il6_day84_checkpointed_efast_20260504" / "susilo_fig5_il6_day84_checkpointed_efast_20260504"
DIGITIZED_ROOT = REPO_ROOT / "assets" / "digitization_hosseini_2020_vpop"
OUT_ROOT = Path(
    os.environ.get(
        "SUSILO_EFAST_HOSSEINI_VPOP_OUT_DIR",
        REPO_ROOT / "generated" / "figures" / "vpop_pruning" / "susilo_efast_hosseini2020_fig5_vpop250_20260505",
    )
)

N_RESIM = int(os.environ.get("SUSILO_EFAST_VPOP_N_RESIM", "20000"))
SEED = int(os.environ.get("SUSILO_EFAST_VPOP_SEED", "20260505"))
N_PREFILTER_POOL = int(os.environ.get("SUSILO_EFAST_VPOP_PREFILTER_POOL", "150000"))
N_BEST_ALWAYS = int(os.environ.get("SUSILO_EFAST_VPOP_BEST_ALWAYS", "5000"))

KEY_COLS = ["seed", "efast_eval_idx", "parameter_block", "sample_in_block"]
REGIMEN_ORDER = [
    "c1d1_20_only_then_20_q3w",
    "c1d1_d8_d15_6p7_then_20_q3w",
    "c1d1_d8_d15_1p6_10_10_then_20_q3w",
    "c1d1_d8_1p6_20_then_20_q3w",
]
REGIMEN_LABELS = {
    "c1d1_20_only_then_20_q3w": "20 mg on C1D1, then 20 mg q3w",
    "c1d1_d8_d15_6p7_then_20_q3w": "6.7/6.7/6.7 mg, then 20 mg q3w",
    "c1d1_d8_d15_1p6_10_10_then_20_q3w": "1.6/10/10 mg, then 20 mg q3w",
    "c1d1_d8_1p6_20_then_20_q3w": "1.6/20 mg, then 20 mg q3w",
}
DIGITIZED_REGIMEN_TO_CODE = {
    "20 mg C1D1; 20 mg Day 1 subsequent cycles": "c1d1_20_only_then_20_q3w",
    "6.7 mg C1D1/C1D8/C1D15; 20 mg Day 1 subsequent cycles": "c1d1_d8_d15_6p7_then_20_q3w",
    "1.6/10/10 mg C1D1/C1D8/C1D15; 20 mg Day 1 subsequent cycles": "c1d1_d8_d15_1p6_10_10_then_20_q3w",
    "1.6/20 mg C1D1/C1D8; 20 mg Day 1 subsequent cycles": "c1d1_d8_1p6_20_then_20_q3w",
}
FORCED_PARAMS = {
    "tumor_burden_factor",
    "BT_ratio_tumor_init",
    "kBtumorprolif",
    "kIL6prod",
    "thalfIL6",
    "IL6_tiss_contribution",
}


def log(msg: str) -> None:
    print(msg, flush=True)


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def candidate_id(df: pd.DataFrame) -> pd.Series:
    return (
        df["seed"].astype(str)
        + "__"
        + df["efast_eval_idx"].astype(str)
        + "__"
        + df["parameter_block"].astype(str)
        + "__"
        + df["sample_in_block"].astype(str)
    )


def high_importance_parameters() -> pd.DataFrame:
    full_sig = pd.read_csv(FULL_ROOT / "significant_parameter_screen.csv")
    ckpt_sig = pd.read_csv(CKPT_ROOT / "significant_parameter_screen.csv")
    full_sig = full_sig[full_sig["significant"].astype(bool)].copy()
    ckpt_sig = ckpt_sig[ckpt_sig["significant"].astype(bool)].copy()

    full_b = (
        full_sig.groupby(["parameter", "parameter_family"], as_index=False)
        .agg(full_n_endpoints=("endpoint", "nunique"), full_max_ST=("median_ST", "max"), full_median_ST=("median_ST", "median"))
    )
    ckpt_b = (
        ckpt_sig.groupby(["parameter", "parameter_family"], as_index=False)
        .agg(ckpt_n_endpoints=("endpoint", "nunique"), ckpt_max_ST=("median_ST", "max"), ckpt_median_ST=("median_ST", "median"))
    )
    imp = full_b.merge(ckpt_b, on=["parameter", "parameter_family"], how="outer")
    for col in ["full_n_endpoints", "full_max_ST", "full_median_ST", "ckpt_n_endpoints", "ckpt_max_ST", "ckpt_median_ST"]:
        imp[col] = imp[col].fillna(0.0)
    imp["forced_susilo_il6"] = imp["parameter"].isin(FORCED_PARAMS)
    imp["selected_high_importance"] = (
        (imp["ckpt_n_endpoints"] >= 2)
        | (imp["full_n_endpoints"] >= 30)
        | (imp["forced_susilo_il6"] & ((imp["ckpt_n_endpoints"] >= 1) | (imp["full_n_endpoints"] >= 1)))
    )
    imp = imp[imp["selected_high_importance"]].copy()
    imp = imp.sort_values(
        ["ckpt_n_endpoints", "full_n_endpoints", "ckpt_max_ST", "full_max_ST", "parameter"],
        ascending=[False, False, False, False, True],
    )
    imp.insert(0, "importance_rank", np.arange(1, len(imp) + 1))
    if len(imp) != 33:
        raise RuntimeError(f"Expected 33 high-importance parameters, found {len(imp)}")
    return imp


def digitized_scalar_targets() -> pd.DataFrame:
    traj = pd.read_csv(DIGITIZED_ROOT / "hosseini2020_fig5_digitized_IL6_Tcell_green_gray_trajectories_anchor_scaled.csv")
    wf = pd.read_csv(DIGITIZED_ROOT / "hosseini2020_fig5_digitized_tumor_waterfall_shapes.csv")
    rows = []
    for digitized_regimen, regimen in DIGITIZED_REGIMEN_TO_CODE.items():
        for biomarker in ["IL6", "Tcell"]:
            sub = traj[(traj["regimen"] == digitized_regimen) & (traj["biomarker"] == biomarker)].copy()
            for curve in ["median_green", "p5_green", "p95_green"]:
                c = sub[sub["curve"] == curve]
                if c.empty:
                    raise RuntimeError(f"Missing target {digitized_regimen}/{biomarker}/{curve}")
                if biomarker == "IL6":
                    value = float(c[c["time_day"] <= 7.0]["value"].max())
                    metric = "first_peak_il6_0_7d"
                else:
                    value = float(c[c["time_day"] <= 42.0]["value"].max())
                    metric = "tafraction_peak_0_42d"
                rows.append((regimen, biomarker, metric, curve, value))
        w = wf[wf["regimen"] == digitized_regimen].sort_values("patient_rank_fraction")
        rows.append((regimen, "tumor", "waterfall_min_day84", "digitized", float(w["tumor_change_percent_day84"].min())))
        rows.append((regimen, "tumor", "waterfall_max_day84", "digitized", float(w["tumor_change_percent_day84"].max())))
        rows.append((regimen, "tumor", "responder_gt50_day84_fraction", "digitized", float((w["tumor_change_percent_day84"] <= -50.0).mean())))
    return pd.DataFrame(rows, columns=["regimen", "biomarker", "metric", "curve", "target_value"])


def load_metric_features() -> pd.DataFrame:
    usecols = KEY_COLS + [
        "regimen",
        "status",
        "first_peak_il6_0_7d",
        "tafraction_peak_0_42d",
        "day84_tumor_size_change_pct",
        "responder_gt50_day84",
    ]
    path = FULL_ROOT / "simulation_metrics_long.csv"
    log(f"Reading scalar metrics from {path}")
    df = pd.read_csv(path, usecols=usecols)
    expected_rows = 3_809_280
    if len(df) != expected_rows:
        raise RuntimeError(f"Expected {expected_rows} scalar metric rows, found {len(df)}")
    df = df[df["regimen"].isin(["no_dose"] + REGIMEN_ORDER)].copy()
    df["candidate_id"] = candidate_id(df)
    status = df.pivot(index="candidate_id", columns="regimen", values="status")
    ok_ids = status.index[(status == "success").all(axis=1)]
    df = df[df["candidate_id"].isin(ok_ids)].copy()
    if df.empty:
        raise RuntimeError("No successful candidate rows after status filtering.")

    id_frame = df[["candidate_id"] + KEY_COLS].drop_duplicates("candidate_id").set_index("candidate_id")
    features = id_frame.copy()
    for regimen in REGIMEN_ORDER:
        sub = df[df["regimen"] == regimen].set_index("candidate_id")
        features[f"{regimen}__first_peak_il6_0_7d"] = sub["first_peak_il6_0_7d"].reindex(features.index)
        features[f"{regimen}__tafraction_peak_0_42d"] = sub["tafraction_peak_0_42d"].reindex(features.index)
        features[f"{regimen}__day84_tumor_size_change_pct"] = sub["day84_tumor_size_change_pct"].reindex(features.index)
        features[f"{regimen}__responder_gt50_day84"] = sub["responder_gt50_day84"].reindex(features.index)
    features = features.reset_index()
    return features


def add_scalar_prefilter_score(features: pd.DataFrame, targets: pd.DataFrame) -> pd.DataFrame:
    target_lookup: Dict[tuple, float] = {
        (r.regimen, r.metric, r.curve): float(r.target_value) for r in targets.itertuples()
    }
    penalties: List[np.ndarray] = []
    for regimen in REGIMEN_ORDER:
        il6 = np.log10(np.maximum(features[f"{regimen}__first_peak_il6_0_7d"].to_numpy(float), 1e-12))
        il6_med = np.log10(max(target_lookup[(regimen, "first_peak_il6_0_7d", "median_green")], 1e-12))
        il6_p5 = np.log10(max(target_lookup[(regimen, "first_peak_il6_0_7d", "p5_green")], 1e-12))
        il6_p95 = np.log10(max(target_lookup[(regimen, "first_peak_il6_0_7d", "p95_green")], 1e-12))
        il6_scale = max(abs(il6_p95 - il6_p5), 0.25)
        penalties.append(np.abs(il6 - il6_med) / il6_scale)

        tc = features[f"{regimen}__tafraction_peak_0_42d"].to_numpy(float)
        tc_med = target_lookup[(regimen, "tafraction_peak_0_42d", "median_green")]
        tc_p5 = target_lookup[(regimen, "tafraction_peak_0_42d", "p5_green")]
        tc_p95 = target_lookup[(regimen, "tafraction_peak_0_42d", "p95_green")]
        tc_scale = max(abs(tc_p95 - tc_p5), 5.0)
        penalties.append(np.abs(tc - tc_med) / tc_scale)

        tumor = features[f"{regimen}__day84_tumor_size_change_pct"].to_numpy(float)
        tmin = target_lookup[(regimen, "waterfall_min_day84", "digitized")]
        tmax = target_lookup[(regimen, "waterfall_max_day84", "digitized")]
        below = np.maximum(tmin - tumor, 0.0)
        above = np.maximum(tumor - tmax, 0.0)
        penalties.append((below + above) / 100.0)
    features["scalar_prefilter_score"] = np.nanmean(np.column_stack(penalties), axis=1)
    features = features.replace([np.inf, -np.inf], np.nan)
    features = features.dropna(subset=["scalar_prefilter_score"]).copy()
    features["mean_day84_tumor_change_pct"] = np.nanmean(
        np.column_stack([features[f"{r}__day84_tumor_size_change_pct"].to_numpy(float) for r in REGIMEN_ORDER]),
        axis=1,
    )
    features["mean_log10_il6_peak"] = np.nanmean(
        np.column_stack([np.log10(np.maximum(features[f"{r}__first_peak_il6_0_7d"].to_numpy(float), 1e-12)) for r in REGIMEN_ORDER]),
        axis=1,
    )
    return features


def load_candidate_params(parameter_names: Iterable[str]) -> pd.DataFrame:
    cols = KEY_COLS + list(dict.fromkeys(parameter_names))
    path = FULL_ROOT / "efast_sample_matrix.csv"
    log(f"Reading parameter matrix columns from {path}")
    params = pd.read_csv(path, usecols=cols)
    params["candidate_id"] = candidate_id(params)
    return params


def choose_resim_candidates(features: pd.DataFrame, params: pd.DataFrame, high_params: List[str]) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    merged = features.merge(params, on=["candidate_id"] + KEY_COLS, how="inner", validate="one_to_one")
    if len(merged) < N_RESIM:
        raise RuntimeError(f"Only {len(merged)} eligible candidates; need {N_RESIM}.")
    pool_n = min(max(N_PREFILTER_POOL, N_RESIM), len(merged))
    pool = merged.nsmallest(pool_n, "scalar_prefilter_score").copy()

    selected_ids: set[str] = set()
    best_n = min(N_BEST_ALWAYS, N_RESIM, len(pool))
    selected_ids.update(pool.nsmallest(best_n, "scalar_prefilter_score")["candidate_id"].tolist())

    # Add high-importance parameter tails to avoid selecting only median-like parameter values.
    tail_per_param = max(10, min(120, N_RESIM // (len(high_params) * 2)))
    for param in high_params:
        selected_ids.update(pool.nsmallest(tail_per_param, param)["candidate_id"].tolist())
        selected_ids.update(pool.nlargest(tail_per_param, param)["candidate_id"].tolist())

    remaining = pool[~pool["candidate_id"].isin(selected_ids)].copy()
    n_need = N_RESIM - len(selected_ids)
    if n_need > 0:
        score = remaining["scalar_prefilter_score"].to_numpy(float)
        z = (score - np.nanmin(score)) / max(np.nanstd(score), 1e-9)
        weights = np.exp(-0.75 * z)
        weights = weights / weights.sum()
        chosen_idx = rng.choice(remaining.index.to_numpy(), size=n_need, replace=False, p=weights)
        selected_ids.update(remaining.loc[chosen_idx, "candidate_id"].tolist())

    selected = merged[merged["candidate_id"].isin(selected_ids)].copy()
    if len(selected) > N_RESIM:
        selected = selected.sort_values(["scalar_prefilter_score", "candidate_id"]).head(N_RESIM).copy()
    if len(selected) != N_RESIM:
        raise RuntimeError(f"Selected {len(selected)} resim candidates, expected {N_RESIM}.")
    selected = selected.sort_values(["scalar_prefilter_score", "candidate_id"]).reset_index(drop=True)
    selected.insert(0, "resim_candidate_rank", np.arange(1, len(selected) + 1))
    return selected


def main() -> None:
    for path in [
        FULL_ROOT / "simulation_metrics_long.csv",
        FULL_ROOT / "efast_sample_matrix.csv",
        FULL_ROOT / "significant_parameter_screen.csv",
        CKPT_ROOT / "significant_parameter_screen.csv",
        DIGITIZED_ROOT / "hosseini2020_fig5_digitized_IL6_Tcell_green_gray_trajectories_anchor_scaled.csv",
        DIGITIZED_ROOT / "hosseini2020_fig5_digitized_tumor_waterfall_shapes.csv",
    ]:
        require(path)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    high = high_importance_parameters()
    high.to_csv(OUT_ROOT / "high_importance_parameters.csv", index=False)
    targets = digitized_scalar_targets()
    targets.to_csv(OUT_ROOT / "digitized_scalar_targets.csv", index=False)

    features = add_scalar_prefilter_score(load_metric_features(), targets)
    features.to_csv(OUT_ROOT / "efast_candidate_scalar_features.csv", index=False)
    universe = pd.read_csv(FULL_ROOT / "parameter_universe.csv")
    all_parameters = universe["parameter"].astype(str).tolist()
    params = load_candidate_params(all_parameters)
    selected = choose_resim_candidates(features, params, high["parameter"].tolist())

    id_cols = ["resim_candidate_rank", "candidate_id"] + KEY_COLS + [
        "scalar_prefilter_score",
        "mean_day84_tumor_change_pct",
        "mean_log10_il6_peak",
    ]
    selected[id_cols].to_csv(OUT_ROOT / "trajectory_resim_candidate_ids.csv", index=False)
    param_cols = ["resim_candidate_rank", "candidate_id"] + KEY_COLS + all_parameters
    selected[param_cols].to_csv(OUT_ROOT / "trajectory_resim_candidate_parameters.csv", index=False)

    meta = {
        "created_by": "scripts/build_susilo_efast_hosseini2020_vpop250_candidate_set.py",
        "full_efast_root": str(FULL_ROOT),
        "checkpointed_efast_root": str(CKPT_ROOT),
        "digitized_root": str(DIGITIZED_ROOT),
        "output_root": str(OUT_ROOT),
        "n_resim": N_RESIM,
        "seed": SEED,
        "n_high_importance_parameters": int(len(high)),
        "n_eligible_candidates": int(len(features)),
        "n_selected_resim_candidates": int(len(selected)),
        "prefilter": {
            "pool_size": int(min(max(N_PREFILTER_POOL, N_RESIM), len(features))),
            "n_best_always": N_BEST_ALWAYS,
            "tail_per_high_importance_parameter": max(10, min(120, N_RESIM // (len(high) * 2))),
        },
    }
    (OUT_ROOT / "vpop250_prefilter_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    log(f"Wrote eFAST-derived resimulation candidate set to {OUT_ROOT}")


if __name__ == "__main__":
    main()
