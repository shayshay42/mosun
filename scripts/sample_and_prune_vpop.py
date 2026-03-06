#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class CalibrationTarget:
    target_id: str
    evaluation_mode: str
    feature: str
    regimen: str | None = None
    target_value: float | None = None
    tolerance: float = 0.05
    weight: float = 1.0
    regimen_low: str | None = None
    regimen_high: str | None = None
    min_difference: float | None = None
    max_difference: float | None = None


def normalize_regimen_name(name: str) -> str:
    return str(name).replace(".", "_")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def merge_plan_with_standard_regimens(
    plan: dict[str, Any],
    standard_regimens_csv: Path,
    include_standard: bool,
) -> dict[str, Any]:
    if not include_standard:
        return plan
    if not standard_regimens_csv.is_file():
        return plan

    merged = json.loads(json.dumps(plan))
    regimens = dict(merged.get("regimens", {}))
    std = pd.read_csv(standard_regimens_csv).copy()
    std["regimen_norm"] = std["regimen"].map(normalize_regimen_name)

    for reg, sub in std.groupby("regimen_norm", sort=False):
        if reg in regimens:
            continue
        sub = sub.sort_values("event_idx")
        regimens[reg] = {
            "regimen_type": str(sub["regimen_type"].iloc[0]),
            "time_day": [float(x) for x in sub["time_day"].tolist()],
            "dose_mg": [float(x) for x in sub["dose_mg"].tolist()],
            "clinical_context": "phase1_standard_regimen",
        }
    merged["regimens"] = regimens
    return merged


def regimen_events_from_plan(plan: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for reg_name_raw, reg in plan["regimens"].items():
        reg_name = normalize_regimen_name(reg_name_raw)
        times = reg["time_day"]
        doses = reg["dose_mg"]
        if len(times) != len(doses):
            raise ValueError(f"Regimen {reg_name} has mismatched time/day and dose lengths.")
        for i, (t, d) in enumerate(zip(times, doses), start=1):
            rows.append(
                {
                    "regimen": reg_name,
                    "regimen_type": reg.get("regimen_type", "custom"),
                    "event_idx": i,
                    "time_day": float(t),
                    "dose_mg": float(d),
                }
            )
    return pd.DataFrame(rows)


def load_param_defaults(params_xlsx: Path) -> tuple[pd.DataFrame, dict[str, float], float, float, float, float, float]:
    ptab = pd.read_excel(params_xlsx, sheet_name="Sheet1")
    dlbcl = ptab[ptab["DLBCL"].notna()][["NAME", "DLBCL"]].copy()
    dlbcl.columns = ["name", "value"]
    dmap = {str(r["name"]): float(r["value"]) for _, r in dlbcl.iterrows()}

    # Some parameter sheets only provide *_ref values for DLBCL; use them as bo fallbacks.
    base_bpbo = float(dmap.get("Bpbo_perml", dmap.get("Bpbref_perml", 250_000.0)))
    base_bpbref = float(dmap.get("Bpbref_perml", base_bpbo))
    base_trpbo = float(dmap.get("Trpbo_perml", dmap.get("Trpbref_perml", 500_000.0)))
    base_trpbref = float(dmap.get("Trpbref_perml", base_trpbo))
    base_btumor_perml = float(dmap.get("Btumor_perml", 3.25e9))
    return dlbcl, dmap, base_bpbo, base_bpbref, base_trpbo, base_trpbref, base_btumor_perml


def get_hard_prior_bounds(targets_json: dict[str, Any], spread_mode: str) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    k_bounds = None
    bt_bounds = None
    spread_primary = None
    spread_sens = None
    for t in targets_json.get("hard_priors", []):
        tid = t.get("id", "")
        if tid == "tumor_bcell_proliferation_rate":
            b = t["bounds"]
            k_bounds = (float(b["lower"]), float(b["upper"]))
        elif tid == "tumor_b_to_t_ratio":
            b = t["bounds"]
            bt_bounds = (float(b["lower"]), float(b["upper"]))
        elif tid == "baseline_tumor_bcell_from_spd_primary":
            b = t["bounds_relative_factor"]
            spread_primary = (float(b["lower"]), float(b["upper"]))
        elif tid == "baseline_tumor_bcell_from_spd_alternative":
            b = t["bounds_relative_factor"]
            spread_sens = (float(b["lower"]), float(b["upper"]))

    if k_bounds is None or bt_bounds is None:
        raise ValueError("Missing required hard prior bounds in targets JSON.")
    if spread_primary is None:
        spread_primary = (0.9, 1.1)
    if spread_sens is None:
        spread_sens = (0.5, 1.5)

    spread_bounds = spread_primary if spread_mode == "primary" else spread_sens
    return k_bounds, bt_bounds, spread_bounds


def parse_extra_param_list(arg_value: str, default_names: list[str]) -> list[str]:
    txt = str(arg_value or "").strip()
    if not txt:
        return default_names
    return [x.strip() for x in txt.split(",") if x.strip()]


def load_extra_param_priors(priors_csv: Path, names: list[str]) -> dict[str, tuple[float, float]]:
    if not names or not priors_csv.is_file():
        return {}
    tab = pd.read_csv(priors_csv).copy()
    tab["name"] = tab["name"].astype(str)
    out: dict[str, tuple[float, float]] = {}
    for n in names:
        sub = tab[tab["name"] == n]
        if sub.empty:
            continue
        lb = float(sub["cyno_lb"].iloc[0])
        ub = float(sub["cyno_ub"].iloc[0])
        if np.isfinite(lb) and np.isfinite(ub) and ub > lb:
            out[n] = (lb, ub)
    return out


def sample_from_bounds(rng: np.random.Generator, lb: float, ub: float) -> float:
    if lb > 0 and ub / lb >= 1.5:
        return float(np.exp(rng.uniform(np.log(lb), np.log(ub))))
    return float(rng.uniform(lb, ub))


def sample_candidates(
    n_candidates: int,
    seed: int,
    base_bpbo: float,
    base_bpbref: float,
    base_trpbo: float,
    base_trpbref: float,
    base_btumor_perml: float,
    k_bounds: tuple[float, float],
    bt_bounds: tuple[float, float],
    spread_bounds: tuple[float, float],
    extra_priors: dict[str, tuple[float, float]],
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, float | int]] = []

    for i in range(1, n_candidates + 1):
        bt_ratio = float(np.exp(rng.uniform(np.log(bt_bounds[0]), np.log(bt_bounds[1]))))
        k_btumor = float(rng.uniform(k_bounds[0], k_bounds[1]))
        burden_factor = float(rng.uniform(spread_bounds[0], spread_bounds[1]))
        btumor_perml = float(base_btumor_perml * burden_factor)

        kbptumor = float(btumor_perml / base_bpbo)
        ktrptumor = float(btumor_perml / (bt_ratio * base_trpbo))

        row: dict[str, float | int] = {
            "patient_id": i,
            "Bpbo_perml": base_bpbo,
            "Bpbref_perml": base_bpbref,
            "Trpbo_perml": base_trpbo,
            "Trpbref_perml": base_trpbref,
            "kBtumorprolif": k_btumor,
            "KBptumor": kbptumor,
            "KTrptumor": ktrptumor,
            "BT_ratio_tumor_init": bt_ratio,
            "Btumor_perml_init": btumor_perml,
            "tumor_burden_factor": burden_factor,
        }
        for pname, (lb, ub) in extra_priors.items():
            row[pname] = sample_from_bounds(rng, lb, ub)
        rows.append(row)

    return pd.DataFrame(rows)


def run_julia_sweep(
    repo: Path,
    design_dir: Path,
    sim_out_dir: Path,
    plan: dict[str, Any],
    julia_project: str | None,
    julia_bin: str | None,
) -> Path:
    sim_cfg = plan["simulation_defaults"]
    env = os.environ.copy()
    env["PHASE1_DESIGN_DIR"] = str(design_dir)
    env["PHASE1_JULIA_OUT_DIR"] = str(sim_out_dir)
    env["TCE_ENGINE"] = str(sim_cfg.get("engine", "canonical"))
    env["TCE_SOLVER"] = str(sim_cfg.get("solver", "cvode_bdf"))
    env["TCE_ABSTOL"] = str(sim_cfg.get("abstol", 1e-8))
    env["TCE_RELTOL"] = str(sim_cfg.get("reltol", 1e-5))
    env["PHASE1_VARIANTS_MODE"] = str(sim_cfg.get("variants_mode", "matlab_empty"))

    if julia_bin:
        julia_exec = julia_bin
    else:
        julia_exec = shutil.which("julia")
        if not julia_exec:
            fallback = Path("/Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia")
            julia_exec = str(fallback) if fallback.is_file() else None
    if not julia_exec:
        raise FileNotFoundError("Could not find Julia executable. Use --julia-bin /path/to/julia.")

    cmd = [julia_exec]
    if julia_project:
        cmd.append(f"--project={julia_project}")
    cmd.append(str(repo / "julia" / "run_phase1_callback_sweep.jl"))
    subprocess.run(cmd, check=True, cwd=repo, env=env)

    out_path = sim_out_dir / "phase1_metrics_julia.csv"
    if not out_path.is_file():
        raise FileNotFoundError(f"Expected Julia output not found: {out_path}")
    return out_path


def build_feature_table(metrics: pd.DataFrame, regimen_names: list[str]) -> pd.DataFrame:
    ok = metrics.copy()
    if "status" in ok.columns:
        ok = ok[ok["status"] == "ok"].copy()
    ok["regimen_norm"] = ok["regimen"].map(normalize_regimen_name)

    metric_cols = [c for c in ok.columns if c not in {"regimen", "regimen_type", "patient_id", "status", "regimen_norm"}]
    by_reg: dict[str, pd.DataFrame] = {}
    for reg in regimen_names:
        sub = ok[ok["regimen_norm"] == reg].copy()
        if sub.empty:
            raise ValueError(f"No successful simulation rows found for regimen {reg}")
        sub = sub.drop_duplicates(subset=["patient_id"], keep="first")
        if "tumor_resid_day42" in sub.columns:
            sub["responder_day42_proxy"] = (sub["tumor_resid_day42"] <= 0.5).astype(float)
        else:
            sub["responder_day42_proxy"] = np.nan
        by_reg[reg] = sub.set_index("patient_id")

    common_ids = None
    for reg in regimen_names:
        ids = set(by_reg[reg].index.tolist())
        common_ids = ids if common_ids is None else (common_ids & ids)
    if not common_ids:
        raise ValueError("No patient IDs have successful runs across all required regimens.")

    ids_sorted = sorted(common_ids)
    out_cols: dict[str, np.ndarray] = {"patient_id": np.asarray(ids_sorted, dtype=int)}
    include_features = ["responder_day42_proxy"] + metric_cols
    include_features = list(dict.fromkeys(include_features))
    for reg in regimen_names:
        sub = by_reg[reg].loc[ids_sorted]
        for feat in include_features:
            if feat in sub.columns:
                out_cols[f"{feat}__{reg}"] = sub[feat].to_numpy(dtype=float)
    return pd.DataFrame(out_cols)


def parse_targets(plan: dict[str, Any]) -> list[CalibrationTarget]:
    out: list[CalibrationTarget] = []
    for t in plan.get("target_bindings", []):
        mode = str(t.get("evaluation_mode", "population_fraction"))
        if mode == "population_fraction":
            out.append(
                CalibrationTarget(
                    target_id=str(t["target_id"]),
                    evaluation_mode=mode,
                    regimen=normalize_regimen_name(str(t["regimen"])),
                    feature=str(t["feature"]),
                    target_value=float(t["target_fraction"]),
                    tolerance=float(t.get("tolerance", 0.05)),
                    weight=float(t.get("weight", 1.0)),
                )
            )
        elif mode == "population_mean":
            out.append(
                CalibrationTarget(
                    target_id=str(t["target_id"]),
                    evaluation_mode=mode,
                    regimen=normalize_regimen_name(str(t["regimen"])),
                    feature=str(t["feature"]),
                    target_value=float(t["target_mean"]),
                    tolerance=float(t.get("tolerance", 0.05)),
                    weight=float(t.get("weight", 1.0)),
                )
            )
        elif mode == "pairwise_difference":
            out.append(
                CalibrationTarget(
                    target_id=str(t["target_id"]),
                    evaluation_mode=mode,
                    feature=str(t["feature"]),
                    tolerance=float(t.get("tolerance", 0.05)),
                    weight=float(t.get("weight", 1.0)),
                    regimen_low=normalize_regimen_name(str(t["regimen_low"])),
                    regimen_high=normalize_regimen_name(str(t["regimen_high"])),
                    min_difference=float(t["min_difference"]) if "min_difference" in t else None,
                    max_difference=float(t["max_difference"]) if "max_difference" in t else None,
                )
            )
    return out


def crs_feature_to_il6_metric(feature: str) -> str:
    mapping = {
        "crs_proxy_day0_2": "il6_peak_0_2",
        "crs_proxy_day0_21": "il6_peak_0_21",
    }
    return mapping.get(feature, "il6_peak_0_2")


def choose_crs_thresholds(feature_table: pd.DataFrame, targets: list[CalibrationTarget]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    crs_features = sorted({t.feature for t in targets if t.feature.startswith("crs_proxy_")})
    for feat in crs_features:
        feat_targets = [t for t in targets if t.feature == feat and t.target_value is not None and t.regimen]
        if not feat_targets:
            continue
        ref = max(feat_targets, key=lambda t: float(t.target_value))
        il6_metric = crs_feature_to_il6_metric(feat)
        col = f"{il6_metric}__{ref.regimen}"
        if col not in feature_table.columns:
            raise ValueError(f"Missing required feature column for {feat}: {col}")
        q = max(0.0, min(1.0, 1.0 - float(ref.target_value)))
        threshold = float(np.quantile(feature_table[col].to_numpy(dtype=float), q))
        out[feat] = {
            "threshold": threshold,
            "reference_regimen": ref.regimen,
            "reference_target": float(ref.target_value),
            "source_column": col,
        }
    return out


def build_regimen_metadata(reg_events: pd.DataFrame) -> pd.DataFrame:
    g = reg_events.groupby("regimen", as_index=False).agg(
        regimen_type=("regimen_type", "first"),
        total_dose_mg=("dose_mg", "sum"),
        max_dose_mg=("dose_mg", "max"),
    )
    d14 = (
        reg_events[np.isclose(reg_events["time_day"].to_numpy(dtype=float), 14.0)]
        .groupby("regimen", as_index=False)["dose_mg"]
        .sum()
        .rename(columns={"dose_mg": "dose_day14_mg"})
    )
    g = g.merge(d14, on="regimen", how="left")
    g["dose_day14_mg"] = g["dose_day14_mg"].fillna(0.0)
    return g


def build_shape_config(
    feature_table: pd.DataFrame,
    regimen_meta: pd.DataFrame,
    objective_defaults: dict[str, Any],
) -> dict[str, Any]:
    enabled = bool(int(objective_defaults.get("enable_multiregimen_shape", 1)))
    if not enabled:
        return {"enabled": False}

    within_type = bool(int(objective_defaults.get("shape_within_regimen_type", 1)))
    resp_feature = str(objective_defaults.get("shape_response_feature", "responder_day42_proxy"))
    tox_feature = str(objective_defaults.get("shape_toxicity_feature", "il6_peak_0_21"))
    if not any(c.startswith(f"{tox_feature}__") for c in feature_table.columns):
        tox_feature = "il6_peak_0_2"
    resp_weight = float(objective_defaults.get("shape_response_weight", 0.25))
    tox_weight = float(objective_defaults.get("shape_toxicity_weight", 0.25))
    resp_tol = float(objective_defaults.get("shape_response_tolerance", 0.04))
    tox_tol = float(objective_defaults.get("shape_toxicity_tolerance", 0.0))

    pairs: list[tuple[str, str, str]] = []
    groups = regimen_meta.groupby("regimen_type", sort=False) if within_type else [("all", regimen_meta)]
    for gname, sub in groups:
        ss = sub.sort_values(["total_dose_mg", "dose_day14_mg", "max_dose_mg", "regimen"]).copy()
        regs = ss["regimen"].tolist()
        for i in range(len(regs) - 1):
            pairs.append((regs[i], regs[i + 1], str(gname)))

    if tox_tol <= 0.0:
        means = []
        for reg in regimen_meta["regimen"].tolist():
            col = f"{tox_feature}__{reg}"
            if col in feature_table.columns:
                means.append(float(np.nanmean(feature_table[col].to_numpy(dtype=float))))
        if len(means) >= 2:
            tox_tol = max(1e-12, 0.1 * (max(means) - min(means)))
        else:
            tox_tol = 50.0

    return {
        "enabled": True,
        "within_type": within_type,
        "response_feature": resp_feature,
        "toxicity_feature": tox_feature,
        "response_weight": resp_weight,
        "toxicity_weight": tox_weight,
        "response_tolerance": max(resp_tol, 1e-12),
        "toxicity_tolerance": max(tox_tol, 1e-12),
        "pairs": pairs,
    }


def feature_values(
    rows: pd.DataFrame,
    feature: str,
    regimen: str,
    crs_thresholds: dict[str, dict[str, Any]],
) -> np.ndarray:
    if feature.startswith("crs_proxy_"):
        if feature not in crs_thresholds:
            raise ValueError(f"Missing threshold for CRS proxy feature: {feature}")
        il6_metric = crs_feature_to_il6_metric(feature)
        col = f"{il6_metric}__{regimen}"
        if col not in rows.columns:
            raise ValueError(f"Missing feature column: {col}")
        thr = float(crs_thresholds[feature]["threshold"])
        x = rows[col].to_numpy(dtype=float)
        return (x >= thr).astype(float)
    col = f"{feature}__{regimen}"
    if col not in rows.columns:
        raise ValueError(f"Missing feature column: {col}")
    return rows[col].to_numpy(dtype=float)


def safe_mean(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan")
    return float(np.mean(x))


def score_subset(
    feature_table: pd.DataFrame,
    subset_idx: np.ndarray,
    targets: list[CalibrationTarget],
    crs_thresholds: dict[str, dict[str, Any]],
    shape_cfg: dict[str, Any],
) -> tuple[float, list[dict[str, Any]]]:
    rows = feature_table.iloc[subset_idx]
    total = 0.0
    details: list[dict[str, Any]] = []

    for t in targets:
        if t.evaluation_mode in {"population_fraction", "population_mean"}:
            if not t.regimen:
                continue
            achieved = safe_mean(feature_values(rows, t.feature, t.regimen, crs_thresholds))
            target = float(t.target_value) if t.target_value is not None else np.nan
            if not np.isfinite(achieved) or not np.isfinite(target):
                term = 1e6
            else:
                err = (achieved - target) / max(t.tolerance, 1e-12)
                term = t.weight * (err ** 2)
            total += term
            details.append(
                {
                    "target_id": t.target_id,
                    "evaluation_mode": t.evaluation_mode,
                    "regimen": t.regimen,
                    "feature": t.feature,
                    "achieved_value": achieved,
                    "target_value": target,
                    "tolerance": t.tolerance,
                    "weighted_term": term,
                }
            )
        elif t.evaluation_mode == "pairwise_difference":
            if not t.regimen_low or not t.regimen_high:
                continue
            lo = safe_mean(feature_values(rows, t.feature, t.regimen_low, crs_thresholds))
            hi = safe_mean(feature_values(rows, t.feature, t.regimen_high, crs_thresholds))
            diff = hi - lo
            violation = 0.0
            if t.min_difference is not None:
                violation = max(violation, float(t.min_difference) - diff)
            if t.max_difference is not None:
                violation = max(violation, diff - float(t.max_difference))
            if not np.isfinite(diff):
                term = 1e6
            else:
                term = t.weight * (violation / max(t.tolerance, 1e-12)) ** 2
            total += term
            details.append(
                {
                    "target_id": t.target_id,
                    "evaluation_mode": t.evaluation_mode,
                    "feature": t.feature,
                    "regimen_low": t.regimen_low,
                    "regimen_high": t.regimen_high,
                    "achieved_difference": diff,
                    "min_difference": t.min_difference,
                    "max_difference": t.max_difference,
                    "tolerance": t.tolerance,
                    "weighted_term": term,
                }
            )

    if shape_cfg.get("enabled", False):
        resp_feat = str(shape_cfg["response_feature"])
        tox_feat = str(shape_cfg["toxicity_feature"])
        resp_w = float(shape_cfg["response_weight"])
        tox_w = float(shape_cfg["toxicity_weight"])
        resp_tol = float(shape_cfg["response_tolerance"])
        tox_tol = float(shape_cfg["toxicity_tolerance"])
        for reg_lo, reg_hi, group_name in shape_cfg.get("pairs", []):
            resp_lo = safe_mean(feature_values(rows, resp_feat, reg_lo, crs_thresholds))
            resp_hi = safe_mean(feature_values(rows, resp_feat, reg_hi, crs_thresholds))
            tox_lo = safe_mean(feature_values(rows, tox_feat, reg_lo, crs_thresholds))
            tox_hi = safe_mean(feature_values(rows, tox_feat, reg_hi, crs_thresholds))

            resp_violation = max(0.0, resp_lo - resp_hi)
            tox_violation = max(0.0, tox_lo - tox_hi)
            term_resp = resp_w * (resp_violation / max(resp_tol, 1e-12)) ** 2
            term_tox = tox_w * (tox_violation / max(tox_tol, 1e-12)) ** 2
            total += term_resp + term_tox

            details.append(
                {
                    "target_id": f"shape_response_{reg_lo}_to_{reg_hi}",
                    "evaluation_mode": "shape_order",
                    "group": group_name,
                    "feature": resp_feat,
                    "regimen_low": reg_lo,
                    "regimen_high": reg_hi,
                    "achieved_difference": resp_hi - resp_lo,
                    "tolerance": resp_tol,
                    "weighted_term": term_resp,
                }
            )
            details.append(
                {
                    "target_id": f"shape_toxicity_{reg_lo}_to_{reg_hi}",
                    "evaluation_mode": "shape_order",
                    "group": group_name,
                    "feature": tox_feat,
                    "regimen_low": reg_lo,
                    "regimen_high": reg_hi,
                    "achieved_difference": tox_hi - tox_lo,
                    "tolerance": tox_tol,
                    "weighted_term": term_tox,
                }
            )

    return total, details


def random_subset_search(
    feature_table: pd.DataFrame,
    n_select: int,
    n_draws: int,
    seed: int,
    targets: list[CalibrationTarget],
    crs_thresholds: dict[str, dict[str, Any]],
    shape_cfg: dict[str, Any],
) -> tuple[np.ndarray, float, list[dict[str, Any]], pd.DataFrame]:
    rng = np.random.default_rng(seed)
    n = len(feature_table)
    if n_select > n:
        raise ValueError(f"n_select={n_select} cannot exceed available patients with complete simulations ({n}).")

    best_idx = None
    best_score = np.inf
    best_details: list[dict[str, Any]] = []
    trace_rows: list[dict[str, Any]] = []

    all_idx = np.arange(n)
    for i in range(n_draws):
        idx = np.sort(rng.choice(all_idx, size=n_select, replace=False))
        score, details = score_subset(feature_table, idx, targets, crs_thresholds, shape_cfg)
        trace_rows.append({"draw": i + 1, "objective": score})
        if score < best_score:
            best_score = score
            best_idx = idx
            best_details = details

    assert best_idx is not None
    return best_idx, float(best_score), best_details, pd.DataFrame(trace_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Sample candidate VPs, simulate in Julia, and prune against calibration targets.")
    parser.add_argument("--targets-json", type=Path, default=Path("generated/vpop_targets/vpop_calibration_targets.json"))
    parser.add_argument("--eval-plan-json", type=Path, default=Path("generated/vpop_targets/vpop_target_evaluation_plan.json"))
    parser.add_argument("--params-xlsx", type=Path, default=Path("assets/params_41540_2020_145_MOESM2_ESM.xlsx"))
    parser.add_argument("--extra-priors-csv", type=Path, default=Path("generated/vpop_targets/parameter_prior_bounds_from_params_sheet.csv"))
    parser.add_argument("--sample-extra-params", type=str, default="")
    parser.add_argument("--include-standard-regimens", type=int, default=-1, help="-1 use plan default, 0 disable, 1 enable")
    parser.add_argument("--run-tag", type=str, default="")
    parser.add_argument("--n-candidates", type=int, default=0)
    parser.add_argument("--n-select", type=int, default=0)
    parser.add_argument("--n-random-subsets", type=int, default=0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--subset-seed", type=int, default=0, help="RNG seed for random subset search (default: seed+17)")
    parser.add_argument("--spread-mode", choices=["primary", "sensitivity"], default="primary")
    parser.add_argument("--julia-project", type=str, default="")
    parser.add_argument("--julia-bin", type=str, default="")
    parser.add_argument("--skip-sim", action="store_true")
    parser.add_argument("--metrics-csv", type=str, default="")
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[1]
    targets_json = load_json(repo / args.targets_json)
    plan = load_json(repo / args.eval_plan_json)

    sampling_defaults = plan.get("sampling_defaults", {})
    objective_defaults = plan.get("objective_defaults", {})
    n_candidates = args.n_candidates or int(sampling_defaults.get("n_candidates", 400))
    n_select = args.n_select or int(sampling_defaults.get("n_select", 140))
    seed = args.seed or int(sampling_defaults.get("seed", 20260228))
    n_draws = args.n_random_subsets or int(objective_defaults.get("n_random_subsets", 3000))

    include_default = int(objective_defaults.get("include_standard_phase1_regimens", 1))
    include_standard = bool(include_default if args.include_standard_regimens < 0 else args.include_standard_regimens)
    standard_regimens_csv = repo / "generated" / "phase1_design_standard" / "regimen_events.csv"
    plan = merge_plan_with_standard_regimens(plan, standard_regimens_csv, include_standard=include_standard)

    if not args.run_tag:
        tag = datetime.now(timezone.utc).strftime("run_%Y%m%dT%H%M%SZ")
    else:
        tag = args.run_tag

    run_dir = repo / "generated" / "vpop_pruning" / tag
    design_dir = run_dir / "design"
    sim_dir = run_dir / "sim"
    run_dir.mkdir(parents=True, exist_ok=True)
    design_dir.mkdir(parents=True, exist_ok=True)
    sim_dir.mkdir(parents=True, exist_ok=True)

    dlbcl_overrides, _dmap, base_bpbo, base_bpbref, base_trpbo, base_trpbref, base_btumor = load_param_defaults(repo / args.params_xlsx)
    k_bounds, bt_bounds, spread_bounds = get_hard_prior_bounds(targets_json, args.spread_mode)

    default_extra_names = [
        "kBapop",
        "kBkill",
        "fBkill",
        "kBprolif",
        "kTaexit",
        "kTact",
        "fTadeact",
        "fTap",
        "kTaapop",
        "fTaprolif",
        "fTrapop",
    ]
    extra_names = parse_extra_param_list(args.sample_extra_params, default_extra_names)
    extra_priors = load_extra_param_priors(repo / args.extra_priors_csv, extra_names)

    patients = sample_candidates(
        n_candidates=n_candidates,
        seed=seed,
        base_bpbo=base_bpbo,
        base_bpbref=base_bpbref,
        base_trpbo=base_trpbo,
        base_trpbref=base_trpbref,
        base_btumor_perml=base_btumor,
        k_bounds=k_bounds,
        bt_bounds=bt_bounds,
        spread_bounds=spread_bounds,
        extra_priors=extra_priors,
    )
    reg_events = regimen_events_from_plan(plan)

    patients_path = design_dir / "patients.csv"
    reg_path = design_dir / "regimen_events.csv"
    dlbcl_path = design_dir / "dlbcl_param_overrides.csv"
    patients.to_csv(patients_path, index=False)
    reg_events.to_csv(reg_path, index=False)
    dlbcl_overrides.to_csv(dlbcl_path, index=False)

    if args.skip_sim:
        if not args.metrics_csv.strip():
            raise ValueError("--skip-sim requires --metrics-csv")
        metrics_path_raw = Path(args.metrics_csv)
        metrics_path = (repo / metrics_path_raw) if not metrics_path_raw.is_absolute() else metrics_path_raw
    else:
        julia_project = args.julia_project.strip() or None
        julia_bin = args.julia_bin.strip() or None
        metrics_path = run_julia_sweep(repo, design_dir, sim_dir, plan, julia_project, julia_bin)

    metrics = pd.read_csv(metrics_path)
    regimen_names = sorted(reg_events["regimen"].astype(str).unique().tolist())
    feature_table = build_feature_table(metrics, regimen_names)

    targets = parse_targets(plan)
    crs_thresholds = choose_crs_thresholds(feature_table, targets)
    regimen_meta = build_regimen_metadata(reg_events)
    shape_cfg = build_shape_config(feature_table, regimen_meta, objective_defaults)

    subset_seed = args.subset_seed or (seed + 17)

    best_idx, best_score, best_details, trace = random_subset_search(
        feature_table=feature_table,
        n_select=n_select,
        n_draws=n_draws,
        seed=subset_seed,
        targets=targets,
        crs_thresholds=crs_thresholds,
        shape_cfg=shape_cfg,
    )

    selected_ids = feature_table.iloc[best_idx]["patient_id"].tolist()
    selected_patients = patients[patients["patient_id"].isin(selected_ids)].copy().sort_values("patient_id")

    feature_path = run_dir / "candidate_features.csv"
    trace_path = run_dir / "objective_trace.csv"
    selected_path = run_dir / "selected_patients.csv"
    summary_path = run_dir / "pruning_summary.json"

    feature_table.to_csv(feature_path, index=False)
    trace.to_csv(trace_path, index=False)
    selected_patients.to_csv(selected_path, index=False)

    summary = {
        "schema_version": "2.0.0",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "run_tag": tag,
        "inputs": {
            "targets_json": str(args.targets_json),
            "eval_plan_json": str(args.eval_plan_json),
            "params_xlsx": str(args.params_xlsx),
            "extra_priors_csv": str(args.extra_priors_csv),
            "sampled_extra_params": sorted(extra_priors.keys()),
            "include_standard_regimens": include_standard,
            "n_candidates": n_candidates,
            "n_select": n_select,
            "seed": seed,
            "subset_seed": subset_seed,
            "n_random_subsets": n_draws,
            "spread_mode": args.spread_mode,
            "julia_project": args.julia_project,
            "julia_bin": args.julia_bin,
            "metrics_csv": str(metrics_path),
        },
        "simulation": {
            "n_regimens": int(reg_events["regimen"].nunique()),
            "regimens": regimen_names,
            "n_metrics_rows": int(len(metrics)),
            "n_eligible_patients_complete": int(len(feature_table)),
            "metric_columns_available": sorted([c for c in metrics.columns if c not in {"status"}]),
        },
        "toxicity_proxy_thresholds": crs_thresholds,
        "objective": {
            "best_score": best_score,
            "best_details": best_details,
            "shape_config": shape_cfg,
        },
        "not_directly_evaluable_with_current_model_outputs": plan.get(
            "not_directly_evaluable_with_current_model_outputs", []
        ),
        "outputs": {
            "patients_candidates_csv": str(patients_path),
            "regimen_events_csv": str(reg_path),
            "selected_patients_csv": str(selected_path),
            "candidate_features_csv": str(feature_path),
            "objective_trace_csv": str(trace_path),
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2))

    print(run_dir)
    print(selected_path)
    print(summary_path)
    print(f"best_objective={best_score:.6g}")
    n_shape = 0
    for d in best_details:
        if d.get("evaluation_mode") == "shape_order":
            n_shape += 1
            continue
        tid = d.get("target_id", "")
        ach = d.get("achieved_value", d.get("achieved_difference", np.nan))
        tgt = d.get("target_value", d.get("min_difference", np.nan))
        tol = d.get("tolerance", np.nan)
        print(f"{tid}: achieved={ach:.4f} target={tgt:.4f} tol={tol:.4f}")
    if n_shape:
        print(f"shape_terms_suppressed={n_shape}")


if __name__ == "__main__":
    main()
