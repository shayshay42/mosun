#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import pandas as pd


def build_parameter_bounds(params_xlsx: Path, out_csv: Path) -> int:
    df = pd.read_excel(params_xlsx, sheet_name="Sheet1")
    keep = df[df["Cyno_LB"].notna() | df["Cyno_UB"].notna()].copy()
    keep["group"] = keep["GROUP"].fillna("").astype(str)
    keep["name"] = keep["NAME"].astype(str)
    keep["comment"] = keep["Comments"].fillna("").astype(str)
    keep["reference"] = keep["Reference"].fillna("").astype(str)
    keep["prior_type"] = keep["comment"].str.contains("Estimated", case=False, na=False).map(
        {True: "estimated", False: "other"}
    )
    cols = [
        "group",
        "name",
        "prior_type",
        "Cyno_LB",
        "Cyno_UB",
        "Cyno",
        "Human",
        "ALL",
        "DLBCL",
        "comment",
        "reference",
    ]
    out = keep[cols].rename(
        columns={
            "Cyno_LB": "cyno_lb",
            "Cyno_UB": "cyno_ub",
            "Cyno": "cyno_value",
            "Human": "human_value",
            "ALL": "all_value",
            "DLBCL": "dlbcl_value",
        }
    )
    out.to_csv(out_csv, index=False)
    return len(out)


def build_flat_targets(target_json: Path, out_csv: Path) -> int:
    obj = json.loads(target_json.read_text())
    rows: list[dict[str, object]] = []

    def push(
        id_: str,
        category: str,
        population: str,
        regimen: str,
        metric: str,
        value: object,
        units: str,
        priority: str,
        source: str,
    ) -> None:
        rows.append(
            {
                "id": id_,
                "category": category,
                "population": population,
                "regimen": regimen,
                "metric": metric,
                "target_value": value,
                "units": units,
                "priority": priority,
                "source": source,
            }
        )

    for target in obj.get("hard_priors", []):
        sid = target["id"]
        priority = target.get("priority", "")
        source = target.get("source", {}).get("location") or target.get("source", {}).get("url", "")
        if "bounds" in target:
            bounds = target["bounds"]
            push(sid, "hard_prior", "all", "na", "lower_bound", bounds.get("lower"), target.get("units", ""), priority, source)
            push(sid, "hard_prior", "all", "na", "upper_bound", bounds.get("upper"), target.get("units", ""), priority, source)
        if "bounds_relative_factor" in target:
            bounds = target["bounds_relative_factor"]
            push(sid, "hard_prior", "all", "na", "relative_lower_factor", bounds.get("lower"), "fold", priority, source)
            push(sid, "hard_prior", "all", "na", "relative_upper_factor", bounds.get("upper"), "fold", priority, source)
        if "target_fraction" in target:
            push(sid, "hard_prior", "all", "na", "target_fraction", target["target_fraction"], "fraction", priority, source)

    for target in obj.get("clinical_endpoint_targets", []):
        sid = target["id"]
        priority = target.get("priority", "")
        population = target.get("population", "all")
        regimen = target.get("regimen", "na")
        source = target.get("source", {}).get("location") or target.get("source", {}).get("url", "")
        if "metrics" in target:
            for metric, value in target["metrics"].items():
                push(sid, "clinical_endpoint", population, regimen, metric, value, "percent_or_months", priority, source)

    for target in obj.get("pk_and_baseline_distribution_targets", []):
        sid = target["id"]
        priority = target.get("priority", "")
        population = target.get("population", "all")
        source = target.get("source", {}).get("location") or target.get("source", {}).get("url", "")
        if "target_summary" in target:
            for metric, value in target["target_summary"].items():
                unit = "mm2" if "SPD" in target.get("metric", "") else "na"
                push(sid, "pk_or_baseline", population, "na", metric, value, unit, priority, source)

    if not rows:
        out_csv.write_text("id,category,population,regimen,metric,target_value,units,priority,source\n")
        return 0

    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build VPop target bundle files.")
    parser.add_argument(
        "--params-xlsx",
        type=Path,
        default=Path("assets/params_41540_2020_145_MOESM2_ESM.xlsx"),
    )
    parser.add_argument(
        "--target-json",
        type=Path,
        default=Path("generated/vpop_targets/vpop_calibration_targets.json"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("generated/vpop_targets"),
    )
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    bounds_csv = args.out_dir / "parameter_prior_bounds_from_params_sheet.csv"
    flat_csv = args.out_dir / "vpop_calibration_targets_flat.csv"

    n_bounds = build_parameter_bounds(args.params_xlsx, bounds_csv)
    n_flat = build_flat_targets(args.target_json, flat_csv)

    print(bounds_csv)
    print(flat_csv)
    print(f"n_parameter_bounds={n_bounds}")
    print(f"n_flat_targets={n_flat}")


if __name__ == "__main__":
    main()
