#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def regimen_events() -> pd.DataFrame:
    fixed = [0.05, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 2.8]
    step = [
        (0.4, 1.0, 2.8),
        (0.8, 2.0, 4.2),
        (1.0, 1.0, 3.0),
        (1.0, 2.0, 6.0),
        (0.8, 2.0, 6.0),
        (1.0, 2.0, 9.0),
        (1.0, 2.0, 13.5),
        (1.0, 2.0, 20.0),
        (1.0, 2.0, 27.0),
    ]

    rows: list[dict] = []

    for d in fixed:
        reg = f"fixed_{d:g}mg"
        times = [0.0, 21.0, 42.0, 63.0]
        amounts = [d, d, d, d]
        for i, (t, a) in enumerate(zip(times, amounts), start=1):
            rows.append(
                {
                    "regimen": reg,
                    "regimen_type": "fixed_q3w",
                    "event_idx": i,
                    "time_day": t,
                    "dose_mg": a,
                }
            )

    for d1, d2, dt in step:
        reg = f"step_{d1:g}_{d2:g}_{dt:g}mg"
        times = [0.0, 7.0, 14.0, 21.0, 42.0, 63.0]
        amounts = [d1, d2, dt, dt, dt, dt]
        for i, (t, a) in enumerate(zip(times, amounts), start=1):
            rows.append(
                {
                    "regimen": reg,
                    "regimen_type": "stepup_q3w",
                    "event_idx": i,
                    "time_day": t,
                    "dose_mg": a,
                }
            )

    return pd.DataFrame(rows)


def sample_patients(
    n_patients: int,
    seed: int,
    tumor_burden_spread: float,
    base_bpbo: float,
    base_bpbref: float,
    base_trpbo: float,
    base_trpbref: float,
    base_btumor_perml: float,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows: list[dict] = []

    for i in range(1, n_patients + 1):
        # Follow-up paper bounds:
        # - tumor B:T in [4, 80]
        # - tumor proliferation in [0, 0.15] day^-1
        bt_ratio = float(np.exp(rng.uniform(np.log(4.0), np.log(80.0))))
        k_btumor = float(rng.uniform(0.0, 0.15))
        bfactor = float(rng.uniform(1.0 - tumor_burden_spread, 1.0 + tumor_burden_spread))
        btumor_perml = float(base_btumor_perml * bfactor)

        kbptumor = float(btumor_perml / base_bpbo)
        ktrptumor = float(btumor_perml / (bt_ratio * base_trpbo))

        rows.append(
            {
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
            }
        )

    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate sampled phase-1 design for MATLAB/Julia sweeps.")
    parser.add_argument("--n-patients", type=int, default=24)
    parser.add_argument("--seed", type=int, default=20260226)
    parser.add_argument(
        "--tumor-burden-spread",
        type=float,
        default=0.10,
        help="Fractional spread around reference tumor B-cell burden (0.10 = +/-10%%).",
    )
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[1]
    out_dir = repo / "generated" / "phase1_design"
    out_dir.mkdir(parents=True, exist_ok=True)

    params_path = repo / "assets" / "params_41540_2020_145_MOESM2_ESM.xlsx"
    ptab = pd.read_excel(params_path, sheet_name="Sheet1")
    dlbcl = {str(r["NAME"]): r["DLBCL"] for _, r in ptab.iterrows() if pd.notna(r.get("DLBCL"))}

    # Follow-up table may provide *_ref values only; use them as bo fallbacks when needed.
    base_bpbo = float(dlbcl.get("Bpbo_perml", dlbcl.get("Bpbref_perml", 250_000.0)))
    base_bpbref = float(dlbcl.get("Bpbref_perml", base_bpbo))
    base_trpbo = float(dlbcl.get("Trpbo_perml", dlbcl.get("Trpbref_perml", 500_000.0)))
    base_trpbref = float(dlbcl.get("Trpbref_perml", base_trpbo))
    base_btumor_perml = float(dlbcl.get("Btumor_perml", 3.25e9))

    reg_df = regimen_events()
    pat_df = sample_patients(
        n_patients=args.n_patients,
        seed=args.seed,
        tumor_burden_spread=args.tumor_burden_spread,
        base_bpbo=base_bpbo,
        base_bpbref=base_bpbref,
        base_trpbo=base_trpbo,
        base_trpbref=base_trpbref,
        base_btumor_perml=base_btumor_perml,
    )

    reg_path = out_dir / "regimen_events.csv"
    pat_path = out_dir / "patients.csv"
    reg_df.to_csv(reg_path, index=False)
    pat_df.to_csv(pat_path, index=False)

    print(reg_path)
    print(pat_path)
    print(f"n_regimens={reg_df['regimen'].nunique()} n_patients={len(pat_df)} n_runs={reg_df['regimen'].nunique() * len(pat_df)}")


if __name__ == "__main__":
    main()
