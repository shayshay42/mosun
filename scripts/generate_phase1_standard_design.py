#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

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

    rows: list[dict[str, object]] = []
    for d in fixed:
        reg = f"fixed_{d:g}mg"
        for i, t in enumerate([0.0, 21.0, 42.0, 63.0], start=1):
            rows.append(
                {
                    "regimen": reg,
                    "regimen_type": "fixed_q3w",
                    "event_idx": i,
                    "time_day": t,
                    "dose_mg": d,
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


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    out_dir = repo / "generated" / "phase1_design_standard"
    out_dir.mkdir(parents=True, exist_ok=True)

    ptab = pd.read_excel(repo / "assets" / "params_41540_2020_145_MOESM2_ESM.xlsx", sheet_name="Sheet1")
    dlbcl = ptab[ptab["DLBCL"].notna()][["NAME", "DLBCL"]].copy()
    dlbcl.columns = ["name", "value"]
    dlbcl.to_csv(out_dir / "dlbcl_param_overrides.csv", index=False)

    dmap = {str(r["name"]): float(r["value"]) for _, r in dlbcl.iterrows()}

    bpbo = 250_000.0
    trpbo = 500_000.0
    kbptumor = float(dmap["KBptumor"])
    ktrptumor = float(dmap["KTrptumor"])
    kbtumorprolif = float(dmap["kBtumorprolif"])
    btumor = float(dmap["Btumor_perml"])
    bt_ratio = (kbptumor * bpbo) / max(ktrptumor * trpbo, 1e-12)

    patient = pd.DataFrame(
        [
            {
                "patient_id": 1,
                "Bpbo_perml": bpbo,
                "Bpbref_perml": bpbo,
                "Trpbo_perml": trpbo,
                "Trpbref_perml": trpbo,
                "kBtumorprolif": kbtumorprolif,
                "KBptumor": kbptumor,
                "KTrptumor": ktrptumor,
                "BT_ratio_tumor_init": bt_ratio,
                "Btumor_perml_init": btumor,
            }
        ]
    )

    reg = regimen_events()
    patient.to_csv(out_dir / "patients.csv", index=False)
    reg.to_csv(out_dir / "regimen_events.csv", index=False)

    print(out_dir / "patients.csv")
    print(out_dir / "regimen_events.csv")
    print(out_dir / "dlbcl_param_overrides.csv")
    print(f"n_patients={len(patient)} n_regimens={reg['regimen'].nunique()} n_runs={len(patient) * reg['regimen'].nunique()}")


if __name__ == "__main__":
    main()

