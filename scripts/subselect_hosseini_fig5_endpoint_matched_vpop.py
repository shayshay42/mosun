#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE = REPO_ROOT / "generated/figures/vpop_pruning/susilo_efast_hosseini2020_fig5_vpop1000_20260516/vpop250_full1000_hybrid_sinkhorn1000_hosseini_fig5_endpoint_violin_source.csv"
REFERENCE = REPO_ROOT / "generated/figures/vpop_pruning/susilo_efast_hosseini2020_fig5_vpop1000_20260516/hosseini2020_fig5_digitized_endpoint_reference_for_violins.csv"
DIGITIZED_WATERFALL = REPO_ROOT / "assets/digitization_hosseini_2020_vpop/hosseini2020_fig5_digitized_tumor_waterfall_shapes.csv"
OUT_DIR = REPO_ROOT / "generated/figures/vpop_pruning/endpoint_matched_subvpop_from_all_fig5_patients_20260517"

REGIMEN_ORDER = [
    "c1d1_20_only_then_20_q3w",
    "c1d1_d8_d15_6p7_then_20_q3w",
    "c1d1_d8_d15_1p6_10_10_then_20_q3w",
    "c1d1_d8_1p6_20_then_20_q3w",
]
REGIMEN_LABELS = {
    "c1d1_20_only_then_20_q3w": "20 mg C1D1",
    "c1d1_d8_d15_6p7_then_20_q3w": "6.7/6.7/6.7 mg",
    "c1d1_d8_d15_1p6_10_10_then_20_q3w": "1.6/10/10 mg",
    "c1d1_d8_1p6_20_then_20_q3w": "1.6/20 mg",
}
DIGITIZED_REGIMEN_TO_CODE = {
    "20 mg C1D1; 20 mg Day 1 subsequent cycles": "c1d1_20_only_then_20_q3w",
    "6.7 mg C1D1/C1D8/C1D15; 20 mg Day 1 subsequent cycles": "c1d1_d8_d15_6p7_then_20_q3w",
    "1.6/10/10 mg C1D1/C1D8/C1D15; 20 mg Day 1 subsequent cycles": "c1d1_d8_d15_1p6_10_10_then_20_q3w",
    "1.6/20 mg C1D1/C1D8; 20 mg Day 1 subsequent cycles": "c1d1_d8_1p6_20_then_20_q3w",
}
REGIMEN_TO_DIGITIZED = {v: k for k, v in DIGITIZED_REGIMEN_TO_CODE.items()}
METRICS = [
    "profile_peak_il6_pg_ml",
    "profile_peak_cd69_cd8_pct",
    "day84_tumor_change_pct",
]
Q = np.array([5, 50, 95], dtype=float)
SOURCE_ORDER = ["VPop250", "Old full1000", "Hybrid", "Sinkhorn VPop1000", "Endpoint-matched VPop1000", "Endpoint-matched VPop250"]
COLORS = {
    "VPop250": "#4C78A8",
    "Old full1000": "#6A3D9A",
    "Hybrid": "#009E73",
    "Sinkhorn VPop1000": "#D55E00",
    "Endpoint-matched VPop1000": "#C61A09",
    "Endpoint-matched VPop250": "#111111",
}


def fs_path(path: Path) -> str:
    resolved = str(path.resolve())
    if not resolved.startswith("\\\\?\\"):
        return "\\\\?\\" + resolved
    return resolved


def canonical_key(vpop: str, candidate_id: str) -> str:
    cid = str(candidate_id)
    if "__" in cid and cid.startswith("2026"):
        return f"efast::{cid}"
    return f"{vpop}::{cid}"


def load_pool() -> tuple[pd.DataFrame, pd.DataFrame]:
    long = pd.read_csv(SOURCE)
    long = long[long["regimen"].isin(REGIMEN_ORDER) & long["metric"].isin(METRICS)].copy()
    long["canonical_id"] = [canonical_key(v, c) for v, c in zip(long["vpop"], long["candidate_id"])]
    source_members = (
        long[["canonical_id", "vpop", "candidate_id"]]
        .drop_duplicates()
        .groupby("canonical_id", as_index=False)
        .agg(source_vpops=("vpop", lambda s: ";".join(sorted(set(map(str, s))))), source_candidate_ids=("candidate_id", lambda s: ";".join(sorted(set(map(str, s))))))
    )
    long = long.sort_values(["canonical_id", "regimen", "metric", "vpop"]).drop_duplicates(["canonical_id", "regimen", "metric"], keep="first")
    wide = long.pivot_table(index="canonical_id", columns=["regimen", "metric"], values="value", aggfunc="first")
    missing = int(wide.isna().sum().sum())
    if missing:
        raise RuntimeError(f"Missing endpoint values after pivot: {missing}")
    return wide, source_members


def load_reference() -> tuple[pd.DataFrame, dict[str, tuple[np.ndarray, np.ndarray]]]:
    ref = pd.read_csv(REFERENCE)
    wf = pd.read_csv(DIGITIZED_WATERFALL)
    waterfall_targets: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for regimen in REGIMEN_ORDER:
        sub = wf[wf["regimen"] == REGIMEN_TO_DIGITIZED[regimen]].sort_values("patient_rank_fraction")
        waterfall_targets[regimen] = (sub["patient_rank_fraction"].to_numpy(float), sub["tumor_change_percent_day84"].to_numpy(float))
    return ref, waterfall_targets


def arrays_from_wide(wide: pd.DataFrame) -> dict[tuple[str, str], np.ndarray]:
    return {(regimen, metric): wide[(regimen, metric)].to_numpy(float) for regimen in REGIMEN_ORDER for metric in METRICS}


def score_subset(arrays: dict[tuple[str, str], np.ndarray], subset_idx: np.ndarray, ref: pd.DataFrame, waterfall_targets: dict[str, tuple[np.ndarray, np.ndarray]]) -> dict[str, float]:
    il6_scores = []
    tcell_scores = []
    waterfall_scores = []
    responder_scores = []
    for regimen in REGIMEN_ORDER:
        for metric, bucket in [("profile_peak_il6_pg_ml", il6_scores), ("profile_peak_cd69_cd8_pct", tcell_scores)]:
            values = arrays[(regimen, metric)][subset_idx]
            target = ref[(ref["regimen"] == regimen) & (ref["metric"] == metric)].iloc[0]
            target_values = np.array([target["lower"], target["median"], target["upper"]], dtype=float)
            if metric == "profile_peak_il6_pg_ml":
                sim = np.log10(np.maximum(np.nanpercentile(values, Q), 1e-12))
                tar = np.log10(np.maximum(target_values, 1e-12))
                bucket.append(float(np.sqrt(np.mean((sim - tar) ** 2))))
            else:
                sim = np.nanpercentile(values, Q)
                bucket.append(float(np.sqrt(np.mean((sim - target_values) ** 2))) / 100.0)

        tumor = arrays[(regimen, "day84_tumor_change_pct")][subset_idx]
        tumor_sorted = np.sort(tumor)[::-1]
        rank = np.linspace(0.0, 1.0, len(tumor_sorted))
        target_rank, target_y = waterfall_targets[regimen]
        interp = np.interp(target_rank, rank, tumor_sorted)
        waterfall_scores.append(float(np.sqrt(np.nanmean((interp - target_y) ** 2))) / 100.0)
        responder_scores.append(abs(float(np.mean(tumor <= -50.0)) - float(np.mean(target_y <= -50.0))))

    il6 = float(np.mean(il6_scores))
    tcell = float(np.mean(tcell_scores))
    waterfall = float(np.mean(waterfall_scores))
    responder = float(np.mean(responder_scores))
    tumor = 0.8 * waterfall + 0.2 * responder
    return {"IL6": il6, "Tcell": tcell, "waterfall": waterfall, "responder": responder, "tumor": tumor, "total": (il6 + tcell + tumor) / 3.0}


def source_score_table(wide: pd.DataFrame, arrays: dict[tuple[str, str], np.ndarray], source_members: pd.DataFrame, ref: pd.DataFrame, waterfall_targets: dict[str, tuple[np.ndarray, np.ndarray]]) -> pd.DataFrame:
    rows = []
    for source in ["VPop250", "Full1000", "Hybrid", "Sinkhorn VPop1000"]:
        label = "Old full1000" if source == "Full1000" else source
        ids = source_members[source_members["source_vpops"].str.contains(source, regex=False)]["canonical_id"]
        idx = wide.index.get_indexer(ids)
        idx = idx[idx >= 0]
        if len(idx) == 0:
            continue
        scores = score_subset(arrays, idx, ref, waterfall_targets)
        rows.append({"vpop": label, "n": int(len(idx)), **scores})
    return pd.DataFrame(rows)


def scalar_candidate_weights(arrays: dict[tuple[str, str], np.ndarray], n_pool: int, ref: pd.DataFrame) -> np.ndarray:
    penalties = np.zeros(n_pool, dtype=float)
    n_terms = 0
    for regimen in REGIMEN_ORDER:
        for metric in ["profile_peak_il6_pg_ml", "profile_peak_cd69_cd8_pct"]:
            target = ref[(ref["regimen"] == regimen) & (ref["metric"] == metric)].iloc[0]
            lo, med, hi = float(target["lower"]), float(target["median"]), float(target["upper"])
            values = arrays[(regimen, metric)]
            if metric == "profile_peak_il6_pg_ml":
                values = np.log10(np.maximum(values, 1e-12))
                lo, med, hi = np.log10(max(lo, 1e-12)), np.log10(max(med, 1e-12)), np.log10(max(hi, 1e-12))
            scale = max(abs(hi - lo), 1e-6)
            anchors = np.array([lo, med, hi], dtype=float)
            penalties += np.min(np.abs(values[:, None] - anchors[None, :]), axis=1) / scale
            n_terms += 1
    penalties = penalties / max(n_terms, 1)
    weights = np.exp(-2.0 * penalties)
    weights = np.maximum(weights, 1e-12)
    return weights / weights.sum()


def make_initial_subset(
    rng: np.random.Generator,
    n_pool: int,
    n_select: int,
    source_members: pd.DataFrame,
    wide: pd.DataFrame,
    weights: np.ndarray,
    mode: int,
) -> np.ndarray:
    if mode == 0:
        return rng.choice(n_pool, size=n_select, replace=False)
    if mode == 1:
        return rng.choice(n_pool, size=n_select, replace=False, p=weights)
    source_name = ["VPop250", "Sinkhorn VPop1000", "Full1000", "Hybrid"][(mode - 2) % 4]
    source_ids = source_members[source_members["source_vpops"].str.contains(source_name, regex=False)]["canonical_id"]
    source_idx = wide.index.get_indexer(source_ids)
    source_idx = source_idx[source_idx >= 0]
    if len(source_idx) >= n_select:
        return rng.choice(source_idx, size=n_select, replace=False)
    chosen = list(source_idx)
    remaining = np.setdiff1d(np.arange(n_pool), np.array(chosen, dtype=int), assume_unique=False)
    fill = rng.choice(remaining, size=n_select - len(chosen), replace=False, p=(weights[remaining] / weights[remaining].sum()))
    return np.array(chosen + fill.tolist(), dtype=int)


def optimize_subset(
    wide: pd.DataFrame,
    arrays: dict[tuple[str, str], np.ndarray],
    source_members: pd.DataFrame,
    ref: pd.DataFrame,
    waterfall_targets: dict[str, tuple[np.ndarray, np.ndarray]],
    n_select: int,
    seed: int,
    starts: int,
    swaps: int,
) -> tuple[np.ndarray, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    n_pool = len(wide)
    weights = scalar_candidate_weights(arrays, n_pool, ref)
    best_idx: np.ndarray | None = None
    best_score = np.inf
    trace_rows = []
    all_indices = np.arange(n_pool)
    for start in range(starts):
        current = make_initial_subset(rng, n_pool, n_select, source_members, wide, weights, start % 6)
        selected = np.zeros(n_pool, dtype=bool)
        selected[current] = True
        current_score = score_subset(arrays, current, ref, waterfall_targets)["total"]
        start_best = current_score
        temp0 = 0.030 if n_select >= 1000 else 0.045
        for it in range(swaps):
            selected_idx = np.flatnonzero(selected)
            unselected_idx = all_indices[~selected]
            rem = int(rng.choice(selected_idx))
            if rng.random() < 0.7:
                p = weights[unselected_idx] / weights[unselected_idx].sum()
                add = int(rng.choice(unselected_idx, p=p))
            else:
                add = int(rng.choice(unselected_idx))
            selected[rem] = False
            selected[add] = True
            proposal = np.flatnonzero(selected)
            proposal_score = score_subset(arrays, proposal, ref, waterfall_targets)["total"]
            temp = temp0 * (1.0 - it / max(swaps - 1, 1)) + 0.002
            accept = proposal_score <= current_score or rng.random() < np.exp((current_score - proposal_score) / temp)
            if accept:
                current_score = proposal_score
                current = proposal
                if proposal_score < best_score:
                    best_score = proposal_score
                    best_idx = proposal.copy()
            else:
                selected[add] = False
                selected[rem] = True
        trace_rows.append({"n_select": n_select, "start": start + 1, "initial_score": start_best, "final_score": current_score, "global_best_score": best_score})
        print(f"N={n_select} start {start+1}/{starts}: final={current_score:.4f} global_best={best_score:.4f}", flush=True)
    if best_idx is None:
        raise RuntimeError("No subset selected.")
    return best_idx, pd.DataFrame(trace_rows)


def selected_long_table(wide: pd.DataFrame, selected_idx: np.ndarray, source_members: pd.DataFrame, label: str) -> pd.DataFrame:
    ids = wide.index[selected_idx].to_numpy()
    records = []
    member_lookup = source_members.set_index("canonical_id").to_dict("index")
    for cid in ids:
        row = member_lookup.get(cid, {})
        for regimen in REGIMEN_ORDER:
            for metric in METRICS:
                records.append(
                    {
                        "vpop": label,
                        "canonical_id": cid,
                        "source_vpops": row.get("source_vpops", ""),
                        "source_candidate_ids": row.get("source_candidate_ids", ""),
                        "regimen": regimen,
                        "regimen_label": REGIMEN_LABELS[regimen],
                        "metric": metric,
                        "value": float(wide.loc[cid, (regimen, metric)]),
                    }
                )
    return pd.DataFrame(records)


def write_membership(wide: pd.DataFrame, selected_idx: np.ndarray, source_members: pd.DataFrame, label: str, scores: dict[str, float]) -> pd.DataFrame:
    ids = wide.index[selected_idx]
    membership = pd.DataFrame({"selected_rank": np.arange(1, len(ids) + 1), "canonical_id": ids})
    membership = membership.merge(source_members, on="canonical_id", how="left")
    for key, value in scores.items():
        membership[f"subset_score_{key}"] = value
    return membership


def plot_score_comparison(score_df: pd.DataFrame) -> list[Path]:
    outputs = []
    plot = score_df.melt(id_vars=["vpop", "n"], value_vars=["IL6", "Tcell", "tumor", "total"], var_name="term", value_name="score")
    order = [x for x in SOURCE_ORDER if x in set(score_df["vpop"])]
    fig, ax = plt.subplots(figsize=(9.5, 4.6))
    terms = ["IL6", "Tcell", "tumor", "total"]
    x = np.arange(len(terms))
    width = 0.13
    offsets = np.linspace(-width * (len(order) - 1) / 2, width * (len(order) - 1) / 2, len(order))
    for off, vpop in zip(offsets, order):
        vals = plot[plot["vpop"] == vpop].set_index("term").reindex(terms)["score"].to_numpy(float)
        ax.bar(x + off, np.maximum(vals, 1e-4), width=width, color=COLORS[vpop], alpha=0.82, label=vpop)
    ax.set_yscale("log")
    ax.set_xticks(x, terms)
    ax.set_ylabel("Endpoint-target score, log scale; lower is better")
    ax.set_title("Endpoint-matched subVPop score comparison")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, ncol=3, fontsize=8)
    fig.tight_layout()
    for ext in [".png", ".pdf"]:
        path = OUT_DIR / f"endpoint_matched_subvpop_score_comparison{ext}"
        fig.savefig(fs_path(path), dpi=240, bbox_inches="tight")
        outputs.append(path)
    plt.close(fig)
    return outputs


def plot_violins(combined_long: pd.DataFrame, ref: pd.DataFrame) -> list[Path]:
    outputs = []
    plot_sources = ["VPop250", "Old full1000", "Sinkhorn VPop1000", "Endpoint-matched VPop1000"]
    positions = {name: 0.7 + i * 0.25 for i, name in enumerate(plot_sources)}
    ref_pos = 0.7 + len(plot_sources) * 0.25
    metric_info = [
        ("profile_peak_il6_pg_ml", "IL6 peak (pg/mL)", "log", (1e-3, 3e4)),
        ("profile_peak_cd69_cd8_pct", "Peak CD69+ CD8+ T cells (%)", "linear", (0, 105)),
        ("day84_tumor_change_pct", "Day-84 tumor change from baseline (%)", "linear", (-120, 320)),
    ]
    fig, axes = plt.subplots(3, 4, figsize=(17.2, 8.8), sharex="col")
    for row, (metric, ylabel, yscale, ylim) in enumerate(metric_info):
        for col, regimen in enumerate(REGIMEN_ORDER):
            ax = axes[row, col]
            for source in plot_sources:
                vals = combined_long[(combined_long["vpop"] == source) & (combined_long["regimen"] == regimen) & (combined_long["metric"] == metric)]["value"].to_numpy(float)
                vals = vals[np.isfinite(vals)]
                if yscale == "log":
                    vals = np.maximum(vals, 1e-12)
                parts = ax.violinplot([vals], positions=[positions[source]], widths=0.18, showmedians=True, showextrema=False)
                for body in parts["bodies"]:
                    body.set_facecolor(COLORS[source])
                    body.set_edgecolor(COLORS[source])
                    body.set_alpha(0.30)
                if "cmedians" in parts:
                    parts["cmedians"].set_color(COLORS[source])
                q05, q50, q95 = np.nanpercentile(vals, [5, 50, 95])
                ax.vlines(positions[source], q05, q95, color=COLORS[source], linewidth=1.1)
                ax.scatter([positions[source]], [q50], s=16, color=COLORS[source], zorder=5)
            target = ref[(ref["regimen"] == regimen) & (ref["metric"] == metric)]
            if not target.empty:
                t = target.iloc[0]
                lo, med, hi = float(t["lower"]), float(t["median"]), float(t["upper"])
                if yscale == "log":
                    lo, med, hi = max(lo, 1e-12), max(med, 1e-12), max(hi, 1e-12)
                ax.vlines(ref_pos, lo, hi, color="#111111", linewidth=1.3)
                ax.hlines(med, ref_pos - 0.08, ref_pos + 0.08, color="#111111", linewidth=2.0)
            ax.set_yscale(yscale)
            ax.set_ylim(ylim)
            ax.set_xlim(0.48, ref_pos + 0.18)
            ax.set_xticks([positions[s] for s in plot_sources] + [ref_pos], ["250", "old\n1000", "Sinkhorn\n1000", "matched\n1000", "Hosseini\nref"], fontsize=8)
            if row == 0:
                ax.set_title(REGIMEN_LABELS[regimen], fontsize=10, fontweight="bold")
            if col == 0:
                ax.set_ylabel(ylabel, fontsize=9)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.grid(False)
    handles = [plt.Line2D([0], [0], color=COLORS[s], marker="s", linestyle="", markersize=8, alpha=0.7, label=s) for s in plot_sources]
    handles.append(plt.Line2D([0], [0], color="#111111", marker="_", linestyle="-", markersize=10, label="Digitized Hosseini reference"))
    fig.legend(handles=handles, frameon=False, ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.015), fontsize=8)
    fig.suptitle("Endpoint-matched VPop1000 selected from pooled Fig. 5 patients", y=1.055, fontsize=15, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    for ext in [".png", ".pdf", ".svg"]:
        path = OUT_DIR / f"endpoint_matched_vpop1000_side_by_side_violins{ext}"
        fig.savefig(fs_path(path), dpi=260, bbox_inches="tight")
        outputs.append(path)
    plt.close(fig)
    return outputs


def main() -> None:
    for path in [SOURCE, REFERENCE, DIGITIZED_WATERFALL]:
        if not path.exists():
            raise FileNotFoundError(path)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    wide, source_members = load_pool()
    arrays = arrays_from_wide(wide)
    ref, waterfall_targets = load_reference()
    print(f"Loaded endpoint pool: {len(wide)} unique patients", flush=True)

    source_scores = source_score_table(wide, arrays, source_members, ref, waterfall_targets)
    outputs = []
    selected_longs = []
    selected_memberships = []
    trace_all = []
    score_rows = [source_scores]
    for n_select, starts, swaps in [(1000, 16, 3500), (250, 12, 3000)]:
        selected_idx, trace = optimize_subset(
            wide,
            arrays,
            source_members,
            ref,
            waterfall_targets,
            n_select=n_select,
            seed=20260517 + n_select,
            starts=starts,
            swaps=swaps,
        )
        label = f"Endpoint-matched VPop{n_select}"
        scores = score_subset(arrays, selected_idx, ref, waterfall_targets)
        score_rows.append(pd.DataFrame([{"vpop": label, "n": n_select, **scores}]))
        selected_long = selected_long_table(wide, selected_idx, source_members, label)
        selected_long.to_csv(OUT_DIR / f"endpoint_matched_vpop{n_select}_endpoint_long.csv", index=False)
        selected_longs.append(selected_long)
        membership = write_membership(wide, selected_idx, source_members, label, scores)
        membership.to_csv(OUT_DIR / f"endpoint_matched_vpop{n_select}_membership.csv", index=False)
        selected_memberships.append(membership)
        trace["label"] = label
        trace_all.append(trace)

    scores_df = pd.concat(score_rows, ignore_index=True)
    scores_df.to_csv(OUT_DIR / "endpoint_matched_score_comparison.csv", index=False)
    pd.concat(trace_all, ignore_index=True).to_csv(OUT_DIR / "endpoint_matched_selection_trace.csv", index=False)
    source_members.to_csv(OUT_DIR / "endpoint_pool_source_membership.csv", index=False)
    outputs.extend(plot_score_comparison(scores_df))

    original_long = pd.read_csv(SOURCE)
    original_long = original_long[original_long["vpop"].isin(["VPop250", "Full1000", "Sinkhorn VPop1000"])].copy()
    original_long["vpop"] = original_long["vpop"].replace({"Full1000": "Old full1000"})
    combined_long = pd.concat([original_long, *selected_longs], ignore_index=True)
    combined_long.to_csv(OUT_DIR / "endpoint_matched_combined_violin_source.csv", index=False)
    outputs.extend(plot_violins(combined_long, ref))

    meta = {
        "created_by": "scripts/subselect_hosseini_fig5_endpoint_matched_vpop.py",
        "source": str(SOURCE),
        "reference": str(REFERENCE),
        "digitized_waterfall": str(DIGITIZED_WATERFALL),
        "pool_unique_patients": int(len(wide)),
        "selection_objective": "Match endpoint IL6 and T-cell p5/median/p95 plus tumor waterfall contour/responder fraction. No trajectory-curve terms beyond endpoint peaks.",
        "selected_sizes": [1000, 250],
        "outputs": [str(p) for p in outputs],
        "caveat": "This is endpoint-metric sub-selection from pooled simulated patients. It is not yet a new coherent parameter sampling distribution unless membership is restricted to parameter-complete sources.",
    }
    (OUT_DIR / "endpoint_matched_subvpop_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(scores_df.to_string(index=False))
    print(f"Wrote endpoint-matched outputs to {OUT_DIR}")


if __name__ == "__main__":
    main()
