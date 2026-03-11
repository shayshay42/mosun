using CSV
using DataFrames
using DiffEqCallbacks
using DifferentialEquations
using JSON3
using SciMLBase
using Statistics

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

const MMC = TCellEngagerQSP.MosunModelCore
const REPO_ROOT = TCellEngagerQSP.REPO_ROOT

function parse_float_list(txt::AbstractString)
    out = Float64[]
    for tok in split(String(txt), ",")
        s = strip(tok)
        isempty(s) && continue
        push!(out, parse(Float64, s))
    end
    return out
end

function max_norm_abs_err(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real})
    @assert size(X) == size(Y)
    mask = isfinite.(X) .& isfinite.(Y)
    any(mask) || return Inf
    amp = max(maximum(abs.(Y[mask])), 1e-12)
    return maximum(abs.(X[mask] .- Y[mask])) / amp
end

function build_base_param_map()
    p = joinpath(REPO_ROOT, "generated", "phase1_design_standard", "dlbcl_param_overrides.csv")
    df = DataFrame(CSV.File(p))
    out = Dict{String, Float64}()
    for r in eachrow(df)
        out[String(r.name)] = Float64(r.value)
    end
    return out
end

function row_to_param_overrides(row, cols)
    out = Dict{String, Float64}()
    for cname in cols
        cname == "patient_id" && continue
        raw = row[Symbol(cname)]
        if !ismissing(raw)
            v = try
                Float64(raw)
            catch
                NaN
            end
            if isfinite(v)
                out[cname] = v
            end
        end
    end
    return out
end

function merged_param_map(base::Dict{String,Float64}, overrides::Dict{String,Float64}, end_time::Float64)
    p = copy(base)
    merge!(p, overrides)
    p["PKflag"] = 1.0
    p["fvalidation"] = 0.0
    p["VPid"] = 1.0
    p["end_time"] = end_time
    return p
end

function solve_realdata_regimen_direct(pmap::Dict{String,Float64}, dose_map::Dict{Float64,Float64}, horizon_days::Float64)
    params = MMC.params_from_dict(pmap; ignore_unknown = true)
    regimen = MMC.bolus_regimen(:TDBc_ugperkg, dose_map)
    _, sol = MMC.solve_regimen(
        regimen,
        params,
        TCellEngagerQSP.make_solver_alg("cvode_bdf");
        tspan = (0.0, horizon_days),
        callback_mode = :callback,
        abstol = TCellEngagerQSP.SOLVER_ABSTOL,
        reltol = TCellEngagerQSP.SOLVER_RELTOL,
    )
    sol.retcode == SciMLBase.ReturnCode.Success || error("direct solve failed with retcode=$(sol.retcode)")
    return params, sol
end

function solve_realdata_regimen_legacy(pmap::Dict{String,Float64}, dose_map::Dict{Float64,Float64}, horizon_days::Float64)
    pnames = sort(collect(keys(pmap)))
    pvals = [pmap[n] for n in pnames]
    mdl = TCellEngagerQSP.build_model_with_variant_ids(Int[], pnames, pvals; core_mode_override = :legacy_reference)
    target_idx = TCellEngagerQSP.canonical_state_to_idx(mdl)["TDBc_ugperkg"]
    local_dose_map = copy(dose_map)
    dose_at_t0 = get(local_dose_map, 0.0, 0.0)
    if dose_at_t0 != 0.0
        local_dose_map[0.0] = 0.0
    end
    cb_times = sort([t for t in keys(local_dose_map) if t > 0.0])

    function affect!(integrator)
        t = Float64(integrator.t)
        amt = get(local_dose_map, t, NaN)
        if isnan(amt)
            for (tt, vv) in local_dose_map
                if isapprox(t, tt; atol = 1e-8, rtol = 0.0)
                    amt = vv
                    break
                end
            end
        end
        if !isnan(amt)
            integrator.u[target_idx] += amt
        end
    end

    u0 = copy(TCellEngagerQSP.canonical_u0(mdl))
    if dose_at_t0 != 0.0
        u0[target_idx] += dose_at_t0
    end
    ctx = TCellEngagerQSP.CanonicalSimContext(copy(mdl.z), mdl.pvals, TCellEngagerQSP.InfusionEvent[])
    prob = ODEProblem(mdl.canonical_rhs, u0, (0.0, horizon_days), ctx)
    cb = PresetTimeCallback(cb_times, affect!; save_positions = (true, true))
    sol = solve(
        prob,
        TCellEngagerQSP.make_solver_alg("cvode_bdf");
        abstol = TCellEngagerQSP.SOLVER_ABSTOL,
        reltol = TCellEngagerQSP.SOLVER_RELTOL,
        callback = cb,
        tstops = cb_times,
        d_discontinuities = cb_times,
    )
    sol.retcode == SciMLBase.ReturnCode.Success || error("legacy solve failed with retcode=$(sol.retcode)")
    return mdl, sol
end

function sample_direct(params::MMC.MosunParams, sol, tgrid::Vector{Float64}, names::Vector{Symbol})
    cache = MMC.zero_observables_cache()
    X = Matrix{Float64}(undef, length(tgrid), length(names))
    for (i, t) in enumerate(tgrid)
        u = Float64.(sol(t))
        for (j, nm) in enumerate(names)
            X[i, j] = MMC.value_at(u, params, t, nm, cache)
        end
    end
    return X
end

function sample_legacy(mdl, sol, tgrid::Vector{Float64}, names::Vector{String})
    X = Matrix{Float64}(undef, length(tgrid), length(names))
    for (i, t) in enumerate(tgrid)
        u = Float64.(sol(t))
        for (j, nm) in enumerate(names)
            X[i, j] = TCellEngagerQSP.eval_symbol(mdl, u, t, nm)
        end
    end
    return X
end

function best_spd_from_bt(bt::AbstractVector{<:Real})
    bt0 = max(Float64(bt[1]), 1e-12)
    spd_pct = @. 100.0 * (Float64(bt) / bt0 - 1.0)
    spd_eval = length(spd_pct) > 1 ? spd_pct[2:end] : spd_pct
    return minimum(spd_eval)
end

function realdata_regimen_spec()
    ref = DataFrame(CSV.File(joinpath(REPO_ROOT, "generated", "reference", "musun_2022_88patients_SPD_from_vpop_generation_3a50cd1.csv")))
    row = ref[1, :]
    dose_days = parse_float_list(String(row.dose_days))
    dose_mg = parse_float_list(String(row.dose_mg))
    bw_kg = 70.0
    dose_map = Dict{Float64, Float64}()
    for (t, mg) in zip(dose_days, dose_mg)
        dose_map[t] = get(dose_map, t, 0.0) + mg * 1000.0 / bw_kg
    end
    return (
        scenario = String(row.scenario),
        regimen_label = String(row.regimen_label),
        dose_map = dose_map,
        horizon_days = Float64(row.horizon_days),
    )
end

function build_realdata_cases()
    base = build_base_param_map()
    phase1_patients = DataFrame(CSV.File(joinpath(REPO_ROOT, "generated", "phase1_design", "patients.csv")))
    selected = DataFrame(CSV.File(joinpath(REPO_ROOT, "generated", "vpop_pruning", "vpop_wasserstein_spd88_postopt_20260306", "selected_patients.csv")))
    spec = realdata_regimen_spec()

    cases = NamedTuple[]
    push!(cases, (
        case_key = "realdata_regimen_dlbcl_baseline",
        param_source = "phase1_design_standard/dlbcl_param_overrides",
        pmap = merged_param_map(base, Dict{String, Float64}(), spec.horizon_days),
        dose_map = spec.dose_map,
        horizon_days = spec.horizon_days,
        regimen_label = spec.regimen_label,
    ))

    prow = phase1_patients[1, :]
    push!(cases, (
        case_key = "realdata_regimen_phase1_patient1",
        param_source = "generated/phase1_design/patients.csv:1",
        pmap = merged_param_map(base, row_to_param_overrides(prow, names(phase1_patients)), spec.horizon_days),
        dose_map = spec.dose_map,
        horizon_days = spec.horizon_days,
        regimen_label = spec.regimen_label,
    ))

    sel_rows = [1, cld(nrow(selected), 2), nrow(selected)]
    for idx in sel_rows
        srow = selected[idx, :]
        pid = Int(round(Float64(srow.patient_id)))
        push!(cases, (
            case_key = "realdata_regimen_selected_patient_$(pid)",
            param_source = "generated/vpop_pruning/vpop_wasserstein_spd88_postopt_20260306/selected_patients.csv:$(idx + 1)",
            pmap = merged_param_map(base, row_to_param_overrides(srow, names(selected)), spec.horizon_days),
            dose_map = spec.dose_map,
            horizon_days = spec.horizon_days,
            regimen_label = spec.regimen_label,
        ))
    end

    return cases
end

function main()
    out_dir_env = get(ENV, "PROD_LEGACY_COMPARE_OUT_DIR", joinpath("generated", "benchmarks", "production_vs_legacy_detailed_20260310"))
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(REPO_ROOT, out_dir_env)
    mkpath(out_dir)

    summary = DataFrame(
        case_key = String[],
        case_group = String[],
        regimen_label = String[],
        param_source = String[],
        note = String[],
        n_times = Int[],
        max_rel_err = Float64[],
        mean_rel_err = Float64[],
        best_spd_direct = Float64[],
        best_spd_legacy = Float64[],
        best_spd_abs_diff = Float64[],
    )
    trajectories = DataFrame(
        case_key = String[],
        time = Float64[],
        output = String[],
        direct = Float64[],
        legacy = Float64[],
        abs_diff = Float64[],
    )

    manifest = DataFrame(CSV.File(joinpath(REPO_ROOT, "generated", "matlab_reference", "manifest_three_cases.csv")))
    for row in eachrow(manifest)
        save_times, outvec, X_direct, ref_df = TCellEngagerQSP.run_manifest_row(row; engine = :canonical, core_mode_override = :production)
        _, _, X_legacy, _ = TCellEngagerQSP.run_manifest_row(row; engine = :canonical, core_mode_override = :legacy_reference)
        rel_err = max_norm_abs_err(X_direct, X_legacy)
        mean_rel = mean(abs.(X_direct .- X_legacy)) / max(maximum(abs.(X_legacy)), 1e-12)
        push!(summary, (
            String(row.sim_key),
            "matlab_reference",
            String(row.schedule),
            "manifest_case_$(row.case_no)",
            Float64(row.case_no) in (10.0, 20.0) ? "Expected mismatch: legacy PK lookup branch (PKflag=0) is intentionally removed from production core." : "Production-relevant manifest case (PKflag=1).",
            length(save_times),
            rel_err,
            mean_rel,
            NaN,
            NaN,
            NaN,
        ))
        for (j, nm) in enumerate(outvec)
            for i in eachindex(save_times)
                push!(trajectories, (
                    String(row.sim_key),
                    Float64(save_times[i]),
                    String(nm),
                    Float64(X_direct[i, j]),
                    Float64(X_legacy[i, j]),
                    abs(Float64(X_direct[i, j] - X_legacy[i, j])),
                ))
            end
        end
    end

    output_syms = [:Btumor, :IL6combo, :TDBc_ugperml, :Bpb_perml, :totTpb_perml, :Tafraction_pb]
    output_names = String.(output_syms)
    for case in build_realdata_cases()
        params, sol_direct = solve_realdata_regimen_direct(case.pmap, case.dose_map, case.horizon_days)
        mdl, sol_legacy = solve_realdata_regimen_legacy(case.pmap, case.dose_map, case.horizon_days)
        tgrid = collect(0.0:1.0:case.horizon_days)
        X_direct = sample_direct(params, sol_direct, tgrid, output_syms)
        X_legacy = sample_legacy(mdl, sol_legacy, tgrid, output_names)
        bt_direct = X_direct[:, 1]
        bt_legacy = X_legacy[:, 1]
        spd_direct = best_spd_from_bt(bt_direct)
        spd_legacy = best_spd_from_bt(bt_legacy)
        rel_err = max_norm_abs_err(X_direct, X_legacy)
        mean_rel = mean(abs.(X_direct .- X_legacy)) / max(maximum(abs.(X_legacy)), 1e-12)
        push!(summary, (
            case.case_key,
            "realdata_regimen",
            case.regimen_label,
            case.param_source,
            "Production-relevant real-data regimen (PKflag=1).",
            length(tgrid),
            rel_err,
            mean_rel,
            spd_direct,
            spd_legacy,
            abs(spd_direct - spd_legacy),
        ))
        for (j, nm) in enumerate(output_names)
            for i in eachindex(tgrid)
                push!(trajectories, (
                    case.case_key,
                    tgrid[i],
                    nm,
                    Float64(X_direct[i, j]),
                    Float64(X_legacy[i, j]),
                    abs(Float64(X_direct[i, j] - X_legacy[i, j])),
                ))
            end
        end
    end

    sort!(summary, [:case_group, :case_key])
    CSV.write(joinpath(out_dir, "summary.csv"), summary)
    CSV.write(joinpath(out_dir, "trajectories.csv"), trajectories)

    open(joinpath(out_dir, "run_manifest.json"), "w") do io
        JSON3.pretty(
            io,
            Dict(
                "out_dir" => out_dir,
                "realdata_regimen_source" => joinpath(REPO_ROOT, "generated", "reference", "musun_2022_88patients_SPD_from_vpop_generation_3a50cd1.csv"),
                "selected_patients_source" => joinpath(REPO_ROOT, "generated", "vpop_pruning", "vpop_wasserstein_spd88_postopt_20260306", "selected_patients.csv"),
                "phase1_patients_source" => joinpath(REPO_ROOT, "generated", "phase1_design", "patients.csv"),
                "matlab_manifest_source" => joinpath(REPO_ROOT, "generated", "matlab_reference", "manifest_three_cases.csv"),
            ),
        )
    end

    println("Wrote $(joinpath(out_dir, "summary.csv"))")
    println("Wrote $(joinpath(out_dir, "trajectories.csv"))")
end

main()
