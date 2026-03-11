using CSV
using DataFrames
using DifferentialEquations
using DiffEqCallbacks
using Sundials
using SciMLBase

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

const MMC = TCellEngagerQSP.MosunModelCore

repo_root = TCellEngagerQSP.REPO_ROOT
design_dir_default = joinpath(repo_root, "generated", "phase1_design")
design_dir_env = get(ENV, "PHASE1_DESIGN_DIR", design_dir_default)
design_dir = isabspath(design_dir_env) ? design_dir_env : joinpath(repo_root, design_dir_env)
patients_path = joinpath(design_dir, "patients.csv")
reg_path = joinpath(design_dir, "regimen_events.csv")

if !isfile(patients_path) || !isfile(reg_path)
    error("Design files missing. Run scripts/generate_phase1_design.py first.")
end

patients = DataFrame(CSV.File(patients_path))
reg_events = DataFrame(CSV.File(reg_path))

out_dir_default = joinpath(repo_root, "generated", "phase1_julia")
out_dir_env = get(ENV, "PHASE1_JULIA_OUT_DIR", out_dir_default)
out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(repo_root, out_dir_env)
mkpath(out_dir)
save_trajectories = get(ENV, "PHASE1_SAVE_TRAJECTORIES", "0") == "1"
trace_out_dir_default = joinpath(out_dir, "traces")
trace_out_dir_env = get(ENV, "PHASE1_TRACE_OUT_DIR", trace_out_dir_default)
trace_out_dir = isabspath(trace_out_dir_env) ? trace_out_dir_env : joinpath(repo_root, trace_out_dir_env)
trace_dt_days = parse(Float64, get(ENV, "PHASE1_TRACE_DT_DAYS", "0.1"))
if save_trajectories
    mkpath(trace_out_dir)
end
horizon_days = parse(Float64, get(ENV, "PHASE1_HORIZON_DAYS", "84.0"))
if horizon_days < 42.0
    error("PHASE1_HORIZON_DAYS must be >= 42.0 to support day-42 endpoint extraction.")
end

phase1_active_variant_ids = [5, 9, 14, 20, 24, 25, 27, 28]
variants_mode = lowercase(get(ENV, "PHASE1_VARIANTS_MODE", "matlab_empty"))
variant_ids =
    if variants_mode == "matlab_empty"
        # Match MATLAB script behavior in run_phase1_mosu_matlab_reference.m:
        # sbiosimulate(model, cs, [], d) -> no variants applied.
        Int[]
    elseif variants_mode == "active_phase1"
        phase1_active_variant_ids
    else
        error("Unsupported PHASE1_VARIANTS_MODE=$variants_mode (use matlab_empty or active_phase1)")
    end
bw_kg_for_dose_conversion = 70.0
dlbcl_overrides_path = joinpath(design_dir, "dlbcl_param_overrides.csv")
dlbcl_overrides = isfile(dlbcl_overrides_path) ? DataFrame(CSV.File(dlbcl_overrides_path)) : DataFrame(name = String[], value = Float64[])

solver_name_full = lowercase(get(ENV, "TCE_SOLVER", "cvode_bdf"))
use_mtk_jac = endswith(solver_name_full, "_mtkjac") || get(ENV, "TCE_MTK_JAC", "0") == "1"
use_mtk_sparse_jac = get(ENV, "TCE_MTK_JAC_SPARSE", "0") == "1"
solver_name = replace(solver_name_full, "_mtkjac" => "")
hybrid_mode = startswith(solver_name, "hybrid")
alg = hybrid_mode ? nothing : TCellEngagerQSP.make_solver_alg(solver_name)
hybrid_stiff_alg = hybrid_mode ? TCellEngagerQSP.make_solver_alg(get(ENV, "TCE_HYBRID_STIFF_SOLVER", "cvode_bdf")) : nothing
hybrid_nonstiff_alg = hybrid_mode ? TCellEngagerQSP.make_solver_alg(get(ENV, "TCE_HYBRID_NONSTIFF_SOLVER", "tsit5")) : nothing
hybrid_stiff_window_days = parse(Float64, get(ENV, "TCE_HYBRID_STIFF_WINDOW_DAYS", "2.0"))
hybrid_nonstiff_fixed_dt_days = parse(Float64, get(ENV, "TCE_HYBRID_NONSTIFF_FIXED_DT_DAYS", "0.0"))
post_dose_proposed_dt_days = parse(Float64, get(ENV, "TCE_POST_DOSE_PROPOSED_DT_DAYS", "0.0"))
engine = Symbol(lowercase(get(ENV, "TCE_ENGINE", "canonical")))
if !(engine in (:legacy, :canonical))
    error("Unsupported TCE_ENGINE=$engine. Use legacy or canonical.")
end

function interp1_linear(x::Vector{Float64}, y::Vector{Float64}, tq::Float64)
    if tq <= x[1]
        return y[1]
    elseif tq >= x[end]
        return y[end]
    end
    i = searchsortedlast(x, tq)
    x1, x2 = x[i], x[i + 1]
    y1, y2 = y[i], y[i + 1]
    if x2 == x1
        return y1
    end
    return y1 + (y2 - y1) * (tq - x1) / (x2 - x1)
end

function window_peak(t::Vector{Float64}, y::Vector{Float64}, t0::Float64, t1::Float64)
    if isempty(t) || isempty(y)
        return NaN
    end
    if t1 < t0
        t0, t1 = t1, t0
    end
    vals = Float64[]
    push!(vals, interp1_linear(t, y, t0))
    push!(vals, interp1_linear(t, y, t1))
    for i in eachindex(t)
        ti = t[i]
        if t0 <= ti <= t1
            push!(vals, y[i])
        end
    end
    isempty(vals) && return NaN
    return maximum(vals)
end

function window_auc(t::Vector{Float64}, y::Vector{Float64}, t0::Float64, t1::Float64)
    if isempty(t) || isempty(y)
        return NaN
    end
    if t1 < t0
        t0, t1 = t1, t0
    end
    if t1 == t0
        return 0.0
    end

    tw = Float64[t0]
    yw = Float64[interp1_linear(t, y, t0)]
    for i in eachindex(t)
        ti = t[i]
        if t0 < ti < t1
            push!(tw, ti)
            push!(yw, y[i])
        end
    end
    push!(tw, t1)
    push!(yw, interp1_linear(t, y, t1))

    ord = sortperm(tw)
    tw = tw[ord]
    yw = yw[ord]

    auc = 0.0
    for i in 1:(length(tw) - 1)
        dt = tw[i + 1] - tw[i]
        if dt > 0.0
            auc += 0.5 * dt * (yw[i] + yw[i + 1])
        end
    end
    return auc
end

function append_solution!(sol_t::Vector{Float64}, sol_u::Vector{Vector{Float64}}, ts::Vector{Float64}, us::Vector)
    isempty(ts) && return
    if isempty(sol_t)
        append!(sol_t, Float64.(ts))
        for u in us
            push!(sol_u, Float64.(u))
        end
        return
    end
    start_idx = isapprox(ts[1], sol_t[end]; atol = 1e-12, rtol = 0.0) ? 2 : 1
    for i in start_idx:length(ts)
        push!(sol_t, Float64(ts[i]))
        push!(sol_u, Float64.(us[i]))
    end
end

function run_segment!(
        sol_t::Vector{Float64},
        sol_u::Vector{Vector{Float64}},
        ode_rhs,
        ctx,
        u_start::Vector{Float64},
        t0::Float64,
        t1::Float64,
        seg_alg;
        fixed_dt_days::Float64 = 0.0)
    if t1 <= t0
        return u_start
    end
    prob = ODEProblem(ode_rhs, u_start, (t0, t1), ctx)
    if fixed_dt_days > 0.0
        sol = solve(
            prob,
            seg_alg;
            abstol = TCellEngagerQSP.SOLVER_ABSTOL,
            reltol = TCellEngagerQSP.SOLVER_RELTOL,
            tstops = [t1],
            adaptive = false,
            dt = fixed_dt_days,
        )
    else
        sol = solve(
            prob,
            seg_alg;
            abstol = TCellEngagerQSP.SOLVER_ABSTOL,
            reltol = TCellEngagerQSP.SOLVER_RELTOL,
            tstops = [t1],
        )
    end
    if sol.retcode != SciMLBase.ReturnCode.Success
        error("segment solve failed with retcode=$(sol.retcode)")
    end
    append_solution!(sol_t, sol_u, sol.t, sol.u)
    return Float64.(sol.u[end])
end

function run_hybrid_segmented(
        ode_rhs,
        ctx,
        u0::Vector{Float64},
        dose_map::Dict{Float64,Float64},
        target_idx::Int,
        tf::Float64,
        stiff_alg,
        nonstiff_alg,
        stiff_window_days::Float64,
        nonstiff_fixed_dt_days::Float64)
    local_dose_map = copy(dose_map)
    dose_at_t0 = get(local_dose_map, 0.0, 0.0)
    if dose_at_t0 != 0.0
        local_dose_map[0.0] = 0.0
    end

    u_curr = copy(u0)
    if dose_at_t0 != 0.0
        u_curr[target_idx] += dose_at_t0
    end

    dose_times = sort([t for t in keys(local_dose_map) if t > 0.0])
    segment_ends = sort(unique(vcat(dose_times, [tf])))
    stiff_starts = Set(vcat([0.0], dose_times))

    sol_t = Float64[]
    sol_u = Vector{Vector{Float64}}()
    t_curr = 0.0

    for te in segment_ends
        if te > t_curr
            stiff_here = any(isapprox(t_curr, ts; atol = 1e-10, rtol = 0.0) for ts in stiff_starts)
            if stiff_here && stiff_window_days > 0.0
                t_mid = min(te, t_curr + stiff_window_days)
                u_curr = run_segment!(sol_t, sol_u, ode_rhs, ctx, u_curr, t_curr, t_mid, stiff_alg)
                if te > t_mid
                    u_curr = run_segment!(
                        sol_t,
                        sol_u,
                        ode_rhs,
                        ctx,
                        u_curr,
                        t_mid,
                        te,
                        nonstiff_alg;
                        fixed_dt_days = nonstiff_fixed_dt_days,
                    )
                end
            else
                u_curr = run_segment!(
                    sol_t,
                    sol_u,
                    ode_rhs,
                    ctx,
                    u_curr,
                    t_curr,
                    te,
                    nonstiff_alg;
                    fixed_dt_days = nonstiff_fixed_dt_days,
                )
            end
        end
        if haskey(local_dose_map, te)
            u_curr[target_idx] += local_dose_map[te]
        end
        t_curr = te
    end
    return sol_t, sol_u
end

rows = DataFrame(
    regimen = String[],
    regimen_type = String[],
    patient_id = Int[],
    il6_peak_0_2 = Float64[],
    il6_peak_0_21 = Float64[],
    il6_auc_0_2 = Float64[],
    il6_auc_0_21 = Float64[],
    btumor_t0 = Float64[],
    btumor_day42 = Float64[],
    tumor_resid_day42 = Float64[],
    tumor_cfbl_day42 = Float64[],
    best_spd_pct = Float64[],
    best_spd_le_minus90 = Float64[],
    best_spd_le_minus50 = Float64[],
    best_spd_gt_zero = Float64[],
    best_spd_floor_hit = Float64[],
    loss_raw = Float64[],
    bt_ratio_tumor_init = Float64[],
    status = String[],
)
trace_manifest = DataFrame(
    regimen = String[],
    regimen_type = String[],
    patient_id = Int[],
    trace_csv = String[],
    n_timepoints = Int[],
    status = String[],
)

regimens = unique(reg_events.regimen)
max_regimens = parse(Int, get(ENV, "PHASE1_MAX_REGIMENS", string(length(regimens))))
max_patients = parse(Int, get(ENV, "PHASE1_MAX_PATIENTS", string(nrow(patients))))
regimens = regimens[1:min(end, max_regimens)]
for reg in regimens
    sub = reg_events[reg_events.regimen .== reg, :]
    sort!(sub, :event_idx)

    dose_map = Dict{Float64, Float64}()
    for rr in eachrow(sub)
        t = Float64(rr.time_day)
        amt = Float64(rr.dose_mg) * 1000.0 / bw_kg_for_dose_conversion
        dose_map[t] = get(dose_map, t, 0.0) + amt
    end
    dose_times = sort(collect(keys(dose_map)))
    reg_type = String(sub.regimen_type[1])

    for prow in eachrow(patients[1:min(end, max_patients), :])
        pmap = Dict{String, Float64}()
        for rr in eachrow(dlbcl_overrides)
            pmap[String(rr.name)] = Float64(rr.value)
        end
        for cname in names(patients)
            if cname == "patient_id"
                continue
            end
            raw = prow[cname]
            if !ismissing(raw)
                v = try
                    Float64(raw)
                catch
                    NaN
                end
                if isfinite(v)
                    pmap[cname] = v
                end
            end
        end
        pmap["PKflag"] = 1.0
        pmap["fvalidation"] = 0.0
        pmap["VPid"] = 1.0
        pmap["end_time"] = horizon_days
        pnames = sort(collect(keys(pmap)))
        pvals = [pmap[n] for n in pnames]

        try
            il6_series = Float64[]
            bt_series = Float64[]
            t_series = Float64[]
            use_direct_core = engine == :canonical && isempty(variant_ids) && !use_mtk_jac && !hybrid_mode && post_dose_proposed_dt_days == 0.0

            if use_direct_core
                params = MMC.params_from_dict(pmap; ignore_unknown = true)
                regimen = MMC.bolus_regimen(:TDBc_ugperkg, dose_map)
                _, sol = MMC.solve_regimen(
                    regimen,
                    params,
                    alg;
                    tspan = (0.0, horizon_days),
                    callback_mode = :callback,
                    abstol = TCellEngagerQSP.SOLVER_ABSTOL,
                    reltol = TCellEngagerQSP.SOLVER_RELTOL,
                )
                if sol.retcode != SciMLBase.ReturnCode.Success
                    error("solve failed with retcode=$(sol.retcode)")
                end
                cache = MMC.zero_observables_cache()
                for i in eachindex(sol.t)
                    t = Float64(sol.t[i])
                    u = Float64.(sol.u[i])
                    push!(t_series, t)
                    push!(il6_series, MMC.value_at(u, params, t, :IL6combo, cache))
                    push!(bt_series, MMC.value_at(u, params, t, :Btumor, cache))
                end
            else
                mdl = TCellEngagerQSP.build_model_with_variant_ids(variant_ids, pnames, pvals)
                rhs_fun = TCellEngagerQSP.rhs!
                ctx = nothing
                if engine == :legacy
                    ctx = TCellEngagerQSP.make_legacy_context(mdl, TCellEngagerQSP.InfusionEvent[])
                    rhs_fun = TCellEngagerQSP.rhs!
                else
                    ctx = TCellEngagerQSP.CanonicalSimContext(
                        copy(mdl.z),
                        mdl.pvals,
                        TCellEngagerQSP.InfusionEvent[],
                    )
                    rhs_fun = mdl.canonical_rhs
                end
                jac_fun = nothing
                jac_prototype = nothing
                if use_mtk_jac
                    if use_mtk_sparse_jac
                        jac_fun = TCellEngagerQSP.make_mtk_sparse_jacobian(mdl)
                        jac_prototype = TCellEngagerQSP.make_mtk_sparse_jacobian_prototype(mdl)
                    else
                        jac_fun = TCellEngagerQSP.make_mtk_dense_jacobian(mdl)
                    end
                end
                ode_rhs =
                    if isnothing(jac_fun)
                        rhs_fun
                    elseif isnothing(jac_prototype)
                        ODEFunction(rhs_fun; jac = jac_fun)
                    else
                        ODEFunction(rhs_fun; jac = jac_fun, jac_prototype = jac_prototype)
                    end

                target_idx = TCellEngagerQSP.canonical_state_to_idx(mdl)["TDBc_ugperkg"]

                local_dose_map = copy(dose_map)
                dose_at_t0 = get(local_dose_map, 0.0, 0.0)
                if dose_at_t0 != 0.0
                    local_dose_map[0.0] = 0.0
                end
                cb_times = filter(t -> t > 0.0, dose_times)

                function affect!(integrator)
                    t = integrator.t
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
                        if post_dose_proposed_dt_days > 0.0
                            SciMLBase.set_proposed_dt!(integrator, post_dose_proposed_dt_days)
                        end
                    end
                end

                u0 = copy(TCellEngagerQSP.canonical_u0(mdl))
                if dose_at_t0 != 0.0
                    u0[target_idx] += dose_at_t0
                end
                if hybrid_mode
                    sol_t, sol_u = run_hybrid_segmented(
                        ode_rhs,
                        ctx,
                        u0,
                        local_dose_map,
                        target_idx,
                        horizon_days,
                        hybrid_stiff_alg,
                        hybrid_nonstiff_alg,
                        hybrid_stiff_window_days,
                        hybrid_nonstiff_fixed_dt_days,
                    )
                    for i in eachindex(sol_t)
                        t = Float64(sol_t[i])
                        u = Float64.(sol_u[i])
                        push!(t_series, t)
                        push!(il6_series, TCellEngagerQSP.eval_symbol(mdl, u, t, "IL6combo"))
                        push!(bt_series, TCellEngagerQSP.eval_symbol(mdl, u, t, "Btumor"))
                    end
                else
                    cb = PresetTimeCallback(cb_times, affect!; save_positions = (true, true))
                    prob = ODEProblem(ode_rhs, u0, (0.0, horizon_days), ctx)
                    sol = solve(
                        prob,
                        alg;
                        abstol = TCellEngagerQSP.SOLVER_ABSTOL,
                        reltol = TCellEngagerQSP.SOLVER_RELTOL,
                        callback = cb,
                        tstops = cb_times,
                        d_discontinuities = cb_times,
                    )
                    if sol.retcode != SciMLBase.ReturnCode.Success
                        error("solve failed with retcode=$(sol.retcode)")
                    end
                    for i in eachindex(sol.t)
                        t = Float64(sol.t[i])
                        u = Float64.(sol.u[i])
                        push!(t_series, t)
                        push!(il6_series, TCellEngagerQSP.eval_symbol(mdl, u, t, "IL6combo"))
                        push!(bt_series, TCellEngagerQSP.eval_symbol(mdl, u, t, "Btumor"))
                    end
                end
            end

            bt0 = bt_series[1]
            tuniq = Float64[]
            btuniq = Float64[]
            il6uniq = Float64[]
            seen = Set{Float64}()
            for (t, b, il6v) in zip(t_series, bt_series, il6_series)
                if !(t in seen)
                    push!(seen, t)
                    push!(tuniq, t)
                    push!(btuniq, b)
                    push!(il6uniq, il6v)
                end
            end
            il6_peak_0_2 = window_peak(tuniq, il6uniq, 0.0, 2.0)
            il6_peak_0_21 = window_peak(tuniq, il6uniq, 0.0, 21.0)
            il6_auc_0_2 = window_auc(tuniq, il6uniq, 0.0, 2.0)
            il6_auc_0_21 = window_auc(tuniq, il6uniq, 0.0, 21.0)
            bt42 = interp1_linear(tuniq, btuniq, 42.0)
            tumor_resid_day42 = bt42 / max(bt0, 1e-12)
            tumor_cfbl_day42 = (bt42 - bt0) / max(bt0, 1e-12)
            spd_pct = @. 100.0 * (btuniq / max(bt0, 1e-12) - 1.0)
            spd_eval = length(spd_pct) > 1 ? spd_pct[2:end] : spd_pct
            best_spd_pct = minimum(spd_eval)
            best_spd_le_minus90 = best_spd_pct <= -90.0 ? 1.0 : 0.0
            best_spd_le_minus50 = best_spd_pct <= -50.0 ? 1.0 : 0.0
            best_spd_gt_zero = best_spd_pct > 0.0 ? 1.0 : 0.0
            best_spd_floor_hit = best_spd_pct <= -99.9 ? 1.0 : 0.0
            loss_raw = 0.5 * log10(1 + max(il6_peak_0_2, 0.0)) + 0.5 * tumor_resid_day42

            if save_trajectories
                patient_id = Int(prow.patient_id)
                safe_reg = replace(String(reg), r"[^A-Za-z0-9_]" => "_")
                trace_csv = joinpath(trace_out_dir, "regimen_$(safe_reg)_patient_$(patient_id).csv")
                if trace_dt_days > 0.0
                    ttrace = collect(0.0:trace_dt_days:horizon_days)
                    if isempty(ttrace) || !isapprox(ttrace[end], horizon_days; atol = 1e-10, rtol = 0.0)
                        push!(ttrace, horizon_days)
                    end
                    il6trace = [interp1_linear(tuniq, il6uniq, tq) for tq in ttrace]
                    bttrace = [interp1_linear(tuniq, btuniq, tq) for tq in ttrace]
                    trace_df = DataFrame(
                        time = ttrace,
                        IL6combo = il6trace,
                        Btumor = bttrace,
                        regimen = fill(String(reg), length(ttrace)),
                        patient_id = fill(patient_id, length(ttrace)),
                    )
                    CSV.write(trace_csv, trace_df)
                    push!(trace_manifest, (String(reg), reg_type, patient_id, trace_csv, length(ttrace), "ok"))
                else
                    trace_df = DataFrame(
                        time = tuniq,
                        IL6combo = il6uniq,
                        Btumor = btuniq,
                        regimen = fill(String(reg), length(tuniq)),
                        patient_id = fill(patient_id, length(tuniq)),
                    )
                    CSV.write(trace_csv, trace_df)
                    push!(trace_manifest, (String(reg), reg_type, patient_id, trace_csv, length(tuniq), "ok"))
                end
            end

            push!(
                rows,
                (
                    String(reg),
                    reg_type,
                    Int(prow.patient_id),
                    il6_peak_0_2,
                    il6_peak_0_21,
                    il6_auc_0_2,
                    il6_auc_0_21,
                    bt0,
                    bt42,
                    tumor_resid_day42,
                    tumor_cfbl_day42,
                    best_spd_pct,
                    best_spd_le_minus90,
                    best_spd_le_minus50,
                    best_spd_gt_zero,
                    best_spd_floor_hit,
                    loss_raw,
                    Float64(prow.BT_ratio_tumor_init),
                    "ok",
                ),
            )
        catch err
            if save_trajectories
                patient_id = Int(prow.patient_id)
                safe_reg = replace(String(reg), r"[^A-Za-z0-9_]" => "_")
                trace_csv = joinpath(trace_out_dir, "regimen_$(safe_reg)_patient_$(patient_id).csv")
                push!(trace_manifest, (String(reg), reg_type, patient_id, trace_csv, 0, "error"))
            end
            push!(
                rows,
                (
                    String(reg),
                    reg_type,
                    Int(prow.patient_id),
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    Float64(prow.BT_ratio_tumor_init),
                    "error",
                ),
            )
            @warn "Simulation failed" regimen = reg patient = Int(prow.patient_id) err
        end
    end
end

out_path = joinpath(out_dir, "phase1_metrics_julia.csv")
CSV.write(out_path, rows)
println("Wrote " * out_path)
if save_trajectories
    trace_manifest_path = joinpath(out_dir, "phase1_trace_manifest.csv")
    CSV.write(trace_manifest_path, trace_manifest)
    println("Wrote " * trace_manifest_path)
end
