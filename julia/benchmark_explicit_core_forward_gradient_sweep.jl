using CSV
using DataFrames
using DifferentialEquations
using DiffEqCallbacks
using ForwardDiff
using SciMLBase
using Sundials
using Statistics

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

const MMC = TCellEngagerQSP.MosunModelCore
const REPO_ROOT = TCellEngagerQSP.REPO_ROOT

elapsed_s(f::Function) = (t0 = time_ns(); val = f(); ((time_ns() - t0) / 1e9, val))

function parse_float_list(txt::AbstractString)
    out = Float64[]
    for tok in split(String(txt), ",")
        s = strip(tok)
        isempty(s) && continue
        push!(out, parse(Float64, s))
    end
    return out
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
        dose_days = dose_days,
        dose_mg = dose_mg,
        dose_map = dose_map,
        horizon_days = Float64(row.horizon_days),
    )
end

function build_parameter_cases(horizon_days::Float64)
    base = build_base_param_map()
    phase1_patients = DataFrame(CSV.File(joinpath(REPO_ROOT, "generated", "phase1_design", "patients.csv")))
    cases = NamedTuple[]
    push!(cases, (
        case_key = "dlbcl_baseline",
        source = "generated/phase1_design_standard/dlbcl_param_overrides.csv",
        pmap = merged_param_map(base, Dict{String, Float64}(), horizon_days),
    ))
    prow = phase1_patients[1, :]
    push!(cases, (
        case_key = "phase1_patient1",
        source = "generated/phase1_design/patients.csv:1",
        pmap = merged_param_map(base, row_to_param_overrides(prow, names(phase1_patients)), horizon_days),
    ))
    return cases
end

function build_saveat(dose_times::Vector{Float64}, horizon_days::Float64)
    ts = Float64[]
    append!(ts, 0.0:7.0:horizon_days)
    push!(ts, horizon_days)
    for td in dose_times
        for dt in 0.0:0.25:3.0
            t = td + dt
            if t <= horizon_days + 1e-8
                push!(ts, round(t; digits = 8))
            end
        end
    end
    sort!(unique!(ts))
    return ts
end

function smoothmax(vals::AbstractVector, tau)
    m = maximum(vals)
    return m + tau * log(sum(exp.((vals .- m) ./ tau)))
end

function trapz(vals::AbstractVector, times::AbstractVector{<:Real})
    total = zero(eltype(vals))
    for i in 1:(length(vals) - 1)
        dt = times[i + 1] - times[i]
        total += dt * (vals[i] + vals[i + 1]) / 2
    end
    return total
end

function build_vector_callback(target_idx::Int, dose_map::Dict{Float64,Float64}; post_event_proposed_dt = nothing)
    local_dose_map = copy(dose_map)
    t0_dose = get(local_dose_map, 0.0, 0.0)
    if t0_dose != 0.0
        delete!(local_dose_map, 0.0)
    end
    cb_times = sort(collect(keys(local_dose_map)))
    callback = nothing
    if !isempty(cb_times)
        function affect!(integrator)
            integrator.u[target_idx] += local_dose_map[Float64(integrator.t)]
            if !isnothing(post_event_proposed_dt)
                SciMLBase.set_proposed_dt!(integrator, post_event_proposed_dt)
            end
        end
        callback = PresetTimeCallback(cb_times, affect!; save_positions = (false, false))
    end
    return t0_dose, callback, cb_times
end

function solver_configurations()
    return [
        (label = "cvode_bdf_postdt_1e-2", alg = CVODE_BDF(), post_event_dt = 1e-2),
        (label = "qndf", alg = QNDF(autodiff = false), post_event_dt = nothing),
        (label = "qndf_postdt_1e-4", alg = QNDF(autodiff = false), post_event_dt = 1e-4),
        (label = "qndf_postdt_1e-3", alg = QNDF(autodiff = false), post_event_dt = 1e-3),
        (label = "kencarp4", alg = KenCarp4(autodiff = false), post_event_dt = nothing),
        (label = "kencarp4_postdt_1e-3", alg = KenCarp4(autodiff = false), post_event_dt = 1e-3),
        (label = "trbdf2", alg = TRBDF2(autodiff = false), post_event_dt = nothing),
        (label = "rodas4p", alg = Rodas4P(autodiff = false), post_event_dt = nothing),
        (label = "tsit5", alg = Tsit5(), post_event_dt = nothing),
    ]
end

function build_case_context(case, regimen_spec, saveat)
    params = MMC.params_from_dict(case.pmap; ignore_unknown = true)
    base_u0 = MMC.pack_state(MMC.initial_state(params))
    p_base = MMC.pack_params(params)
    p_idx = Dict(sym => i for (i, sym) in enumerate(MMC.PARAMETER_NAMES))
    u_idx = Dict(sym => i for (i, sym) in enumerate(MMC.DYNAMIC_STATE_NAMES))
    fit_syms = [:fTadeact, :kBtumorprolif, :kTaapop, :fTaprolif, :fBkill, :fBprolif, :fBexit, :fTact]
    fit_idxs = [p_idx[s] for s in fit_syms]
    fit_base = p_base[fit_idxs]
    return (
        case_key = case.case_key,
        source = case.source,
        p_base = p_base,
        fit_syms = fit_syms,
        fit_idxs = fit_idxs,
        fit_base = fit_base,
        base_u0 = base_u0,
        saveat = saveat,
        dose_map = regimen_spec.dose_map,
        horizon_days = regimen_spec.horizon_days,
        u_idx = u_idx,
        p_idx = p_idx,
    )
end

function make_dual_safe_objective(ctx, cfg; tau = 50.0)
    target_idx = ctx.u_idx[:TDBc_ugperkg]
    il6pb_idx = ctx.u_idx[:IL6pb]
    il6tiss_idx = ctx.u_idx[:IL6tiss]
    il6tiss2_idx = ctx.u_idx[:IL6tiss2]
    il6tiss3_idx = ctx.u_idx[:IL6tiss3]
    il6tumor_idx = ctx.u_idx[:IL6tumor]
    btumor_idx = ctx.u_idx[:Btumor]
    auc_idx = ctx.u_idx[:TDBc_ugperml_AUC]

    idx_vpb = ctx.p_idx[:Vpb]
    idx_vtissue = ctx.p_idx[:Vtissue]
    idx_vtissue2 = ctx.p_idx[:Vtissue2]
    idx_vtissue3 = ctx.p_idx[:Vtissue3]
    idx_vtumor = ctx.p_idx[:Vtumor]
    idx_il6_tiss_contrib = ctx.p_idx[:IL6_tiss_contribution]
    idx_vc_tdb = ctx.p_idx[:Vc_tdb]
    dose_times = sort(collect(keys(ctx.dose_map)))
    t0_dose = get(ctx.dose_map, 0.0, 0.0)
    zero_x = zeros(Float64, length(ctx.fit_idxs))

    function solve_with_x(x::AbstractVector)
        T = promote_type(eltype(ctx.p_base), eltype(x))
        p = Vector{T}(undef, length(ctx.p_base))
        copyto!(p, ctx.p_base)
        @inbounds for i in eachindex(ctx.fit_idxs)
            p[ctx.fit_idxs[i]] = ctx.fit_base[i] * exp(x[i])
        end
        u0 = Vector{T}(undef, length(ctx.base_u0))
        copyto!(u0, ctx.base_u0)
        t0_dose != 0.0 && (u0[target_idx] += t0_dose)
        t_hist = Float64[]
        u_hist = Vector{typeof(u0)}()
        if 0.0 in ctx.saveat
            push!(t_hist, 0.0)
            push!(u_hist, copy(u0))
        end

        t_curr = 0.0
        u_curr = copy(u0)
        for td in dose_times
            td <= t_curr + 1e-8 && continue
            td > ctx.horizon_days + 1e-8 && continue
            seg_requested = [t for t in ctx.saveat if t_curr < t <= td + 1e-8]
            seg_save = sort(unique(vcat(seg_requested, [td])))
            prob = ODEProblem(MMC.mosun_rhs_vector!, u_curr, (t_curr, td), p)
            sol = solve(
                prob,
                cfg.alg;
                abstol = TCellEngagerQSP.SOLVER_ABSTOL,
                reltol = TCellEngagerQSP.SOLVER_RELTOL,
                saveat = seg_save,
                save_start = false,
                save_everystep = false,
                tstops = [td],
            )
            sol.retcode == SciMLBase.ReturnCode.Success || error("solve failed with retcode=$(sol.retcode)")
            for i in eachindex(sol.t)
                if sol.t[i] in seg_requested
                    push!(t_hist, Float64(sol.t[i]))
                    push!(u_hist, copy(sol.u[i]))
                end
            end
            u_curr = copy(sol.u[end])
            u_curr[target_idx] += ctx.dose_map[td]
            t_curr = td
        end

        if t_curr < ctx.horizon_days - 1e-8
            seg_requested = [t for t in ctx.saveat if t_curr < t <= ctx.horizon_days + 1e-8]
            seg_save = sort(unique(vcat(seg_requested, [ctx.horizon_days])))
            prob = ODEProblem(MMC.mosun_rhs_vector!, u_curr, (t_curr, ctx.horizon_days), p)
            sol = solve(
                prob,
                cfg.alg;
                abstol = TCellEngagerQSP.SOLVER_ABSTOL,
                reltol = TCellEngagerQSP.SOLVER_RELTOL,
                saveat = seg_save,
                save_start = false,
                save_everystep = false,
                tstops = [ctx.horizon_days],
            )
            sol.retcode == SciMLBase.ReturnCode.Success || error("solve failed with retcode=$(sol.retcode)")
            for i in eachindex(sol.t)
                if sol.t[i] in seg_requested
                    push!(t_hist, Float64(sol.t[i]))
                    push!(u_hist, copy(sol.u[i]))
                end
            end
        end

        return p, (t = t_hist, u = u_hist)
    end

    p0, sol0 = solve_with_x(zero_x)
    bt0 = sol0.u[1][btumor_idx]

    function il6combo(u, p)
        return u[il6pb_idx] +
            p[idx_il6_tiss_contrib] * (
                u[il6tiss_idx] * p[idx_vtissue] +
                u[il6tiss2_idx] * p[idx_vtissue2] +
                u[il6tiss3_idx] * p[idx_vtissue3] +
                u[il6tumor_idx] * p[idx_vtumor]
            ) / p[idx_vpb]
    end

    function objective_components(x::AbstractVector)
        p, sol = solve_with_x(x)
        tumor_ratio = similar(sol.t, eltype(sol.u[1]))
        il6_vals = similar(sol.t, eltype(sol.u[1]))
        for i in eachindex(sol.t)
            u = sol.u[i]
            tumor_ratio[i] = u[btumor_idx] / bt0
            il6_vals[i] = il6combo(u, p)
        end
        peak_il6 = smoothmax(il6_vals, tau)
        tumor_auc_ratio = trapz(tumor_ratio, sol.t) / (sol.t[end] - sol.t[1] + eps(Float64))
        final_drug_auc = sol.u[end][auc_idx]
        final_tdbc = begin
            val = sol.u[end][target_idx] / p[idx_vc_tdb]
            ifelse(val > 1e-5, val, zero(val))
        end
        return peak_il6, tumor_auc_ratio, final_drug_auc, final_tdbc
    end

    c0 = objective_components(zero_x)
    peak_scale = max(abs(Float64(c0[1])), 1e-12)
    tumor_scale = max(abs(Float64(c0[2])), 1e-12)
    auc_scale = max(abs(Float64(c0[3])), 1e-12)
    drug_scale = max(abs(Float64(c0[4])), 1e-12)

    function objective(x::AbstractVector)
        peak_il6, tumor_auc_ratio, final_drug_auc, final_tdbc = objective_components(x)
        return 0.45 * (peak_il6 / peak_scale) +
            0.30 * (tumor_auc_ratio / tumor_scale) +
            0.15 * (final_drug_auc / auc_scale) +
            0.10 * (final_tdbc / drug_scale)
    end

    return objective, c0
end

function finite_difference_gradient(f, x::Vector{Float64}; rel_step::Float64 = 1e-4, abs_step::Float64 = 1e-6)
    g = zeros(Float64, length(x))
    for i in eachindex(x)
        h = max(abs_step, rel_step * max(abs(x[i]), 1.0))
        xp = copy(x)
        xm = copy(x)
        xp[i] += h
        xm[i] -= h
        g[i] = (f(xp) - f(xm)) / (2h)
    end
    return g
end

function make_chunk(chunk_size::Int)
    if chunk_size == 1
        return ForwardDiff.Chunk{1}()
    elseif chunk_size == 2
        return ForwardDiff.Chunk{2}()
    elseif chunk_size == 4
        return ForwardDiff.Chunk{4}()
    elseif chunk_size == 8
        return ForwardDiff.Chunk{8}()
    end
    throw(ArgumentError("Unsupported EXPLICIT_CORE_FWD_GRAD_CHUNK=$chunk_size. Use 1, 2, 4, or 8."))
end

function gradient_metrics(g::Vector{Float64}, g_fd::Vector{Float64})
    abs_err = abs.(g .- g_fd)
    rel_err = abs_err ./ (abs.(g_fd) .+ 1e-12)
    cos_sim = dot(g, g_fd) / (norm(g) * norm(g_fd) + 1e-12)
    return (
        max_abs_err = maximum(abs_err),
        median_abs_err = median(abs_err),
        max_rel_err = maximum(rel_err),
        median_rel_err = median(rel_err),
        cosine_similarity = cos_sim,
    )
end

function main()
    out_dir_env = get(ENV, "EXPLICIT_CORE_FWD_GRAD_OUT_DIR", joinpath("generated", "benchmarks", "explicit_core_forward_gradient_sweep_20260310"))
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(REPO_ROOT, out_dir_env)
    mkpath(out_dir)

    repeats = parse(Int, get(ENV, "EXPLICIT_CORE_FWD_GRAD_REPEATS", "2"))
    regimen_spec = realdata_regimen_spec()
    horizon_override = parse(Float64, get(ENV, "EXPLICIT_CORE_FWD_GRAD_HORIZON_DAYS", "42.0"))
    bench_horizon_days = min(regimen_spec.horizon_days, horizon_override)
    bench_regimen_spec = (
        scenario = regimen_spec.scenario,
        regimen_label = regimen_spec.regimen_label,
        dose_days = [t for t in regimen_spec.dose_days if t <= bench_horizon_days + 1e-8],
        dose_mg = regimen_spec.dose_mg[1:length([t for t in regimen_spec.dose_days if t <= bench_horizon_days + 1e-8])],
        dose_map = Dict(t => amt for (t, amt) in regimen_spec.dose_map if t <= bench_horizon_days + 1e-8),
        horizon_days = bench_horizon_days,
    )
    saveat = build_saveat(bench_regimen_spec.dose_days, bench_regimen_spec.horizon_days)
    cases = build_parameter_cases(bench_regimen_spec.horizon_days)
    x0 = zeros(Float64, 8)
    chunk_size = parse(Int, get(ENV, "EXPLICIT_CORE_FWD_GRAD_CHUNK", "1"))
    chunk = make_chunk(chunk_size)

    runs = DataFrame(
        config = String[],
        case_key = String[],
        parameter_source = String[],
        rep = Int[],
        is_warmup = Bool[],
        objective_time_s = Float64[],
        gradient_time_s = Float64[],
        status = String[],
        objective_value = Float64[],
        grad_norm = Float64[],
        max_abs_err_vs_fd = Float64[],
        median_abs_err_vs_fd = Float64[],
        max_rel_err_vs_fd = Float64[],
        median_rel_err_vs_fd = Float64[],
        cosine_similarity_vs_fd = Float64[],
        note = String[],
    )

    summary = DataFrame(
        config = String[],
        case_key = String[],
        parameter_source = String[],
        n_success = Int[],
        warmup_gradient_time_s = Float64[],
        mean_gradient_time_s = Float64[],
        mean_objective_time_s = Float64[],
        gradient_over_objective = Float64[],
        mean_objective_value = Float64[],
        mean_grad_norm = Float64[],
        mean_max_abs_err_vs_fd = Float64[],
        mean_max_rel_err_vs_fd = Float64[],
        mean_cosine_similarity_vs_fd = Float64[],
    )

    meta = DataFrame(
        key = ["scenario", "regimen_label", "dose_days", "dose_mg", "horizon_days", "fit_parameters", "saveat", "repeats", "chunk_size"],
        value = [
            bench_regimen_spec.scenario,
            bench_regimen_spec.regimen_label,
            join(bench_regimen_spec.dose_days, ","),
            join(bench_regimen_spec.dose_mg, ","),
            string(bench_regimen_spec.horizon_days),
            join(string.([:fTadeact, :kBtumorprolif, :kTaapop, :fTaprolif, :fBkill, :fBprolif, :fBexit, :fTact]), ","),
            join(saveat, ","),
            string(repeats),
            string(chunk_size),
        ],
    )
    CSV.write(joinpath(out_dir, "meta.csv"), meta)

    configs = solver_configurations()
    for case in cases
        ctx = build_case_context(case, bench_regimen_spec, saveat)
        for cfg in configs
            note = ""
            try
                obj, base_components = make_dual_safe_objective(ctx, cfg)
                fd_time, g_fd = elapsed_s() do
                    finite_difference_gradient(obj, x0)
                end
                grad_cfg = ForwardDiff.GradientConfig(obj, x0, chunk)

                warm_obj_time, warm_obj = elapsed_s(() -> obj(x0))
                warm_grad_time, warm_g = elapsed_s(() -> ForwardDiff.gradient(obj, x0, grad_cfg))
                warm_metrics = gradient_metrics(Float64.(warm_g), g_fd)
                push!(runs, (
                    cfg.label,
                    ctx.case_key,
                    ctx.source,
                    0,
                    true,
                    warm_obj_time,
                    warm_grad_time,
                    "ok",
                    Float64(warm_obj),
                    norm(Float64.(warm_g)),
                    warm_metrics.max_abs_err,
                    warm_metrics.median_abs_err,
                    warm_metrics.max_rel_err,
                    warm_metrics.median_rel_err,
                    warm_metrics.cosine_similarity,
                    "fd_time_s=$(fd_time); base_components=$(base_components)",
                ))

                for rep in 1:repeats
                    obj_time, obj_val = elapsed_s(() -> obj(x0))
                    grad_time, g = elapsed_s(() -> ForwardDiff.gradient(obj, x0, grad_cfg))
                    gm = gradient_metrics(Float64.(g), g_fd)
                    push!(runs, (
                        cfg.label,
                        ctx.case_key,
                        ctx.source,
                        rep,
                        false,
                        obj_time,
                        grad_time,
                        "ok",
                        Float64(obj_val),
                        norm(Float64.(g)),
                        gm.max_abs_err,
                        gm.median_abs_err,
                        gm.max_rel_err,
                        gm.median_rel_err,
                        gm.cosine_similarity,
                        note,
                    ))
                end
            catch err
                push!(runs, (
                    cfg.label,
                    ctx.case_key,
                    ctx.source,
                    0,
                    true,
                    NaN,
                    NaN,
                    "error",
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    sprint(showerror, err),
                ))
            end
        end
    end

    for sub in groupby(runs, [:config, :case_key, :parameter_source])
        warm = sub[sub.is_warmup .== true, :]
        timed = sub[(sub.is_warmup .== false) .& (sub.status .== "ok"), :]
        push!(summary, (
            String(sub.config[1]),
            String(sub.case_key[1]),
            String(sub.parameter_source[1]),
            nrow(timed),
            isempty(warm) ? NaN : Float64(warm.gradient_time_s[1]),
            isempty(timed) ? NaN : mean(timed.gradient_time_s),
            isempty(timed) ? NaN : mean(timed.objective_time_s),
            isempty(timed) ? NaN : mean(timed.gradient_time_s) / max(mean(timed.objective_time_s), 1e-12),
            isempty(timed) ? NaN : mean(timed.objective_value),
            isempty(timed) ? NaN : mean(timed.grad_norm),
            isempty(timed) ? NaN : mean(timed.max_abs_err_vs_fd),
            isempty(timed) ? NaN : mean(timed.max_rel_err_vs_fd),
            isempty(timed) ? NaN : mean(timed.cosine_similarity_vs_fd),
        ))
    end

    aggregate = combine(
        groupby(summary, :config),
        :n_success => sum => :n_success_total,
        :mean_gradient_time_s => mean => :mean_gradient_time_s,
        :mean_objective_time_s => mean => :mean_objective_time_s,
        :gradient_over_objective => mean => :mean_gradient_over_objective,
        :warmup_gradient_time_s => mean => :mean_warmup_gradient_time_s,
        :mean_max_abs_err_vs_fd => mean => :mean_max_abs_err_vs_fd,
        :mean_max_rel_err_vs_fd => mean => :mean_max_rel_err_vs_fd,
        :mean_cosine_similarity_vs_fd => mean => :mean_cosine_similarity_vs_fd,
    )
    sort!(aggregate, :mean_gradient_time_s)

    CSV.write(joinpath(out_dir, "runs.csv"), runs)
    CSV.write(joinpath(out_dir, "summary_by_case.csv"), summary)
    CSV.write(joinpath(out_dir, "summary_aggregate.csv"), aggregate)
    println("Wrote forward-gradient benchmark results to $(out_dir)")
end

main()
