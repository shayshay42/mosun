using CSV
using DataFrames
using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using SciMLBase
using SciMLSensitivity
using Statistics

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

const MMC = TCellEngagerQSP.MosunModelCore
const REPO_ROOT = TCellEngagerQSP.REPO_ROOT

elapsed_s(f::Function) = (t0 = time_ns(); val = f(); ((time_ns() - t0) / 1e9, val))

function build_case()
    patients = DataFrame(CSV.File(joinpath(REPO_ROOT, "generated", "phase1_design", "patients.csv")))
    overrides = DataFrame(CSV.File(joinpath(REPO_ROOT, "generated", "phase1_design", "dlbcl_param_overrides.csv")))
    prow = patients[findfirst(patients.patient_id .== 1), :]

    pmap = Dict{String, Float64}()
    for rr in eachrow(overrides)
        pmap[String(rr.name)] = Float64(rr.value)
    end
    for nm in ("Bpbo_perml", "Bpbref_perml", "Trpbo_perml", "Trpbref_perml", "KBptumor", "KTrptumor", "kBtumorprolif")
        pmap[nm] = Float64(prow[Symbol(nm)])
    end
    pmap["PKflag"] = 1.0
    pmap["fvalidation"] = 0.0
    pmap["VPid"] = 1.0
    pmap["end_time"] = 7.0

    params = MMC.params_from_dict(pmap; ignore_unknown = true)
    p_base = MMC.pack_params(params)
    p_idx = Dict(sym => i for (i, sym) in enumerate(MMC.PARAMETER_NAMES))
    u_idx = Dict(sym => i for (i, sym) in enumerate(MMC.DYNAMIC_STATE_NAMES))
    fit_syms = [:fTadeact, :kBtumorprolif, :kTaapop, :fTaprolif, :fBkill, :fBprolif, :fBexit, :fTact]
    fit_idxs = [p_idx[s] for s in fit_syms]
    fit_base = p_base[fit_idxs]
    regimen = MMC.MosunRegimen(events = [
        MMC.MosunRegimenEvent(target = :TDBc_ugperkg, time = 0.0, amount = 0.8 * 1000.0 / 70.0, rate = 0.0),
    ])
    saveat = collect(0.0:1.0:7.0)
    return (; p_base, p_idx, u_idx, fit_syms, fit_idxs, fit_base, regimen, saveat)
end

function il6combo(u::AbstractVector, p::AbstractVector, u_idx, p_idx)
    return u[u_idx[:IL6pb]] +
        p[p_idx[:IL6_tiss_contribution]] * (
            u[u_idx[:IL6tiss]] * p[p_idx[:Vtissue]] +
            u[u_idx[:IL6tiss2]] * p[p_idx[:Vtissue2]] +
            u[u_idx[:IL6tiss3]] * p[p_idx[:Vtissue3]] +
            u[u_idx[:IL6tumor]] * p[p_idx[:Vtumor]]
        ) / p[p_idx[:Vpb]]
end

function solver_configs()
    return [
        (label = "tsit5", alg = Tsit5()),
        (label = "qndf", alg = QNDF(autodiff = false)),
        (label = "rodas4p", alg = Rodas4P(autodiff = false)),
    ]
end

function build_solution(ctx, x::AbstractVector, alg; sensealg = nothing)
    T = promote_type(eltype(ctx.p_base), eltype(x))
    p = T.(ctx.p_base)
    for i in eachindex(ctx.fit_idxs)
        p[ctx.fit_idxs[i]] = ctx.fit_base[i] * exp(x[i])
    end
    built = MMC.build_problem_vector(ctx.regimen, p; tspan = (0.0, 7.0), saveat = ctx.saveat, callback_mode = :callback)
    solve_kwargs = (
        abstol = TCellEngagerQSP.SOLVER_ABSTOL,
        reltol = TCellEngagerQSP.SOLVER_RELTOL,
        callback = built.callback,
        tstops = built.tstops,
        d_discontinuities = built.d_discontinuities,
        saveat = built.saveat,
        save_everystep = false,
    )
    sol = isnothing(sensealg) ?
        solve(built.prob, alg; solve_kwargs...) :
        solve(built.prob, alg; solve_kwargs..., sensealg = sensealg)
    sol.retcode == SciMLBase.ReturnCode.Success || error("solve failed with retcode=$(sol.retcode)")
    return p, built, sol
end

function objective_from_solution(ctx, p::AbstractVector, built, sol, scales; tau = 50.0)
    il6_vals = [il6combo(u, p, ctx.u_idx, ctx.p_idx) for u in sol.u]
    m = maximum(il6_vals)
    w = exp.((il6_vals .- m) ./ tau)
    softmax = w ./ sum(w)
    peak_proxy = m + tau * log(sum(w))
    bt0 = built.initial_u[ctx.u_idx[:Btumor]]
    final_tumor_ratio = sol.u[end][ctx.u_idx[:Btumor]] / (bt0 + 1e-12)
    final_auc = sol.u[end][ctx.u_idx[:TDBc_ugperml_AUC]]
    loss = 0.45 * (peak_proxy / scales.peak) +
        0.35 * (final_tumor_ratio / scales.tumor) +
        0.20 * (final_auc / scales.auc)
    return loss, (; softmax, bt0, peak_proxy, final_tumor_ratio, final_auc)
end

function make_objective(ctx, alg, scales)
    return function (x::AbstractVector; sensealg = nothing)
        p, built, sol = build_solution(ctx, x, alg; sensealg = sensealg)
        loss, _ = objective_from_solution(ctx, p, built, sol, scales)
        return loss
    end
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

function gradient_metrics(g::AbstractVector, g_fd::Vector{Float64})
    g_vec = Float64.(g)
    abs_err = abs.(g_vec .- g_fd)
    rel_err = abs_err ./ (abs.(g_fd) .+ 1e-12)
    cos_sim = dot(g_vec, g_fd) / (norm(g_vec) * norm(g_fd) + 1e-12)
    return (
        grad = g_vec,
        max_abs_err = maximum(abs_err),
        max_rel_err = maximum(rel_err),
        median_rel_err = median(rel_err),
        cosine_similarity = cos_sim,
    )
end

function reverse_adjoint_gradient(ctx, alg, scales, x::Vector{Float64})
    p, built, sol = build_solution(ctx, x, alg)
    loss, aux = objective_from_solution(ctx, p, built, sol, scales)
    n_times = length(sol.t)
    peak_weight = 0.45 / scales.peak
    tumor_weight = 0.35 / scales.tumor
    auc_weight = 0.20 / scales.auc
    il6pb_idx = ctx.u_idx[:IL6pb]
    il6tiss_idx = ctx.u_idx[:IL6tiss]
    il6tiss2_idx = ctx.u_idx[:IL6tiss2]
    il6tiss3_idx = ctx.u_idx[:IL6tiss3]
    il6tumor_idx = ctx.u_idx[:IL6tumor]
    bt_idx = ctx.u_idx[:Btumor]
    auc_idx = ctx.u_idx[:TDBc_ugperml_AUC]
    il6_tiss_coeff = p[ctx.p_idx[:IL6_tiss_contribution]]
    vpb = p[ctx.p_idx[:Vpb]]
    coeff_tiss = il6_tiss_coeff * p[ctx.p_idx[:Vtissue]] / vpb
    coeff_tiss2 = il6_tiss_coeff * p[ctx.p_idx[:Vtissue2]] / vpb
    coeff_tiss3 = il6_tiss_coeff * p[ctx.p_idx[:Vtissue3]] / vpb
    coeff_tumor = il6_tiss_coeff * p[ctx.p_idx[:Vtumor]] / vpb

    function dg(out, u, p_local, t, i)
        out .= 0.0
        w = peak_weight * aux.softmax[i]
        out[il6pb_idx] += w
        out[il6tiss_idx] += w * coeff_tiss
        out[il6tiss2_idx] += w * coeff_tiss2
        out[il6tiss3_idx] += w * coeff_tiss3
        out[il6tumor_idx] += w * coeff_tumor
        if i == n_times
            out[bt_idx] += tumor_weight / (aux.bt0 + 1e-12)
            out[auc_idx] += auc_weight
        end
        return nothing
    end

    _, dp = adjoint_sensitivities(
        sol,
        alg;
        t = sol.t,
        dgdu_discrete = dg,
        sensealg = QuadratureAdjoint(autojacvec = EnzymeVJP()),
        callback = built.callback,
        abstol = TCellEngagerQSP.SOLVER_ABSTOL,
        reltol = TCellEngagerQSP.SOLVER_RELTOL,
    )
    gx = similar(x)
    for i in eachindex(ctx.fit_idxs)
        gx[i] = dp[ctx.fit_idxs[i]] * (ctx.fit_base[i] * exp(x[i]))
    end
    return loss, gx
end

function main()
    out_dir_env = get(ENV, "MOSUN_CORE_AD_BENCH_OUT_DIR", joinpath("generated", "benchmarks", "mosun_core_gradient_modes_20260310"))
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(REPO_ROOT, out_dir_env)
    mkpath(out_dir)

    ctx = build_case()
    x0 = zeros(Float64, length(ctx.fit_idxs))
    repeats = parse(Int, get(ENV, "MOSUN_CORE_AD_BENCH_REPEATS", "1"))

    base_p, base_built, base_sol = build_solution(ctx, x0, Tsit5())
    base_loss, base_aux = objective_from_solution(ctx, base_p, base_built, base_sol, (peak = 1.0, tumor = 1.0, auc = 1.0))
    scales = (
        peak = abs(base_aux.peak_proxy) > 0 ? abs(base_aux.peak_proxy) : 1.0,
        tumor = abs(base_aux.final_tumor_ratio) > 0 ? abs(base_aux.final_tumor_ratio) : 1.0,
        auc = abs(base_aux.final_auc) > 0 ? abs(base_aux.final_auc) : 1.0,
    )

    meta = DataFrame(
        key = ["fit_parameters", "saveat_days", "reference_loss"],
        value = [join(string.(ctx.fit_syms), ","), join(string.(ctx.saveat), ","), string(base_loss)],
    )
    CSV.write(joinpath(out_dir, "meta.csv"), meta)

    rows = DataFrame(
        solver = String[],
        method = String[],
        status = String[],
        objective_value = Float64[],
        objective_time_s = Float64[],
        gradient_time_s = Float64[],
        grad_norm = Float64[],
        max_abs_err_vs_fd = Float64[],
        max_rel_err_vs_fd = Float64[],
        median_rel_err_vs_fd = Float64[],
        cosine_similarity_vs_fd = Float64[],
        note = String[],
    )

    for cfg in solver_configs()
        obj = make_objective(ctx, cfg.alg, scales)
        fd_obj = x -> obj(x; sensealg = nothing)
        fd_time, g_fd = elapsed_s(() -> finite_difference_gradient(fd_obj, x0))
        for rep in 1:repeats
            obj_time, obj_val = elapsed_s(() -> fd_obj(x0))
            push!(rows, (
                cfg.label,
                "objective",
                "ok",
                Float64(obj_val),
                obj_time,
                NaN,
                NaN,
                0.0,
                0.0,
                0.0,
                1.0,
                "fd_time_s=$(fd_time); rep=$(rep)",
            ))

            for (method, gradfun) in (
                ("dto_forward", () -> ForwardDiff.gradient(x -> obj(x; sensealg = nothing), x0)),
                ("otd_forward", () -> ForwardDiff.gradient(x -> obj(x; sensealg = ForwardSensitivity()), x0)),
                ("adjoint_reverse_enzyme", () -> reverse_adjoint_gradient(ctx, cfg.alg, scales, x0)[2]),
            )
                try
                    grad_time, g = elapsed_s(gradfun)
                    gm = gradient_metrics(g, g_fd)
                    push!(rows, (
                        cfg.label,
                        method,
                        "ok",
                        Float64(obj_val),
                        obj_time,
                        grad_time,
                        norm(gm.grad),
                        gm.max_abs_err,
                        gm.max_rel_err,
                        gm.median_rel_err,
                        gm.cosine_similarity,
                        "",
                    ))
                catch err
                    push!(rows, (
                        cfg.label,
                        method,
                        "error",
                        Float64(obj_val),
                        obj_time,
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
    end

    summary = combine(
        groupby(rows[rows.method .!= "objective", :], [:solver, :method]),
        :status => (x -> count(==("ok"), x)) => :n_success,
        :gradient_time_s => mean => :mean_gradient_time_s,
        :grad_norm => mean => :mean_grad_norm,
        :max_abs_err_vs_fd => mean => :mean_max_abs_err_vs_fd,
        :max_rel_err_vs_fd => mean => :mean_max_rel_err_vs_fd,
        :cosine_similarity_vs_fd => mean => :mean_cosine_similarity_vs_fd,
    )
    sort!(summary, [:method, :mean_gradient_time_s])

    CSV.write(joinpath(out_dir, "runs.csv"), rows)
    CSV.write(joinpath(out_dir, "summary.csv"), summary)
    println("Wrote gradient benchmark results to $(out_dir)")
end

main()
