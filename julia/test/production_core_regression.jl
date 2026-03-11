using CSV
using DataFrames
using DifferentialEquations
using DiffEqCallbacks
using JSON3
using SciMLBase
using Test
using ForwardDiff

const MMC = TCellEngagerQSP.MosunModelCore

function max_norm_abs_err(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real})
    @assert size(X) == size(Y)
    mask = isfinite.(X) .& isfinite.(Y)
    any(mask) || return Inf
    amp = max(maximum(abs.(Y[mask])), 1e-12)
    return maximum(abs.(X[mask] .- Y[mask])) / amp
end

function legacy_u_from_core(mdl_legacy, mdl_prod, u_core::Vector{Float64})
    u_legacy = copy(mdl_legacy.u0)
    for (nm, idx_core) in mdl_prod.core_state_to_idx
        u_legacy[mdl_legacy.state_to_idx[nm]] = u_core[idx_core]
    end
    return u_legacy
end

function build_phase1_patient_case(; horizon_days::Float64 = 42.0)
    design_dir = joinpath(TCellEngagerQSP.REPO_ROOT, "generated", "phase1_design")
    patients = DataFrame(CSV.File(joinpath(design_dir, "patients.csv")))
    reg_events = DataFrame(CSV.File(joinpath(design_dir, "regimen_events.csv")))
    dlbcl_path = joinpath(design_dir, "dlbcl_param_overrides.csv")
    dlbcl_overrides = isfile(dlbcl_path) ? DataFrame(CSV.File(dlbcl_path)) : DataFrame(name = String[], value = Float64[])

    regimen_names = unique(reg_events.regimen)
    reg_name = "fixed_0.05mg" in regimen_names ? "fixed_0.05mg" : String(regimen_names[1])
    reg_sub = reg_events[reg_events.regimen .== reg_name, :]
    sort!(reg_sub, :event_idx)

    prow = patients[1, :]
    pmap = Dict{String, Float64}()
    for rr in eachrow(dlbcl_overrides)
        pmap[String(rr.name)] = Float64(rr.value)
    end
    for cname in names(patients)
        cname == "patient_id" && continue
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

    bw_kg_for_dose_conversion = 70.0
    dose_map = Dict{Float64, Float64}()
    for rr in eachrow(reg_sub)
        t = Float64(rr.time_day)
        amt = Float64(rr.dose_mg) * 1000.0 / bw_kg_for_dose_conversion
        dose_map[t] = get(dose_map, t, 0.0) + amt
    end
    return pnames, pvals, dose_map, reg_name
end

function solve_canonical_with_dose_map(mdl, dose_map::Dict{Float64,Float64}, horizon_days::Float64)
    state_to_idx = TCellEngagerQSP.canonical_state_to_idx(mdl)
    target_idx = state_to_idx["TDBc_ugperkg"]
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
    return solve(
        prob,
        TCellEngagerQSP.make_solver_alg("cvode_bdf");
        abstol = TCellEngagerQSP.SOLVER_ABSTOL,
        reltol = TCellEngagerQSP.SOLVER_RELTOL,
        callback = cb,
        tstops = cb_times,
        d_discontinuities = cb_times,
    )
end

function solve_production_with_dose_map(pnames::Vector{String}, pvals::Vector{Float64}, dose_map::Dict{Float64,Float64}, horizon_days::Float64)
    params = MMC.params_from_named_values(pnames, pvals; ignore_unknown = true)
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
    return params, sol
end

function sampled_series(mdl, sol, tgrid::Vector{Float64}, name::String)
    vals = Float64[]
    for t in tgrid
        u = Float64.(sol(t))
        push!(vals, TCellEngagerQSP.eval_symbol(mdl, u, t, name))
    end
    return vals
end

function sampled_series(params::MMC.MosunParams, sol, tgrid::Vector{Float64}, name::Symbol)
    vals = Float64[]
    cache = MMC.zero_observables_cache()
    for t in tgrid
        u = Float64.(sol(t))
        push!(vals, MMC.value_at(u, params, t, name, cache))
    end
    return vals
end

function best_spd_from_bt(bt::Vector{Float64})
    bt0 = bt[1]
    spd_pct = @. 100.0 * (bt / max(bt0, 1e-12) - 1.0)
    spd_eval = length(spd_pct) > 1 ? spd_pct[2:end] : spd_pct
    return minimum(spd_eval)
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

@testset "Production Core Structure" begin
    src = read(joinpath(TCellEngagerQSP.REPO_ROOT, "julia", "src", "MosunModelCore.jl"), String)
    for forbidden in ("using CSV", "using DataFrames", "using JSON3", "using ModelingToolkit", "using RuntimeGeneratedFunctions")
        @test !occursin(forbidden, src)
    end
    @test MMC.DYNAMIC_STATE_COUNT == 36
    @test length(MMC.OBSERVABLE_NAMES) == 89
    @test MMC.RHS_OBSERVABLE_COUNT == 58
    @test length(MMC.DEAD_LEGACY_NAMES) == 6
    @test :Tafraction_pb_init in MMC.DEAD_LEGACY_NAMES
end

@testset "Production Core Initial State And Observables" begin
    pnames, pvals, _, _ = build_phase1_patient_case()
    mdl_prod = TCellEngagerQSP.build_model_with_variant_ids(Int[], pnames, pvals; core_mode_override = :production)
    mdl_legacy = TCellEngagerQSP.build_model_with_variant_ids(Int[], pnames, pvals; core_mode_override = :legacy_reference)

    u_prod = MMC.pack_state(MMC.initial_state(mdl_prod.core_params))
    for (nm, idx_core) in mdl_prod.core_state_to_idx
        idx_legacy = mdl_legacy.state_to_idx[nm]
        @test isapprox(u_prod[idx_core], mdl_legacy.u0[idx_legacy]; atol = 1e-12, rtol = 0.0)
    end

    u_probe = copy(u_prod)
    u_probe[mdl_prod.core_state_to_idx["TDBc_ugperkg"]] += 0.7
    u_probe[mdl_prod.core_state_to_idx["actTpb"]] += 10.0
    u_probe[mdl_prod.core_state_to_idx["Btumor"]] *= 0.85
    u_probe[mdl_prod.core_state_to_idx["IL6pb"]] += 3.0
    t_probe = 3.25

    legacy_u = legacy_u_from_core(mdl_legacy, mdl_prod, u_probe)
    legacy_z = copy(mdl_legacy.z)
    TCellEngagerQSP.apply_repeated_rules_direct!(legacy_z, legacy_u, mdl_legacy.pvals, t_probe)
    cache = MMC.zero_observables_cache()
    MMC.update_observables!(cache, u_probe, mdl_prod.core_params, t_probe)

    check_names = (
        "TDBc_ugperml",
        "restTpb_perml",
        "Tafraction_tiss2",
        "B1920Trratio_tiss3",
        "Tafraction_tumor",
        "IL6combo",
    )
    for nm in check_names
        @test isapprox(MMC.observable(cache, Symbol(nm)), legacy_z[mdl_legacy.name_to_idx[nm]]; atol = 1e-10, rtol = 1e-10)
    end

    @test_throws ArgumentError MMC.state_or_observable(u_probe, cache, :Tafraction_pb_init)
end

@testset "Production Core Case30 Parity" begin
    manifest = DataFrame(CSV.File(joinpath(TCellEngagerQSP.REPO_ROOT, "generated", "matlab_reference", "manifest_three_cases.csv")))
    row = manifest[manifest.sim_key .== "case30_id1", :][1, :]
    _, outvec, X_prod, ref_df = TCellEngagerQSP.run_manifest_row(row; engine = :canonical, core_mode_override = :production)
    _, _, X_legacy, _ = TCellEngagerQSP.run_manifest_row(row; engine = :canonical, core_mode_override = :legacy_reference)
    X_ref = hcat([Float64.(ref_df[!, Symbol(name)]) for name in outvec]...)

    @test max_norm_abs_err(X_prod, X_legacy) <= 2e-4
    @test max_norm_abs_err(X_prod, X_ref) <= 2e-4
end

@testset "Production Core Phase1 Workflow Regression" begin
    pnames, pvals, dose_map, _ = build_phase1_patient_case()
    mdl_legacy = TCellEngagerQSP.build_model_with_variant_ids(Int[], pnames, pvals; core_mode_override = :legacy_reference)

    horizon_days = 42.0
    params_prod, sol_prod = solve_production_with_dose_map(pnames, pvals, dose_map, horizon_days)
    sol_legacy = solve_canonical_with_dose_map(mdl_legacy, dose_map, horizon_days)
    @test sol_prod.retcode == SciMLBase.ReturnCode.Success
    @test sol_legacy.retcode == SciMLBase.ReturnCode.Success

    tgrid = collect(0.0:1.0:horizon_days)
    il6_prod = sampled_series(params_prod, sol_prod, tgrid, :IL6combo)
    il6_legacy = sampled_series(mdl_legacy, sol_legacy, tgrid, "IL6combo")
    bt_prod = sampled_series(params_prod, sol_prod, tgrid, :Btumor)
    bt_legacy = sampled_series(mdl_legacy, sol_legacy, tgrid, "Btumor")

    @test maximum(abs.(il6_prod .- il6_legacy)) / max(maximum(abs.(il6_legacy)), 1e-12) <= 2e-5
    @test maximum(abs.(bt_prod .- bt_legacy)) / max(maximum(abs.(bt_legacy)), 1e-12) <= 1e-5
    @test isapprox(best_spd_from_bt(bt_prod), best_spd_from_bt(bt_legacy); atol = 5e-5, rtol = 5e-5)
end

@testset "Production Core ForwardDiff Path" begin
    pnames, pvals, _, _ = build_phase1_patient_case(; horizon_days = 7.0)
    params = MMC.params_from_named_values(pnames, pvals; ignore_unknown = true)
    p_base = MMC.pack_params(params)
    p_idx = Dict(sym => i for (i, sym) in enumerate(MMC.PARAMETER_NAMES))
    u_idx = Dict(sym => i for (i, sym) in enumerate(MMC.DYNAMIC_STATE_NAMES))
    fit_syms = [:fTadeact, :kBtumorprolif]
    fit_idxs = [p_idx[s] for s in fit_syms]
    fit_base = p_base[fit_idxs]
    regimen = MMC.MosunRegimen(events = [
        MMC.MosunRegimenEvent(target = :TDBc_ugperkg, time = 0.0, amount = 0.8 * 1000.0 / 70.0, rate = 0.0),
    ])

    function objective(x::AbstractVector)
        T = promote_type(eltype(x), Float64)
        p = T.(p_base)
        for i in eachindex(fit_idxs)
            p[fit_idxs[i]] = fit_base[i] * exp(x[i])
        end
        built = MMC.build_problem_vector(regimen, p; tspan = (0.0, 7.0), saveat = [0.0, 7.0], callback_mode = :callback)
        sol = solve(
            built.prob,
            Tsit5();
            abstol = TCellEngagerQSP.SOLVER_ABSTOL,
            reltol = TCellEngagerQSP.SOLVER_RELTOL,
            callback = built.callback,
            tstops = built.tstops,
            d_discontinuities = built.d_discontinuities,
            saveat = built.saveat,
        )
        @test sol.retcode == SciMLBase.ReturnCode.Success
        bt0 = built.initial_u[u_idx[:Btumor]]
        return sol.u[end][u_idx[:Btumor]] / (bt0 + T(1e-12))
    end

    x0 = zeros(Float64, length(fit_idxs))
    g_fd = finite_difference_gradient(objective, x0)
    g_ad = ForwardDiff.gradient(objective, x0)
    @test all(isfinite, g_ad)
    @test maximum(abs.(g_ad .- g_fd)) <= 1e-5
end

@testset "Production Core Callback Dose AD Path" begin
    pnames, pvals, _, _ = build_phase1_patient_case(; horizon_days = 7.0)
    params = MMC.params_from_named_values(pnames, pvals; ignore_unknown = true)
    pvec = MMC.pack_params(params)
    u_idx = Dict(sym => i for (i, sym) in enumerate(MMC.DYNAMIC_STATE_NAMES))

    function objective(x::AbstractVector)
        regimen = MMC.MosunRegimen(events = [
            MMC.MosunRegimenEvent(target = :TDBc_ugperkg, time = 0.0, amount = x[1] * (1000.0 / 70.0), rate = zero(eltype(x))),
            MMC.MosunRegimenEvent(target = :TDBc_ugperkg, time = 7.0, amount = x[2] * (1000.0 / 70.0), rate = zero(eltype(x))),
        ])
        built = MMC.build_problem_vector(regimen, pvec; tspan = (0.0, 7.0), saveat = [0.0, 7.0], callback_mode = :callback)
        sol = solve(
            built.prob,
            Tsit5();
            abstol = TCellEngagerQSP.SOLVER_ABSTOL,
            reltol = TCellEngagerQSP.SOLVER_RELTOL,
            callback = built.callback,
            tstops = built.tstops,
            d_discontinuities = built.d_discontinuities,
            saveat = built.saveat,
        )
        @test sol.retcode == SciMLBase.ReturnCode.Success
        bt0 = built.initial_u[u_idx[:Btumor]]
        return sol.u[end][u_idx[:Btumor]] / (bt0 + eltype(sol.u[end])(1e-12))
    end

    x0 = [0.8, 2.0]
    g_fd = finite_difference_gradient(objective, x0)
    g_ad = ForwardDiff.gradient(objective, x0)
    @test all(isfinite, g_ad)
    @test maximum(abs.(g_ad .- g_fd)) <= 1e-4
end
