using CSV
using DataFrames
using DifferentialEquations
using DiffEqCallbacks
using FFTW
using JSON3
using SciMLBase
using Sundials
using Random

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

const REPO_ROOT = TCellEngagerQSP.REPO_ROOT

function load_json(path::String)
    return JSON3.read(read(path, String))
end

function read_param_bounds_csv(path::String)
    df = DataFrame(CSV.File(path))
    out = Dict{String, Tuple{Float64, Float64}}()
    for r in eachrow(df)
        name = String(r.name)
        lb = try
            Float64(r.cyno_lb)
        catch
            NaN
        end
        ub = try
            Float64(r.cyno_ub)
        catch
            NaN
        end
        if isfinite(lb) && isfinite(ub) && ub > lb
            out[name] = (lb, ub)
        end
    end
    return out
end

function get_hard_bounds(targets_json, spread_mode::String)
    k_bounds = nothing
    bt_bounds = nothing
    spread_primary = nothing
    spread_sensitivity = nothing
    for t in targets_json["hard_priors"]
        tid = String(t["id"])
        if tid == "tumor_bcell_proliferation_rate"
            b = t["bounds"]
            k_bounds = (Float64(b["lower"]), Float64(b["upper"]))
        elseif tid == "tumor_b_to_t_ratio"
            b = t["bounds"]
            bt_bounds = (Float64(b["lower"]), Float64(b["upper"]))
        elseif tid == "baseline_tumor_bcell_from_spd_primary"
            b = t["bounds_relative_factor"]
            spread_primary = (Float64(b["lower"]), Float64(b["upper"]))
        elseif tid == "baseline_tumor_bcell_from_spd_alternative"
            b = t["bounds_relative_factor"]
            spread_sensitivity = (Float64(b["lower"]), Float64(b["upper"]))
        end
    end
    if k_bounds === nothing || bt_bounds === nothing
        error("Missing hard prior bounds in generated/vpop_targets/vpop_calibration_targets.json")
    end
    if spread_primary === nothing
        spread_primary = (0.9, 1.1)
    end
    if spread_sensitivity === nothing
        spread_sensitivity = (0.5, 1.5)
    end
    spread_bounds = lowercase(spread_mode) == "sensitivity" ? spread_sensitivity : spread_primary
    return k_bounds, bt_bounds, spread_bounds
end

function load_base_param_map()
    p = joinpath(REPO_ROOT, "generated", "phase1_design_standard", "dlbcl_param_overrides.csv")
    if !isfile(p)
        error("Missing DLBCL overrides at $p")
    end
    df = DataFrame(CSV.File(p))
    pmap = Dict{String, Float64}()
    for r in eachrow(df)
        pmap[String(r.name)] = Float64(r.value)
    end
    # Consistent baseline bo/ref values for generated phase1 standard design.
    pmap["Bpbo_perml"] = get(pmap, "Bpbo_perml", get(pmap, "Bpbref_perml", 500_000.0))
    pmap["Trpbo_perml"] = get(pmap, "Trpbo_perml", get(pmap, "Trpbref_perml", 500_000.0))
    pmap["Bpbref_perml"] = get(pmap, "Bpbref_perml", pmap["Bpbo_perml"])
    pmap["Trpbref_perml"] = get(pmap, "Trpbref_perml", pmap["Trpbo_perml"])
    return pmap
end

function salib_fast_sample(
        bounds::Vector{Tuple{Float64, Float64}},
        names::Vector{String},
        N::Int,
        M::Int,
        rng::AbstractRNG)
    D = length(names)
    if N <= 4 * M^2
        error("Sample size N must satisfy N > 4*M^2 for FAST/eFAST (got N=$N, M=$M).")
    end
    omega = zeros(Int, D)
    omega[1] = floor(Int, (N - 1) / (2 * M))
    m = floor(Int, omega[1] / (2 * M))
    if m < 1
        error("Computed m=$m; increase N or reduce M.")
    end
    if m >= (D - 1)
        omega[2:end] = floor.(Int, range(1, m; length = D - 1))
    else
        omega[2:end] = [mod(i - 1, m) + 1 for i in 1:(D - 1)]
    end

    s = (2 * π / N) .* collect(0:(N - 1))
    X = zeros(Float64, N * D, D)
    omega2 = zeros(Int, D)

    for i in 1:D
        fill!(omega2, 0)
        omega2[i] = omega[1]
        rem_idx = [j for j in 1:D if j != i]
        omega2[rem_idx] = omega[2:end]
        z0 = (i - 1) * N
        ϕ = 2 * π * rand(rng)
        for j in 1:D
            g = @. 0.5 + (1 / π) * asin(sin(omega2[j] * s + ϕ))
            X[(z0 + 1):(z0 + N), j] = g
        end
    end

    for j in 1:D
        lb, ub = bounds[j]
        X[:, j] .= lb .+ X[:, j] .* (ub - lb)
    end
    return X, omega[1]
end

function salib_fast_analyze(Y::Vector{Float64}, D::Int, N::Int, M::Int)
    if length(Y) != N * D
        error("Expected length(Y)=N*D=$(N * D), got $(length(Y))")
    end
    omega0 = floor(Int, (N - 1) / (2 * M))
    S1 = fill(NaN, D)
    ST = fill(NaN, D)
    for i in 1:D
        y = Y[((i - 1) * N + 1):(i * N)]
        f = fft(y)
        sp = (abs.(f[2:ceil(Int, N / 2)]) ./ N) .^ 2
        V = 2.0 * sum(sp)
        if !(isfinite(V) && V > 0)
            continue
        end
        idx_d1 = (collect(1:M) .* omega0)
        idx_d1 = idx_d1[idx_d1 .<= length(sp)]
        d1 = 2.0 * sum(sp[idx_d1])
        idx_dt = 1:floor(Int, omega0 / 2)
        idx_dt = idx_dt[idx_dt .<= length(sp)]
        dt = 2.0 * sum(sp[idx_dt])
        S1[i] = d1 / V
        ST[i] = 1.0 - dt / V
    end
    return S1, ST
end

function build_regimen_8cycle_stepup()
    dose_days = Float64[0.0, 7.0, 14.0, 21.0, 42.0, 63.0, 84.0, 105.0, 126.0, 147.0]
    dose_mg = Float64[1.0, 2.0, 60.0, 60.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0]
    return dose_days, dose_mg
end

function eval_best_spd_for_sample(
        sample_vals::Dict{String, Float64},
        base_pmap::Dict{String, Float64},
        extra_names::Vector{String},
        variant_ids::Vector{Int},
        dose_days::Vector{Float64},
        dose_mg::Vector{Float64},
        horizon_days::Float64,
        save_dt::Float64,
        bw_kg::Float64,
        alg,
        use_mtk_jac::Bool)
    pmap = copy(base_pmap)

    base_bpbo = pmap["Bpbo_perml"]
    base_trpbo = pmap["Trpbo_perml"]
    # Base DLBCL tumor burden in cells/ml from KBptumor * Bpbo_perml.
    base_btumor_perml = pmap["KBptumor"] * base_bpbo

    k_btumor = sample_vals["kBtumorprolif"]
    bt_ratio = sample_vals["BT_ratio_tumor_init"]
    burden_factor = sample_vals["tumor_burden_factor"]
    btumor_perml = base_btumor_perml * burden_factor

    pmap["kBtumorprolif"] = k_btumor
    pmap["KBptumor"] = btumor_perml / base_bpbo
    pmap["KTrptumor"] = btumor_perml / (max(bt_ratio, 1e-12) * base_trpbo)
    for nm in extra_names
        pmap[nm] = sample_vals[nm]
    end

    pmap["PKflag"] = 1.0
    pmap["fvalidation"] = 0.0
    pmap["VPid"] = 1.0
    pmap["end_time"] = horizon_days

    pnames = sort(collect(keys(pmap)))
    pvals = [pmap[n] for n in pnames]

    mdl = TCellEngagerQSP.build_model_with_variant_ids(variant_ids, pnames, pvals)

    ctx = TCellEngagerQSP.CanonicalSimContext(
        copy(mdl.z),
        mdl.pvals,
        TCellEngagerQSP.InfusionEvent[],
    )
    rhs_fun = mdl.canonical_rhs
    jac_fun = use_mtk_jac ? TCellEngagerQSP.make_mtk_dense_jacobian(mdl) : nothing
    ode_rhs = isnothing(jac_fun) ? rhs_fun : ODEFunction(rhs_fun; jac = jac_fun)
    target_idx = mdl.state_to_idx["TDBc_ugperkg"]
    bt_idx = mdl.state_to_idx["Btumor"]

    dose_map = Dict{Float64, Float64}()
    for (t, d) in zip(dose_days, dose_mg)
        amt = d * 1000.0 / bw_kg
        dose_map[t] = get(dose_map, t, 0.0) + amt
    end
    dose_at_t0 = get(dose_map, 0.0, 0.0)
    dose_map[0.0] = 0.0
    cb_times = sort([t for t in keys(dose_map) if t > 0.0])

    function affect!(integrator)
        t = integrator.t
        amt = get(dose_map, t, NaN)
        if isnan(amt)
            for (tt, vv) in dose_map
                if isapprox(t, tt; atol = 1e-8, rtol = 0.0)
                    amt = vv
                    break
                end
            end
        end
        if !isnan(amt) && amt != 0.0
            integrator.u[target_idx] += amt
        end
    end

    u0 = copy(mdl.u0)
    if dose_at_t0 != 0.0
        u0[target_idx] += dose_at_t0
    end
    cb = PresetTimeCallback(cb_times, affect!; save_positions = (true, true))
    prob = ODEProblem(ode_rhs, u0, (0.0, horizon_days), ctx)
    saveat = collect(0.0:save_dt:horizon_days)
    sol = solve(
        prob,
        alg;
        abstol = TCellEngagerQSP.SOLVER_ABSTOL,
        reltol = TCellEngagerQSP.SOLVER_RELTOL,
        callback = cb,
        tstops = cb_times,
        d_discontinuities = cb_times,
        saveat = saveat,
    )

    if sol.retcode != SciMLBase.ReturnCode.Success
        return NaN
    end
    bt = [Float64(u[bt_idx]) for u in sol.u]
    bt0 = max(bt[1], 1e-12)
    spd_pct = @. 100.0 * (bt / bt0 - 1.0)
    return minimum(spd_pct)
end

function main()
    out_dir_env = get(ENV, "EFAST_OUT_DIR", joinpath(REPO_ROOT, "generated", "figures", "sensitivity", "efast_best_spd"))
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(REPO_ROOT, out_dir_env)
    mkpath(out_dir)

    N = parse(Int, get(ENV, "EFAST_N", "65"))
    M = parse(Int, get(ENV, "EFAST_M", "4"))
    seed = parse(Int, get(ENV, "EFAST_SEED", "20260305"))
    bw_kg = parse(Float64, get(ENV, "EFAST_BW_KG", "70.0"))
    horizon_days = parse(Float64, get(ENV, "EFAST_HORIZON_DAYS", "168.0"))
    save_dt = parse(Float64, get(ENV, "EFAST_SAVE_DT", "0.25"))
    spread_mode = lowercase(get(ENV, "EFAST_SPREAD_MODE", "primary"))
    max_evals = parse(Int, get(ENV, "EFAST_MAX_EVALS", "0"))
    fail_fill = parse(Float64, get(ENV, "EFAST_FAIL_FILL_BEST_SPD", "100.0"))

    solver_name_full = lowercase(get(ENV, "TCE_SOLVER", "cvode_bdf"))
    solver_name = replace(solver_name_full, "_mtkjac" => "")
    use_mtk_jac = endswith(solver_name_full, "_mtkjac") || get(ENV, "TCE_MTK_JAC", "0") == "1"
    alg = TCellEngagerQSP.make_solver_alg(solver_name)

    phase1_active_variant_ids = [5, 9, 14, 20, 24, 25, 27, 28]
    variants_mode = lowercase(get(ENV, "EFAST_VARIANTS_MODE", "matlab_empty"))
    variant_ids = if variants_mode == "matlab_empty"
        Int[]
    elseif variants_mode == "active_phase1"
        phase1_active_variant_ids
    else
        error("Unsupported EFAST_VARIANTS_MODE=$variants_mode (use matlab_empty or active_phase1)")
    end

    targets_json = load_json(joinpath(REPO_ROOT, "generated", "vpop_targets", "vpop_calibration_targets.json"))
    k_bounds, bt_bounds, spread_bounds = get_hard_bounds(targets_json, spread_mode)
    priors = read_param_bounds_csv(joinpath(REPO_ROOT, "generated", "vpop_targets", "parameter_prior_bounds_from_params_sheet.csv"))

    extra_names = String[
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

    missing_priors = [nm for nm in extra_names if !haskey(priors, nm)]
    if !isempty(missing_priors)
        error("Missing prior bounds for: $(join(missing_priors, ", "))")
    end

    param_names = String["kBtumorprolif", "BT_ratio_tumor_init", "tumor_burden_factor"]
    bounds = Tuple{Float64, Float64}[k_bounds, bt_bounds, spread_bounds]
    for nm in extra_names
        push!(param_names, nm)
        push!(bounds, priors[nm])
    end
    D = length(param_names)

    rng = MersenneTwister(seed)
    X, omega0 = salib_fast_sample(bounds, param_names, N, M, rng)
    n_total = size(X, 1)
    n_eval = max_evals > 0 ? min(max_evals, n_total) : n_total

    base_pmap = load_base_param_map()
    dose_days, dose_mg = build_regimen_8cycle_stepup()

    Y = fill(NaN, n_total)
    status = fill("not_run", n_total)
    failures = 0

    println("eFAST best%SPD run:")
    println("  D=$D, N=$N, total_evals=$n_total, evals_to_run=$n_eval, omega0=$omega0")
    println("  solver=$(solver_name_full), variants_mode=$variants_mode, spread_mode=$spread_mode")
    println("  horizon_days=$horizon_days, save_dt=$save_dt")

    for i in 1:n_eval
        vals = Dict{String, Float64}()
        for (j, nm) in enumerate(param_names)
            vals[nm] = X[i, j]
        end
        y = try
            eval_best_spd_for_sample(
                vals,
                base_pmap,
                extra_names,
                variant_ids,
                dose_days,
                dose_mg,
                horizon_days,
                save_dt,
                bw_kg,
                alg,
                use_mtk_jac,
            )
        catch err
            @warn "Sample failed" idx = i err
            NaN
        end

        if isfinite(y)
            Y[i] = y
            status[i] = "ok"
        else
            Y[i] = fail_fill
            status[i] = "failed_fill"
            failures += 1
        end
        if i % 25 == 0 || i == n_eval
            println("  progress $i / $n_eval")
        end
    end

    if n_eval < n_total
        Y[(n_eval + 1):end] .= fail_fill
        status[(n_eval + 1):end] .= "not_run_fill"
    end

    S1, ST = salib_fast_analyze(Y, D, N, M)
    ord = sortperm(ST; rev = true)
    rank = zeros(Int, D)
    for (r, idx) in enumerate(ord)
        rank[idx] = r
    end

    idx_rows = DataFrame(
        parameter = param_names,
        lower_bound = [b[1] for b in bounds],
        upper_bound = [b[2] for b in bounds],
        S1 = S1,
        ST = ST,
        rank_ST = rank,
    )
    sort!(idx_rows, :rank_ST)

    sample_rows = DataFrame()
    sample_rows[!, "sample_idx"] = collect(1:n_total)
    for (j, nm) in enumerate(param_names)
        sample_rows[!, nm] = X[:, j]
    end
    sample_rows[!, "best_spd_pct"] = Y
    sample_rows[!, "status"] = status

    idx_path = joinpath(out_dir, "efast_best_spd_indices.csv")
    samples_path = joinpath(out_dir, "efast_best_spd_samples.csv")
    meta_path = joinpath(out_dir, "efast_best_spd_meta.json")

    CSV.write(idx_path, idx_rows)
    CSV.write(samples_path, sample_rows)

    meta = Dict(
        "method" => "eFAST (SALib-compatible sampling/analysis implementation)",
        "output" => "best_spd_pct = min_t 100*(Btumor(t)/Btumor(0)-1), t∈[0,horizon]",
        "horizon_days" => horizon_days,
        "dose_days" => dose_days,
        "dose_mg" => dose_mg,
        "n_params" => D,
        "n_base_samples" => N,
        "n_total_evals" => n_total,
        "n_eval_executed" => n_eval,
        "M" => M,
        "omega0" => omega0,
        "solver" => solver_name_full,
        "variants_mode" => variants_mode,
        "spread_mode" => spread_mode,
        "failures" => failures,
        "fail_fill_value" => fail_fill,
        "param_names" => param_names,
        "bounds" => [Dict("name" => nm, "lb" => b[1], "ub" => b[2]) for (nm, b) in zip(param_names, bounds)],
    )
    open(meta_path, "w") do io
        JSON3.pretty(io, meta)
    end

    println("Wrote:")
    println("  " * idx_path)
    println("  " * samples_path)
    println("  " * meta_path)

    println("Top ST parameters:")
    for r in eachrow(idx_rows[1:min(8, nrow(idx_rows)), :])
        println("  ", r.parameter, " : ST=", round(r.ST, digits = 4), " (S1=", round(r.S1, digits = 4), ")")
    end
end

main()
