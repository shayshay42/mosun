using CSV
using DataFrames
using Dates
using DifferentialEquations
using ForwardDiff
using JSON3
using Optim
using Printf
using Random
using SciMLBase
using Statistics
using Base.Threads

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

const REPO_ROOT = TCellEngagerQSP.REPO_ROOT

struct MatchConfig
    param_names::Vector{String}
    bounds::Vector{Tuple{Float64, Float64}}
    param_index::Dict{String, Int}
    base_pvals::Vector{Float64}
    base_bpbo::Float64
    base_trpbo::Float64
    base_btumor_perml::Float64
    base_bt_ratio::Float64
    idx_kbtumorprolif::Int
    idx_kbptumor::Int
    idx_ktrptumor::Int
    rhs_fun::Function
    u0::Vector{Float64}
    bt_idx::Int
    tdbc_idx::Int
    dose_days::Vector{Float64}
    dose_mg::Vector{Float64}
    bw_kg::Float64
    horizon_days::Float64
    save_dt_opt::Float64
    save_dt_report::Float64
    smooth_tau::Float64
    target_sigma::Float64
    l2_seed_weight::Float64
    alg
    abstol::Float64
    reltol::Float64
    maxiters::Int
end

function parse_float_list(txt::String)
    vals = Float64[]
    for tok in split(txt, ",")
        s = strip(tok)
        isempty(s) && continue
        push!(vals, parse(Float64, s))
    end
    return vals
end

function parse_name_list(txt::String)
    vals = String[]
    for tok in split(txt, ",")
        s = strip(tok)
        isempty(s) && continue
        push!(vals, s)
    end
    return vals
end

function sigmoid(x)
    return inv(one(x) + exp(-x))
end

function logit(x::Float64)
    xc = min(max(x, 1e-9), 1.0 - 1e-9)
    return log(xc / (1.0 - xc))
end

function smoothmin(v::AbstractVector{T}, tau::Float64) where {T}
    m = minimum(v)
    ex = exp.((m .- v) ./ tau)
    return m - T(tau) * log(sum(ex))
end

function row_float(row, nm::String, fallback::Real)
    sym = Symbol(nm)
    if !(sym in propertynames(row))
        return Float64(fallback)
    end
    v = row[sym]
    if ismissing(v)
        return Float64(fallback)
    end
    fv = try
        Float64(v)
    catch
        NaN
    end
    return isfinite(fv) ? fv : Float64(fallback)
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
    pmap["Bpbo_perml"] = get(pmap, "Bpbo_perml", get(pmap, "Bpbref_perml", 500_000.0))
    pmap["Trpbo_perml"] = get(pmap, "Trpbo_perml", get(pmap, "Trpbref_perml", 500_000.0))
    pmap["Bpbref_perml"] = get(pmap, "Bpbref_perml", pmap["Bpbo_perml"])
    pmap["Trpbref_perml"] = get(pmap, "Trpbref_perml", pmap["Trpbo_perml"])
    return pmap
end

function load_bounds(param_names::Vector{String}, spread_mode::String)
    priors_csv = joinpath(REPO_ROOT, "generated", "vpop_targets", "parameter_prior_bounds_from_params_sheet.csv")
    targets_json = joinpath(REPO_ROOT, "generated", "vpop_targets", "vpop_calibration_targets.json")
    if !isfile(priors_csv)
        error("Missing priors CSV at $priors_csv")
    end
    if !isfile(targets_json)
        error("Missing targets JSON at $targets_json")
    end
    pri = DataFrame(CSV.File(priors_csv))
    pri.name = string.(pri.name)
    tj = JSON3.read(read(targets_json, String))

    k_bounds = nothing
    bt_bounds = nothing
    spread_primary = (0.9, 1.1)
    spread_sens = (0.5, 1.5)
    for t in tj["hard_priors"]
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
            spread_sens = (Float64(b["lower"]), Float64(b["upper"]))
        end
    end
    if k_bounds === nothing || bt_bounds === nothing
        error("Missing hard prior bounds in vpop_calibration_targets.json")
    end
    spread_bounds = lowercase(spread_mode) == "sensitivity" ? spread_sens : spread_primary

    out = Tuple{Float64, Float64}[]
    for nm in param_names
        if nm == "kBtumorprolif"
            push!(out, k_bounds)
        elseif nm == "BT_ratio_tumor_init"
            push!(out, bt_bounds)
        elseif nm == "tumor_burden_factor"
            push!(out, spread_bounds)
        else
            sub = pri[pri.name .== nm, :]
            if nrow(sub) < 1
                error("Missing prior bounds for parameter '$nm' in $priors_csv")
            end
            lb = Float64(sub.cyno_lb[1])
            ub = Float64(sub.cyno_ub[1])
            if !(isfinite(lb) && isfinite(ub) && ub > lb)
                error("Invalid bounds for '$nm': lb=$lb ub=$ub")
            end
            push!(out, (lb, ub))
        end
    end
    return out
end

function load_target_spd_values(path::String, col::String)
    if !isfile(path)
        error("Missing SPD data CSV at $path")
    end
    df = DataFrame(CSV.File(path))
    nms = string.(names(df))
    if !(col in nms)
        error("Column '$col' not present in $path")
    end
    cname = nms[findfirst(==(col), nms)]
    vals = Float64[]
    for v in df[!, cname]
        if ismissing(v)
            continue
        end
        fv = try
            Float64(v)
        catch
            NaN
        end
        if isfinite(fv)
            push!(vals, fv)
        end
    end
    isempty(vals) && error("No finite values found in $path for column '$col'")
    return vals
end

function sample_from_bounds(rng::AbstractRNG, lb::Float64, ub::Float64)
    if lb > 0.0 && ub / lb >= 1.5
        return exp(rand(rng) * (log(ub) - log(lb)) + log(lb))
    end
    return lb + rand(rng) * (ub - lb)
end

function theta_to_bounded(theta::AbstractVector, bounds::Vector{Tuple{Float64, Float64}})
    x = similar(theta)
    for i in eachindex(theta)
        lb, ub = bounds[i]
        s = sigmoid(theta[i])
        x[i] = lb + (ub - lb) * s
    end
    return x
end

function build_pfull(cfg::MatchConfig, x::AbstractVector{T}) where {T}
    pfull = T.(cfg.base_pvals)
    xmap = Dict{String, T}()
    for (nm, xv) in zip(cfg.param_names, x)
        xmap[nm] = xv
    end

    kbt = get(xmap, "kBtumorprolif", pfull[cfg.idx_kbtumorprolif])
    bt_ratio = get(xmap, "BT_ratio_tumor_init", T(cfg.base_bt_ratio))
    burden = get(xmap, "tumor_burden_factor", one(T))
    btumor_perml = T(cfg.base_btumor_perml) * burden
    pfull[cfg.idx_kbtumorprolif] = kbt
    pfull[cfg.idx_kbptumor] = btumor_perml / T(cfg.base_bpbo)
    pfull[cfg.idx_ktrptumor] = btumor_perml / (max(bt_ratio, T(1e-12)) * T(cfg.base_trpbo))

    for (nm, xv) in xmap
        if nm in ("kBtumorprolif", "BT_ratio_tumor_init", "tumor_burden_factor")
            continue
        end
        if haskey(cfg.param_index, nm)
            pfull[cfg.param_index[nm]] = xv
        end
    end
    return pfull
end

function simulate_best_spd(cfg::MatchConfig, pfull::AbstractVector{T}, save_dt::Float64, do_smooth::Bool) where {T}
    u = T.(cfg.u0)
    dose_map = Dict{Float64, T}()
    for (t, d) in zip(cfg.dose_days, cfg.dose_mg)
        amt = T(d * 1000.0 / cfg.bw_kg)
        dose_map[t] = get(dose_map, t, zero(T)) + amt
    end
    if haskey(dose_map, 0.0)
        u[cfg.tdbc_idx] += dose_map[0.0]
        dose_map[0.0] = zero(T)
    end

    segment_ends = sort(unique(vcat([t for t in keys(dose_map) if t > 0.0], [cfg.horizon_days])))
    bt_trace = T[]
    t_curr = 0.0

    for te in segment_ends
        if te > t_curr
            tvec = collect(t_curr:save_dt:te)
            if isempty(tvec) || !isapprox(tvec[end], te; atol = 1e-10, rtol = 0.0)
                push!(tvec, te)
            end
            prob = ODEProblem(cfg.rhs_fun, u, (t_curr, te), pfull)
            sol = solve(
                prob,
                cfg.alg;
                abstol = cfg.abstol,
                reltol = cfg.reltol,
                saveat = tvec,
                tstops = [te],
                maxiters = cfg.maxiters,
            )
            if sol.retcode != SciMLBase.ReturnCode.Success
                return T(100.0)
            end
            start_idx = isempty(bt_trace) ? 1 : 2
            for i in start_idx:length(sol.u)
                push!(bt_trace, sol.u[i][cfg.bt_idx])
            end
            u = sol.u[end]
        end
        if haskey(dose_map, te)
            u[cfg.tdbc_idx] += dose_map[te]
        end
        t_curr = te
    end

    isempty(bt_trace) && return T(100.0)
    bt0 = max(bt_trace[1], T(1e-12))
    spd = @. T(100.0) * (bt_trace / bt0 - one(T))
    spd_eval = length(spd) > 1 ? spd[2:end] : spd
    return do_smooth ? smoothmin(spd_eval, cfg.smooth_tau) : minimum(spd_eval)
end

function row_to_theta(row, cfg::MatchConfig)
    theta = zeros(Float64, length(cfg.param_names))
    for (i, nm) in enumerate(cfg.param_names)
        lb, ub = cfg.bounds[i]
        x0 = if nm == "kBtumorprolif"
            row_float(row, nm, cfg.base_pvals[cfg.idx_kbtumorprolif])
        elseif nm == "BT_ratio_tumor_init"
            row_float(row, nm, cfg.base_bt_ratio)
        elseif nm == "tumor_burden_factor"
            row_float(row, nm, 1.0)
        else
            row_float(row, nm, (lb + ub) / 2)
        end
        frac = (x0 - lb) / max(ub - lb, 1e-12)
        theta[i] = logit(frac)
    end
    return theta
end

function theta_to_row_dict(theta::AbstractVector{Float64}, cfg::MatchConfig)
    x = Float64.(theta_to_bounded(theta, cfg.bounds))
    xmap = Dict{String, Float64}(nm => xv for (nm, xv) in zip(cfg.param_names, x))
    bt_ratio = get(xmap, "BT_ratio_tumor_init", cfg.base_bt_ratio)
    burden = get(xmap, "tumor_burden_factor", 1.0)
    btumor_perml = cfg.base_btumor_perml * burden
    kbptumor = btumor_perml / cfg.base_bpbo
    ktrptumor = btumor_perml / (max(bt_ratio, 1e-12) * cfg.base_trpbo)

    row = Dict{String, Any}()
    row["Bpbo_perml"] = cfg.base_bpbo
    row["Bpbref_perml"] = cfg.base_bpbo
    row["Trpbo_perml"] = cfg.base_trpbo
    row["Trpbref_perml"] = cfg.base_trpbo
    row["kBtumorprolif"] = xmap["kBtumorprolif"]
    row["KBptumor"] = kbptumor
    row["KTrptumor"] = ktrptumor
    row["BT_ratio_tumor_init"] = bt_ratio
    row["Btumor_perml_init"] = btumor_perml
    row["tumor_burden_factor"] = burden
    for nm in cfg.param_names
        row[nm] = xmap[nm]
    end
    return row
end

function build_random_seed_tab(cfg::MatchConfig, n_seed::Int, rng_seed::Int)
    rng = MersenneTwister(rng_seed)
    rows = Vector{Dict{String, Any}}(undef, n_seed)
    for i in 1:n_seed
        theta = zeros(Float64, length(cfg.param_names))
        for j in eachindex(cfg.param_names)
            lb, ub = cfg.bounds[j]
            x = sample_from_bounds(rng, lb, ub)
            frac = (x - lb) / max(ub - lb, 1e-12)
            theta[j] = logit(frac)
        end
        row = theta_to_row_dict(theta, cfg)
        row["patient_id"] = i
        rows[i] = row
    end
    return DataFrame(rows)
end

function wasserstein_distance_1d(x::Vector{Float64}, y::Vector{Float64})
    xs = x[isfinite.(x)]
    ys = y[isfinite.(y)]
    if isempty(xs) || isempty(ys)
        return Inf
    end
    sort!(xs)
    sort!(ys)
    z = sort(vcat(xs, ys))
    if length(z) < 2
        return 0.0
    end
    dz = diff(z)
    xcdf = searchsortedlast.(Ref(xs), z[1:(end - 1)]) ./ length(xs)
    ycdf = searchsortedlast.(Ref(ys), z[1:(end - 1)]) ./ length(ys)
    return sum(abs.(xcdf .- ycdf) .* dz)
end

function quantile_pick_indices(n_pool::Int, n_target::Int)
    if n_target == 1
        return [clamp(div(n_pool + 1, 2), 1, n_pool)]
    end
    out = Int[]
    for j in 1:n_target
        idx = floor(Int, (j - 1) * (n_pool - 1) / (n_target - 1)) + 1
        push!(out, idx)
    end
    return out
end

function optimal_monotone_seed_pick(target_vals::Vector{Float64}, seed_vals_sorted::Vector{Float64})
    m = length(target_vals)
    n = length(seed_vals_sorted)
    n < m && error("Need at least as many seeds as targets for monotone assignment.")
    dp = fill(Inf, m + 1, n + 1)
    take = falses(m + 1, n + 1)
    dp[1, :] .= 0.0

    for k in 2:(m + 1)
        for j in 2:(n + 1)
            best = dp[k, j - 1]
            use_take = false
            cand = dp[k - 1, j - 1] + abs(seed_vals_sorted[j - 1] - target_vals[k - 1])
            if cand < best
                best = cand
                use_take = true
            end
            dp[k, j] = best
            take[k, j] = use_take
        end
    end

    idx = zeros(Int, m)
    k = m + 1
    j = n + 1
    while k > 1 && j > 1
        if take[k, j]
            idx[k - 1] = j - 1
            k -= 1
            j -= 1
        else
            j -= 1
        end
    end
    any(==(0), idx) && error("Failed to reconstruct monotone seed assignment.")
    return idx
end

function fit_target(
    target_spd::Float64,
    seed_row,
    cfg::MatchConfig,
    opt_iters::Int,
)
    theta_seed = row_to_theta(seed_row, cfg)
    f_scalar = function(th)
        x = theta_to_bounded(th, cfg.bounds)
        pfull = build_pfull(cfg, x)
        best_spd = simulate_best_spd(cfg, pfull, cfg.save_dt_opt, true)
        data_term = 0.5 * ((best_spd - target_spd) / cfg.target_sigma)^2
        reg_term = cfg.l2_seed_weight * sum((th .- theta_seed) .* (th .- theta_seed))
        return data_term + reg_term
    end

    theta_best = copy(theta_seed)
    obj_best = Float64(f_scalar(theta_seed))
    status = "seed"
    if opt_iters > 0
        f_opt = x -> Float64(f_scalar(x))
        g_opt! = (G, x) -> ForwardDiff.gradient!(G, f_scalar, x)
        opt = Optim.Options(iterations = opt_iters, show_trace = false)
        res = try
            optimize(f_opt, g_opt!, theta_seed, LBFGS(), opt)
        catch
            nothing
        end
        if res !== nothing
            theta_best = Optim.minimizer(res)
            obj_best = Float64(Optim.minimum(res))
            status = string(Optim.converged(res) ? "optimized_converged" : "optimized")
        else
            status = "opt_failed_seed_kept"
        end
    end

    xbest = theta_to_bounded(theta_best, cfg.bounds)
    pfull = build_pfull(cfg, xbest)
    achieved_spd = Float64(simulate_best_spd(cfg, pfull, cfg.save_dt_report, false))
    return (; theta_best, achieved_spd, objective = obj_best, status)
end

function curve_df(vals::Vector{Float64}, dataset::String)
    xs = sort(copy(vals))
    n = length(xs)
    rank = n <= 1 ? [0.0] : collect(range(0.0, 1.0; length = n))
    return DataFrame(dataset = fill(dataset, n), rank_frac = rank, best_spd_pct = xs)
end

function main()
    param_names = parse_name_list(get(
        ENV,
        "SPD_MATCH_PARAM_NAMES",
        "kBtumorprolif,BT_ratio_tumor_init,tumor_burden_factor,fTadeact,kTaapop,fTaprolif,fBkill",
    ))
    spread_mode = lowercase(get(ENV, "SPD_MATCH_SPREAD_MODE", "sensitivity"))
    target_sigma = parse(Float64, get(ENV, "SPD_MATCH_TARGET_SIGMA", "4.0"))
    l2_seed_weight = parse(Float64, get(ENV, "SPD_MATCH_L2_SEED_WEIGHT", "1e-3"))
    smooth_tau = parse(Float64, get(ENV, "SPD_MATCH_SMOOTH_TAU", "2.0"))
    opt_iters = parse(Int, get(ENV, "SPD_MATCH_OPT_ITERS", "8"))
    n_targets_cap = parse(Int, get(ENV, "SPD_MATCH_N_TARGETS", "0"))
    n_seed_cap = parse(Int, get(ENV, "SPD_MATCH_N_SEEDS", "0"))
    seed_source = lowercase(get(ENV, "SPD_MATCH_SEED_SOURCE", "csv"))
    seed_rng_seed = parse(Int, get(ENV, "SPD_MATCH_RANDOM_SEED", "20260307"))

    dose_days = parse_float_list(get(ENV, "SPD_MATCH_DOSE_DAYS", "0,7,14,21,42,63,84,105,126,147"))
    dose_mg = parse_float_list(get(ENV, "SPD_MATCH_DOSE_MG", "1,2,60,60,30,30,30,30,30,30"))
    bw_kg = parse(Float64, get(ENV, "SPD_MATCH_BW_KG", "70.0"))
    horizon_days = parse(Float64, get(ENV, "SPD_MATCH_HORIZON_DAYS", "168"))
    save_dt_opt = parse(Float64, get(ENV, "SPD_MATCH_SAVE_DT_OPT", "7.0"))
    save_dt_report = parse(Float64, get(ENV, "SPD_MATCH_SAVE_DT_REPORT", "1.0"))
    solver_name = lowercase(get(ENV, "SPD_MATCH_SOLVER", "tsit5"))
    abstol = parse(Float64, get(ENV, "SPD_MATCH_ABSTOL", "1e-8"))
    reltol = parse(Float64, get(ENV, "SPD_MATCH_RELTOL", "1e-6"))
    maxiters = parse(Int, get(ENV, "SPD_MATCH_MAXITERS", "100000000"))

    target_csv_env = get(
        ENV,
        "SPD_MATCH_TARGET_CSV",
        joinpath("generated", "reference", "musun_2022_88patients_SPD_from_vpop_generation_3a50cd1.csv"),
    )
    target_col = get(ENV, "SPD_MATCH_TARGET_COL", "SPD")
    seed_csv_env = get(
        ENV,
        "SPD_MATCH_SEED_PATIENTS_CSV",
        joinpath("generated", "vpop_best200_all_regimens_design", "patients.csv"),
    )
    out_dir_env = get(
        ENV,
        "SPD_MATCH_OUT_DIR",
        joinpath(REPO_ROOT, "generated", "figures", "optimization", "spd88_quantile_match"),
    )
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(REPO_ROOT, out_dir_env)
    mkpath(out_dir)

    target_csv = isabspath(target_csv_env) ? target_csv_env : joinpath(REPO_ROOT, target_csv_env)
    seed_csv = isabspath(seed_csv_env) ? seed_csv_env : joinpath(REPO_ROOT, seed_csv_env)

    target_vals_all = sort(load_target_spd_values(target_csv, target_col))
    target_vals = copy(target_vals_all)
    if n_targets_cap > 0 && n_targets_cap < length(target_vals_all)
        pick = quantile_pick_indices(length(target_vals_all), n_targets_cap)
        target_vals = target_vals_all[pick]
    end
    base_map = load_base_param_map()
    base_map["PKflag"] = 1.0
    base_map["fvalidation"] = 0.0
    base_map["VPid"] = 1.0
    base_map["end_time"] = horizon_days
    pnames = sort(collect(keys(base_map)))
    pvals = [base_map[n] for n in pnames]
    mdl = TCellEngagerQSP.build_model_with_variant_ids(Int[], pnames, pvals)
    rhs_fun = TCellEngagerQSP.build_symbolic_probe_rhs_function(
        mdl.name_to_idx,
        mdl.repeated_rule_defs,
        mdl.rate_exprs,
        mdl.stoich,
        length(mdl.u0),
        length(mdl.pvals),
    )
    param_index = Dict{String, Int}(nm => i for (i, nm) in enumerate(mdl.param_names))
    for req in ("kBtumorprolif", "KBptumor", "KTrptumor", "Bpbo_perml", "Trpbo_perml")
        haskey(param_index, req) || error("Missing model parameter index for $req")
    end
    base_bpbo = mdl.pvals[param_index["Bpbo_perml"]]
    base_trpbo = mdl.pvals[param_index["Trpbo_perml"]]
    base_btumor_perml = mdl.pvals[param_index["KBptumor"]] * base_bpbo
    base_bt_ratio = (mdl.pvals[param_index["KBptumor"]] * base_bpbo) / max(mdl.pvals[param_index["KTrptumor"]] * base_trpbo, 1e-12)

    cfg = MatchConfig(
        param_names,
        load_bounds(param_names, spread_mode),
        param_index,
        mdl.pvals,
        base_bpbo,
        base_trpbo,
        base_btumor_perml,
        base_bt_ratio,
        param_index["kBtumorprolif"],
        param_index["KBptumor"],
        param_index["KTrptumor"],
        rhs_fun,
        mdl.u0,
        mdl.state_to_idx["Btumor"],
        mdl.state_to_idx["TDBc_ugperkg"],
        dose_days,
        dose_mg,
        bw_kg,
        horizon_days,
        save_dt_opt,
        save_dt_report,
        smooth_tau,
        target_sigma,
        l2_seed_weight,
        TCellEngagerQSP.make_solver_alg(solver_name),
        abstol,
        reltol,
        maxiters,
    )

    seed_tab = if seed_source == "random"
        n_seed = n_seed_cap > 0 ? n_seed_cap : max(200, length(target_vals))
        build_random_seed_tab(cfg, n_seed, seed_rng_seed)
    else
        isfile(seed_csv) || error("Missing seed patients CSV at $seed_csv")
        tab = DataFrame(CSV.File(seed_csv))
        if n_seed_cap > 0
            tab = tab[1:min(n_seed_cap, nrow(tab)), :]
        end
        tab
    end
    if nrow(seed_tab) < length(target_vals)
        error("Need at least as many seed rows as target rows. Have $(nrow(seed_tab)), need $(length(target_vals)).")
    end

    n_seed = nrow(seed_tab)
    n_target = length(target_vals)
    println("Evaluating seed pool for target regimen ...")
    println("  n_seed=$n_seed n_target=$n_target threads=$(Threads.nthreads()) opt_iters=$opt_iters")
    seed_best = fill(NaN, n_seed)
    seed_lock = ReentrantLock()
    progress = Atomic{Int}(0)
    Threads.@threads for i in 1:n_seed
        row = seed_tab[i, :]
        theta0 = row_to_theta(row, cfg)
        pfull = build_pfull(cfg, theta_to_bounded(theta0, cfg.bounds))
        seed_best[i] = Float64(simulate_best_spd(cfg, pfull, cfg.save_dt_opt, false))
        done = atomic_add!(progress, 1) + 1
        if done % 20 == 0 || done == n_seed
            lock(seed_lock) do
                println("  seed progress $done / $n_seed")
            end
        end
    end
    seed_eval = copy(seed_tab)
    seed_eval[!, :seed_best_spd_pct] = seed_best
    sort!(seed_eval, :seed_best_spd_pct)
    CSV.write(joinpath(out_dir, "seed_pool_best_spd.csv"), seed_eval)

    seed_pick = optimal_monotone_seed_pick(target_vals, Float64.(seed_eval.seed_best_spd_pct))

    # Compile the optimization path once before threading.
    println("Compiling single-target fit path ...")
    _ = fit_target(target_vals[1], seed_eval[seed_pick[1], :], cfg, min(opt_iters, 1))

    println("Fitting quantile-matched parameterizations ...")
    fit_progress = Atomic{Int}(0)
    results = Vector{Any}(undef, n_target)
    Threads.@threads for i in 1:n_target
        srow = seed_eval[seed_pick[i], :]
        results[i] = fit_target(target_vals[i], srow, cfg, opt_iters)
        done = atomic_add!(fit_progress, 1) + 1
        if done % 10 == 0 || done == n_target
            lock(seed_lock) do
                println("  fit progress $done / $n_target")
            end
        end
    end

    matched_rows = Vector{Dict{String, Any}}(undef, n_target)
    patient_rows = Vector{Dict{String, Any}}(undef, n_target)
    achieved_vals = Float64[]
    seed_match_vals = Float64[]
    for i in 1:n_target
        res = results[i]
        srow = seed_eval[seed_pick[i], :]
        pmap = theta_to_row_dict(res.theta_best, cfg)
        pid = i
        prow = Dict{String, Any}(
            "patient_id" => pid,
            "Bpbo_perml" => pmap["Bpbo_perml"],
            "Bpbref_perml" => pmap["Bpbref_perml"],
            "Trpbo_perml" => pmap["Trpbo_perml"],
            "Trpbref_perml" => pmap["Trpbref_perml"],
            "kBtumorprolif" => pmap["kBtumorprolif"],
            "KBptumor" => pmap["KBptumor"],
            "KTrptumor" => pmap["KTrptumor"],
            "BT_ratio_tumor_init" => pmap["BT_ratio_tumor_init"],
            "Btumor_perml_init" => pmap["Btumor_perml_init"],
            "tumor_burden_factor" => pmap["tumor_burden_factor"],
        )
        for nm in cfg.param_names
            prow[nm] = pmap[nm]
        end
        patient_rows[i] = prow

        seed_spd = row_float(srow, "seed_best_spd_pct", NaN)
        push!(achieved_vals, res.achieved_spd)
        push!(seed_match_vals, seed_spd)
        mrow = copy(prow)
        mrow["target_rank"] = i
        mrow["target_best_spd_pct"] = target_vals[i]
        mrow["achieved_best_spd_pct"] = res.achieved_spd
        mrow["seed_patient_id"] = Int(round(row_float(srow, "patient_id", i)))
        mrow["seed_best_spd_pct"] = seed_spd
        mrow["seed_abs_error"] = abs(seed_spd - target_vals[i])
        mrow["fit_abs_error"] = abs(res.achieved_spd - target_vals[i])
        mrow["fit_objective"] = res.objective
        mrow["status"] = res.status
        matched_rows[i] = mrow
    end

    matched_df = DataFrame(matched_rows)
    patient_df = DataFrame(patient_rows)
    CSV.write(joinpath(out_dir, "matched_parameterizations_results.csv"), matched_df)
    CSV.write(joinpath(out_dir, "matched_patients.csv"), patient_df)

    curve = vcat(
        curve_df(target_vals, "target_real"),
        curve_df(seed_match_vals, "seed_quantile_match"),
        curve_df(achieved_vals, "optimized_quantile_match"),
    )
    CSV.write(joinpath(out_dir, "spd_distribution_curves.csv"), curve)

    target_sorted = sort(copy(target_vals))
    seed_sorted = sort(copy(seed_match_vals))
    achieved_sorted = sort(copy(achieved_vals))
    summary = Dict(
        "created_utc" => string(Dates.now(Dates.UTC)),
        "description" => "Quantile-matched parameterizations fit independently to real-data best%SPD order statistics.",
        "regimen_dose_days" => dose_days,
        "regimen_dose_mg" => dose_mg,
        "target_csv" => target_csv,
        "target_col" => target_col,
        "seed_source" => seed_source,
        "seed_patients_csv" => seed_csv,
        "param_names" => cfg.param_names,
        "bounds" => [Dict("name" => n, "lb" => b[1], "ub" => b[2]) for (n, b) in zip(cfg.param_names, cfg.bounds)],
        "n_seed" => n_seed,
        "n_target" => n_target,
        "spread_mode" => spread_mode,
        "target_sigma" => target_sigma,
        "l2_seed_weight" => l2_seed_weight,
        "smooth_tau" => smooth_tau,
        "opt_iters" => opt_iters,
        "solver" => solver_name,
        "save_dt_opt" => save_dt_opt,
        "save_dt_report" => save_dt_report,
        "seed_quantile_match_metrics" => Dict(
            "mae" => mean(abs.(seed_sorted .- target_sorted)),
            "rmse" => sqrt(mean((seed_sorted .- target_sorted) .^ 2)),
            "wasserstein_w1" => wasserstein_distance_1d(copy(seed_sorted), copy(target_sorted)),
        ),
        "optimized_match_metrics" => Dict(
            "mae" => mean(abs.(achieved_sorted .- target_sorted)),
            "rmse" => sqrt(mean((achieved_sorted .- target_sorted) .^ 2)),
            "wasserstein_w1" => wasserstein_distance_1d(copy(achieved_sorted), copy(target_sorted)),
        ),
        "outputs" => Dict(
            "matched_results_csv" => joinpath(out_dir, "matched_parameterizations_results.csv"),
            "matched_patients_csv" => joinpath(out_dir, "matched_patients.csv"),
            "seed_pool_csv" => joinpath(out_dir, "seed_pool_best_spd.csv"),
            "distribution_curve_csv" => joinpath(out_dir, "spd_distribution_curves.csv"),
        ),
    )
    open(joinpath(out_dir, "matched_summary.json"), "w") do io
        JSON3.pretty(io, summary)
    end

    println("Wrote $(joinpath(out_dir, "matched_parameterizations_results.csv"))")
    println("Wrote $(joinpath(out_dir, "matched_patients.csv"))")
    println("Wrote $(joinpath(out_dir, "seed_pool_best_spd.csv"))")
    println("Wrote $(joinpath(out_dir, "spd_distribution_curves.csv"))")
    println("Wrote $(joinpath(out_dir, "matched_summary.json"))")
    println(@sprintf("Seed W1 = %.6g", summary["seed_quantile_match_metrics"]["wasserstein_w1"]))
    println(@sprintf("Fit  W1 = %.6g", summary["optimized_match_metrics"]["wasserstein_w1"]))
end

main()
