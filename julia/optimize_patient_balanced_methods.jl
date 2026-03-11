using CSV
using DataFrames
using DifferentialEquations
using FiniteDiff
using ForwardDiff
using JSON3
using Optim
using Printf
using Random
using RuntimeGeneratedFunctions
using SciMLBase

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

RuntimeGeneratedFunctions.init(@__MODULE__)

function smoothmax(v::AbstractVector, tau::Real)
    m = maximum(v)
    return m + tau * log(sum(exp.((v .- m) ./ tau)))
end

function trapz_generic(t::Vector{Float64}, y::Vector{T}) where {T}
    if length(t) < 2 || length(y) < 2
        return zero(T)
    end
    s = zero(T)
    @inbounds for i in 1:(length(t) - 1)
        dt = t[i + 1] - t[i]
        if dt > 0.0
            s += T(0.5 * dt) * (y[i + 1] + y[i])
        end
    end
    return s
end

function collect_window_vals(t::Vector{Float64}, y::Vector{T}, t0::Float64, t1::Float64) where {T}
    vals = T[]
    @inbounds for i in eachindex(t)
        if t[i] >= t0 - 1e-10 && t[i] <= t1 + 1e-10
            push!(vals, y[i])
        end
    end
    if isempty(vals)
        idx = searchsortedfirst(t, t0)
        idx = clamp(idx, 1, length(t))
        return T[y[idx]]
    end
    return vals
end

function parse_float_list(s::AbstractString)
    vals = Float64[]
    for tok in split(String(s), ",")
        st = strip(tok)
        isempty(st) && continue
        push!(vals, parse(Float64, st))
    end
    return vals
end

function load_loss_config(path::String)
    cfg = JSON3.read(read(path, String))
    weights = Dict{String, Float64}()
    scales = Dict{String, Float64}()
    for (k, v) in pairs(cfg["loss_weights"])
        weights[String(k)] = Float64(v)
    end
    for (k, v) in pairs(cfg["loss_scales"])
        scales[String(k)] = Float64(v)
    end
    return (
        weights = weights,
        scales = scales,
        horizon_days = Float64(cfg["horizon_days"]),
        tox_window_days = Float64(cfg["tox_window_days"]),
        tox_tau = Float64(cfg["tox_tau"]),
    )
end

function build_ad_rhs_no_cast(mdl)
    name_to_idx = mdl.name_to_idx
    repeated_rule_defs = mdl.repeated_rule_defs
    rate_exprs = mdl.rate_exprs
    stoich = mdl.stoich

    repeated_assign_exprs = Any[]
    for (lhs_idx, rhs0) in repeated_rule_defs
        rhs =
            if occursin("PK_v26(", rhs0)
                "ifelse(TDBc_ugperkg / Vc_tdb > 1e-5, TDBc_ugperkg / Vc_tdb, 0.0)"
            else
                rhs0
            end
        rhs_expr = TCellEngagerQSP.compile_formula_expr(rhs, name_to_idx)
        push!(repeated_assign_exprs, :(z[$lhs_idx] = real($rhs_expr)))
    end

    reaction_blocks = Any[]
    for j in eachindex(rate_exprs)
        rvar = gensym(:rate)
        rexpr = TCellEngagerQSP.compile_formula_expr(rate_exprs[j], name_to_idx)
        block = Expr(:block, :($rvar = real($rexpr)))
        for (i, coeff) in stoich[j]
            push!(block.args, :(du[$i] += $(Float64(coeff)) * $rvar))
        end
        push!(reaction_blocks, block)
    end

    np = length(mdl.pvals)
    zlen = length(mdl.name_to_idx)
    fexpr = quote
        (du, u, pvals, t) -> begin
            ns = length(u)
            T = eltype(u)
            z = Vector{T}(undef, $zlen)
            @inbounds begin
                z[1:ns] .= u
                z[ns+1:ns+$np] .= pvals
                $(repeated_assign_exprs...)

                fill!(du, zero(T))
                $(reaction_blocks...)
            end
            return nothing
        end
    end

    rhs_rgf = RuntimeGeneratedFunction(
        TCellEngagerQSP,
        TCellEngagerQSP,
        TCellEngagerQSP.normalize_function_expr(fexpr),
    )
    pvals0 = copy(mdl.pvals)
    return (du, u, p, t) -> rhs_rgf(du, u, pvals0, t)
end

function build_symbol_observer_no_cast(mdl, symbol_name::String)
    name_to_idx = mdl.name_to_idx
    repeated_rule_defs = mdl.repeated_rule_defs
    zlen = length(name_to_idx)
    target_idx = name_to_idx[symbol_name]
    np = length(mdl.pvals)
    repeated_assign_exprs = Any[]
    for (lhs_idx, rhs0) in repeated_rule_defs
        rhs =
            if occursin("PK_v26(", rhs0)
                "ifelse(TDBc_ugperkg / Vc_tdb > 1e-5, TDBc_ugperkg / Vc_tdb, 0.0)"
            else
                rhs0
            end
        rhs_expr = TCellEngagerQSP.compile_formula_expr(rhs, name_to_idx)
        push!(repeated_assign_exprs, :(z[$lhs_idx] = real($rhs_expr)))
    end

    fexpr = quote
        (u, t, pvals) -> begin
            ns = length(u)
            T = eltype(u)
            z = Vector{T}(undef, $zlen)
            @inbounds begin
                z[1:ns] .= u
                z[ns+1:ns+$np] .= pvals
                $(repeated_assign_exprs...)
            end
            return z[$target_idx]
        end
    end
    obs_rgf = RuntimeGeneratedFunction(
        TCellEngagerQSP,
        TCellEngagerQSP,
        TCellEngagerQSP.normalize_function_expr(fexpr),
    )
    pvals0 = copy(mdl.pvals)
    return (u, t) -> obs_rgf(u, t, pvals0)
end

struct BalancedOptProblem
    patient_id::Int
    mdl
    rhs_fun::Function
    il6_obs_fun::Function
    alg
    e_tdbc::Vector{Float64}
    bt_idx::Int
    dose_times_days::Vector{Float64}
    x0_mg::Vector{Float64}
    lower_mg::Vector{Float64}
    upper_mg::Vector{Float64}
    mg_to_ugkg::Float64
    horizon_days::Float64
    tox_window_days::Float64
    tox_tau::Float64
    weights::Dict{String, Float64}
    scales::Dict{String, Float64}
    save_dt::Float64
    abstol::Float64
    reltol::Float64
    maxiters::Int
end

mutable struct ObjectiveProgress
    evals::Int
    best_loss::Float64
    best_x::Vector{Float64}
    start_t::Float64
    last_log_t::Float64
    log_every_evals::Int
    log_every_seconds::Float64
    max_evals::Int
    max_seconds::Float64
    method_name::String
    hist_eval::Vector{Int}
    hist_elapsed_s::Vector{Float64}
    hist_current_loss::Vector{Float64}
    hist_best_loss::Vector{Float64}
end

struct EvalBudgetReached <: Exception
    msg::String
end

function solve_segment(prob::BalancedOptProblem, u0, tspan::Tuple{Float64, Float64}, saveat::Vector{Float64})
    ode_prob = ODEProblem(prob.rhs_fun, u0, tspan, [1.0, 1.0, 1.0])
    sol = solve(
        ode_prob,
        prob.alg;
        abstol = prob.abstol,
        reltol = prob.reltol,
        saveat = saveat,
        tstops = [tspan[2]],
        maxiters = prob.maxiters,
    )
    if sol.retcode != SciMLBase.ReturnCode.Success
        error("segment solve failed with retcode=$(sol.retcode) tspan=$tspan")
    end
    return sol
end

function push_point!(t_all::Vector{Float64}, bt_all, il6_all, t::Float64, btv, il6v)
    if !isempty(t_all) && isapprox(t, t_all[end]; atol = 1e-12, rtol = 0.0)
        bt_all[end] = btv
        il6_all[end] = il6v
    else
        push!(t_all, t)
        push!(bt_all, btv)
        push!(il6_all, il6v)
    end
end

function make_saveat(t0::Float64, t1::Float64, dt::Float64)
    if t1 <= t0 + 1e-12
        return Float64[t1]
    end
    ts = collect(t0:dt:t1)
    if isempty(ts) || !isapprox(ts[end], t1; atol = 1e-12, rtol = 0.0)
        push!(ts, t1)
    end
    return ts
end

function simulate_trajectory(prob::BalancedOptProblem, doses_mg::AbstractVector{T}) where {T<:Real}
    nd = length(prob.dose_times_days)
    length(doses_mg) == nd || error("Expected $nd doses, got $(length(doses_mg))")
    d_ugkg = doses_mg .* T(prob.mg_to_ugkg)
    e_tdbc_t = T.(prob.e_tdbc)

    u_curr = T.(prob.mdl.u0)
    t_curr = 0.0

    t_all = Float64[]
    bt_all = T[]
    il6_all = T[]

    for i in 1:nd
        tdose = prob.dose_times_days[i]
        if tdose > t_curr + 1e-12
            saveat = make_saveat(t_curr, tdose, prob.save_dt)
            sol = solve_segment(prob, u_curr, (t_curr, tdose), saveat)
            for j in eachindex(sol.t)
                tj = Float64(sol.t[j])
                uj = sol.u[j]
                il6v = prob.il6_obs_fun(uj, sol.t[j])
                push_point!(t_all, bt_all, il6_all, tj, uj[prob.bt_idx], il6v)
            end
            u_curr = sol.u[end]
            t_curr = tdose
        end
        u_curr = u_curr .+ d_ugkg[i] .* e_tdbc_t
        il6_now = prob.il6_obs_fun(u_curr, t_curr)
        push_point!(t_all, bt_all, il6_all, t_curr, u_curr[prob.bt_idx], il6_now)
    end

    if t_curr < prob.horizon_days - 1e-12
        saveat = make_saveat(t_curr, prob.horizon_days, prob.save_dt)
        sol = solve_segment(prob, u_curr, (t_curr, prob.horizon_days), saveat)
        for j in eachindex(sol.t)
            tj = Float64(sol.t[j])
            uj = sol.u[j]
            il6v = prob.il6_obs_fun(uj, sol.t[j])
            push_point!(t_all, bt_all, il6_all, tj, uj[prob.bt_idx], il6v)
        end
        u_curr = sol.u[end]
    end

    return (t = t_all, bt = bt_all, il6 = il6_all)
end

function raw_metrics(prob::BalancedOptProblem, t::Vector{Float64}, bt, il6)
    T = eltype(bt)
    bt0 = bt[1]
    tox_peaks = T[]
    for td in prob.dose_times_days
        vals = collect_window_vals(t, il6, td, td + prob.tox_window_days)
        push!(tox_peaks, smoothmax(vals, prob.tox_tau))
    end
    return (
        tox_peak_mean = sum(tox_peaks) / T(length(tox_peaks)),
        tox_peak_max = maximum(tox_peaks),
        tox_auc = trapz_generic(t, max.(il6, zero(T))) / T(prob.horizon_days),
        tumor_terminal = bt[end] / (bt0 + T(1e-12)),
        tumor_auc = trapz_generic(t, max.(bt, zero(T))) / (T(prob.horizon_days) * (bt0 + T(1e-12))),
    )
end

function scaled_loss(prob::BalancedOptProblem, rm)
    T = typeof(rm.tox_peak_mean)
    w = prob.weights
    s = prob.scales
    loss =
        T(w["tox_peak_mean"]) * (rm.tox_peak_mean / T(max(s["tox_peak_mean"], 1e-12))) +
        T(w["tox_peak_max"]) * (rm.tox_peak_max / T(max(s["tox_peak_max"], 1e-12))) +
        T(w["tox_auc"]) * (rm.tox_auc / T(max(s["tox_auc"], 1e-12))) +
        T(w["tumor_terminal"]) * (rm.tumor_terminal / T(max(s["tumor_terminal"], 1e-12))) +
        T(w["tumor_auc"]) * (rm.tumor_auc / T(max(s["tumor_auc"], 1e-12)))
    return loss
end

function evaluate_loss(prob::BalancedOptProblem, doses_mg::AbstractVector{T}) where {T<:Real}
    sim = simulate_trajectory(prob, doses_mg)
    rm = raw_metrics(prob, sim.t, sim.bt, sim.il6)
    return scaled_loss(prob, rm)
end

function make_loss_objective(prob::BalancedOptProblem; method_name::String = "method", method_tag::String = "")
    big_penalty = 1.0e9
    base_max_evals = parse_int_env("OPTB_MAX_EVALS", 0)
    method_max_evals = isempty(method_tag) ? base_max_evals : parse_int_env("OPTB_$(method_tag)_MAX_EVALS", base_max_evals)
    base_max_seconds = parse_float_env("OPTB_MAX_SECONDS", 0.0)
    method_max_seconds = isempty(method_tag) ? base_max_seconds : parse_float_env("OPTB_$(method_tag)_MAX_SECONDS", base_max_seconds)
    t_now = time()
    prog = ObjectiveProgress(
        0,
        Inf,
        Float64[],
        t_now,
        t_now,
        parse_int_env("OPTB_PROGRESS_EVERY_EVALS", 25),
        parse_float_env("OPTB_PROGRESS_EVERY_SECONDS", 20.0),
        max(method_max_evals, 0),
        max(method_max_seconds, 0.0),
        method_name,
        Int[],
        Float64[],
        Float64[],
        Float64[],
    )
    obj = function (x::AbstractVector)
        if prog.max_seconds > 0.0 && (time() - prog.start_t) >= prog.max_seconds
            throw(EvalBudgetReached("max seconds reached for $(prog.method_name): $(prog.max_seconds)"))
        end
        is_dual = x[1] isa ForwardDiff.Dual
        if !is_dual
            xv = Float64.(x)
            if any(!isfinite, xv)
                return Float64(big_penalty)
            end
            viol = max.(prob.lower_mg .- xv, 0.0) .+ max.(xv .- prob.upper_mg, 0.0)
            if any(viol .> 0.0)
                return Float64(big_penalty + 1.0e6 * sum(abs2, viol))
            end
        end
        val =
        try
            evaluate_loss(prob, x)
        catch err
            if err isa EvalBudgetReached
                rethrow(err)
            end
            if is_dual
                one(x[1]) * big_penalty
            else
                Float64(big_penalty)
            end
        end
        if !is_dual
            prog.evals += 1
            v = Float64(val)
            if isfinite(v) && v < prog.best_loss
                prog.best_loss = v
                empty!(prog.best_x)
                append!(prog.best_x, Float64.(x))
            end
            if prog.max_evals > 0 && prog.evals >= prog.max_evals
                throw(EvalBudgetReached("max evals reached for $(prog.method_name): $(prog.max_evals)"))
            end
            now_t = time()
            push!(prog.hist_eval, prog.evals)
            push!(prog.hist_elapsed_s, now_t - prog.start_t)
            push!(prog.hist_current_loss, v)
            push!(prog.hist_best_loss, prog.best_loss)
            if (prog.evals % prog.log_every_evals == 0) || ((now_t - prog.last_log_t) >= prog.log_every_seconds)
                @printf("[%s] eval=%d best_loss=%.6g current_loss=%.6g\n", prog.method_name, prog.evals, prog.best_loss, v)
                flush(stdout)
                prog.last_log_t = now_t
            end
        end
        return val
    end
    return obj, prog
end

function parse_methods_env()
    default_methods = "simulated_annealing,nelder_mead,finite_diff_lbfgs,forward_ad_dto_lbfgs"
    raw = lowercase(strip(get(ENV, "OPTB_METHODS", default_methods)))
    methods = String[]
    for tok in split(raw, ",")
        m = strip(tok)
        isempty(m) && continue
        push!(methods, m)
    end
    isempty(methods) && error("No methods selected from OPTB_METHODS='$raw'")
    valid = Set(["simulated_annealing", "nelder_mead", "finite_diff_lbfgs", "forward_ad_dto_lbfgs"])
    for m in methods
        m in valid || error("Unknown method '$m'. Valid: $(collect(valid))")
    end
    return unique(methods)
end

function build_problem()
    repo_root = TCellEngagerQSP.REPO_ROOT
    design_dir_env = get(ENV, "OPTB_DESIGN_DIR", joinpath(repo_root, "generated", "random_regimen_84_design"))
    design_dir = isabspath(design_dir_env) ? design_dir_env : joinpath(repo_root, design_dir_env)
    reg_path = joinpath(design_dir, "regimen_events.csv")
    patients_path = joinpath(design_dir, "patients.csv")
    overrides_path = joinpath(design_dir, "dlbcl_param_overrides.csv")
    reg_df = DataFrame(CSV.File(reg_path))
    patients = DataFrame(CSV.File(patients_path))
    overrides = isfile(overrides_path) ? DataFrame(CSV.File(overrides_path)) : DataFrame(name = String[], value = Float64[])

    regimen_name = get(ENV, "OPTB_REGIMEN_NAME", String(first(reg_df.regimen)))
    sub = reg_df[reg_df.regimen .== regimen_name, :]
    sort!(sub, :event_idx)
    dose_times_days = Float64.(sub.time_day)
    x0_mg = Float64.(sub.dose_mg)

    patient_id = parse(Int, get(ENV, "OPTB_PATIENT_ID", "1"))
    prow_idx = findfirst(patients.patient_id .== patient_id)
    prow_idx === nothing && error("patient_id=$patient_id not found in $patients_path")
    prow = patients[prow_idx, :]

    pmap = Dict{String, Float64}()
    for rr in eachrow(overrides)
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

    loss_cfg_path_env = get(
        ENV,
        "OPTB_LOSS_CONFIG_JSON",
        joinpath(repo_root, "generated", "random_regimen_84_run_20260304", "balanced_loss", "balanced_multicycle_loss_config_random_regimen.json"),
    )
    loss_cfg_path = isabspath(loss_cfg_path_env) ? loss_cfg_path_env : joinpath(repo_root, loss_cfg_path_env)
    cfg = load_loss_config(loss_cfg_path)
    pmap["end_time"] = cfg.horizon_days

    pnames = sort(collect(keys(pmap)))
    pvals = [pmap[n] for n in pnames]

    mdl = TCellEngagerQSP.build_model_with_variant_ids(Int[], pnames, pvals)
    rhs_fun = build_ad_rhs_no_cast(mdl)
    il6_obs_fun = build_symbol_observer_no_cast(mdl, "IL6combo")
    bt_idx = mdl.state_to_idx["Btumor"]
    tdbc_idx = mdl.state_to_idx["TDBc_ugperkg"]
    e_tdbc = zeros(Float64, length(mdl.u0))
    e_tdbc[tdbc_idx] = 1.0

    lower_scalar = parse(Float64, get(ENV, "OPTB_DOSE_LOWER_MG", "0.0"))
    upper_scalar = parse(Float64, get(ENV, "OPTB_DOSE_UPPER_MG", "30.0"))
    lower_mg = fill(lower_scalar, length(x0_mg))
    upper_mg = fill(upper_scalar, length(x0_mg))

    alg = TCellEngagerQSP.make_solver_alg(lowercase(get(ENV, "OPTB_SOLVER", "tsit5")))
    save_dt = parse(Float64, get(ENV, "OPTB_SAVE_DT_DAYS", "0.1"))
    abstol = parse(Float64, get(ENV, "OPTB_ABSTOL", "1e-8"))
    reltol = parse(Float64, get(ENV, "OPTB_RELTOL", "1e-6"))
    maxiters = parse(Int, get(ENV, "OPTB_MAXITERS_SOLVE", "1000000"))

    return BalancedOptProblem(
        patient_id,
        mdl,
        rhs_fun,
        il6_obs_fun,
        alg,
        e_tdbc,
        bt_idx,
        dose_times_days,
        x0_mg,
        lower_mg,
        upper_mg,
        1000.0 / 70.0,
        cfg.horizon_days,
        cfg.tox_window_days,
        cfg.tox_tau,
        cfg.weights,
        cfg.scales,
        save_dt,
        abstol,
        reltol,
        maxiters,
    )
end

function parse_int_env(name::String, default::Int)
    v = strip(get(ENV, name, string(default)))
    isempty(v) && return default
    try
        return parse(Int, v)
    catch
        return default
    end
end

function parse_float_env(name::String, default::Float64)
    v = strip(get(ENV, name, string(default)))
    isempty(v) && return default
    try
        return parse(Float64, v)
    catch
        return default
    end
end

function make_optim_options(method_tag::String, iterations::Int)
    base_f = parse_int_env("OPTB_F_CALLS_LIMIT", 0)
    base_g = parse_int_env("OPTB_G_CALLS_LIMIT", 0)
    base_t = parse_float_env("OPTB_TIME_LIMIT_S", NaN)
    f_calls = parse_int_env("OPTB_$(method_tag)_F_CALLS_LIMIT", base_f)
    g_calls = parse_int_env("OPTB_$(method_tag)_G_CALLS_LIMIT", base_g)
    t_limit = parse_float_env("OPTB_$(method_tag)_TIME_LIMIT_S", base_t)
    return Optim.Options(
        iterations = iterations,
        show_trace = false,
        store_trace = false,
        f_calls_limit = max(f_calls, 0),
        g_calls_limit = max(g_calls, 0),
        time_limit = t_limit,
    )
end

function write_trajectory_csv(path::String, sim)
    df = DataFrame(
        time_day = sim.t,
        Btumor = Float64.(sim.bt),
        IL6combo = Float64.(sim.il6),
    )
    CSV.write(path, df)
end

function run_one_method(prob::BalancedOptProblem, method::String, method_out_dir::String)
    method_tag = if method == "simulated_annealing"
        "SA"
    elseif method == "nelder_mead"
        "NM"
    elseif method == "finite_diff_lbfgs"
        "FD"
    elseif method == "forward_ad_dto_lbfgs"
        "AD"
    else
        ""
    end
    obj, prog = make_loss_objective(prob; method_name = method, method_tag = method_tag)
    lower = copy(prob.lower_mg)
    upper = copy(prob.upper_mg)
    x0 = clamp.(copy(prob.x0_mg), lower, upper)

    sa_iters = parse(Int, get(ENV, "OPTB_SA_ITERS", "120"))
    nm_iters = parse(Int, get(ENV, "OPTB_NM_ITERS", "120"))
    lbfgs_default = get(ENV, "OPTB_LBFGS_ITERS", "80")
    fd_lbfgs_iters = parse(Int, get(ENV, "OPTB_FD_LBFGS_ITERS", lbfgs_default))
    ad_lbfgs_iters = parse(Int, get(ENV, "OPTB_AD_LBFGS_ITERS", lbfgs_default))
    sa_nt = parse(Int, get(ENV, "OPTB_SA_NT", "5"))
    sa_ns = parse(Int, get(ENV, "OPTB_SA_NS", "5"))
    sa_t0 = parse(Float64, get(ENV, "OPTB_SA_T0", "1.0"))
    sa_rt = parse(Float64, get(ENV, "OPTB_SA_RT", "0.9"))

    t0 = time()
    x_best = copy(x0)
    initial_loss = Float64(obj(x_best))
    loss_best = initial_loss
    status = "ok"
    println("[$method] start loss=$(initial_loss)")
    flush(stdout)
    try
        if method == "simulated_annealing"
            Random.seed!(parse(Int, get(ENV, "OPTB_RANDOM_SEED", "20260304")))
            sa_opt = make_optim_options("SA", sa_iters)
            sa_alg = SAMIN(nt = sa_nt, ns = sa_ns, t0 = sa_t0, rt = sa_rt, neps = 5, coverage_ok = true, verbosity = 0)
            res = optimize(obj, lower, upper, copy(x0), sa_alg, sa_opt)
            x_best = clamp.(Optim.minimizer(res), lower, upper)
            loss_best = Float64(obj(x_best))
        elseif method == "nelder_mead"
            nm_opt = make_optim_options("NM", nm_iters)
            res = optimize(obj, copy(x0), NelderMead(), nm_opt)
            x_best = clamp.(Optim.minimizer(res), lower, upper)
            loss_best = Float64(obj(x_best))
        elseif method == "finite_diff_lbfgs"
            g_fd! = (G, x) -> begin
                copyto!(G, FiniteDiff.finite_difference_gradient(obj, x))
                return G
            end
            lb_opt = make_optim_options("FD", fd_lbfgs_iters)
            res = optimize(obj, g_fd!, copy(x0), LBFGS(), lb_opt)
            x_best = clamp.(Optim.minimizer(res), lower, upper)
            loss_best = Float64(obj(x_best))
        elseif method == "forward_ad_dto_lbfgs"
            g_ad! = (G, x) -> begin
                ForwardDiff.gradient!(G, obj, x)
                return G
            end
            lb_opt = make_optim_options("AD", ad_lbfgs_iters)
            res = optimize(obj, g_ad!, copy(x0), LBFGS(), lb_opt)
            x_best = clamp.(Optim.minimizer(res), lower, upper)
            loss_best = Float64(obj(x_best))
        else
            error("Unknown method $method")
        end
    catch err
        status = "error: " * sprint(showerror, err)
        if (err isa EvalBudgetReached) && !isempty(prog.best_x)
            x_best = clamp.(copy(prog.best_x), lower, upper)
            loss_best = prog.best_loss
            status = "max_budget_reached"
        end
    end

    if !isempty(prog.best_x) && prog.best_loss < loss_best
        x_best = clamp.(copy(prog.best_x), lower, upper)
        loss_best = prog.best_loss
    end

    runtime_s = time() - t0
    println("[$method] done status=$status evals=$(prog.evals) runtime_s=$(round(runtime_s, digits=2)) best_loss=$(loss_best)")
    flush(stdout)

    mkpath(method_out_dir)
    hist_df = DataFrame(
        eval = prog.hist_eval,
        elapsed_s = prog.hist_elapsed_s,
        current_loss = prog.hist_current_loss,
        best_loss = prog.hist_best_loss,
    )
    CSV.write(joinpath(method_out_dir, "loss_history.csv"), hist_df)

    sim_initial = simulate_trajectory(prob, x0)
    sim_best = simulate_trajectory(prob, x_best)
    write_trajectory_csv(joinpath(method_out_dir, "trajectory_initial.csv"), sim_initial)
    write_trajectory_csv(joinpath(method_out_dir, "trajectory_best.csv"), sim_best)

    rm_initial = raw_metrics(prob, sim_initial.t, sim_initial.bt, sim_initial.il6)
    rm_best = raw_metrics(prob, sim_best.t, sim_best.bt, sim_best.il6)
    method_meta = Dict(
        "method" => method,
        "status" => status,
        "runtime_s" => runtime_s,
        "evals" => prog.evals,
        "initial_loss" => initial_loss,
        "best_loss" => loss_best,
        "improvement_abs" => initial_loss - loss_best,
        "improvement_pct" => initial_loss > 0 ? 100.0 * (initial_loss - loss_best) / initial_loss : 0.0,
        "initial_doses_mg" => x0,
        "best_doses_mg" => x_best,
        "initial_raw_terms" => Dict(
            "tox_peak_mean" => Float64(rm_initial.tox_peak_mean),
            "tox_peak_max" => Float64(rm_initial.tox_peak_max),
            "tox_auc" => Float64(rm_initial.tox_auc),
            "tumor_terminal" => Float64(rm_initial.tumor_terminal),
            "tumor_auc" => Float64(rm_initial.tumor_auc),
        ),
        "best_raw_terms" => Dict(
            "tox_peak_mean" => Float64(rm_best.tox_peak_mean),
            "tox_peak_max" => Float64(rm_best.tox_peak_max),
            "tox_auc" => Float64(rm_best.tox_auc),
            "tumor_terminal" => Float64(rm_best.tumor_terminal),
            "tumor_auc" => Float64(rm_best.tumor_auc),
        ),
    )
    open(joinpath(method_out_dir, "method_diagnostics.json"), "w") do io
        JSON3.pretty(io, method_meta)
    end

    return (
        method = method,
        initial_loss = initial_loss,
        best_loss = loss_best,
        runtime_s = runtime_s,
        status = status,
        evals = prog.evals,
        doses_mg = join(string.(round.(x_best, digits = 4)), ","),
    )
end

function run_methods(prob::BalancedOptProblem, methods::Vector{String}, out_dir::String)
    results = DataFrame(
        method = String[],
        initial_loss = Float64[],
        best_loss = Float64[],
        improvement_abs = Float64[],
        improvement_pct = Float64[],
        runtime_s = Float64[],
        status = String[],
        evals = Int[],
        doses_mg = String[],
    )
    single_method_mode = length(methods) == 1
    out_dir_name = lowercase(basename(normpath(out_dir)))
    for method in methods
        method_out_dir = (single_method_mode && out_dir_name == lowercase(method)) ? out_dir : joinpath(out_dir, method)
        row = run_one_method(prob, method, method_out_dir)
        improvement_abs = row.initial_loss - row.best_loss
        improvement_pct = row.initial_loss > 0 ? 100.0 * improvement_abs / row.initial_loss : 0.0
        push!(results, (row.method, row.initial_loss, row.best_loss, improvement_abs, improvement_pct, row.runtime_s, row.status, row.evals, row.doses_mg))
    end
    sort!(results, :best_loss)
    return results
end

function main()
    prob = build_problem()
    methods = parse_methods_env()
    out_dir_default = joinpath(TCellEngagerQSP.REPO_ROOT, "generated", "random_regimen_84_run_20260304", "optimization_patient_balanced")
    out_dir_env = get(ENV, "OPTB_OUT_DIR", out_dir_default)
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(TCellEngagerQSP.REPO_ROOT, out_dir_env)
    mkpath(out_dir)

    println("patient_id=$(prob.patient_id)")
    println("n_doses=$(length(prob.dose_times_days)), horizon_days=$(prob.horizon_days)")
    println("start_doses_mg=$(prob.x0_mg)")
    println("methods=$(methods)")
    flush(stdout)

    results = run_methods(prob, methods, out_dir)

    summary_path = joinpath(out_dir, "optimization_methods_summary.csv")
    CSV.write(summary_path, results)

    meta = Dict(
        "patient_id" => prob.patient_id,
        "dose_times_days" => prob.dose_times_days,
        "horizon_days" => prob.horizon_days,
        "tox_window_days" => prob.tox_window_days,
        "tox_tau" => prob.tox_tau,
        "loss_formula" => "L=w1*(mean_dose_window_IL6_peak/s1)+w2*(max_dose_window_IL6_peak/s2)+w3*(IL6_AUC/s3)+w4*(Btumor_end/Btumor0/s4)+w5*(Btumor_AUC_ratio/s5)",
        "weights" => prob.weights,
        "scales" => prob.scales,
        "solver" => string(prob.alg),
        "start_doses_mg" => prob.x0_mg,
        "bounds_mg" => Dict("lower" => prob.lower_mg[1], "upper" => prob.upper_mg[1]),
    )
    meta_path = joinpath(out_dir, "optimization_methods_meta.json")
    open(meta_path, "w") do io
        JSON3.pretty(io, meta)
    end

    println("Wrote " * summary_path)
    println("Wrote " * meta_path)
    println("Best:")
    println(first(results, 1))
end

main()
