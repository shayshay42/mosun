using CSV
using DataFrames
using DifferentialEquations
using FiniteDiff
using ForwardDiff
using JSON3
using LinearAlgebra
using Optim
using Printf
using Random
using ReverseDiff
using SciMLBase
using SciMLSensitivity
using Statistics

include(joinpath(@__DIR__, "optimize_clinical_lhs_cohort_dosing.jl"))

const DTO_PASSTHROUGH = SensitivityADPassThrough()
const OTD_FORWARD = ForwardSensitivity()
const OTD_REVERSE = ReverseDiffAdjoint()

mutable struct TraceLogger
    objective_name::String
    method::String
    t0::Float64
    eval::Int
    best_loss::Float64
    trace_path::String
    progress_path::String
    print_every::Int
    rows::Vector{NamedTuple}
end

Base.@kwdef struct MethodBenchmarkCase
    sample_id::Int
    params::MMC.MosunParams
    pvec::Vector{Float64}
    p_idx::Dict{Symbol, Int}
    u_idx::Dict{Symbol, Int}
    dose_times_days::Vector{Float64}
    decision_labels::Vector{String}
    decision_groups::Vector{Vector{Int}}
    x_clinical::Vector{Float64}
    x0::Vector{Float64}
    lower_mg::Vector{Float64}
    upper_mg::Vector{Float64}
    bw_kg::Float64
    horizon_days::Float64
    saveat::Vector{Float64}
    tox_window_days::Float64
    tox_tau::Float64
    alg
    abstol::Float64
    reltol::Float64
    maxiters::Int
    post_event_proposed_dt
    target_trajectory
    simple_scales
    clinical_scales
    clinical_weights
    objective_name::String
    objective_label::String
    scenario_label::String
    solver_label::String
    random_seed::Int
end

function parse_solver(name::AbstractString)
    lname = lowercase(String(name))
    if lname == "tsit5"
        return Tsit5(), "Tsit5"
    elseif lname == "qndf"
        return QNDF(autodiff = false), "QNDF"
    elseif lname == "rodas4p"
        return Rodas4P(autodiff = false), "Rodas4P"
    else
        error("Unsupported SINGLE_BENCH_SOLVER=$name")
    end
end

function append_row!(path::AbstractString, row::NamedTuple)
    tbl = DataFrame([row])
    write_header = !isfile(path) || filesize(path) == 0
    CSV.write(path, tbl; append = !write_header, writeheader = write_header)
    return nothing
end

function write_progress!(path::AbstractString, payload)
    open(path, "w") do io
        JSON3.pretty(io, payload)
    end
    return nothing
end

function progress_payload(logger::TraceLogger; status::String = "running")
    return Dict(
        "objective" => logger.objective_name,
        "method" => logger.method,
        "status" => status,
        "eval" => logger.eval,
        "elapsed_s" => time() - logger.t0,
        "best_loss" => logger.best_loss,
        "updated_at_unix_s" => time(),
    )
end

function print_progress!(logger::TraceLogger, loss::Float64)
    println(
        @sprintf(
            "[%s] method=%s eval=%d loss=%.6g best=%.6g elapsed=%.1fs",
            logger.objective_name,
            logger.method,
            logger.eval,
            loss,
            logger.best_loss,
            time() - logger.t0,
        ),
    )
    flush(stdout)
    return nothing
end

function make_trace_logger(objective_name::String, method::String, out_dir::AbstractString)
    print_every = max(1, parse(Int, get(ENV, "SINGLE_BENCH_PRINT_EVERY", "5")))
    return TraceLogger(
        objective_name,
        method,
        time(),
        0,
        Inf,
        joinpath(out_dir, "optimization_traces.csv"),
        joinpath(out_dir, "optimization_progress.json"),
        print_every,
        NamedTuple[],
    )
end

function log_eval!(logger::TraceLogger, x::AbstractVector, loss::Real)
    logger.eval += 1
    loss64 = Float64(loss)
    logger.best_loss = min(logger.best_loss, loss64)
    d = Dict{Symbol, Any}(
        :objective => logger.objective_name,
        :method => logger.method,
        :eval => logger.eval,
        :elapsed_s => time() - logger.t0,
        :loss => loss64,
        :best_loss_so_far => logger.best_loss,
    )
    for i in eachindex(x)
        d[Symbol("dose$(i)_mg")] = Float64(x[i])
    end
    row = (; d...)
    push!(logger.rows, row)
    append_row!(logger.trace_path, row)
    write_progress!(logger.progress_path, progress_payload(logger))
    if logger.eval == 1 || logger.eval % logger.print_every == 0
        print_progress!(logger, loss64)
    end
    return nothing
end

function make_logged_objective(base_obj::Function, logger::TraceLogger)
    return function (x::AbstractVector)
        val = base_obj(x)
        if eltype(x) <: Real && !(x[1] isa ForwardDiff.Dual)
            log_eval!(logger, x, val)
        end
        return val
    end
end

function smoothmax_generic(v::AbstractVector, tau::Real)
    isempty(v) && return 0.0
    m = maximum(v)
    return m + tau * log(sum(exp.((v .- m) ./ tau)))
end

function trapz_generic(x::AbstractVector{<:Real}, y::AbstractVector)
    n = length(x)
    n == length(y) || throw(ArgumentError("x/y length mismatch"))
    n <= 1 && return zero(eltype(y))
    s = zero(eltype(y))
    @inbounds for i in 1:(n - 1)
        dt = x[i + 1] - x[i]
        s += (y[i + 1] + y[i]) * (dt / 2)
    end
    return s
end

function expand_decision_doses_generic(decision_doses::AbstractVector, decision_groups::Vector{Vector{Int}}, n_events::Int)
    T = eltype(decision_doses)
    full = zeros(T, n_events)
    for (j, grp) in enumerate(decision_groups)
        for idx in grp
            full[idx] = decision_doses[j]
        end
    end
    return full
end

function regimen_from_decision_doses(case::MethodBenchmarkCase, decision_doses::AbstractVector)
    full = expand_decision_doses_generic(decision_doses, case.decision_groups, length(case.dose_times_days))
    events = map(zip(case.dose_times_days, full)) do (t, dose_mg)
        amt = dose_mg * (1000.0 / case.bw_kg)
        MMC.MosunRegimenEvent(target = :TDBc_ugperkg, time = t, amount = amt, rate = zero(amt))
    end
    return MMC.MosunRegimen(events = collect(events))
end

@inline function tdbc_ugperml(u::AbstractVector, p::AbstractVector, u_idx::Dict{Symbol, Int}, p_idx::Dict{Symbol, Int})
    val = u[u_idx[:TDBc_ugperkg]] / p[p_idx[:Vc_tdb]]
    return ifelse(val > 1e-5, val, zero(val))
end

@inline function il6combo(u::AbstractVector, p::AbstractVector, u_idx::Dict{Symbol, Int}, p_idx::Dict{Symbol, Int})
    return u[u_idx[:IL6pb]] +
        p[p_idx[:IL6_tiss_contribution]] * (
            u[u_idx[:IL6tiss]] * p[p_idx[:Vtissue]] +
            u[u_idx[:IL6tiss2]] * p[p_idx[:Vtissue2]] +
            u[u_idx[:IL6tiss3]] * p[p_idx[:Vtissue3]] +
            u[u_idx[:IL6tumor]] * p[p_idx[:Vtumor]]
        ) / p[p_idx[:Vpb]]
end

function solve_with_callback(case::MethodBenchmarkCase, decision_doses::AbstractVector; sensealg = nothing)
    regimen = regimen_from_decision_doses(case, decision_doses)
    built = MMC.build_problem_vector(
        regimen,
        case.pvec;
        tspan = (0.0, case.horizon_days),
        saveat = case.saveat,
        callback_mode = :callback,
        post_event_proposed_dt = case.post_event_proposed_dt,
    )
    kwargs = (
        abstol = case.abstol,
        reltol = case.reltol,
        callback = built.callback,
        tstops = built.tstops,
        d_discontinuities = built.d_discontinuities,
        saveat = built.saveat,
        save_everystep = false,
        maxiters = case.maxiters,
    )
    sol = isnothing(sensealg) ?
        solve(built.prob, case.alg; kwargs...) :
        solve(built.prob, case.alg; kwargs..., sensealg = sensealg)
    sol.retcode == SciMLBase.ReturnCode.Success || error("solve failed with retcode=$(sol.retcode)")
    return built, sol
end

function trajectory_metrics(case::MethodBenchmarkCase, decision_doses::AbstractVector; sensealg = nothing)
    built, sol = solve_with_callback(case, decision_doses; sensealg = sensealg)
    p = case.pvec
    bt_idx = case.u_idx[:Btumor]
    auc_idx = case.u_idx[:TDBc_ugperml_AUC]
    t = Float64.(sol.t)
    bt = [u[bt_idx] for u in sol.u]
    il6 = [il6combo(u, p, case.u_idx, case.p_idx) for u in sol.u]
    tdbc = [tdbc_ugperml(u, p, case.u_idx, case.p_idx) for u in sol.u]
    bt0 = built.initial_u[bt_idx]
    best_spd = minimum((bt ./ (bt0 + eltype(bt0)(1e-12)) .- 1) .* 100)
    window_peaks = [smoothmax_generic([il6[i] for i in eachindex(t) if t[i] >= td && t[i] <= td + case.tox_window_days + 1e-10], case.tox_tau) for td in case.dose_times_days]
    return (
        t = t,
        bt = bt,
        il6 = il6,
        tdbc = tdbc,
        bt0 = bt0,
        best_spd = best_spd,
        tox_peak_simple = smoothmax_generic(il6, case.tox_tau),
        tox_peak_mean = sum(window_peaks) / length(window_peaks),
        tox_peak_max = maximum(window_peaks),
        tox_auc = trapz_generic(t, max.(il6, zero(eltype(il6[1])))) / case.horizon_days,
        tumor_terminal = bt[end] / (bt0 + eltype(bt0)(1e-12)),
        tumor_auc = trapz_generic(t, bt ./ (bt0 + eltype(bt0)(1e-12))) / case.horizon_days,
        auc_tdbc = sol.u[end][auc_idx],
    )
end

function simple_objective(case::MethodBenchmarkCase, m)
    return 0.5 * (m.tox_peak_simple / case.simple_scales.tox) +
        0.5 * (m.tumor_terminal / case.simple_scales.tumor)
end

function clinical_objective(case::MethodBenchmarkCase, m)
    w = case.clinical_weights
    s = case.clinical_scales
    return w.tox_peak_mean * (m.tox_peak_mean / s.tox_peak_mean) +
        w.tox_peak_max * (m.tox_peak_max / s.tox_peak_max) +
        w.tox_auc * (m.tox_auc / s.tox_auc) +
        w.tumor_terminal * (m.tumor_terminal / s.tumor_terminal) +
        w.tumor_auc * (m.tumor_auc / s.tumor_auc)
end

function tracking_objective(case::MethodBenchmarkCase, m)
    target = case.target_trajectory
    bt_rel = log.(m.bt ./ (m.bt0 + eltype(m.bt0)(1e-12)) .+ eltype(m.bt0)(1e-12))
    tbtc = log.(target.bt ./ (target.bt0 + 1e-12) .+ 1e-12)
    il6_cur = log1p.(max.(m.il6, zero(eltype(m.il6[1]))) ./ target.il6_scale)
    il6_tgt = target.il6_log
    tdbc_cur = log1p.(max.(m.tdbc, zero(eltype(m.tdbc[1]))) ./ target.tdbc_scale)
    tdbc_tgt = target.tdbc_log
    bt_loss = trapz_generic(target.t, (bt_rel .- tbtc) .^ 2) / case.horizon_days
    il6_loss = trapz_generic(target.t, (il6_cur .- il6_tgt) .^ 2) / case.horizon_days
    tdbc_loss = trapz_generic(target.t, (tdbc_cur .- tdbc_tgt) .^ 2) / case.horizon_days
    return 0.50 * bt_loss + 0.30 * il6_loss + 0.20 * tdbc_loss
end

function objective_value(case::MethodBenchmarkCase, x::AbstractVector, sensealg, allow_catch::Bool)
    big_penalty = 1.0e9
    if allow_catch && !(x[1] isa ForwardDiff.Dual)
        xv = Float64.(x)
        if any(!isfinite, xv)
            return big_penalty
        end
        viol = max.(case.lower_mg .- xv, 0.0) .+ max.(xv .- case.upper_mg, 0.0)
        if any(viol .> 0.0)
            return big_penalty + 1.0e6 * sum(abs2, viol)
        end
    end
    if allow_catch
        try
                m = trajectory_metrics(case, x; sensealg = sensealg)
            if case.objective_name == "simple"
                return simple_objective(case, m)
            elseif case.objective_name == "clinical"
                return clinical_objective(case, m)
            elseif case.objective_name == "tracking"
                return tracking_objective(case, m)
            else
                error("Unknown objective $(case.objective_name)")
            end
        catch
            if x[1] isa ForwardDiff.Dual
                return one(x[1]) * big_penalty
            end
            return big_penalty
        end
    end
    m = trajectory_metrics(case, x; sensealg = sensealg)
    if case.objective_name == "simple"
        return simple_objective(case, m)
    elseif case.objective_name == "clinical"
        return clinical_objective(case, m)
    elseif case.objective_name == "tracking"
        return tracking_objective(case, m)
    end
    error("Unknown objective $(case.objective_name)")
end

function make_objective(case::MethodBenchmarkCase; sensealg = nothing, allow_catch::Bool = true)
    return x -> objective_value(case, x, sensealg, allow_catch)
end

function do_random_search(
    obj::Function,
    x0::Vector{Float64},
    lower::Vector{Float64},
    upper::Vector{Float64};
    n_evals::Int,
    local_prob::Float64,
    local_sigma_frac::Float64,
    rng_seed::Int,
)
    rng = MersenneTwister(rng_seed)
    nd = length(x0)
    x_best = clamp.(copy(x0), lower, upper)
    f_best = Float64(obj(x_best))
    for _ in 2:max(n_evals, 1)
        x = similar(x_best)
        if rand(rng) < local_prob
            for j in 1:nd
                span = upper[j] - lower[j]
                x[j] = clamp(x_best[j] + randn(rng) * local_sigma_frac * span, lower[j], upper[j])
            end
        else
            for j in 1:nd
                x[j] = lower[j] + rand(rng) * (upper[j] - lower[j])
            end
        end
        f = Float64(obj(x))
        if f < f_best
            x_best = copy(x)
            f_best = f
        end
    end
    return x_best, f_best
end

function do_metropolis_anneal(
    obj::Function,
    x0::Vector{Float64},
    lower::Vector{Float64},
    upper::Vector{Float64};
    n_steps::Int,
    temp0::Float64,
    tempf::Float64,
    proposal_sigma_frac::Float64,
    rng_seed::Int,
)
    rng = MersenneTwister(rng_seed)
    nd = length(x0)
    x = clamp.(copy(x0), lower, upper)
    f = Float64(obj(x))
    x_best = copy(x)
    f_best = f
    accepted = 0
    denom = max(n_steps - 1, 1)
    for k in 1:max(n_steps, 1)
        α = (k - 1) / denom
        T = exp(log(max(temp0, 1e-12)) * (1 - α) + log(max(tempf, 1e-12)) * α)
        x_prop = similar(x)
        for j in 1:nd
            span = upper[j] - lower[j]
            x_prop[j] = clamp(x[j] + randn(rng) * proposal_sigma_frac * span, lower[j], upper[j])
        end
        f_prop = Float64(obj(x_prop))
        Δ = f_prop - f
        accept = (Δ <= 0.0) || (rand(rng) < exp(-Δ / max(T, 1e-12)))
        if accept
            x = x_prop
            f = f_prop
            accepted += 1
            if f < f_best
                x_best = copy(x)
                f_best = f
            end
        end
    end
    return x_best, f_best, accepted / max(n_steps, 1)
end

function push_summary_row!(rows::Vector{NamedTuple}, case::MethodBenchmarkCase, method::String, x::Vector{Float64}, runtime_s::Float64, status::String)
    m = trajectory_metrics(case, x; sensealg = nothing)
    push!(rows, (
        objective = case.objective_name,
        method = method,
        best_loss = Float64(objective_value(case, x, nothing, true)),
        runtime_s = runtime_s,
        status = status,
        dose1_mg = x[1],
        dose2_mg = x[2],
        dose3_mg = x[3],
        dose4_mg = x[4],
        best_spd_pct = Float64(m.best_spd),
        peak_il6 = Float64(maximum(Float64.(m.il6))),
        auc_tdbc = Float64(m.auc_tdbc),
    ))
end

function run_methods(case::MethodBenchmarkCase, out_dir::AbstractString)
    lower = copy(case.lower_mg)
    upper = copy(case.upper_mg)
    x0 = clamp.(copy(case.x0), lower, upper)
    loss_primal = make_objective(case; sensealg = nothing)
    loss_dto = make_objective(case; sensealg = DTO_PASSTHROUGH)
    loss_dto_reverse = make_objective(case; sensealg = DTO_PASSTHROUGH, allow_catch = false)
    loss_otd_fwd = make_objective(case; sensealg = OTD_FORWARD)
    loss_otd_rev = make_objective(case; sensealg = OTD_REVERSE, allow_catch = false)

    nm_iters = parse(Int, get(ENV, "SINGLE_BENCH_NM_ITERS", "60"))
    lbfgs_iters = parse(Int, get(ENV, "SINGLE_BENCH_LBFGS_ITERS", "25"))
    rs_evals = parse(Int, get(ENV, "SINGLE_BENCH_RS_EVALS", "180"))
    mc_steps = parse(Int, get(ENV, "SINGLE_BENCH_MC_STEPS", "260"))
    filter_env = strip(get(ENV, "SINGLE_BENCH_METHOD_FILTER", ""))
    method_filter = isempty(filter_env) ? nothing : Set(lowercase.(split(filter_env, ",")))

    trace_path = joinpath(out_dir, "optimization_traces.csv")
    summary_path = joinpath(out_dir, "optimization_summary.csv")
    progress_path = joinpath(out_dir, "optimization_progress.json")
    rm(trace_path; force = true)
    rm(summary_path; force = true)
    rm(progress_path; force = true)

    summary_rows = NamedTuple[]
    trace_rows = NamedTuple[]

    function finish_method(method::String, logger::TraceLogger, x_best::Vector{Float64}, status::String)
        append!(trace_rows, logger.rows)
        runtime_s = time() - logger.t0
        pre_len = length(summary_rows)
        push_summary_row!(summary_rows, case, method, x_best, runtime_s, status)
        row = summary_rows[pre_len + 1]
        append_row!(summary_path, row)
        write_progress!(progress_path, progress_payload(logger; status = status))
        println(
            @sprintf(
                "[%s] method=%s finished status=%s best_loss=%.6g runtime=%.1fs",
                case.objective_name,
                method,
                status,
                row.best_loss,
                runtime_s,
            ),
        )
        flush(stdout)
        return nothing
    end

    function wants(method::String)
        isnothing(method_filter) && return true
        return lowercase(method) in method_filter
    end

    # 1) Metropolis annealing
    if wants("metropolis_anneal")
        println("running method=metropolis_anneal objective=$(case.objective_name)")
        flush(stdout)
        logger = make_trace_logger(case.objective_name, "metropolis_anneal", out_dir)
        f_logged = make_logged_objective(loss_primal, logger)
        x_best = copy(x0)
        status = "ok"
        f_logged(x0)
        try
            x_best, _, acc = do_metropolis_anneal(
                f_logged,
                x0,
                lower,
                upper;
                n_steps = mc_steps,
                temp0 = 0.03,
                tempf = 2e-4,
                proposal_sigma_frac = 0.08,
                rng_seed = case.random_seed + 11,
            )
            status = @sprintf("ok_accept=%.3f", acc)
        catch err
            status = "error: " * sprint(showerror, err)
        end
        finish_method("metropolis_anneal", logger, x_best, status)
    end

    # 2) Random search
    if wants("random_search")
        println("running method=random_search objective=$(case.objective_name)")
        flush(stdout)
        logger = make_trace_logger(case.objective_name, "random_search", out_dir)
        f_logged = make_logged_objective(loss_primal, logger)
        x_best = copy(x0)
        status = "ok"
        f_logged(x0)
        try
            x_best, _ = do_random_search(
                f_logged,
                x0,
                lower,
                upper;
                n_evals = rs_evals,
                local_prob = 0.75,
                local_sigma_frac = 0.15,
                rng_seed = case.random_seed + 17,
            )
        catch err
            status = "error: " * sprint(showerror, err)
        end
        finish_method("random_search", logger, x_best, status)
    end

    # 3) Nelder-Mead
    if wants("nelder_mead")
        println("running method=nelder_mead objective=$(case.objective_name)")
        flush(stdout)
        logger = make_trace_logger(case.objective_name, "nelder_mead", out_dir)
        f_logged = make_logged_objective(loss_primal, logger)
        x_best = copy(x0)
        status = "ok"
        f_logged(x0)
        try
            nm_opts = Optim.Options(iterations = nm_iters, show_trace = false, store_trace = false)
            res = optimize(f_logged, lower, upper, copy(x0), Fminbox(NelderMead()), nm_opts)
            x_best = clamp.(Optim.minimizer(res), lower, upper)
            status = string(Optim.converged(res) ? "ok" : "not_converged")
        catch err
            status = "error: " * sprint(showerror, err)
        end
        finish_method("nelder_mead", logger, x_best, status)
    end

    function run_lbfgs(method::String, grad_mode::Symbol, grad_fun)
        wants(method) || return
        println("running method=$(method) objective=$(case.objective_name)")
        flush(stdout)
        logger = make_trace_logger(case.objective_name, method, out_dir)
        f_logged = make_logged_objective(loss_primal, logger)
        x_best = copy(x0)
        status = "ok"
        f_logged(x0)
        try
            g! = (G, x) -> begin
                grad_fun(G, x)
                return G
            end
            opts = Optim.Options(iterations = lbfgs_iters, show_trace = false, store_trace = false)
            res = optimize(f_logged, g!, lower, upper, copy(x0), Fminbox(LBFGS()), opts)
            x_best = clamp.(Optim.minimizer(res), lower, upper)
            status = string(Optim.converged(res) ? "ok" : "not_converged")
        catch err
            status = "error: " * sprint(showerror, err)
        end
        finish_method(method, logger, x_best, status)
    end

    run_lbfgs("finite_diff_lbfgs", :fd, (G, x) -> FiniteDiff.finite_difference_gradient!(G, loss_primal, x))
    run_lbfgs("dto_forward_ad_lbfgs", :dtof, (G, x) -> ForwardDiff.gradient!(G, loss_dto, x))
    run_lbfgs("dto_reverse_ad_lbfgs", :dtor, (G, x) -> copyto!(G, ReverseDiff.gradient(z -> objective_value(case, z, DTO_PASSTHROUGH, false), x)))
    run_lbfgs("otd_forward_ad_lbfgs", :otdf, (G, x) -> ForwardDiff.gradient!(G, loss_otd_fwd, x))
    run_lbfgs("otd_reverse_ad_lbfgs", :otdr, (G, x) -> copyto!(G, ReverseDiff.gradient(z -> objective_value(case, z, OTD_REVERSE, false), x)))

    return DataFrame(summary_rows), DataFrame(trace_rows)
end

function make_case(objective_name::String)
    regimen = load_realdata_regimen()
    samples = load_lhs_cohort(regimen.horizon_days)
    seed = parse(Int, get(ENV, "SINGLE_BENCH_RANDOM_SEED", "20260311"))
    rng = MersenneTwister(seed)
    sample = samples[rand(rng, eachindex(samples))]

    summary = JSON3.read(read(joinpath(REPO_ROOT, "generated", "figures", "optimization", "clinical_lhs_100_cohort_dose_lbfgs_20260310", "optimization_summary.json"), String))
    x0 = Float64.(summary["initial_random_decision_doses_mg"])
    lower = Float64.(summary["lower_mg"])
    upper = Float64.(summary["upper_mg"])
    decision_labels, decision_groups, clinical_decision = build_decision_scheme(regimen.dose_times_days, regimen.dose_mg)

    alg, solver_label = parse_solver(get(ENV, "SINGLE_BENCH_SOLVER", "rodas4p"))
    saveat = build_saveat(regimen.horizon_days, parse(Float64, get(ENV, "SINGLE_BENCH_SAVEAT_DT", "1.0")))
    tox_tau = parse(Float64, get(ENV, "SINGLE_BENCH_TOX_TAU", "50.0"))
    tox_window_days = parse(Float64, get(ENV, "SINGLE_BENCH_TOX_WINDOW_DAYS", "2.0"))
    post_event_proposed_dt = parse(Float64, get(ENV, "SINGLE_BENCH_POST_EVENT_DT", "0.01"))
    bw_kg = parse(Float64, get(ENV, "SINGLE_BENCH_BW_KG", "70.0"))

    pvec = MMC.pack_params(sample.params)
    p_idx = Dict(sym => i for (i, sym) in enumerate(MMC.PARAMETER_NAMES))
    u_idx = Dict(sym => i for (i, sym) in enumerate(MMC.DYNAMIC_STATE_NAMES))

    case_stub = MethodBenchmarkCase(
        sample_id = sample.sample_id,
        params = sample.params,
        pvec = pvec,
        p_idx = p_idx,
        u_idx = u_idx,
        dose_times_days = regimen.dose_times_days,
        decision_labels = decision_labels,
        decision_groups = decision_groups,
        x_clinical = Float64.(clinical_decision),
        x0 = x0,
        lower_mg = lower,
        upper_mg = upper,
        bw_kg = bw_kg,
        horizon_days = regimen.horizon_days,
        saveat = saveat,
        tox_window_days = tox_window_days,
        tox_tau = tox_tau,
        alg = alg,
        abstol = TCellEngagerQSP.SOLVER_ABSTOL,
        reltol = TCellEngagerQSP.SOLVER_RELTOL,
        maxiters = 1_000_000,
        post_event_proposed_dt = post_event_proposed_dt,
        target_trajectory = nothing,
        simple_scales = nothing,
        clinical_scales = nothing,
        clinical_weights = normalize_weights(LossWeights(0.2, 0.25, 0.15, 0.2, 0.2)),
        objective_name = objective_name,
        objective_label = objective_name,
        scenario_label = "DLBCL-like virtual patient under 8-cycle mosunetuzumab step-up dosing",
        solver_label = solver_label,
        random_seed = seed,
    )

    m0 = trajectory_metrics(case_stub, x0; sensealg = nothing)
    mclin = trajectory_metrics(case_stub, clinical_decision; sensealg = nothing)
    simple_scales = (
        tox = max(abs(Float64(m0.tox_peak_simple)), 1e-6),
        tumor = max(abs(Float64(m0.tumor_terminal)), 1e-8),
    )
    clinical_scales = LossScales(
        max(abs(Float64(m0.tox_peak_mean)), 1e-6),
        max(abs(Float64(m0.tox_peak_max)), 1e-6),
        max(abs(Float64(m0.tox_auc)), 1e-6),
        max(abs(Float64(m0.tumor_terminal)), 1e-8),
        max(abs(Float64(m0.tumor_auc)), 1e-8),
    )
    target_trajectory = (
        t = mclin.t,
        bt = Float64.(mclin.bt),
        bt0 = Float64(mclin.bt0),
        il6_log = log1p.(max.(Float64.(mclin.il6), 0.0) ./ max(maximum(Float64.(mclin.il6)), 1.0)),
        tdbc_log = log1p.(max.(Float64.(mclin.tdbc), 0.0) ./ max(maximum(Float64.(mclin.tdbc)), 1e-6)),
        il6_scale = max(maximum(Float64.(mclin.il6)), 1.0),
        tdbc_scale = max(maximum(Float64.(mclin.tdbc)), 1e-6),
    )
    labels = Dict(
        "simple" => "Simple toxicity + tumor proxy",
        "clinical" => "Clinical 5-term proxy loss",
        "tracking" => "Quadratic tracking loss to nominal clinical trajectory",
    )

    return MethodBenchmarkCase(
        sample_id = sample.sample_id,
        params = sample.params,
        pvec = pvec,
        p_idx = p_idx,
        u_idx = u_idx,
        dose_times_days = regimen.dose_times_days,
        decision_labels = decision_labels,
        decision_groups = decision_groups,
        x_clinical = Float64.(clinical_decision),
        x0 = x0,
        lower_mg = lower,
        upper_mg = upper,
        bw_kg = bw_kg,
        horizon_days = regimen.horizon_days,
        saveat = saveat,
        tox_window_days = tox_window_days,
        tox_tau = tox_tau,
        alg = alg,
        abstol = TCellEngagerQSP.SOLVER_ABSTOL,
        reltol = TCellEngagerQSP.SOLVER_RELTOL,
        maxiters = 1_000_000,
        post_event_proposed_dt = post_event_proposed_dt,
        target_trajectory = target_trajectory,
        simple_scales = simple_scales,
        clinical_scales = clinical_scales,
        clinical_weights = normalize_weights(LossWeights(0.2, 0.25, 0.15, 0.2, 0.2)),
        objective_name = objective_name,
        objective_label = labels[objective_name],
        scenario_label = "DLBCL-like virtual patient under 8-cycle mosunetuzumab step-up dosing",
        solver_label = solver_label,
        random_seed = seed,
    )
end

function main()
    objective_name = lowercase(get(ENV, "SINGLE_BENCH_OBJECTIVE", "simple"))
    objective_name in ("simple", "clinical", "tracking") || error("Unsupported SINGLE_BENCH_OBJECTIVE=$objective_name")
    out_root_env = get(ENV, "SINGLE_BENCH_OUT_ROOT", joinpath("generated", "figures", "optimization", "single_sample_dose_method_benchmark_20260311"))
    out_root = isabspath(out_root_env) ? out_root_env : joinpath(REPO_ROOT, out_root_env)
    out_dir = joinpath(out_root, objective_name)
    mkpath(out_dir)

    case = make_case(objective_name)
    println(@sprintf("benchmarking objective=%s sample_id=%d solver=%s", case.objective_name, case.sample_id, case.solver_label))
    results, traces = run_methods(case, out_dir)

    CSV.write(joinpath(out_dir, "optimization_summary.csv"), results)
    CSV.write(joinpath(out_dir, "optimization_traces.csv"), traces)

    meta = Dict(
        "objective_name" => case.objective_name,
        "objective_label" => case.objective_label,
        "sample_id" => case.sample_id,
        "decision_labels" => case.decision_labels,
        "x0_mg" => case.x0,
        "x_clinical_mg" => case.x_clinical,
        "dose_times_days" => case.dose_times_days,
        "decision_groups" => case.decision_groups,
        "horizon_days" => case.horizon_days,
        "saveat_dt_days" => case.saveat[2] - case.saveat[1],
        "solver" => case.solver_label,
        "callback_mode" => "PresetTimeCallback",
        "callbacks_ad_compatible" => true,
        "post_event_proposed_dt_days" => case.post_event_proposed_dt,
        "scenario_label" => case.scenario_label,
        "clinical_weights" => Dict(
            "tox_peak_mean" => case.clinical_weights.tox_peak_mean,
            "tox_peak_max" => case.clinical_weights.tox_peak_max,
            "tox_auc" => case.clinical_weights.tox_auc,
            "tumor_terminal" => case.clinical_weights.tumor_terminal,
            "tumor_auc" => case.clinical_weights.tumor_auc,
        ),
        "objective_definition" => case.objective_label,
    )
    open(joinpath(out_dir, "optimization_meta.json"), "w") do io
        JSON3.pretty(io, meta)
    end

    println("Wrote ", joinpath(out_dir, "optimization_summary.csv"))
    println("Wrote ", joinpath(out_dir, "optimization_traces.csv"))
    println("Wrote ", joinpath(out_dir, "optimization_meta.json"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
