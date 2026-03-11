using CSV
using DataFrames
using DifferentialEquations
using FiniteDiff
using JSON3
using LinearAlgebra
using Optim
using Printf
using Random
using SciMLBase
using Statistics
using Sundials
using Base.Threads

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

const MMC = TCellEngagerQSP.MosunModelCore
const REPO_ROOT = TCellEngagerQSP.REPO_ROOT

BLAS.set_num_threads(1)

struct LossScales
    tox_peak_mean::Float64
    tox_peak_max::Float64
    tox_auc::Float64
    tumor_terminal::Float64
    tumor_auc::Float64
end

struct LossWeights
    tox_peak_mean::Float64
    tox_peak_max::Float64
    tox_auc::Float64
    tumor_terminal::Float64
    tumor_auc::Float64
end

struct CohortSample
    sample_id::Int
    params::MMC.MosunParams
end

Base.@kwdef struct CohortOptimizationProblem
    cohort::Vector{CohortSample}
    dose_times_days::Vector{Float64}
    clinical_doses_mg::Vector{Float64}
    decision_labels::Vector{String}
    decision_groups::Vector{Vector{Int}}
    clinical_decision_doses_mg::Vector{Float64}
    lower_mg::Vector{Float64}
    upper_mg::Vector{Float64}
    horizon_days::Float64
    bw_kg::Float64
    objective_dt::Float64
    objective_saveat::Vector{Float64}
    sample_dt::Float64
    saveat::Vector{Float64}
    tox_window_days::Float64
    tox_tau::Float64
    alg
    abstol::Float64
    reltol::Float64
    maxiters::Int
    post_event_proposed_dt
    out_dir::String
    seed::Int
end

mutable struct LossEvaluator
    prob::CohortOptimizationProblem
    scales::LossScales
    weights::LossWeights
    memo_x::Vector{Float64}
    memo_result
    has_memo::Bool
end

function parse_float_list(txt::AbstractString)
    vals = Float64[]
    for tok in split(String(txt), ",")
        s = strip(tok)
        isempty(s) && continue
        push!(vals, parse(Float64, s))
    end
    return vals
end

function smoothmax(v::AbstractVector, tau::Real)
    isempty(v) && return 0.0
    m = maximum(v)
    return m + tau * log(sum(exp.((v .- m) ./ tau)))
end

function trapz(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    n == length(y) || throw(ArgumentError("x/y length mismatch"))
    n <= 1 && return 0.0
    acc = 0.0
    @inbounds for i in 1:(n - 1)
        dx = Float64(x[i + 1] - x[i])
        acc += 0.5 * dx * (Float64(y[i]) + Float64(y[i + 1]))
    end
    return acc
end

function robust_scale(x::Vector{Float64}, floor_val::Float64, fallback::Float64)
    vals = filter(v -> isfinite(v) && v > 0.0, x)
    isempty(vals) && return max(floor_val, fallback)
    return max(floor_val, median(vals))
end

function normalize_weights(w::LossWeights)
    s = w.tox_peak_mean + w.tox_peak_max + w.tox_auc + w.tumor_terminal + w.tumor_auc
    s > 0.0 || error("Loss weights must sum to > 0")
    return LossWeights(
        w.tox_peak_mean / s,
        w.tox_peak_max / s,
        w.tox_auc / s,
        w.tumor_terminal / s,
        w.tumor_auc / s,
    )
end

function make_solver_alg(name::AbstractString)
    lname = lowercase(String(name))
    if lname == "qndf"
        return QNDF(autodiff = false)
    elseif lname == "rodas4p"
        return Rodas4P(autodiff = false)
    elseif lname == "trbdf2"
        return TRBDF2(autodiff = false)
    elseif lname == "cvode_bdf"
        return Sundials.CVODE_BDF()
    else
        error("Unsupported COHORT_OPT_SOLVER=$name")
    end
end

function load_realdata_regimen()
    ref = DataFrame(CSV.File(joinpath(REPO_ROOT, "generated", "reference", "musun_2022_88patients_SPD_from_vpop_generation_3a50cd1.csv")))
    row = ref[1, :]
    return (
        dose_times_days = parse_float_list(String(row.dose_days)),
        dose_mg = parse_float_list(String(row.dose_mg)),
        horizon_days = Float64(row.horizon_days),
        regimen_label = String(row.regimen_label),
    )
end

function load_base_params()
    p = MMC.default_params()
    df = DataFrame(CSV.File(joinpath(REPO_ROOT, "generated", "phase1_design_standard", "dlbcl_param_overrides.csv")))
    for r in eachrow(df)
        MMC.has_parameter(r.name) || continue
        MMC.set_param!(p, r.name, Float64(r.value))
    end
    return p
end

function load_lhs_cohort(horizon_days::Float64)
    base = load_base_params()
    df = DataFrame(CSV.File(joinpath(REPO_ROOT, "generated", "figures", "clinical_lhs_regimen_20260310", "lhs_samples_metrics.csv")))
    ok = df[df.status .== "ok", :]
    max_samples = parse(Int, get(ENV, "COHORT_OPT_MAX_SAMPLES", string(nrow(ok))))
    ok = ok[1:min(max_samples, nrow(ok)), :]
    samples = CohortSample[]
    for r in eachrow(ok)
        p = deepcopy(base)
        for cname in names(ok)
            MMC.has_parameter(cname) || continue
            v = try
                Float64(r[Symbol(cname)])
            catch
                NaN
            end
            isfinite(v) || continue
            MMC.set_param!(p, cname, v)
        end
        MMC.set_param!(p, :PKflag, 1.0)
        MMC.set_param!(p, :VPid, 1.0)
        MMC.set_param!(p, :fvalidation, 0.0)
        MMC.set_param!(p, :end_time, horizon_days)
        push!(samples, CohortSample(Int(r.sample_id), p))
    end
    isempty(samples) && error("No successful cohort samples found")
    return samples
end

function dose_vector_to_regimen(dose_times_days::Vector{Float64}, doses_mg::AbstractVector{<:Real}, bw_kg::Float64)
    length(dose_times_days) == length(doses_mg) || error("dose time / amount length mismatch")
    dose_map = Dict{Float64, Float64}()
    for (t, dose_mg) in zip(dose_times_days, doses_mg)
        amt = Float64(dose_mg) * 1000.0 / bw_kg
        dose_map[Float64(t)] = get(dose_map, Float64(t), 0.0) + amt
    end
    return MMC.bolus_regimen(:TDBc_ugperkg, dose_map)
end

function build_decision_scheme(dose_times_days::Vector{Float64}, clinical_doses_mg::Vector{Float64})
    mode = lowercase(get(ENV, "COHORT_OPT_DECISION_MODE", "clinical_stepup_q3w"))
    if mode == "clinical_stepup_q3w"
        length(dose_times_days) == 10 || error("clinical_stepup_q3w decision mode expects 10 dose events")
        groups = [Int[1], Int[2], Int[3, 4], collect(5:10)]
        labels = ["C1D1", "C1D8", "C1D15_C2D1", "C3plus_q3w"]
    elseif mode == "full10"
        groups = [Int[i] for i in eachindex(dose_times_days)]
        labels = ["dose$(i)" for i in eachindex(dose_times_days)]
    else
        error("Unsupported COHORT_OPT_DECISION_MODE=$mode")
    end
    decision_vals = Float64[]
    for grp in groups
        push!(decision_vals, mean(clinical_doses_mg[grp]))
    end
    return labels, groups, decision_vals
end

function expand_decision_doses(decision_doses::AbstractVector{<:Real}, decision_groups::Vector{Vector{Int}}, n_events::Int)
    length(decision_doses) == length(decision_groups) || error("Decision dose length mismatch")
    full = zeros(Float64, n_events)
    for (j, grp) in enumerate(decision_groups)
        for idx in grp
            full[idx] = Float64(decision_doses[j])
        end
    end
    return full
end

function build_saveat(horizon_days::Float64, dt::Float64)
    vals = collect(0.0:dt:horizon_days)
    if isempty(vals) || vals[end] < horizon_days - 1e-12
        push!(vals, horizon_days)
    else
        vals[end] = horizon_days
    end
    return vals
end

function collect_window_vals(t::Vector{Float64}, y::Vector{Float64}, t0::Float64, t1::Float64)
    vals = Float64[]
    @inbounds for i in eachindex(t)
        if t[i] >= t0 - 1e-10 && t[i] <= t1 + 1e-10
            push!(vals, y[i])
        end
    end
    if isempty(vals)
        idx = searchsortedfirst(t, t0)
        if idx <= 1
            push!(vals, y[1])
        elseif idx > length(t)
            push!(vals, y[end])
        else
            push!(vals, y[idx])
        end
    end
    return vals
end

function solution_series(sol, p::MMC.MosunParams, cache::MMC.MosunObservablesCache)
    bt_idx = MMC.dynamic_state_index(:Btumor)
    n = length(sol.t)
    t = Float64.(sol.t)
    bt = Vector{Float64}(undef, n)
    il6 = Vector{Float64}(undef, n)
    tdbc = Vector{Float64}(undef, n)
    @inbounds for j in eachindex(sol.t)
        uj = sol.u[j]
        bt[j] = Float64(uj[bt_idx])
        MMC.update_observables!(cache, uj, p, sol.t[j])
        il6[j] = cache.IL6combo
        tdbc[j] = cache.TDBc_ugperml
    end
    return t, bt, il6, tdbc
end

function raw_metrics_from_series(prob::CohortOptimizationProblem, t::Vector{Float64}, bt::Vector{Float64}, il6::Vector{Float64})
    tox_peaks = Float64[]
    for tdose in prob.dose_times_days
        vals = collect_window_vals(t, il6, tdose, tdose + prob.tox_window_days)
        push!(tox_peaks, smoothmax(vals, prob.tox_tau))
    end
    bt0 = bt[1]
    il6_auc_abs = trapz(t, max.(il6, 0.0))
    tumor_auc_abs = trapz(t, max.(bt, 0.0))
    tox_peak_mean = mean(tox_peaks)
    tox_peak_max = maximum(tox_peaks)
    tox_auc = il6_auc_abs / prob.horizon_days
    tumor_terminal = bt[end] / (bt0 + 1e-12)
    tumor_auc = tumor_auc_abs / (prob.horizon_days * (bt0 + 1e-12))
    best_spd_pct = 100.0 * (minimum(bt) / (bt0 + 1e-12) - 1.0)
    peak_il6 = maximum(il6)
    _, peak_idx = findmax(il6)
    day_of_global_peak_il6 = t[peak_idx]
    return (
        tox_peak_mean = tox_peak_mean,
        tox_peak_max = tox_peak_max,
        tox_auc = tox_auc,
        il6_auc_abs = il6_auc_abs,
        tumor_terminal = tumor_terminal,
        tumor_auc = tumor_auc,
        tumor_auc_abs = tumor_auc_abs,
        best_spd_pct = best_spd_pct,
        peak_il6 = peak_il6,
        day_of_global_peak_il6 = day_of_global_peak_il6,
    )
end

function scaled_loss(raw, scales::LossScales, weights::LossWeights)
    tox_peak_mean_s = raw.tox_peak_mean / scales.tox_peak_mean
    tox_peak_max_s = raw.tox_peak_max / scales.tox_peak_max
    tox_auc_s = raw.tox_auc / scales.tox_auc
    tumor_terminal_s = raw.tumor_terminal / scales.tumor_terminal
    tumor_auc_s = raw.tumor_auc / scales.tumor_auc
    loss = weights.tox_peak_mean * tox_peak_mean_s +
           weights.tox_peak_max * tox_peak_max_s +
           weights.tox_auc * tox_auc_s +
           weights.tumor_terminal * tumor_terminal_s +
           weights.tumor_auc * tumor_auc_s
    return (
        tox_peak_mean_s = tox_peak_mean_s,
        tox_peak_max_s = tox_peak_max_s,
        tox_auc_s = tox_auc_s,
        tumor_terminal_s = tumor_terminal_s,
        tumor_auc_s = tumor_auc_s,
        loss = loss,
    )
end

function same_dose_vector(a::Vector{Float64}, b::AbstractVector)
    length(a) == length(b) || return false
    @inbounds for i in eachindex(a)
        a[i] == Float64(b[i]) || return false
    end
    return true
end

function evaluate_cohort(prob::CohortOptimizationProblem, doses_mg::AbstractVector{<:Real}, scales::LossScales, weights::LossWeights)
    length(doses_mg) == length(prob.decision_groups) || error("Decision dose vector length mismatch")
    xv = Float64.(doses_mg)
    if any(!isfinite, xv)
        return (ok = false, penalty = true, raw = nothing, scaled = (; loss = 1.0e12), n_failed = length(prob.cohort))
    end
    viol = max.(prob.lower_mg .- xv, 0.0) .+ max.(xv .- prob.upper_mg, 0.0)
    if any(viol .> 0.0)
        penalty = 1.0e9 + 1.0e6 * sum(abs2, viol)
        return (ok = false, penalty = true, raw = nothing, scaled = (; loss = penalty), n_failed = length(prob.cohort))
    end

    full_doses = expand_decision_doses(xv, prob.decision_groups, length(prob.dose_times_days))
    regimen = dose_vector_to_regimen(prob.dose_times_days, full_doses, prob.bw_kg)
    n = length(prob.cohort)
    tox_peak_mean = fill(NaN, n)
    tox_peak_max = fill(NaN, n)
    tox_auc = fill(NaN, n)
    tumor_terminal = fill(NaN, n)
    tumor_auc = fill(NaN, n)
    best_spd = fill(NaN, n)
    peak_il6 = fill(NaN, n)
    ok_mask = falses(n)
    caches = [MMC.zero_observables_cache() for _ in 1:Threads.maxthreadid()]

    @threads for i in 1:n
        sample = prob.cohort[i]
        cache = caches[threadid()]
        try
            built = MMC.build_problem(
                regimen,
                sample.params;
                tspan = (0.0, prob.horizon_days),
                saveat = prob.objective_saveat,
                callback_mode = :callback,
                post_event_proposed_dt = prob.post_event_proposed_dt,
            )
            sol = MMC.solve_problem(
                built,
                prob.alg;
                abstol = prob.abstol,
                reltol = prob.reltol,
                maxiters = prob.maxiters,
            )
            sol.retcode == SciMLBase.ReturnCode.Success || error("retcode=$(sol.retcode)")
            t, bt, il6, _ = solution_series(sol, sample.params, cache)
            raw = raw_metrics_from_series(prob, t, bt, il6)
            tox_peak_mean[i] = raw.tox_peak_mean
            tox_peak_max[i] = raw.tox_peak_max
            tox_auc[i] = raw.tox_auc
            tumor_terminal[i] = raw.tumor_terminal
            tumor_auc[i] = raw.tumor_auc
            best_spd[i] = raw.best_spd_pct
            peak_il6[i] = raw.peak_il6
            ok_mask[i] = true
        catch
            ok_mask[i] = false
        end
    end

    n_failed = count(!, ok_mask)
    if n_failed > 0
        penalty = 1.0e9 + 1.0e7 * n_failed
        return (ok = false, penalty = true, raw = nothing, scaled = (; loss = penalty), n_failed = n_failed)
    end

    raw = (
        tox_peak_mean = mean(tox_peak_mean),
        tox_peak_max = mean(tox_peak_max),
        tox_auc = mean(tox_auc),
        tumor_terminal = mean(tumor_terminal),
        tumor_auc = mean(tumor_auc),
        best_spd_mean = mean(best_spd),
        peak_il6_mean = mean(peak_il6),
    )
    scaled = scaled_loss(raw, scales, weights)
    return (ok = true, penalty = false, raw = raw, scaled = scaled, n_failed = 0)
end

function evaluate!(ev::LossEvaluator, x::AbstractVector{<:Real})
    if ev.has_memo && same_dose_vector(ev.memo_x, x)
        return ev.memo_result
    end
    result = evaluate_cohort(ev.prob, x, ev.scales, ev.weights)
    ev.memo_x = Float64.(x)
    ev.memo_result = result
    ev.has_memo = true
    return result
end

function calibrate_scales(prob::CohortOptimizationProblem, weights::LossWeights, x0::Vector{Float64}; n_samples::Int, seed::Int)
    rng = MersenneTwister(seed)
    tox_peak_mean = Float64[]
    tox_peak_max = Float64[]
    tox_auc = Float64[]
    tumor_terminal = Float64[]
    tumor_auc = Float64[]

    candidates = Vector{Vector{Float64}}()
    push!(candidates, copy(x0))
    for _ in 1:n_samples
        x = similar(x0)
        @inbounds for j in eachindex(x)
            x[j] = prob.lower_mg[j] + rand(rng) * (prob.upper_mg[j] - prob.lower_mg[j])
        end
        push!(candidates, x)
    end

    scales_seed = LossScales(1.0, 1.0, 1.0, 1.0, 1.0)
    for x in candidates
        result = evaluate_cohort(prob, x, scales_seed, weights)
        result.ok || continue
        push!(tox_peak_mean, result.raw.tox_peak_mean)
        push!(tox_peak_max, result.raw.tox_peak_max)
        push!(tox_auc, result.raw.tox_auc)
        push!(tumor_terminal, result.raw.tumor_terminal)
        push!(tumor_auc, result.raw.tumor_auc)
    end

    isempty(tox_peak_mean) && error("Scale calibration failed; no valid candidate regimens")
    return LossScales(
        robust_scale(tox_peak_mean, 1e-6, tox_peak_mean[1]),
        robust_scale(tox_peak_max, 1e-6, tox_peak_max[1]),
        robust_scale(tox_auc, 1e-6, tox_auc[1]),
        robust_scale(tumor_terminal, 1e-8, tumor_terminal[1]),
        robust_scale(tumor_auc, 1e-8, tumor_auc[1]),
    )
end

function format_dose_label(doses::AbstractVector{<:Real})
    chunks = String[]
    for i in 1:5:length(doses)
        j = min(i + 4, length(doses))
        push!(chunks, join([@sprintf("%.2f", doses[k]) for k in i:j], ", "))
    end
    return join(chunks, " | ")
end

function simulate_cohort_trajectories(prob::CohortOptimizationProblem, decision_doses::Vector{Float64})
    full_doses = expand_decision_doses(decision_doses, prob.decision_groups, length(prob.dose_times_days))
    regimen = dose_vector_to_regimen(prob.dose_times_days, full_doses, prob.bw_kg)
    n = length(prob.cohort)
    sample_rows = Vector{NamedTuple}(undef, n)
    traj_rows = [NamedTuple[] for _ in 1:Threads.maxthreadid()]
    caches = [MMC.zero_observables_cache() for _ in 1:Threads.maxthreadid()]

    @threads for i in 1:n
        sample = prob.cohort[i]
        cache = caches[threadid()]
        built = MMC.build_problem(
            regimen,
            sample.params;
            tspan = (0.0, prob.horizon_days),
            saveat = prob.saveat,
            callback_mode = :callback,
            post_event_proposed_dt = prob.post_event_proposed_dt,
        )
        sol = MMC.solve_problem(
            built,
            prob.alg;
            abstol = prob.abstol,
            reltol = prob.reltol,
            maxiters = prob.maxiters,
        )
        sol.retcode == SciMLBase.ReturnCode.Success || error("Trajectory solve failed for sample $(sample.sample_id): retcode=$(sol.retcode)")

        t, bt, il6, tdbc = solution_series(sol, sample.params, cache)
        raw = raw_metrics_from_series(prob, t, bt, il6)
        auc_tdbc = trapz(t, max.(tdbc, 0.0))
        sample_rows[i] = (
            sample_id = sample.sample_id,
            best_spd_pct = raw.best_spd_pct,
            peak_il6 = raw.peak_il6,
            day_of_global_peak_il6 = raw.day_of_global_peak_il6,
            auc_tdbc = auc_tdbc,
            tox_peak_mean = raw.tox_peak_mean,
            tox_peak_max = raw.tox_peak_max,
            tox_auc = raw.tox_auc,
            il6_auc_abs = raw.il6_auc_abs,
            tumor_terminal = raw.tumor_terminal,
            tumor_auc = raw.tumor_auc,
            tumor_auc_abs = raw.tumor_auc_abs,
        )
        local_rows = traj_rows[threadid()]
        for j in eachindex(t)
            push!(local_rows, (
                sample_id = sample.sample_id,
                time_day = t[j],
                Btumor = bt[j],
                IL6combo = il6[j],
                TDBc_ugperml = tdbc[j],
            ))
        end
    end

    return DataFrame(sample_rows), DataFrame(vcat(traj_rows...))
end

function to_float_dict(nt)
    out = Dict{String, Float64}()
    for (k, v) in pairs(nt)
        if v isa Real
            out[String(k)] = Float64(v)
        end
    end
    return out
end

function build_problem()
    regimen = load_realdata_regimen()
    cohort = load_lhs_cohort(regimen.horizon_days)
    bw_kg = parse(Float64, get(ENV, "COHORT_OPT_BW_KG", "70.0"))
    lower_scalar = parse(Float64, get(ENV, "COHORT_OPT_DOSE_LOWER_MG", "0.0"))
    upper_scalar = parse(Float64, get(ENV, "COHORT_OPT_DOSE_UPPER_MG", "60.0"))
    decision_labels, decision_groups, clinical_decision_doses_mg = build_decision_scheme(regimen.dose_times_days, regimen.dose_mg)
    lower_mg = fill(lower_scalar, length(decision_groups))
    upper_mg = fill(upper_scalar, length(decision_groups))
    sample_dt = parse(Float64, get(ENV, "COHORT_OPT_SAMPLE_DT", "0.25"))
    objective_dt = parse(Float64, get(ENV, "COHORT_OPT_OBJECTIVE_DT", "1.0"))
    tox_window_days = parse(Float64, get(ENV, "COHORT_OPT_TOX_WINDOW_DAYS", "2.0"))
    tox_tau = parse(Float64, get(ENV, "COHORT_OPT_TOX_SOFTMAX_TAU", "50.0"))
    solver_name = get(ENV, "COHORT_OPT_SOLVER", "qndf")
    alg = make_solver_alg(solver_name)
    abstol = parse(Float64, get(ENV, "COHORT_OPT_ABSTOL", "1e-8"))
    reltol = parse(Float64, get(ENV, "COHORT_OPT_RELTOL", "1e-6"))
    maxiters = parse(Int, get(ENV, "COHORT_OPT_MAXITERS_SOLVE", "1000000"))
    post_event_dt_txt = strip(get(ENV, "COHORT_OPT_POST_EVENT_DT", "1e-3"))
    post_event_proposed_dt = isempty(post_event_dt_txt) ? nothing : parse(Float64, post_event_dt_txt)
    out_dir_default = joinpath(REPO_ROOT, "generated", "figures", "optimization", "clinical_lhs_100_cohort_dose_lbfgs_20260310")
    out_dir_env = get(ENV, "COHORT_OPT_OUT_DIR", out_dir_default)
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(REPO_ROOT, out_dir_env)
    seed = parse(Int, get(ENV, "COHORT_OPT_SEED", "20260310"))
    prob = CohortOptimizationProblem(
        cohort = cohort,
        dose_times_days = regimen.dose_times_days,
        clinical_doses_mg = regimen.dose_mg,
        decision_labels = decision_labels,
        decision_groups = decision_groups,
        clinical_decision_doses_mg = clinical_decision_doses_mg,
        lower_mg = lower_mg,
        upper_mg = upper_mg,
        horizon_days = regimen.horizon_days,
        bw_kg = bw_kg,
        objective_dt = objective_dt,
        objective_saveat = build_saveat(regimen.horizon_days, objective_dt),
        sample_dt = sample_dt,
        saveat = build_saveat(regimen.horizon_days, sample_dt),
        tox_window_days = tox_window_days,
        tox_tau = tox_tau,
        alg = alg,
        abstol = abstol,
        reltol = reltol,
        maxiters = maxiters,
        post_event_proposed_dt = post_event_proposed_dt,
        out_dir = out_dir,
        seed = seed,
    )
    return prob, regimen.regimen_label
end

function main()
    prob, regimen_label = build_problem()
    mkpath(prob.out_dir)
    println(@sprintf("cohort_n=%d threads=%d solver=%s", length(prob.cohort), nthreads(), string(prob.alg)))

    weights = normalize_weights(
        LossWeights(
            parse(Float64, get(ENV, "COHORT_OPT_W_TOX_PEAK_MEAN", "0.20")),
            parse(Float64, get(ENV, "COHORT_OPT_W_TOX_PEAK_MAX", "0.25")),
            parse(Float64, get(ENV, "COHORT_OPT_W_TOX_AUC", "0.15")),
            parse(Float64, get(ENV, "COHORT_OPT_W_TUMOR_TERMINAL", "0.20")),
            parse(Float64, get(ENV, "COHORT_OPT_W_TUMOR_AUC", "0.20")),
        ),
    )

    rng = MersenneTwister(prob.seed)
    x0 = similar(prob.clinical_decision_doses_mg)
    @inbounds for j in eachindex(x0)
        x0[j] = prob.lower_mg[j] + rand(rng) * (prob.upper_mg[j] - prob.lower_mg[j])
    end

    n_scale_samples = parse(Int, get(ENV, "COHORT_OPT_SCALE_SAMPLES", "12"))
    scales = calibrate_scales(prob, weights, x0; n_samples = n_scale_samples, seed = prob.seed + 17)
    println(
        @sprintf(
            "scales: tox_peak_mean=%.4g tox_peak_max=%.4g tox_auc=%.4g tumor_terminal=%.4g tumor_auc=%.4g",
            scales.tox_peak_mean,
            scales.tox_peak_max,
            scales.tox_auc,
            scales.tumor_terminal,
            scales.tumor_auc,
        ),
    )

    evaluator = LossEvaluator(prob, scales, weights, Float64[], nothing, false)
    objective = x -> begin
        result = evaluate!(evaluator, x)
        return result.ok ? Float64(result.scaled.loss) : Float64(result.scaled.loss)
    end
    gradient! = (G, x) -> begin
        FiniteDiff.finite_difference_gradient!(G, objective, x)
        return G
    end

    trace_rows = NamedTuple[]
    function push_trace!(tag::String, iter::Int, x::Vector{Float64})
        result = evaluate!(evaluator, x)
        if result.ok
            row = (
                tag = tag,
                iteration = iter,
                objective = Float64(result.scaled.loss),
                tox_peak_mean = Float64(result.raw.tox_peak_mean),
                tox_peak_max = Float64(result.raw.tox_peak_max),
                tox_auc = Float64(result.raw.tox_auc),
                tumor_terminal = Float64(result.raw.tumor_terminal),
                tumor_auc = Float64(result.raw.tumor_auc),
                tox_peak_mean_s = Float64(result.scaled.tox_peak_mean_s),
                tox_peak_max_s = Float64(result.scaled.tox_peak_max_s),
                tox_auc_s = Float64(result.scaled.tox_auc_s),
                tumor_terminal_s = Float64(result.scaled.tumor_terminal_s),
                tumor_auc_s = Float64(result.scaled.tumor_auc_s),
                n_failed = result.n_failed,
            )
        else
            row = (
                tag = tag,
                iteration = iter,
                objective = Float64(result.scaled.loss),
                tox_peak_mean = NaN,
                tox_peak_max = NaN,
                tox_auc = NaN,
                tumor_terminal = NaN,
                tumor_auc = NaN,
                tox_peak_mean_s = NaN,
                tox_peak_max_s = NaN,
                tox_auc_s = NaN,
                tumor_terminal_s = NaN,
                tumor_auc_s = NaN,
                n_failed = result.n_failed,
            )
        end
        d = Dict{Symbol, Any}(pairs(row))
        for i in eachindex(x)
            d[Symbol("decision$(i)_mg")] = Float64(x[i])
        end
        full = expand_decision_doses(x, prob.decision_groups, length(prob.dose_times_days))
        for i in eachindex(full)
            d[Symbol("dose$(i)_mg")] = Float64(full[i])
        end
        push!(trace_rows, (; d...))
    end

    push_trace!("initial", 0, copy(x0))

    iter_ref = Ref(0)
    callback = state -> begin
        iter_ref[] += 1
        push_trace!("lbfgs", iter_ref[], Float64.(state.x))
        return false
    end

    lbfgs_iters = parse(Int, get(ENV, "COHORT_OPT_LBFGS_ITERS", "35"))
    f_calls_limit = parse(Int, get(ENV, "COHORT_OPT_LBFGS_F_CALLS", "0"))
    g_tol = parse(Float64, get(ENV, "COHORT_OPT_LBFGS_GTOL", "1e-6"))
    f_tol = parse(Float64, get(ENV, "COHORT_OPT_LBFGS_FTOL", "1e-9"))
    opts = Optim.Options(
        iterations = lbfgs_iters,
        f_calls_limit = f_calls_limit > 0 ? f_calls_limit : typemax(Int),
        show_trace = false,
        store_trace = false,
        callback = callback,
        g_tol = g_tol,
        f_tol = f_tol,
    )

    t0 = time()
    res = optimize(objective, gradient!, prob.lower_mg, prob.upper_mg, copy(x0), Fminbox(LBFGS()), opts)
    runtime_s = time() - t0
    xopt = clamp.(Optim.minimizer(res), prob.lower_mg, prob.upper_mg)
    push_trace!("final", iter_ref[] + 1, copy(xopt))

    initial_result = evaluate!(evaluator, x0)
    optimal_result = evaluate!(evaluator, xopt)
    x0_full = expand_decision_doses(x0, prob.decision_groups, length(prob.dose_times_days))
    xopt_full = expand_decision_doses(xopt, prob.decision_groups, length(prob.dose_times_days))

    println(@sprintf("objective initial=%.6g optimal=%.6g runtime_s=%.2f", initial_result.scaled.loss, optimal_result.scaled.loss, runtime_s))
    println("initial decision doses mg = $(round.(x0; digits = 3))")
    println("optimal decision doses mg = $(round.(xopt; digits = 3))")

    init_metrics, init_traj = simulate_cohort_trajectories(prob, x0)
    opt_metrics, opt_traj = simulate_cohort_trajectories(prob, xopt)

    CSV.write(joinpath(prob.out_dir, "loss_trace.csv"), DataFrame(trace_rows))
    CSV.write(joinpath(prob.out_dir, "trajectory_initial.csv"), init_traj)
    CSV.write(joinpath(prob.out_dir, "trajectory_optimal.csv"), opt_traj)
    CSV.write(joinpath(prob.out_dir, "endpoint_metrics_initial.csv"), init_metrics)
    CSV.write(joinpath(prob.out_dir, "endpoint_metrics_optimal.csv"), opt_metrics)
    CSV.write(
        joinpath(prob.out_dir, "dose_schedule.csv"),
        DataFrame(
            dose_idx = collect(1:length(prob.dose_times_days)),
            time_day = prob.dose_times_days,
            clinical_dose_mg = prob.clinical_doses_mg,
            initial_random_dose_mg = x0_full,
            optimized_dose_mg = xopt_full,
        ),
    )
    CSV.write(
        joinpath(prob.out_dir, "decision_schedule.csv"),
        DataFrame(
            decision_idx = collect(1:length(prob.decision_labels)),
            label = prob.decision_labels,
            lower_mg = prob.lower_mg,
            upper_mg = prob.upper_mg,
            clinical_dose_mg = prob.clinical_decision_doses_mg,
            initial_random_dose_mg = x0,
            optimized_dose_mg = xopt,
        ),
    )

    summary = Dict(
        "regimen_label" => regimen_label,
        "threads" => nthreads(),
        "cohort_size" => length(prob.cohort),
        "solver" => string(prob.alg),
        "objective_dt_days" => prob.objective_dt,
        "sample_dt_days" => prob.sample_dt,
        "horizon_days" => prob.horizon_days,
        "dose_times_days" => prob.dose_times_days,
        "decision_labels" => prob.decision_labels,
        "decision_groups" => prob.decision_groups,
        "lower_mg" => prob.lower_mg,
        "upper_mg" => prob.upper_mg,
        "clinical_decision_doses_mg" => prob.clinical_decision_doses_mg,
        "clinical_doses_mg" => prob.clinical_doses_mg,
        "initial_random_decision_doses_mg" => x0,
        "optimized_decision_doses_mg" => xopt,
        "initial_random_doses_mg" => x0_full,
        "optimized_doses_mg" => xopt_full,
        "initial_dose_label" => format_dose_label(x0_full),
        "optimized_dose_label" => format_dose_label(xopt_full),
        "runtime_s" => runtime_s,
        "iterations_logged" => length(trace_rows),
        "scales" => Dict(
            "tox_peak_mean" => scales.tox_peak_mean,
            "tox_peak_max" => scales.tox_peak_max,
            "tox_auc" => scales.tox_auc,
            "tumor_terminal" => scales.tumor_terminal,
            "tumor_auc" => scales.tumor_auc,
        ),
        "weights" => Dict(
            "tox_peak_mean" => weights.tox_peak_mean,
            "tox_peak_max" => weights.tox_peak_max,
            "tox_auc" => weights.tox_auc,
            "tumor_terminal" => weights.tumor_terminal,
            "tumor_auc" => weights.tumor_auc,
        ),
        "objective_initial" => Float64(initial_result.scaled.loss),
        "objective_optimal" => Float64(optimal_result.scaled.loss),
        "raw_initial" => to_float_dict(initial_result.raw),
        "raw_optimal" => to_float_dict(optimal_result.raw),
        "scaled_initial" => to_float_dict(initial_result.scaled),
        "scaled_optimal" => to_float_dict(optimal_result.scaled),
    )
    open(joinpath(prob.out_dir, "optimization_summary.json"), "w") do io
        JSON3.pretty(io, summary)
    end

    println(joinpath(prob.out_dir, "optimization_summary.json"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
