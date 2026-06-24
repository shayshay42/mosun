using CSV
using Base.Threads
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
include(joinpath(@__DIR__, "src", "MosunModelCoreSupport.jl"))

const DTO_PASSTHROUGH = SensitivityADPassThrough()
const OTD_FORWARD = ForwardSensitivity()
const OTD_REVERSE = ReverseDiffAdjoint()
const MMCS = MosunModelCoreSupport
const TOPPARAM10_VARIANT_IDS = [5, 9, 14, 20, 24, 25, 27, 28]
const PAPER84_ACTUAL_RP2D_DOSES_MG = Float64[1.0, 2.0, 60.0, 60.0, 30.0, 30.0]
const HOSSEINI_EFAST_VPOP250_COHORT_NAMES = Set(["hosseini_efast_vpop250", "susilo_efast_hosseini_vpop250", "efast_hosseini_vpop250"])
const HOSSEINI_EFAST_DUMMY_PARAMS = Set(["dummy_null_1", "dummy_null_2", "dummy_null_3"])
const HOSSEINI_EFAST_TUMOR_LATENT_PARAMS = Set(["tumor_burden_factor", "BT_ratio_tumor_init"])
const ACHIEVABLE_IDEAL_OBJECTIVES = Set([
    "patient_achievable_ideal_tracking",
    "global_best_achievable_ideal_tracking",
    "patient_achievable_ideal_tracking_calibrated",
    "achievable_endpoint_cycle20",
    "achievable_endpoint_cycle60",
    "patient_achievable_ideal_tracking_calibrated_cycle20",
    "patient_achievable_ideal_tracking_calibrated_cycle60",
])
const ACHIEVABLE_TARGET_CURVES_CACHE = Ref{Union{Nothing, DataFrame}}(nothing)
const ACHIEVABLE_TARGET_SUMMARY_CACHE = Ref{Union{Nothing, DataFrame}}(nothing)

env_float_bench(name::String, default::Float64) = parse(Float64, get(ENV, name, string(default)))
is_ideal_tracking_objective(objective_name::AbstractString) = lowercase(strip(String(objective_name))) == "patient_ideal_tracking" || lowercase(strip(String(objective_name))) in ACHIEVABLE_IDEAL_OBJECTIVES
is_calibrated_achievable_objective(objective_name::AbstractString) = lowercase(strip(String(objective_name))) in ("patient_achievable_ideal_tracking_calibrated", "achievable_endpoint_cycle20", "achievable_endpoint_cycle60", "patient_achievable_ideal_tracking_calibrated_cycle20", "patient_achievable_ideal_tracking_calibrated_cycle60")
is_cycle_target_achievable_objective(objective_name::AbstractString) = lowercase(strip(String(objective_name))) in ("achievable_endpoint_cycle20", "achievable_endpoint_cycle60", "patient_achievable_ideal_tracking_calibrated_cycle20", "patient_achievable_ideal_tracking_calibrated_cycle60")
is_endpoint_achievable_objective(objective_name::AbstractString) = lowercase(strip(String(objective_name))) in ("achievable_endpoint_cycle20", "achievable_endpoint_cycle60")

mutable struct TraceLogger
    objective_name::String
    method::String
    t0::Float64
    eval::Int
    best_loss::Float64
    best_x::Vector{Float64}
    trace_path::String
    progress_path::String
    print_every::Int
    max_evals::Int
    rows::Vector{NamedTuple}
end

struct MaxEvalStop <: Exception
    max_evals::Int
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
    remission_refs
    paper_stepup_refs = nothing
    objective_name::String
    objective_label::String
    scenario_label::String
    solver_label::String
    random_seed::Int
    sample_selection_mode::String
    x0_mode::String
    cohort_source::String
    cohort_dir::String
end

function clone_case_with_solver(case::MethodBenchmarkCase, alg, solver_label::AbstractString)
    kwargs = Dict{Symbol, Any}()
    for fn in fieldnames(MethodBenchmarkCase)
        if fn == :alg
            kwargs[fn] = alg
        elseif fn == :solver_label
            kwargs[fn] = String(solver_label)
        else
            kwargs[fn] = getfield(case, fn)
        end
    end
    return MethodBenchmarkCase(; kwargs...)
end

function safe_anchor_metrics(case::MethodBenchmarkCase, decision_doses::AbstractVector; fallback_solver_names::Tuple = ("qndf",))
    try
        return trajectory_metrics(case, decision_doses; sensealg = nothing), case.solver_label
    catch err
        last_err = err
        for solver_name in fallback_solver_names
            alg, label = parse_solver(solver_name)
            fallback_case = clone_case_with_solver(case, alg, label)
            try
                return trajectory_metrics(fallback_case, decision_doses; sensealg = nothing), label
            catch err2
                last_err = err2
            end
        end
        throw(last_err)
    end
end

function trapz_until_generic(t::AbstractVector, y::AbstractVector, horizon::Real)
    idxs = [i for i in eachindex(t) if t[i] <= horizon + 1e-8]
    isempty(idxs) && return zero(promote_type(eltype(t), eltype(y), Float64))
    return trapz_generic(t[idxs], y[idxs])
end

function first_window_peak_value(t::AbstractVector, y::AbstractVector; horizon::Real = 7.0)
    idxs = [i for i in eachindex(t) if t[i] <= horizon + 1e-8]
    isempty(idxs) && return zero(promote_type(eltype(t), eltype(y), Float64))
    return maximum(y[idxs])
end

function smooth_first_window_peak(case::MethodBenchmarkCase, t::AbstractVector, y::AbstractVector; horizon::Real = 7.0)
    vals = [y[i] for i in eachindex(t) if t[i] <= horizon + 1e-8]
    isempty(vals) && return zero(promote_type(eltype(t), eltype(y), Float64))
    return smoothmax_generic(vals, case.tox_tau)
end

function btumor_auc_0_84d_value(t::AbstractVector, bt::AbstractVector)
    return trapz_until_generic(t, smooth_pos.(bt), 84.0)
end

function day84_tumor_size_change_pct_value(bt::AbstractVector, bt0)
    T = promote_type(eltype(bt), typeof(bt0), Float64)
    epsT = T(1e-12)
    return T(100.0) * (bt[end] / (bt0 + epsT) - one(T))
end

function make_stepup_anchor_ref(case::MethodBenchmarkCase, m; anchor_regimen::AbstractString, anchor_label::AbstractString, anchor_solver_label::AbstractString, anchor_doses_mg = nothing, eps::Float64 = 1e-12)
    return (
        peak_il6_0_7d_ref = Float64(smooth_first_window_peak(case, m.t, m.il6; horizon = 7.0)),
        il6_auc_0_84d_ref = Float64(trapz_until_generic(m.t, smooth_pos.(m.il6), 84.0)),
        day84_tumor_size_change_pct_ref = Float64(day84_tumor_size_change_pct_value(m.bt, m.bt0)),
        btumor_auc_0_84d_ref = Float64(btumor_auc_0_84d_value(m.t, m.bt)),
        total_dose_mg_ref = isnothing(anchor_doses_mg) ? NaN : Float64(sum(Float64.(anchor_doses_mg))),
        eps = eps,
        anchor_regimen = String(anchor_regimen),
        anchor_label = String(anchor_label),
        anchor_solver_label = String(anchor_solver_label),
    )
end

function resolve_stepup_anchor_refs(case::MethodBenchmarkCase; anchor_key::Union{Nothing, Symbol} = nothing)
    refs = case.paper_stepup_refs
    refs === nothing && error("paper_stepup_refs not configured for objective=$(case.objective_name)")
    if isnothing(anchor_key)
        anchor_key = case.objective_name == "rp2d_stepup" ? :rp2d : :paper
    end
    if :peak_il6_0_7d_ref in keys(refs)
        return refs
    elseif anchor_key in keys(refs)
        return getproperty(refs, anchor_key)
    end
    error("Anchor $(anchor_key) not configured in paper_stepup_refs for objective=$(case.objective_name)")
end

function stepup_anchor_terms(case::MethodBenchmarkCase, m; anchor_key::Union{Nothing, Symbol} = nothing)
    refs = resolve_stepup_anchor_refs(case; anchor_key = anchor_key)
    T = promote_type(eltype(m.t), eltype(m.bt), eltype(m.il6), typeof(m.bt0), Float64)
    epsT = T(get(refs, :eps, 1e-12))
    first_peak = smooth_first_window_peak(case, m.t, m.il6; horizon = 7.0)
    day84_change = day84_tumor_size_change_pct_value(m.bt, m.bt0)
    bt_auc = btumor_auc_0_84d_value(m.t, m.bt)
    safety_term = log((first_peak + epsT) / (T(refs.peak_il6_0_7d_ref) + epsT))
    day84_gap = smooth_pos(day84_change - T(refs.day84_tumor_size_change_pct_ref))
    day84_penalty = (day84_gap / T(10.0))^2
    auc_gap = smooth_pos(log((bt_auc + epsT) / (T(refs.btumor_auc_0_84d_ref) + epsT)))
    auc_penalty = auc_gap^2
    total = safety_term + T(10.0) * day84_penalty + auc_penalty
    return (
        total = total,
        safety_term = safety_term,
        day84_penalty = day84_penalty,
        auc_penalty = auc_penalty,
        first_peak_il6_0_7d = first_peak,
        day84_tumor_size_change_pct = day84_change,
        btumor_auc_0_84d = bt_auc,
    )
end

stepup_paper_terms(case::MethodBenchmarkCase, m) = stepup_anchor_terms(case, m; anchor_key = :paper)
rp2d_stepup_terms(case::MethodBenchmarkCase, m) = stepup_anchor_terms(case, m; anchor_key = :rp2d)

function total_expanded_dose_mg(case::MethodBenchmarkCase, decision_doses::AbstractVector)
    full = expand_decision_doses_generic(decision_doses, case.decision_groups, length(case.dose_times_days))
    return sum(full)
end

function safety_constrained_ti_terms(case::MethodBenchmarkCase, m, decision_doses::AbstractVector; anchor_key::Symbol = :rp2d)
    refs = resolve_stepup_anchor_refs(case; anchor_key = anchor_key)
    T = promote_type(eltype(m.t), eltype(m.bt), eltype(m.il6), typeof(m.bt0), eltype(decision_doses), Float64)
    epsT = T(get(refs, :eps, 1e-12))

    first_peak = smooth_first_window_peak(case, m.t, m.il6; horizon = 7.0)
    il6_auc = trapz_until_generic(m.t, smooth_pos.(m.il6), 84.0)
    day84_change = day84_tumor_size_change_pct_value(m.bt, m.bt0)
    bt_auc = btumor_auc_0_84d_value(m.t, m.bt)
    total_dose = total_expanded_dose_mg(case, decision_doses)

    ref_peak = T(refs.peak_il6_0_7d_ref)
    ref_il6_auc = T(get(refs, :il6_auc_0_84d_ref, 84.0 * ref_peak))
    ref_day84 = T(refs.day84_tumor_size_change_pct_ref)
    ref_bt_auc = T(refs.btumor_auc_0_84d_ref)
    ref_dose_raw = get(refs, :total_dose_mg_ref, NaN)
    ref_dose = isfinite(Float64(ref_dose_raw)) ? T(ref_dose_raw) : T(sum(PAPER84_ACTUAL_RP2D_DOSES_MG))

    peak_log_ratio = log((first_peak + epsT) / (ref_peak + epsT))
    il6_auc_log_ratio = log((il6_auc + epsT) / (ref_il6_auc + epsT))
    day84_gap_pct = day84_change - ref_day84
    bt_auc_log_ratio = log((bt_auc + epsT) / (ref_bt_auc + epsT))
    dose_log_ratio = log((total_dose + epsT) / (ref_dose + epsT))

    peak_safety_excess = smooth_pos(peak_log_ratio)
    auc_safety_excess = smooth_pos(il6_auc_log_ratio)
    safety_reward_gate = exp(-T(25.0) * (peak_safety_excess + auc_safety_excess))

    peak_penalty = peak_safety_excess^2
    il6_auc_penalty = auc_safety_excess^2
    day84_penalty = (smooth_pos(day84_gap_pct) / T(10.0))^2
    tumor_auc_penalty = smooth_pos(bt_auc_log_ratio)^2

    peak_reward = smooth_pos(-peak_log_ratio)
    il6_auc_reward = smooth_pos(-il6_auc_log_ratio)
    day84_reward = smooth_pos(-day84_gap_pct) / T(10.0)
    tumor_auc_reward = smooth_pos(-bt_auc_log_ratio)
    dose_sparing_term = dose_log_ratio

    gated_reward = safety_reward_gate * (
        T(1.0) * peak_reward +
        T(1.0) * il6_auc_reward +
        T(2.0) * day84_reward +
        T(1.0) * tumor_auc_reward -
        T(0.15) * dose_sparing_term
    )

    total = T(20.0) * peak_penalty +
        T(10.0) * il6_auc_penalty +
        T(10.0) * day84_penalty +
        T(5.0) * tumor_auc_penalty -
        gated_reward

    return (
        total = total,
        peak_safety_penalty = peak_penalty,
        il6_auc_safety_penalty = il6_auc_penalty,
        day84_penalty = day84_penalty,
        tumor_auc_penalty = tumor_auc_penalty,
        peak_safety_reward = safety_reward_gate * peak_reward,
        il6_auc_safety_reward = safety_reward_gate * il6_auc_reward,
        day84_reward = safety_reward_gate * day84_reward,
        tumor_auc_reward = safety_reward_gate * tumor_auc_reward,
        dose_sparing_term = dose_sparing_term,
        safety_reward_gate = safety_reward_gate,
        first_peak_il6_0_7d = first_peak,
        first_peak_il6_0_7d_ref = ref_peak,
        il6_auc_0_84d = il6_auc,
        il6_auc_0_84d_ref = ref_il6_auc,
        day84_tumor_size_change_pct = day84_change,
        day84_tumor_size_change_pct_ref = ref_day84,
        btumor_auc_0_84d = bt_auc,
        btumor_auc_0_84d_ref = ref_bt_auc,
        total_dose_mg = total_dose,
        total_dose_mg_ref = ref_dose,
        peak_log_ratio = peak_log_ratio,
        il6_auc_log_ratio = il6_auc_log_ratio,
        day84_gap_pct = day84_gap_pct,
        btumor_auc_log_ratio = bt_auc_log_ratio,
        dose_log_ratio = dose_log_ratio,
    )
end

safety_constrained_ti_objective(case::MethodBenchmarkCase, m, decision_doses::AbstractVector) =
    safety_constrained_ti_terms(case, m, decision_doses; anchor_key = :rp2d).total

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

function top_param_vpop_dir()
    out_dir_env = get(
        ENV,
        "SINGLE_BENCH_TOPPARAM10_DIR",
        joinpath("generated", "figures", "vpop", "paper84_topparam10_loglhs_vpop1000_20260316"),
    )
    return isabspath(out_dir_env) ? out_dir_env : joinpath(REPO_ROOT, out_dir_env)
end

function paper_regimen_catalog()
    return Dict(
        "paper_fractionated_flat" => (
            label = "Paper flat-fractionation (6.7/6.7/6.7->20)",
            decision_doses_mg = Float64[6.7, 6.7, 6.7, 20.0],
        ),
        "paper_single_step" => (
            label = "Paper single-step (1.6/20/20->20)",
            decision_doses_mg = Float64[1.6, 20.0, 20.0, 20.0],
        ),
        "paper_recommended" => (
            label = "Paper double-step (1.6/10/10->20)",
            decision_doses_mg = Float64[1.6, 10.0, 10.0, 20.0],
        ),
        "rp2d" => (
            label = "RP2D 1/2/60/60/30/30",
            decision_doses_mg = copy(PAPER84_ACTUAL_RP2D_DOSES_MG),
        ),
    )
end

function load_benchmark_scenario()
    scenario = lowercase(get(ENV, "SINGLE_BENCH_SCENARIO", "realdata168"))
    if scenario == "realdata168"
        regimen = load_realdata_regimen()
        return (
            scenario_name = scenario,
            scenario_label = "DLBCL-like virtual patient under 8-cycle mosunetuzumab step-up dosing",
            dose_times_days = regimen.dose_times_days,
            dose_mg = regimen.dose_mg,
            horizon_days = regimen.horizon_days,
            lower_mg = nothing,
            upper_mg = nothing,
        )
    elseif scenario == "paper84"
        decision_mode = lowercase(strip(get(ENV, "SINGLE_BENCH_PAPER84_DECISION_MODE", "four_stepup")))
        n_controls = decision_mode in ("six_independent", "six", "all_events") ? 6 : 4
        lower_mg = fill(parse(Float64, get(ENV, "SINGLE_BENCH_PAPER84_MIN_MG", "0.0")), n_controls)
        upper_mg = fill(parse(Float64, get(ENV, "SINGLE_BENCH_PAPER84_MAX_MG", "20.0")), n_controls)
        return (
            scenario_name = scenario,
            scenario_label = "DLBCL-like virtual patient under the 84-day paper step-up dosing scenario",
            dose_times_days = Float64[0.0, 7.0, 14.0, 21.0, 42.0, 63.0],
            dose_mg = Float64[1.6, 10.0, 10.0, 20.0, 20.0, 20.0],
            horizon_days = 84.0,
            lower_mg = lower_mg,
            upper_mg = upper_mg,
        )
    else
        error("Unsupported SINGLE_BENCH_SCENARIO=$scenario")
    end
end

function build_benchmark_decision_scheme(dose_times_days::Vector{Float64}, clinical_doses_mg::Vector{Float64}, scenario_name::AbstractString)
    s = lowercase(String(scenario_name))
    if s == "realdata168"
        return build_decision_scheme(dose_times_days, clinical_doses_mg)
    elseif s == "paper84"
        length(dose_times_days) == 6 || error("paper84 scenario expects 6 dose events")
        decision_mode = lowercase(strip(get(ENV, "SINGLE_BENCH_PAPER84_DECISION_MODE", "four_stepup")))
        groups = if decision_mode in ("six_independent", "six", "all_events")
            [Int[1], Int[2], Int[3], Int[4], Int[5], Int[6]]
        elseif decision_mode in ("four_stepup", "four", "c2plus_tied")
            [Int[1], Int[2], Int[3], Int[4, 5, 6]]
        else
            error("Unsupported SINGLE_BENCH_PAPER84_DECISION_MODE=$decision_mode")
        end
        labels = if length(groups) == 6
            ["C1D1", "C1D8", "C1D15", "C2D1", "C3D1", "C4D1"]
        else
            ["C1D1", "C1D8", "C1D15", "C2plus_q21"]
        end
        decision_vals = Float64[]
        for grp in groups
            push!(decision_vals, mean(clinical_doses_mg[grp]))
        end
        return labels, groups, decision_vals
    else
        error("Unsupported benchmark scenario=$scenario_name")
    end
end

function pick_sample(samples::Vector, seed::Int)
    sample_id_env = strip(get(ENV, "SINGLE_BENCH_SAMPLE_ID", ""))
    if !isempty(sample_id_env)
        wanted = parse(Int, sample_id_env)
        for sample in samples
            sample.sample_id == wanted && return sample
        end
        error("SINGLE_BENCH_SAMPLE_ID=$wanted not found in cohort")
    end
    rng = MersenneTwister(seed)
    return samples[rand(rng, eachindex(samples))]
end

function build_topparam10_base_params(horizon_days::Float64)
    p = deepcopy(MMC.default_params())
    variant_overrides = MMCS.load_variant_overrides()
    for vid in TOPPARAM10_VARIANT_IDS
        for (name, value) in get(variant_overrides, vid, Pair{String, Float64}[])
            MMC.has_parameter(name) || continue
            MMC.set_param!(p, name, value)
        end
    end
    MMC.has_parameter("PKflag") && MMC.set_param!(p, "PKflag", 1.0)
    MMC.has_parameter("VPid") && MMC.set_param!(p, "VPid", 1.0)
    MMC.has_parameter("fvalidation") && MMC.set_param!(p, "fvalidation", 0.0)
    MMC.has_parameter("end_time") && MMC.set_param!(p, "end_time", horizon_days)
    return p
end

function load_topparam10_cohort(horizon_days::Float64)
    cohort_dir = top_param_vpop_dir()
    samples_path = joinpath(cohort_dir, "topparam_vpop_samples.csv")
    isfile(samples_path) || error("topparam10 cohort CSV not found: $samples_path")
    df = DataFrame(CSV.File(samples_path))
    samples = CohortSample[]
    for r in eachrow(df)
        p = build_topparam10_base_params(horizon_days)
        for cname in names(df)
            MMC.has_parameter(cname) || continue
            v = try
                Float64(r[Symbol(cname)])
            catch
                NaN
            end
            isfinite(v) || continue
            MMC.set_param!(p, cname, v)
        end
        push!(samples, CohortSample(Int(r.sample_id), p))
    end
    isempty(samples) && error("No samples found in topparam10 cohort CSV")
    return samples, "topparam10", cohort_dir
end

function resolve_repo_path(path_str::AbstractString)
    s = String(path_str)
    return isabspath(s) ? s : joinpath(REPO_ROOT, s)
end

function first_existing_path(paths::Vector{String}, label::String)
    for path in paths
        isfile(path) && return path
    end
    error("Could not find $(label). Checked: $(join(paths, ", "))")
end

function hosseini_efast_vpop250_dir()
    raw = get(
        ENV,
        "HOSSEINI_EFAST_VPOP250_DIR",
        joinpath("generated", "figures", "vpop_pruning", "susilo_efast_hosseini2020_fig5_vpop250_20260505"),
    )
    return resolve_repo_path(raw)
end

function hosseini_efast_vpop250_params_csv()
    raw = get(ENV, "HOSSEINI_EFAST_VPOP250_PARAMETERS_CSV", joinpath(hosseini_efast_vpop250_dir(), "selected_vpop250_parameters.csv"))
    return resolve_repo_path(raw)
end

function hosseini_efast_vpop250_universe_csv()
    env_path = strip(get(ENV, "HOSSEINI_EFAST_VPOP250_PARAMETER_UNIVERSE_CSV", ""))
    if !isempty(env_path)
        path = resolve_repo_path(env_path)
        isfile(path) || error("HOSSEINI_EFAST_VPOP250_PARAMETER_UNIVERSE_CSV not found: $path")
        return path
    end
    return first_existing_path(
        [
            joinpath(hosseini_efast_vpop250_dir(), "parameter_universe.csv"),
            joinpath(REPO_ROOT, "generated", "server_results", "susilo_fig5_il6_dummy_efast_20260502", "susilo_fig5_il6_dummy_efast_20260502", "parameter_universe.csv"),
            joinpath(REPO_ROOT, "generated", "figures", "sensitivity", "susilo_fig5_il6_dummy_efast_20260502", "parameter_universe.csv"),
        ],
        "Hosseini eFAST parameter_universe.csv",
    )
end

function boolish(x)
    if x isa Bool
        return x
    elseif x isa Number
        return x != 0
    else
        return lowercase(strip(String(x))) in ("1", "true", "t", "yes", "y")
    end
end

function baseline_float(p, name::String, default::Float64)
    sym = Symbol(name)
    hasproperty(p, sym) ? Float64(getproperty(p, sym)) : default
end

function apply_hosseini_efast_row_to_params(base_params, universe::DataFrame, row)
    p = deepcopy(base_params)
    latent = Dict{String, Float64}()
    for urow in eachrow(universe)
        name = String(urow.parameter)
        sym = Symbol(name)
        hasproperty(row, sym) || error("VPop250 parameter row missing column $(name)")
        value = Float64(row[sym])
        if name in HOSSEINI_EFAST_DUMMY_PARAMS
            continue
        elseif name in HOSSEINI_EFAST_TUMOR_LATENT_PARAMS
            latent[name] = value
            continue
        elseif boolish(urow.applied_to_model) && MMC.has_parameter(name)
            MMC.set_param!(p, name, value)
        end
    end

    base_bpbo = baseline_float(base_params, "Bpbo_perml", 1.0)
    base_trpbo = baseline_float(base_params, "Trpbo_perml", 1.0)
    base_kbp = baseline_float(base_params, "KBptumor", 1.0)
    burden_factor = get(latent, "tumor_burden_factor", 1.0)
    bt_ratio = get(latent, "BT_ratio_tumor_init", 20.0)
    base_btumor_perml = base_kbp * max(base_bpbo, 1e-12)
    btumor_perml = base_btumor_perml * burden_factor
    MMC.has_parameter("KBptumor") && MMC.set_param!(p, "KBptumor", btumor_perml / max(base_bpbo, 1e-12))
    MMC.has_parameter("KTrptumor") && MMC.set_param!(p, "KTrptumor", btumor_perml / max(bt_ratio * base_trpbo, 1e-12))
    return p
end

function load_hosseini_efast_vpop250_cohort(horizon_days::Float64)
    params_path = hosseini_efast_vpop250_params_csv()
    universe_path = hosseini_efast_vpop250_universe_csv()
    isfile(params_path) || error("VPop250 parameter CSV not found: $params_path")
    df = DataFrame(CSV.File(params_path))
    universe = DataFrame(CSV.File(universe_path))
    "vpop_id" in names(df) || error("selected_vpop250_parameters.csv must contain vpop_id")
    "parameter" in names(universe) || error("parameter_universe.csv must contain parameter")
    "applied_to_model" in names(universe) || error("parameter_universe.csv must contain applied_to_model")

    base = build_topparam10_base_params(horizon_days)
    samples = CohortSample[]
    for row in eachrow(df)
        p = apply_hosseini_efast_row_to_params(base, universe, row)
        push!(samples, CohortSample(Int(row.vpop_id), p))
    end
    length(samples) == length(unique(s.sample_id for s in samples)) || error("Duplicate vpop_id values in $params_path")
    sort!(samples, by = s -> s.sample_id)
    return samples, "hosseini_efast_vpop250", dirname(params_path)
end

function load_hosseini_efast_vpop250_mapping(ids::AbstractVector{Int})
    params_path = hosseini_efast_vpop250_params_csv()
    df = DataFrame(CSV.File(params_path))
    wanted = Set(ids)
    keep = df[[Int(row.vpop_id) in wanted for row in eachrow(df)], :]
    sort!(keep, :vpop_id)
    return keep
end

function write_hosseini_efast_vpop250_mapping(out_dir::AbstractString, ids::AbstractVector{Int})
    mapping = load_hosseini_efast_vpop250_mapping(ids)
    path = joinpath(out_dir, "vpop250_candidate_mapping.csv")
    CSV.write(path, mapping)
    return path
end

function load_benchmark_samples(horizon_days::Float64)
    cohort = lowercase(strip(get(ENV, "SINGLE_BENCH_COHORT", "legacy_lhs")))
    if cohort in ("legacy_lhs", "legacy", "clinical_lhs")
        return load_lhs_cohort(horizon_days), "legacy_lhs", joinpath(REPO_ROOT, "generated", "figures", "clinical_lhs_regimen_20260310")
    elseif cohort in ("topparam10", "topparam10_vpop")
        return load_topparam10_cohort(horizon_days)
    elseif cohort in HOSSEINI_EFAST_VPOP250_COHORT_NAMES
        return load_hosseini_efast_vpop250_cohort(horizon_days)
    else
        error("Unsupported SINGLE_BENCH_COHORT=$cohort")
    end
end

function load_topparam10_il6_ref()
    metrics_path = joinpath(top_param_vpop_dir(), "topparam_vpop_metrics.csv")
    isfile(metrics_path) || error("topparam10 metrics CSV not found: $metrics_path")
    df = DataFrame(CSV.File(metrics_path))
    keep = df[(df.regimen .== "paper_best") .& (df.status .== "ok"), :]
    vals = [Float64(v) for v in keep.peak_il6 if isfinite(v)]
    isempty(vals) && error("No valid paper_best peak_il6 values found in $metrics_path")
    return max(median(vals), 1e-6)
end

function append_row!(path::AbstractString, row::NamedTuple)
    tbl = DataFrame([row])
    write_header = !isfile(path) || filesize(path) == 0
    CSV.write(path, tbl; append = !write_header, writeheader = write_header)
    return nothing
end

@inline function shard_mode_enabled()
    return lowercase(strip(get(ENV, "SINGLE_BENCH_SHARD_BY_METHOD", "0"))) in ("1", "true", "yes")
end

function method_slug(method::AbstractString)
    return replace(lowercase(String(method)), r"[^a-z0-9]+" => "_")
end

function artifact_path(out_dir::AbstractString, stem::AbstractString; method::Union{Nothing, AbstractString} = nothing, ext::AbstractString = "csv")
    if !isnothing(method) && shard_mode_enabled()
        return joinpath(out_dir, string(stem, "__", method_slug(method), ".", ext))
    end
    return joinpath(out_dir, string(stem, ".", ext))
end

function method_specific_x0(x0_default::Vector{Float64}, lower::Vector{Float64}, upper::Vector{Float64}, clinical_decision::Vector{Float64}, seed::Int, scenario_name::AbstractString)
    x0_mode = lowercase(strip(get(ENV, "SINGLE_BENCH_X0_MODE", lowercase(String(scenario_name)) == "paper84" ? "random" : "baseline")))
    x0 = if x0_mode == "baseline"
        clamp.(x0_default, lower, upper)
    elseif x0_mode == "random"
        rng = MersenneTwister(seed + 101)
        lower .+ rand(rng, length(lower)) .* (upper .- lower)
    elseif x0_mode == "clinical"
        clamp.(clinical_decision, lower, upper)
    else
        error("Unsupported SINGLE_BENCH_X0_MODE=$x0_mode")
    end
    return Float64.(x0), x0_mode
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
    max_evals = if method == "finite_diff_bfgs"
        parse(Int, get(ENV, "SINGLE_BENCH_FD_BFGS_F_CALLS", get(ENV, "SINGLE_BENCH_FD_LBFGS_F_CALLS", "400")))
    elseif method == "dto_forward_ad_bfgs"
        parse(Int, get(ENV, "SINGLE_BENCH_BFGS_F_CALLS", get(ENV, "SINGLE_BENCH_LBFGS_F_CALLS", "400")))
    elseif method == "finite_diff_lbfgs"
        parse(Int, get(ENV, "SINGLE_BENCH_FD_LBFGS_F_CALLS", "400"))
    else
        0
    end
    return TraceLogger(
        objective_name,
        method,
        time(),
        0,
        Inf,
        Float64[],
        artifact_path(out_dir, "optimization_traces"; method = method),
        artifact_path(out_dir, "optimization_progress"; method = method, ext = "json"),
        print_every,
        max_evals,
        NamedTuple[],
    )
end

function log_eval!(logger::TraceLogger, x::AbstractVector, loss::Real)
    logger.eval += 1
    loss64 = Float64(loss)
    if loss64 <= logger.best_loss
        logger.best_loss = loss64
        logger.best_x = Float64.(x)
    end
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
        if logger.max_evals > 0 && logger.eval >= logger.max_evals && !(x[1] isa ForwardDiff.Dual)
            throw(MaxEvalStop(logger.max_evals))
        end
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

function trapz_weights(x::AbstractVector{<:Real})
    n = length(x)
    n <= 1 && return zeros(Float64, n)
    w = zeros(Float64, n)
    w[1] = Float64(x[2] - x[1]) / 2
    @inbounds for i in 2:(n - 1)
        w[i] = Float64(x[i + 1] - x[i - 1]) / 2
    end
    w[n] = Float64(x[n] - x[n - 1]) / 2
    return w
end

function smoothmax_weights(v::AbstractVector, tau::Real)
    isempty(v) && return zero(eltype(v)), zeros(Float64, 0)
    m = maximum(v)
    shifted = (v .- m) ./ tau
    w = exp.(shifted)
    denom = sum(w)
    weights = Float64.(w ./ denom)
    val = m + tau * log(denom)
    return val, weights
end

function expand_decision_doses_generic(decision_doses::AbstractVector, decision_groups::Vector{Vector{Int}}, n_events::Int)
    T = eltype(decision_doses)
    if length(decision_doses) == n_events
        return T.(decision_doses)
    end
    length(decision_doses) == length(decision_groups) || error("Dose vector length $(length(decision_doses)) must match either $(length(decision_groups)) decision groups or $(n_events) dose events")
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
    return MMC.smooth_tdbc_cutoff(val)
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

@inline smooth_pos(x) = MMC.smooth_positive_part(x)
@inline smooth_pos_deriv(x) = MMC.smooth_positive_part_deriv(x)

@inline function smooth_log_ratio(x, scale, eps)
    return log(smooth_pos(x / scale) + eps)
end

function pack_channels(bt::AbstractVector, il6::AbstractVector, tdbc::AbstractVector)
    return vcat(collect(bt), collect(il6), collect(tdbc))
end

function unpack_channels(y::AbstractVector, n_times::Int)
    bt = view(y, 1:n_times)
    il6 = view(y, (n_times + 1):(2 * n_times))
    tdbc = view(y, (2 * n_times + 1):(3 * n_times))
    return bt, il6, tdbc
end

function loss_from_channels(case::MethodBenchmarkCase, bt::AbstractVector, il6::AbstractVector, tdbc::AbstractVector, bt0)
    T = promote_type(eltype(bt), eltype(il6), eltype(tdbc), typeof(bt0), Float64)
    epsT = T(1e-12)
    if case.objective_name == "simple"
        tox_peak_simple = smoothmax_generic(il6, case.tox_tau)
        tumor_terminal = bt[end] / (bt0 + epsT)
        return 0.5 * (tox_peak_simple / case.simple_scales.tox) +
            0.5 * (tumor_terminal / case.simple_scales.tumor)
    elseif case.objective_name == "auc_combo"
        il6_pos = smooth_pos.(il6)
        tox_auc = trapz_generic(case.saveat, il6_pos) / case.horizon_days
        tumor_auc = trapz_generic(case.saveat, bt ./ (bt0 + epsT)) / case.horizon_days
        s = case.clinical_scales
        return (tox_auc / s.tox_auc) + (tumor_auc / s.tumor_auc)
    elseif case.objective_name == "remission_il6"
        lb, li = remission_il6_components(case, case.saveat, bt, il6, bt0)
        return (lb / (case.remission_refs.lb_paper + T(1e-8))) +
            (li / (case.remission_refs.li_paper + T(1e-8)))
    elseif case.objective_name == "clinical"
        window_peaks = [
            smoothmax_generic(
                [il6[i] for i in eachindex(case.saveat) if case.saveat[i] >= td && case.saveat[i] <= td + case.tox_window_days + 1e-10],
                case.tox_tau,
            ) for td in case.dose_times_days
        ]
        tox_peak_mean = sum(window_peaks) / length(window_peaks)
        tox_peak_max = maximum(window_peaks)
        il6_pos = smooth_pos.(il6)
        tox_auc = trapz_generic(case.saveat, il6_pos) / case.horizon_days
        tumor_terminal = bt[end] / (bt0 + epsT)
        tumor_auc = trapz_generic(case.saveat, bt ./ (bt0 + epsT)) / case.horizon_days
        w = case.clinical_weights
        s = case.clinical_scales
        return w.tox_peak_mean * (tox_peak_mean / s.tox_peak_mean) +
            w.tox_peak_max * (tox_peak_max / s.tox_peak_max) +
            w.tox_auc * (tox_auc / s.tox_auc) +
            w.tumor_terminal * (tumor_terminal / s.tumor_terminal) +
            w.tumor_auc * (tumor_auc / s.tumor_auc)
    elseif case.objective_name == "tracking"
        target = case.target_trajectory
        bt_scale = bt0 + epsT
        bt_rel = smooth_log_ratio.(bt, bt_scale, epsT)
        tbtc = smooth_log_ratio.(target.bt, target.bt0 + 1e-12, 1e-12)
        il6_cur = log1p.(smooth_pos.(il6) ./ target.il6_scale)
        il6_tgt = target.il6_log
        tdbc_cur = log1p.(tdbc ./ target.tdbc_scale)
        tdbc_tgt = target.tdbc_log
        bt_loss = trapz_generic(case.saveat, (bt_rel .- tbtc) .^ 2) / case.horizon_days
        il6_loss = trapz_generic(case.saveat, (il6_cur .- il6_tgt) .^ 2) / case.horizon_days
        tdbc_loss = trapz_generic(case.saveat, (tdbc_cur .- tdbc_tgt) .^ 2) / case.horizon_days
        return 0.50 * bt_loss + 0.30 * il6_loss + 0.20 * tdbc_loss
    elseif case.objective_name == "stepup_paper" || case.objective_name == "rp2d_stepup"
        refs = resolve_stepup_anchor_refs(case)
        first_peak = smooth_first_window_peak(case, case.saveat, il6; horizon = 7.0)
        day84_change = day84_tumor_size_change_pct_value(bt, bt0)
        bt_auc = btumor_auc_0_84d_value(case.saveat, bt)
        safety_term = log((first_peak + epsT) / (T(refs.peak_il6_0_7d_ref) + epsT))
        day84_penalty = (smooth_pos(day84_change - T(refs.day84_tumor_size_change_pct_ref)) / T(10.0))^2
        auc_penalty = smooth_pos(log((bt_auc + epsT) / (T(refs.btumor_auc_0_84d_ref) + epsT)))^2
        return safety_term + T(10.0) * day84_penalty + auc_penalty
    end
    error("Unknown objective $(case.objective_name)")
end

function loss_channel_gradient(case::MethodBenchmarkCase, bt::AbstractVector, il6::AbstractVector, tdbc::AbstractVector, bt0)
    T = promote_type(eltype(bt), eltype(il6), eltype(tdbc), typeof(bt0), Float64)
    epsT = T(1e-12)
    n = length(case.saveat)
    g_bt = zeros(Float64, n)
    g_il6 = zeros(Float64, n)
    g_tdbc = zeros(Float64, n)
    trap_w = trapz_weights(case.saveat)

    if case.objective_name == "simple"
        tox_peak_simple, il6_w = smoothmax_weights(il6, case.tox_tau)
        tumor_terminal = bt[end] / (bt0 + epsT)
        loss = 0.5 * (tox_peak_simple / case.simple_scales.tox) +
            0.5 * (tumor_terminal / case.simple_scales.tumor)
        g_il6 .+= 0.5 / case.simple_scales.tox .* il6_w
        g_bt[end] += 0.5 / (case.simple_scales.tumor * (bt0 + epsT))
        return Float64(loss), g_bt, g_il6, g_tdbc
    elseif case.objective_name == "auc_combo"
        il6_pos = smooth_pos.(il6)
        tox_auc = trapz_generic(case.saveat, il6_pos) / case.horizon_days
        tumor_auc = trapz_generic(case.saveat, bt ./ (bt0 + epsT)) / case.horizon_days
        s = case.clinical_scales
        loss = (tox_auc / s.tox_auc) + (tumor_auc / s.tumor_auc)
        for i in eachindex(il6)
            g_il6[i] += trap_w[i] * smooth_pos_deriv(il6[i]) / (case.horizon_days * s.tox_auc)
        end
        for i in eachindex(bt)
            g_bt[i] += trap_w[i] / (case.horizon_days * s.tumor_auc * (bt0 + epsT))
        end
        return Float64(loss), g_bt, g_il6, g_tdbc
    elseif case.objective_name == "remission_il6"
        remission_target = exp.((-log(T(2.0)) / T(42.0)) .* case.saveat)
        bt_scale = bt0 + epsT
        bt_ratio = bt ./ bt_scale
        bt_ratio_pos = smooth_pos.(bt_ratio)
        tumor_gap = log.(bt_ratio_pos .+ epsT) .- log.(remission_target .+ epsT)
        tumor_pp = smooth_pos.(tumor_gap)
        il6_scale = T(case.remission_refs.il6_ref)
        il6_pos = smooth_pos.(il6)
        il6_weight = one(T) .+ T(4.0) .* exp.(-case.saveat ./ T(2.0))
        il6_log = log1p.(il6_pos ./ il6_scale)
        lb = trapz_generic(case.saveat, tumor_pp .^ 2) / case.horizon_days
        li = trapz_generic(case.saveat, il6_weight .* il6_log .^ 2) / case.horizon_days
        loss = (lb / (case.remission_refs.lb_paper + 1e-8)) +
            (li / (case.remission_refs.li_paper + 1e-8))
        lb_scale = case.horizon_days * (case.remission_refs.lb_paper + 1e-8)
        li_scale = case.horizon_days * (case.remission_refs.li_paper + 1e-8)
        for i in eachindex(bt)
            d_gap = smooth_pos_deriv(tumor_gap[i]) * smooth_pos_deriv(bt_ratio[i]) / ((bt_ratio_pos[i] + epsT) * bt_scale)
            g_bt[i] += trap_w[i] * 2 * tumor_pp[i] * d_gap / lb_scale
        end
        for i in eachindex(il6)
            d_log = smooth_pos_deriv(il6[i]) / (il6_scale + il6_pos[i])
            g_il6[i] += trap_w[i] * il6_weight[i] * 2 * il6_log[i] * d_log / li_scale
        end
        return Float64(loss), g_bt, g_il6, g_tdbc
    elseif case.objective_name == "clinical"
        window_peaks = Float64[]
        window_weight_maps = Vector{Vector{Tuple{Int, Float64}}}()
        for td in case.dose_times_days
            idxs = [i for i in eachindex(case.saveat) if case.saveat[i] >= td && case.saveat[i] <= td + case.tox_window_days + 1e-10]
            vals = il6[idxs]
            peak_val, peak_w = smoothmax_weights(vals, case.tox_tau)
            push!(window_peaks, Float64(peak_val))
            push!(window_weight_maps, [(idxs[j], peak_w[j]) for j in eachindex(idxs)])
        end
        tox_peak_mean = sum(window_peaks) / length(window_peaks)
        tox_peak_max, peak_idx = findmax(window_peaks)
        il6_pos = smooth_pos.(il6)
        tox_auc = trapz_generic(case.saveat, il6_pos) / case.horizon_days
        tumor_terminal = bt[end] / (bt0 + epsT)
        tumor_auc = trapz_generic(case.saveat, bt ./ (bt0 + epsT)) / case.horizon_days
        w = case.clinical_weights
        s = case.clinical_scales
        loss = w.tox_peak_mean * (tox_peak_mean / s.tox_peak_mean) +
            w.tox_peak_max * (tox_peak_max / s.tox_peak_max) +
            w.tox_auc * (tox_auc / s.tox_auc) +
            w.tumor_terminal * (tumor_terminal / s.tumor_terminal) +
            w.tumor_auc * (tumor_auc / s.tumor_auc)

        mean_scale = w.tox_peak_mean / (s.tox_peak_mean * length(window_peaks))
        max_scale = w.tox_peak_max / s.tox_peak_max
        for map in window_weight_maps
            for (idx, wt) in map
                g_il6[idx] += mean_scale * wt
            end
        end
        for (idx, wt) in window_weight_maps[peak_idx]
            g_il6[idx] += max_scale * wt
        end
        for i in eachindex(il6)
            g_il6[i] += (w.tox_auc / s.tox_auc) * trap_w[i] * smooth_pos_deriv(il6[i]) / case.horizon_days
        end
        g_bt[end] += w.tumor_terminal / (s.tumor_terminal * (bt0 + epsT))
        for i in eachindex(bt)
            g_bt[i] += (w.tumor_auc / s.tumor_auc) * trap_w[i] / (case.horizon_days * (bt0 + epsT))
        end
        return Float64(loss), g_bt, g_il6, g_tdbc
    elseif case.objective_name == "tracking"
        target = case.target_trajectory
        bt_scale = bt0 + epsT
        bt_ratio = bt ./ bt_scale
        bt_ratio_pos = smooth_pos.(bt_ratio)
        bt_rel = log.(bt_ratio_pos .+ epsT)
        tbtc = target.bt_log
        il6_pos = smooth_pos.(il6)
        il6_cur = log1p.(il6_pos ./ target.il6_scale)
        il6_tgt = target.il6_log
        tdbc_cur = log1p.(tdbc ./ target.tdbc_scale)
        tdbc_tgt = target.tdbc_log
        bt_loss = trapz_generic(case.saveat, (bt_rel .- tbtc) .^ 2) / case.horizon_days
        il6_loss = trapz_generic(case.saveat, (il6_cur .- il6_tgt) .^ 2) / case.horizon_days
        tdbc_loss = trapz_generic(case.saveat, (tdbc_cur .- tdbc_tgt) .^ 2) / case.horizon_days
        loss = 0.50 * bt_loss + 0.30 * il6_loss + 0.20 * tdbc_loss

        for i in eachindex(bt)
            d_bt_rel = smooth_pos_deriv(bt_ratio[i]) / ((bt_ratio_pos[i] + epsT) * bt_scale)
            g_bt[i] += 0.50 * trap_w[i] * 2 * (bt_rel[i] - tbtc[i]) * d_bt_rel / case.horizon_days
        end
        for i in eachindex(il6)
            d_il6 = smooth_pos_deriv(il6[i]) / (target.il6_scale + il6_pos[i])
            g_il6[i] += 0.30 * trap_w[i] * 2 * (il6_cur[i] - il6_tgt[i]) * d_il6 / case.horizon_days
        end
        for i in eachindex(tdbc)
            d_tdbc = 1 / (target.tdbc_scale + tdbc[i])
            g_tdbc[i] += 0.20 * trap_w[i] * 2 * (tdbc_cur[i] - tdbc_tgt[i]) * d_tdbc / case.horizon_days
        end
        return Float64(loss), g_bt, g_il6, g_tdbc
    elseif case.objective_name == "stepup_paper" || case.objective_name == "rp2d_stepup"
        refs = resolve_stepup_anchor_refs(case)
        win_idxs = [i for i in eachindex(case.saveat) if case.saveat[i] <= 7.0 + 1e-8]
        peak_val, peak_w = smoothmax_weights(il6[win_idxs], case.tox_tau)
        safety_term = log((peak_val + epsT) / (T(refs.peak_il6_0_7d_ref) + epsT))
        for (local_idx, idx) in enumerate(win_idxs)
            g_il6[idx] += peak_w[local_idx] / (peak_val + epsT)
        end

        day84_change = day84_tumor_size_change_pct_value(bt, bt0)
        day84_gap = smooth_pos(day84_change - T(refs.day84_tumor_size_change_pct_ref))
        day84_deriv = smooth_pos_deriv(day84_change - T(refs.day84_tumor_size_change_pct_ref))
        day84_penalty = (day84_gap / T(10.0))^2
        g_bt[end] += T(10.0) * (2 * day84_gap * day84_deriv / (bt0 + epsT))

        bt_pos = smooth_pos.(bt)
        bt_auc = trapz_until_generic(case.saveat, bt_pos, 84.0)
        auc_log_ratio = log((bt_auc + epsT) / (T(refs.btumor_auc_0_84d_ref) + epsT))
        auc_gap = smooth_pos(auc_log_ratio)
        auc_gap_deriv = smooth_pos_deriv(auc_log_ratio)
        auc_penalty = auc_gap^2
        for i in eachindex(bt)
            if case.saveat[i] <= 84.0 + 1e-8
                d_auc = trap_w[i] * smooth_pos_deriv(bt[i])
                g_bt[i] += 2 * auc_gap * auc_gap_deriv * d_auc / (bt_auc + epsT)
            end
        end

        loss = safety_term + T(10.0) * day84_penalty + auc_penalty
        return Float64(loss), g_bt, g_il6, g_tdbc
    end
    error("Unknown objective $(case.objective_name)")
end

function loss_from_packed_channels(case::MethodBenchmarkCase, y::AbstractVector, bt0)
    bt, il6, tdbc = unpack_channels(y, length(case.saveat))
    return loss_from_channels(case, bt, il6, tdbc, bt0)
end

function parameterized_dose_solve(case::MethodBenchmarkCase, decision_doses::AbstractVector; sensealg = nothing)
    nbase = length(case.pvec)
    pfull = vcat(case.pvec, decision_doses)
    Tactive = promote_type(Float64, eltype(pfull))
    base_u0 = Tactive.(MMC.initial_state_vector(case.pvec))
    active_rates = zeros(Tactive, MMC.DYNAMIC_STATE_COUNT)
    target_idx = case.u_idx[:TDBc_ugperkg]
    dose_scale = 1000.0 / case.bw_kg
    event_specs = [(case.dose_times_days[ev_idx], nbase + group_idx) for (group_idx, grp) in enumerate(case.decision_groups) for ev_idx in grp]
    event_times = [spec[1] for spec in event_specs]
    event_param_idxs = [spec[2] for spec in event_specs]
    t0_delta = sum((isapprox(t, 0.0; atol = 1e-8, rtol = 0.0) ? pfull[pidx] * dose_scale : zero(Tactive)) for (t, pidx) in event_specs; init = zero(Tactive))
    u0 = [i == target_idx ? base_u0[i] + t0_delta : base_u0[i] for i in eachindex(base_u0)]
    rhs = let nbase_local = nbase, active_rates_ref = active_rates
        (du, u, p, t) -> begin
            MMC.mosun_rhs_vector!(du, u, view(p, 1:nbase_local), t)
            @inbounds for i in eachindex(du, active_rates_ref)
                du[i] += active_rates_ref[i]
            end
            MMC.sanitize_ad_vector!(du)
            return nothing
        end
    end
    prob = ODEProblem(rhs, u0, (0.0, case.horizon_days), pfull)
    cb_times = sort(unique([t for t in event_times if !isapprox(t, 0.0; atol = 1e-8, rtol = 0.0)]))
    callback = nothing
    if !isempty(cb_times)
        function affect!(integrator)
            @inbounds for i in eachindex(event_times)
                if isapprox(integrator.t, event_times[i]; atol = 1e-8, rtol = 0.0)
                    integrator.u[target_idx] += oftype(integrator.u[target_idx], integrator.p[event_param_idxs[i]] * dose_scale)
                end
            end
            if !isnothing(case.post_event_proposed_dt) && applicable(SciMLBase.set_proposed_dt!, integrator, case.post_event_proposed_dt)
                SciMLBase.set_proposed_dt!(integrator, case.post_event_proposed_dt)
            end
            return nothing
        end
        callback = PresetTimeCallback(cb_times, affect!; save_positions = (false, false))
    end
    solve_kwargs = (
        abstol = case.abstol,
        reltol = case.reltol,
        callback = callback,
        tstops = cb_times,
        d_discontinuities = cb_times,
        saveat = case.saveat,
        save_everystep = false,
        maxiters = case.maxiters,
    )
    sol = isnothing(sensealg) ?
        solve(prob, case.alg; solve_kwargs...) :
        solve(prob, case.alg; solve_kwargs..., sensealg = sensealg)
    sol.retcode == SciMLBase.ReturnCode.Success || error("solve failed with retcode=$(sol.retcode)")
    return (; sol, pfull, initial_u = u0, callback)
end

function objective_channels(case::MethodBenchmarkCase, sol, p_base::AbstractVector)
    bt_idx = case.u_idx[:Btumor]
    bt = [u[bt_idx] for u in sol.u]
    il6 = [il6combo(u, p_base, case.u_idx, case.p_idx) for u in sol.u]
    tdbc = [tdbc_ugperml(u, p_base, case.u_idx, case.p_idx) for u in sol.u]
    return bt, il6, tdbc
end

function grouped_bolus_doses_ugkg(case::MethodBenchmarkCase, decision_doses::AbstractVector)
    T = eltype(decision_doses)
    full = expand_decision_doses_generic(decision_doses, case.decision_groups, length(case.dose_times_days))
    unique_times = sort(unique(case.dose_times_days))
    mg_to_ugkg = T(1000.0 / case.bw_kg)
    bolus_ugkg = [
        mg_to_ugkg * sum(
            (isapprox(case.dose_times_days[j], t; atol = 1e-8, rtol = 0.0) ? full[j] : zero(T)) for
            j in eachindex(full);
            init = zero(T),
        ) for t in unique_times
    ]
    return unique_times, bolus_ugkg
end

function solve_segment_vector(case::MethodBenchmarkCase, u0::AbstractVector, tspan::Tuple{Float64, Float64}, local_saveat; sensealg = nothing)
    rhs = (du, u, p, t) -> begin
        MMC.mosun_rhs_vector!(du, u, p, t)
        MMC.sanitize_ad_vector!(du)
        return nothing
    end
    prob = ODEProblem(rhs, u0, tspan, case.pvec)
    solve_kwargs = (
        abstol = case.abstol,
        reltol = case.reltol,
        saveat = local_saveat,
        save_everystep = false,
        tstops = [tspan[2]],
        maxiters = case.maxiters,
    )
    sol = isnothing(sensealg) ?
        solve(prob, case.alg; solve_kwargs...) :
        solve(prob, case.alg; solve_kwargs..., sensealg = sensealg)
    sol.retcode == SciMLBase.ReturnCode.Success || error("segment solve failed with retcode=$(sol.retcode) over tspan=$(tspan)")
    return sol
end

function segmented_channels(case::MethodBenchmarkCase, decision_doses::AbstractVector; sensealg = nothing)
    event_times, bolus_ugkg = grouped_bolus_doses_ugkg(case, decision_doses)
    T = promote_type(Float64, eltype(decision_doses))
    base_u0 = T.(MMC.initial_state_vector(case.pvec))
    target_idx = case.u_idx[:TDBc_ugperkg]
    bt_idx = case.u_idx[:Btumor]
    t0_delta = sum(
        (isapprox(event_times[i], 0.0; atol = 1e-8, rtol = 0.0) ? bolus_ugkg[i] : zero(T)) for
        i in eachindex(event_times);
        init = zero(T),
    )
    u0 = [i == target_idx ? base_u0[i] + t0_delta : base_u0[i] for i in eachindex(base_u0)]
    bt0 = u0[bt_idx]
    positive_event_times = [t for t in event_times if t > 1e-8]
    segment_ends = vcat(positive_event_times, [case.horizon_days])
    dose_at_time = Dict{Float64, Any}(event_times[i] => bolus_ugkg[i] for i in eachindex(event_times))

    function recurse(seg_idx::Int, u_curr::AbstractVector, t_start::Float64)
        t_end = segment_ends[seg_idx]
        local_saveat = [
            t for t in case.saveat if
            ((seg_idx == 1 ? t >= t_start - 1e-8 : t > t_start + 1e-8) && t <= t_end + 1e-8)
        ]
        isempty(local_saveat) && (local_saveat = [t_end])
        sol = solve_segment_vector(case, u_curr, (t_start, t_end), local_saveat; sensealg = sensealg)
        bt = [u[bt_idx] for u in sol.u]
        il6 = [il6combo(u, case.pvec, case.u_idx, case.p_idx) for u in sol.u]
        tdbc = [tdbc_ugperml(u, case.pvec, case.u_idx, case.p_idx) for u in sol.u]
        if seg_idx == length(segment_ends)
            return (; t = sol.t, bt, il6, tdbc)
        end
        delta_next = haskey(dose_at_time, t_end) ? dose_at_time[t_end] : zero(T)
        u_next = [
            i == target_idx ? sol.u[end][i] + delta_next : sol.u[end][i] for
            i in eachindex(sol.u[end])
        ]
        rest = recurse(seg_idx + 1, u_next, t_end)
        return (
            t = vcat(sol.t[1:(end - 1)], rest.t),
            bt = vcat(bt[1:(end - 1)], rest.bt),
            il6 = vcat(il6[1:(end - 1)], rest.il6),
            tdbc = vcat(tdbc[1:(end - 1)], rest.tdbc),
        )
    end

    chans = recurse(1, u0, 0.0)
    length(chans.t) == length(case.saveat) || error("segmented channel length mismatch: got $(length(chans.t)) expected $(length(case.saveat))")
    return (; chans.t, chans.bt, chans.il6, chans.tdbc, bt0)
end

function objective_value_segmented(case::MethodBenchmarkCase, decision_doses::AbstractVector, sensealg)
    chans = segmented_channels(case, decision_doses; sensealg = sensealg)
    return loss_from_channels(case, chans.bt, chans.il6, chans.tdbc, chans.bt0)
end

function segmented_forward_data(case::MethodBenchmarkCase, decision_doses::AbstractVector)
    event_times, bolus_ugkg = grouped_bolus_doses_ugkg(case, decision_doses)
    target_idx = case.u_idx[:TDBc_ugperkg]
    bt_idx = case.u_idx[:Btumor]
    base_u0 = Float64.(MMC.initial_state_vector(case.pvec))
    t0_delta = sum(
        (isapprox(event_times[i], 0.0; atol = 1e-8, rtol = 0.0) ? Float64(bolus_ugkg[i]) : 0.0) for
        i in eachindex(event_times);
        init = 0.0,
    )
    u_curr = copy(base_u0)
    u_curr[target_idx] += t0_delta
    bt0 = u_curr[bt_idx]
    positive_event_times = [Float64(t) for t in event_times if t > 1e-8]
    segment_starts = vcat([0.0], positive_event_times)
    segment_ends = vcat(positive_event_times, [case.horizon_days])
    dose_by_time = Dict{Float64, Float64}(Float64(event_times[i]) => Float64(bolus_ugkg[i]) for i in eachindex(event_times))

    sols = Vector{Any}(undef, length(segment_ends))
    seg_ranges = Vector{UnitRange{Int}}(undef, length(segment_ends))
    bt = Float64[]
    il6 = Float64[]
    tdbc = Float64[]
    cursor = 1

    for seg_idx in eachindex(segment_ends)
        t_start = segment_starts[seg_idx]
        t_end = segment_ends[seg_idx]
        local_saveat = [
            t for t in case.saveat if
            ((seg_idx == 1 ? t >= t_start - 1e-8 : t > t_start + 1e-8) && t <= t_end + 1e-8)
        ]
        isempty(local_saveat) && error("empty local saveat for segment $(seg_idx) over ($(t_start), $(t_end))")
        sol = solve_segment_vector(case, u_curr, (t_start, t_end), local_saveat; sensealg = nothing)
        sols[seg_idx] = sol
        n_local = length(sol.t)
        seg_ranges[seg_idx] = cursor:(cursor + n_local - 1)
        append!(bt, Float64.([u[bt_idx] for u in sol.u]))
        append!(il6, Float64.([il6combo(u, case.pvec, case.u_idx, case.p_idx) for u in sol.u]))
        append!(tdbc, Float64.([tdbc_ugperml(u, case.pvec, case.u_idx, case.p_idx) for u in sol.u]))
        cursor += n_local
        if seg_idx < length(segment_ends)
            u_curr = Float64.(sol.u[end])
            u_curr[target_idx] += dose_by_time[t_end]
        end
    end

    length(bt) == length(case.saveat) || error("segmented forward length mismatch: got $(length(bt)) expected $(length(case.saveat))")
    return (; sols, seg_ranges, segment_starts, bt, il6, tdbc, bt0)
end

function reverse_adjoint_dose_gradient(case::MethodBenchmarkCase, x::Vector{Float64}, sensealg)
    fwd = segmented_forward_data(case, x)
    loss_val, gy_bt, gy_il6, gy_tdbc = loss_channel_gradient(case, fwd.bt, fwd.il6, fwd.tdbc, fwd.bt0)
    p_base = case.pvec
    il6_tiss_coeff = p_base[case.p_idx[:IL6_tiss_contribution]]
    vpb = p_base[case.p_idx[:Vpb]]
    coeff_tiss = il6_tiss_coeff * p_base[case.p_idx[:Vtissue]] / vpb
    coeff_tiss2 = il6_tiss_coeff * p_base[case.p_idx[:Vtissue2]] / vpb
    coeff_tiss3 = il6_tiss_coeff * p_base[case.p_idx[:Vtissue3]] / vpb
    coeff_tumor = il6_tiss_coeff * p_base[case.p_idx[:Vtumor]] / vpb
    il6pb_idx = case.u_idx[:IL6pb]
    il6tiss_idx = case.u_idx[:IL6tiss]
    il6tiss2_idx = case.u_idx[:IL6tiss2]
    il6tiss3_idx = case.u_idx[:IL6tiss3]
    il6tumor_idx = case.u_idx[:IL6tumor]
    bt_idx = case.u_idx[:Btumor]
    tdb_idx = case.u_idx[:TDBc_ugperkg]
    vc_tdb = p_base[case.p_idx[:Vc_tdb]]
    lambda_next = zeros(Float64, MMC.DYNAMIC_STATE_COUNT)
    lambda_by_start = Dict{Float64, Vector{Float64}}()

    for seg_idx in length(fwd.sols):-1:1
        sol = fwd.sols[seg_idx]
        range = fwd.seg_ranges[seg_idx]
        offset = first(range) - 1
        terminal_lambda = copy(lambda_next)
        function dg(out, u, p_local, t, i)
            fill!(out, 0.0)
            gi = offset + i
            out[bt_idx] += gy_bt[gi]
            d_il6 = gy_il6[gi]
            out[il6pb_idx] += d_il6
            out[il6tiss_idx] += d_il6 * coeff_tiss
            out[il6tiss2_idx] += d_il6 * coeff_tiss2
            out[il6tiss3_idx] += d_il6 * coeff_tiss3
            out[il6tumor_idx] += d_il6 * coeff_tumor
            if (u[tdb_idx] / vc_tdb) > 1e-5
                out[tdb_idx] += gy_tdbc[gi] / vc_tdb
            end
            if i == length(sol.t)
                out .+= terminal_lambda
            end
            return nothing
        end

        du0, _ = adjoint_sensitivities(
            sol,
            case.alg;
            t = sol.t,
            dgdu_discrete = dg,
            sensealg = sensealg,
            abstol = case.abstol,
            reltol = case.reltol,
        )
        lambda_next = Float64.(du0)
        lambda_by_start[fwd.segment_starts[seg_idx]] = copy(lambda_next)
    end

    mg_to_ugkg = 1000.0 / case.bw_kg
    event_grads = [
        lambda_by_start[Float64(case.dose_times_days[i])][tdb_idx] * mg_to_ugkg for
        i in eachindex(case.dose_times_days)
    ]
    decision_grads = [sum(event_grads[idx] for idx in grp) for grp in case.decision_groups]
    return loss_val, decision_grads
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
    il6_pos = smooth_pos.(il6)
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
        tox_auc = trapz_generic(t, il6_pos) / case.horizon_days,
        tumor_terminal = bt[end] / (bt0 + eltype(bt0)(1e-12)),
        tumor_auc = trapz_generic(t, bt ./ (bt0 + eltype(bt0)(1e-12))) / case.horizon_days,
        auc_tdbc = sol.u[end][auc_idx],
    )
end

function simple_objective(case::MethodBenchmarkCase, m)
    return 0.5 * (m.tox_peak_simple / case.simple_scales.tox) +
        0.5 * (m.tumor_terminal / case.simple_scales.tumor)
end

function auc_combo_objective(case::MethodBenchmarkCase, m)
    s = case.clinical_scales
    return (m.tox_auc / s.tox_auc) + (m.tumor_auc / s.tumor_auc)
end

function remission_il6_components(case::MethodBenchmarkCase, t::AbstractVector, bt::AbstractVector, il6::AbstractVector, bt0; il6_ref = case.remission_refs.il6_ref)
    T = promote_type(eltype(t), eltype(bt), eltype(il6), typeof(bt0), Float64)
    eps_b = T(1e-12)
    horizon = T(case.horizon_days)
    bt_scale = bt0 + eps_b
    bt_ratio = bt ./ bt_scale
    remission_target = exp.((-log(T(2.0)) / T(42.0)) .* t)
    tumor_gap = smooth_pos.(log.(smooth_pos.(bt_ratio) .+ eps_b) .- log.(remission_target .+ eps_b))
    lb = trapz_generic(t, tumor_gap .^ 2) / horizon

    il6_scale = T(il6_ref)
    il6_pos = smooth_pos.(il6)
    il6_weight = one(T) .+ T(4.0) .* exp.(-t ./ T(2.0))
    li = trapz_generic(t, il6_weight .* log1p.(il6_pos ./ il6_scale) .^ 2) / horizon
    return lb, li
end

function remission_il6_objective(case::MethodBenchmarkCase, m)
    lb, li = remission_il6_components(case, m.t, m.bt, m.il6, m.bt0)
    return (lb / (case.remission_refs.lb_paper + 1e-8)) +
        (li / (case.remission_refs.li_paper + 1e-8))
end

function stepup_paper_objective(case::MethodBenchmarkCase, m)
    return stepup_paper_terms(case, m).total
end

function rp2d_stepup_objective(case::MethodBenchmarkCase, m)
    return rp2d_stepup_terms(case, m).total
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
    bt_rel = smooth_log_ratio.(m.bt, m.bt0 + eltype(m.bt0)(1e-12), eltype(m.bt0)(1e-12))
    tbtc = smooth_log_ratio.(target.bt, target.bt0 + 1e-12, 1e-12)
    il6_cur = log1p.(smooth_pos.(m.il6) ./ target.il6_scale)
    il6_tgt = target.il6_log
    tdbc_cur = log1p.(m.tdbc ./ target.tdbc_scale)
    tdbc_tgt = target.tdbc_log
    bt_loss = trapz_generic(target.t, (bt_rel .- tbtc) .^ 2) / case.horizon_days
    il6_loss = trapz_generic(target.t, (il6_cur .- il6_tgt) .^ 2) / case.horizon_days
    tdbc_loss = trapz_generic(target.t, (tdbc_cur .- tdbc_tgt) .^ 2) / case.horizon_days
    return 0.50 * bt_loss + 0.30 * il6_loss + 0.20 * tdbc_loss
end

function make_patient_ideal_target(case::MethodBenchmarkCase, hosseini_metrics, rp2d_metrics)
    target_shape = lowercase(strip(get(ENV, "PATIENT_IDEAL_TUMOR_TARGET_SHAPE", "sigmoid_log")))
    r_inf = env_float_bench("PATIENT_IDEAL_TUMOR_R_INF", 0.001)
    tau_days = env_float_bench("PATIENT_IDEAL_TUMOR_TAU_DAYS", 14.0)
    sigmoid_t50_days = env_float_bench("PATIENT_IDEAL_TUMOR_SIGMOID_T50_DAYS", 14.0)
    sigmoid_width_days = env_float_bench("PATIENT_IDEAL_TUMOR_SIGMOID_WIDTH_DAYS", 1.25)
    il6_fraction = env_float_bench("PATIENT_IDEAL_IL6_CEILING_FRACTION", 0.75)
    il6_floor = env_float_bench("PATIENT_IDEAL_IL6_FLOOR_PGML", 1.0)
    t = Float64.(hosseini_metrics.t)
    r_floor = max(r_inf, 1e-12)
    target_ratio = if target_shape in ("sigmoid_log", "sigmoidal_log", "logistic_log", "sigmoid")
        width = max(sigmoid_width_days, 1e-6)
        s = 1.0 ./ (1.0 .+ exp.(-((t .- sigmoid_t50_days) ./ width)))
        s0 = 1.0 / (1.0 + exp(sigmoid_t50_days / width))
        progress = clamp.((s .- s0) ./ max(1.0 - s0, 1e-12), 0.0, 1.0)
        10.0 .^ (log10(r_floor) .* progress)
    elseif target_shape in ("exponential", "exp", "exponential_saturation")
        r_floor .+ (1.0 - r_floor) .* exp.(-t ./ max(tau_days, 1e-6))
    else
        error("Unsupported PATIENT_IDEAL_TUMOR_TARGET_SHAPE=$(target_shape)")
    end
    hosseini_il6 = max.(Float64.(hosseini_metrics.il6), 0.0)
    il6_ceiling = max.(il6_floor, il6_fraction .* hosseini_il6)
    first_idxs = [i for i in eachindex(t) if t[i] <= 7.0 + 1e-8]
    first_ceiling = isempty(first_idxs) ? maximum(il6_ceiling) : maximum(il6_ceiling[first_idxs])
    return (
        name = "patient_ideal_tracking",
        label = "Patient-specific ideal tracking target",
        t = t,
        tumor_target_ratio = target_ratio,
        tumor_target_auc = trapz_generic(t, target_ratio) / case.horizon_days,
        tumor_target_shape = target_shape,
        tumor_r_inf = r_inf,
        tumor_tau_days = tau_days,
        tumor_sigmoid_t50_days = sigmoid_t50_days,
        tumor_sigmoid_width_days = sigmoid_width_days,
        il6_ceiling = il6_ceiling,
        il6_ceiling_fraction = il6_fraction,
        il6_floor_pgml = il6_floor,
        first_peak_il6_ceiling = max(first_ceiling, il6_floor),
        global_peak_il6_ceiling = max(maximum(il6_ceiling), il6_floor),
        drug_auc_scale = max(max(Float64(hosseini_metrics.auc_tdbc), Float64(rp2d_metrics.auc_tdbc)), 1e-6),
    )
end

function achievable_ideal_target_cache_dir()
    raw = get(ENV, "ACHIEVABLE_IDEAL_TARGET_CACHE_DIR", joinpath("generated", "server_results", "hosseini_efast_vpop250_achievable_ideal_tracking_20260515", "target_cache"))
    return isabspath(raw) ? raw : joinpath(REPO_ROOT, raw)
end

function load_achievable_target_curves()
    cached = ACHIEVABLE_TARGET_CURVES_CACHE[]
    if cached !== nothing
        return cached
    end
    path = joinpath(achievable_ideal_target_cache_dir(), "achievable_ideal_target_curves.csv")
    isfile(path) || error("Missing achievable ideal target curves: $(path). Run julia/export_hosseini_efast_vpop250_achievable_ideal_targets.jl first.")
    df = CSV.read(path, DataFrame)
    required = ["sample_id", "time_day", "tumor_target_ratio", "il6_target_pgml"]
    missing_cols = setdiff(required, names(df))
    isempty(missing_cols) || error("Target curves file $(path) is missing columns: $(join(missing_cols, ", "))")
    ACHIEVABLE_TARGET_CURVES_CACHE[] = df
    return df
end

function load_achievable_target_summary()
    cached = ACHIEVABLE_TARGET_SUMMARY_CACHE[]
    if cached !== nothing
        return cached
    end
    path = joinpath(achievable_ideal_target_cache_dir(), "achievable_ideal_target_summary.csv")
    isfile(path) || error("Missing achievable ideal target summary: $(path). Run julia/export_hosseini_efast_vpop250_achievable_ideal_targets.jl first.")
    df = CSV.read(path, DataFrame)
    "sample_id" in names(df) || error("Target summary file $(path) is missing sample_id")
    ACHIEVABLE_TARGET_SUMMARY_CACHE[] = df
    return df
end

function achievable_global_best_sample_id()
    summary = load_achievable_target_summary()
    if "global_best_target" in names(summary)
        flagged = summary[Bool.(summary.global_best_target), :]
        nrow(flagged) == 1 || error("Expected exactly one global_best_target row, found $(nrow(flagged))")
        return Int(flagged.sample_id[1])
    end
    required = ["global_rank_sum", "sample_id"]
    missing_cols = setdiff(required, names(summary))
    isempty(missing_cols) || error("Target summary missing global-best columns: $(join(missing_cols, ", "))")
    sort!(summary, [:global_rank_sum, :sample_id])
    return Int(summary.sample_id[1])
end

function source_summary_row(summary::DataFrame, source_sample_id::Int)
    rows = summary[summary.sample_id .== source_sample_id, :]
    nrow(rows) == 1 || error("Expected one target summary row for sample_id=$(source_sample_id), found $(nrow(rows))")
    return rows[1, :]
end

function summary_float(row, col::String, default::Float64)
    if col in names(parent(row))
        val = row[Symbol(col)]
        if !ismissing(val)
            f = Float64(val)
            isfinite(f) && return f
        end
    end
    return default
end

function make_achievable_ideal_target(case::MethodBenchmarkCase, objective_name::String)
    objective = lowercase(strip(objective_name))
    target_mode = objective == "global_best_achievable_ideal_tracking" ? "global_best" : "patient_specific"
    source_sample_id = target_mode == "global_best" ? achievable_global_best_sample_id() : case.sample_id
    curves = load_achievable_target_curves()
    sub = curves[curves.sample_id .== source_sample_id, :]
    nrow(sub) > 0 || error("No achievable target curve rows for source_sample_id=$(source_sample_id)")
    sort!(sub, :time_day)
    t = Float64.(sub.time_day)
    tumor_target_ratio = max.(Float64.(sub.tumor_target_ratio), 1e-12)
    il6_target = max.(Float64.(sub.il6_target_pgml), 0.0)
    first_idxs = [i for i in eachindex(t) if t[i] <= 7.0 + 1e-8]
    first_ceiling = isempty(first_idxs) ? maximum(il6_target) : maximum(il6_target[first_idxs])
    summary = load_achievable_target_summary()
    src = source_summary_row(summary, source_sample_id)
    max_drug_auc = summary_float(src, "max_dose_drug_auc", 1.0)
    min_drug_auc = summary_float(src, "min_1mg_drug_auc", 1.0)
    tumor_curve_scale = summary_float(src, "calibration_tumor_curve_scale", 1.0)
    tumor_auc_scale = summary_float(src, "calibration_tumor_auc_scale", 1.0)
    il6_curve_scale = summary_float(src, "calibration_il6_curve_scale", 1.0)
    first_peak_scale = summary_float(src, "calibration_first_peak_scale", 1.0)
    global_peak_scale = summary_float(src, "calibration_global_peak_scale", 1.0)
    target_day84_tumor_ratio = summary_float(src, "max_dose_day84_tumor_ratio", tumor_target_ratio[end])
    min_day84_tumor_ratio = summary_float(src, "min_1mg_day84_tumor_ratio", 1.0)
    day84_tumor_scale = max(smooth_pos(log((min_day84_tumor_ratio + 1e-12) / (target_day84_tumor_ratio + 1e-12)))^2, 1e-9)
    return (
        name = objective,
        label = target_mode == "global_best" ? "Global best achievable max-tumor/1-mg-IL6 tracking target" : "Patient-specific achievable max-tumor/1-mg-IL6 tracking target",
        target_mode = target_mode,
        target_source_sample_id = source_sample_id,
        t = t,
        tumor_target_ratio = tumor_target_ratio,
        tumor_target_auc = trapz_generic(t, tumor_target_ratio) / case.horizon_days,
        target_day84_tumor_ratio = max(target_day84_tumor_ratio, 1e-12),
        tumor_target_shape = "raw_max_dose_trajectory",
        tumor_r_inf = 0.0,
        tumor_tau_days = 0.0,
        tumor_sigmoid_t50_days = 0.0,
        tumor_sigmoid_width_days = 0.0,
        il6_ceiling = il6_target,
        il6_ceiling_fraction = 1.0,
        il6_floor_pgml = env_float_bench("ACHIEVABLE_IDEAL_IL6_FLOOR_PGML", 1.0),
        first_peak_il6_ceiling = max(first_ceiling, 1e-12),
        global_peak_il6_ceiling = max(maximum(il6_target), 1e-12),
        drug_auc_scale = max(max(max_drug_auc, min_drug_auc), 1e-6),
        calibration_tumor_curve_scale = max(tumor_curve_scale, 1e-9),
        calibration_tumor_auc_scale = max(tumor_auc_scale, 1e-9),
        calibration_il6_curve_scale = max(il6_curve_scale, 1e-9),
        calibration_first_peak_scale = max(first_peak_scale, 1e-9),
        calibration_global_peak_scale = max(global_peak_scale, 1e-9),
        calibration_day84_tumor_scale = max(day84_tumor_scale, 1e-9),
        tumor_target_schedule_mg = "60/60/60/60/60/60",
        il6_target_schedule_mg = "1/1/1/1/1/1",
    )
end

function patient_ideal_tracking_terms(case::MethodBenchmarkCase, m, decision_doses::AbstractVector)
    target = case.target_trajectory
    T = promote_type(eltype(m.t), eltype(m.bt), eltype(m.il6), typeof(m.bt0), eltype(decision_doses), Float64)
    epsT = T(1e-12)
    horizon = T(case.horizon_days)

    bt_ratio = smooth_pos.(m.bt ./ (m.bt0 + epsT))
    target_ratio = T.(target.tumor_target_ratio)
    tumor_gap = smooth_pos.(log.(bt_ratio .+ epsT) .- log.(target_ratio .+ epsT))
    early_weight = one(T) .+ T(2.0) .* exp.(-T.(target.t) ./ T(21.0))
    tumor_curve_loss = trapz_generic(target.t, early_weight .* tumor_gap .^ 2) / max(trapz_generic(target.t, early_weight), epsT)

    tumor_auc_ratio = trapz_generic(m.t, bt_ratio) / horizon
    tumor_auc_excess = smooth_pos(log((tumor_auc_ratio + epsT) / (T(target.tumor_target_auc) + epsT)))
    tumor_auc_loss = tumor_auc_excess^2
    target_day84_tumor_ratio = T(get(target, :target_day84_tumor_ratio, target_ratio[end]))
    day84_tumor_ratio = bt_ratio[end]
    day84_tumor_excess = smooth_pos(log((day84_tumor_ratio + epsT) / (target_day84_tumor_ratio + epsT)))
    day84_tumor_loss = day84_tumor_excess^2

    il6_pos = smooth_pos.(m.il6)
    il6_ceiling = T.(target.il6_ceiling)
    il6_curve_excess = smooth_pos.(log.((il6_pos .+ T(target.il6_floor_pgml)) ./ (il6_ceiling .+ T(target.il6_floor_pgml))))
    il6_curve_loss = trapz_generic(target.t, il6_curve_excess .^ 2) / horizon

    first_peak = first_window_peak_value(m.t, m.il6; horizon = 7.0)
    global_peak = maximum(il6_pos)
    first_peak_excess = smooth_pos(log((first_peak + epsT) / (T(target.first_peak_il6_ceiling) + epsT)))
    global_peak_excess = smooth_pos(log((global_peak + epsT) / (T(target.global_peak_il6_ceiling) + epsT)))
    first_peak_loss = first_peak_excess^2
    global_peak_loss = global_peak_excess^2

    total_dose = total_expanded_dose_mg(case, decision_doses)
    dose_penalty = total_dose / T(360.0)
    drug_auc_penalty = log1p(m.auc_tdbc / T(target.drug_auc_scale))^2
    full_doses = expand_decision_doses_generic(decision_doses, case.decision_groups, length(case.dose_times_days))
    cycle1_total = sum(full_doses[1:min(3, length(full_doses))])
    cycle2_total = length(full_doses) >= 4 ? full_doses[4] : zero(T)
    cycle3_total = length(full_doses) >= 5 ? full_doses[5] : zero(T)
    cycle4_total = length(full_doses) >= 6 ? full_doses[6] : zero(T)
    cycle_target_mg = T(env_float_bench("ACHIEVABLE_IDEAL_CYCLE_DOSE_TARGET_MG", 60.0))
    cycle_scale = max(cycle_target_mg, epsT)
    cycle_totals = T[cycle1_total, cycle2_total, cycle3_total, cycle4_total]
    cycle_dose_target_loss = mean(((cycle_totals .- cycle_target_mg) ./ cycle_scale) .^ 2)
    default_cycle_weight = env_float_bench("ACHIEVABLE_IDEAL_CYCLE60_WEIGHT", 1.0)
    cycle_dose_target_weight = is_cycle_target_achievable_objective(case.objective_name) ?
        T(env_float_bench("ACHIEVABLE_IDEAL_CYCLE_TARGET_WEIGHT", default_cycle_weight)) :
        zero(T)

    tumor_curve_scale = max(T(get(target, :calibration_tumor_curve_scale, 1.0)), epsT)
    tumor_auc_scale = max(T(get(target, :calibration_tumor_auc_scale, 1.0)), epsT)
    il6_curve_scale = max(T(get(target, :calibration_il6_curve_scale, 1.0)), epsT)
    first_peak_scale = max(T(get(target, :calibration_first_peak_scale, 1.0)), epsT)
    global_peak_scale = max(T(get(target, :calibration_global_peak_scale, 1.0)), epsT)
    day84_tumor_scale = max(T(get(target, :calibration_day84_tumor_scale, 1.0)), epsT)
    normalized_tumor_curve_loss = tumor_curve_loss / tumor_curve_scale
    normalized_tumor_auc_loss = tumor_auc_loss / tumor_auc_scale
    normalized_day84_tumor_loss = day84_tumor_loss / day84_tumor_scale
    normalized_il6_curve_loss = il6_curve_loss / il6_curve_scale
    normalized_first_peak_loss = first_peak_loss / first_peak_scale
    normalized_global_peak_loss = global_peak_loss / global_peak_scale

    calibrated = is_calibrated_achievable_objective(case.objective_name)
    endpoint_objective = is_endpoint_achievable_objective(case.objective_name)
    tumor_family_loss = endpoint_objective ?
        normalized_day84_tumor_loss :
        calibrated ?
        T(0.5) * normalized_tumor_curve_loss + T(0.5) * normalized_tumor_auc_loss :
        T(2.0) * tumor_curve_loss + T(1.0) * tumor_auc_loss
    il6_family_loss = endpoint_objective ?
        normalized_global_peak_loss :
        calibrated ?
        (normalized_il6_curve_loss + normalized_first_peak_loss + normalized_global_peak_loss) / T(3.0) :
        T(8.0) * il6_curve_loss + T(8.0) * first_peak_loss + T(4.0) * global_peak_loss
    tumor_family_weight = endpoint_objective ?
        T(env_float_bench("ACHIEVABLE_ENDPOINT_DAY84_WEIGHT", 1.0)) :
        calibrated ?
        T(env_float_bench("ACHIEVABLE_IDEAL_TUMOR_FAMILY_WEIGHT", 1.0)) :
        T(env_float_bench("PATIENT_IDEAL_TUMOR_FAMILY_WEIGHT", 1.0))
    il6_family_weight = endpoint_objective ?
        T(env_float_bench("ACHIEVABLE_ENDPOINT_GLOBAL_IL6_WEIGHT", 1.0)) :
        calibrated ?
        T(env_float_bench("ACHIEVABLE_IDEAL_IL6_FAMILY_WEIGHT", 1.0)) :
        T(env_float_bench("PATIENT_IDEAL_IL6_FAMILY_WEIGHT", 1.0))

    total = endpoint_objective ?
        tumor_family_weight * tumor_family_loss + il6_family_weight * il6_family_loss + cycle_dose_target_weight * cycle_dose_target_loss :
        calibrated ?
        tumor_family_weight * tumor_family_loss + il6_family_weight * il6_family_loss + cycle_dose_target_weight * cycle_dose_target_loss :
        tumor_family_weight * tumor_family_loss + il6_family_weight * il6_family_loss + T(0.02) * dose_penalty + T(0.02) * drug_auc_penalty

    return (
        total = total,
        tumor_curve_loss = tumor_curve_loss,
        tumor_auc_loss = tumor_auc_loss,
        day84_tumor_loss = day84_tumor_loss,
        il6_curve_loss = il6_curve_loss,
        first_peak_loss = first_peak_loss,
        global_peak_loss = global_peak_loss,
        normalized_tumor_curve_loss = normalized_tumor_curve_loss,
        normalized_tumor_auc_loss = normalized_tumor_auc_loss,
        normalized_day84_tumor_loss = normalized_day84_tumor_loss,
        normalized_il6_curve_loss = normalized_il6_curve_loss,
        normalized_first_peak_loss = normalized_first_peak_loss,
        normalized_global_peak_loss = normalized_global_peak_loss,
        tumor_family_loss = tumor_family_loss,
        il6_family_loss = il6_family_loss,
        tumor_family_weight = tumor_family_weight,
        il6_family_weight = il6_family_weight,
        dose_penalty = dose_penalty,
        drug_auc_penalty = drug_auc_penalty,
        cycle_dose_target_loss = cycle_dose_target_loss,
        cycle_dose_target_weight = cycle_dose_target_weight,
        cycle_dose_target_mg = cycle_target_mg,
        cycle1_total_mg = cycle1_total,
        cycle2_total_mg = cycle2_total,
        cycle3_total_mg = cycle3_total,
        cycle4_total_mg = cycle4_total,
        first_peak_il6 = first_peak,
        global_peak_il6 = global_peak,
        first_peak_il6_ceiling = T(target.first_peak_il6_ceiling),
        global_peak_il6_ceiling = T(target.global_peak_il6_ceiling),
        tumor_target_auc = T(target.tumor_target_auc),
        target_day84_tumor_ratio = target_day84_tumor_ratio,
        day84_tumor_ratio = day84_tumor_ratio,
        tumor_auc_ratio = tumor_auc_ratio,
        total_dose_mg = total_dose,
        tumor_r_inf = T(target.tumor_r_inf),
        tumor_tau_days = T(target.tumor_tau_days),
        tumor_sigmoid_t50_days = T(get(target, :tumor_sigmoid_t50_days, NaN)),
        tumor_sigmoid_width_days = T(get(target, :tumor_sigmoid_width_days, NaN)),
        il6_ceiling_fraction = T(target.il6_ceiling_fraction),
        il6_floor_pgml = T(target.il6_floor_pgml),
        calibration_tumor_curve_scale = tumor_curve_scale,
        calibration_tumor_auc_scale = tumor_auc_scale,
        calibration_day84_tumor_scale = day84_tumor_scale,
        calibration_il6_curve_scale = il6_curve_scale,
        calibration_first_peak_scale = first_peak_scale,
        calibration_global_peak_scale = global_peak_scale,
    )
end

function patient_ideal_tracking_objective(case::MethodBenchmarkCase, m, decision_doses::AbstractVector)
    return patient_ideal_tracking_terms(case, m, decision_doses).total
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
            elseif case.objective_name == "auc_combo"
                return auc_combo_objective(case, m)
            elseif case.objective_name == "remission_il6"
                return remission_il6_objective(case, m)
            elseif case.objective_name == "clinical"
                return clinical_objective(case, m)
            elseif case.objective_name == "tracking"
                return tracking_objective(case, m)
            elseif is_ideal_tracking_objective(case.objective_name)
                return patient_ideal_tracking_objective(case, m, x)
            elseif case.objective_name == "stepup_paper"
                return stepup_paper_objective(case, m)
            elseif case.objective_name == "rp2d_stepup"
                return rp2d_stepup_objective(case, m)
            elseif case.objective_name == "safety_constrained_ti"
                return safety_constrained_ti_objective(case, m, x)
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
    elseif case.objective_name == "auc_combo"
        return auc_combo_objective(case, m)
    elseif case.objective_name == "remission_il6"
        return remission_il6_objective(case, m)
    elseif case.objective_name == "clinical"
        return clinical_objective(case, m)
    elseif case.objective_name == "tracking"
        return tracking_objective(case, m)
    elseif is_ideal_tracking_objective(case.objective_name)
        return patient_ideal_tracking_objective(case, m, x)
    elseif case.objective_name == "stepup_paper"
        return stepup_paper_objective(case, m)
    elseif case.objective_name == "rp2d_stepup"
        return rp2d_stepup_objective(case, m)
    elseif case.objective_name == "safety_constrained_ti"
        return safety_constrained_ti_objective(case, m, x)
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
    logger::Union{Nothing, TraceLogger} = nothing,
)
    rng = MersenneTwister(rng_seed)
    nd = length(x0)
    x_best = clamp.(copy(x0), lower, upper)
    f_best = Float64(obj(x_best))
    !isnothing(logger) && log_eval!(logger, x_best, f_best)
    n_remaining = max(n_evals, 1) - 1
    while n_remaining > 0
        batch = min(n_remaining, max(1, nthreads()))
        x_batch = [similar(x_best) for _ in 1:batch]
        for b in 1:batch
            x = x_batch[b]
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
        end
        f_batch = zeros(Float64, batch)
        @threads for b in 1:batch
            f_batch[b] = Float64(obj(x_batch[b]))
        end
        for b in 1:batch
            x = x_batch[b]
            f = f_batch[b]
            !isnothing(logger) && log_eval!(logger, x, f)
            if f < f_best
                x_best = copy(x)
                f_best = f
            end
        end
        n_remaining -= batch
    end
    return x_best, f_best
end

function do_cma_es(
    obj::Function,
    x0::Vector{Float64},
    lower::Vector{Float64},
    upper::Vector{Float64};
    n_evals::Int,
    rng_seed::Int,
    sigma_frac::Float64 = 0.20,
    logger::Union{Nothing, TraceLogger} = nothing,
)
    rng = MersenneTwister(rng_seed)
    n = length(x0)
    λ = max(4, 4 + floor(Int, 3 * log(n)))
    μ = max(2, fld(λ, 2))
    weights = log.(μ .+ 0.5) .- log.(1:μ)
    weights ./= sum(weights)
    μeff = inv(sum(abs2, weights))
    cc = (4 + μeff / n) / (n + 4 + 2 * μeff / n)
    cs = (μeff + 2) / (n + μeff + 5)
    c1 = 2 / ((n + 1.3)^2 + μeff)
    cmu = min(1 - c1, 2 * (μeff - 2 + inv(μeff)) / ((n + 2)^2 + μeff))
    damps = 1 + 2 * max(0, sqrt((μeff - 1) / (n + 1)) - 1) + cs
    chi_n = sqrt(n) * (1 - 1 / (4n) + 1 / (21n^2))

    m = clamp.(copy(x0), lower, upper)
    sigma = max(sigma_frac * mean(upper .- lower), 1e-3)
    pc = zeros(Float64, n)
    ps = zeros(Float64, n)
    C = Matrix{Float64}(I, n, n)
    B = Matrix{Float64}(I, n, n)
    D = ones(Float64, n)
    invsqrtC = Matrix{Float64}(I, n, n)

    x_best = copy(m)
    f_best = Float64(obj(m))
    !isnothing(logger) && log_eval!(logger, m, f_best)
    evals = 1
    generation = 0

    while evals < max(n_evals, 1)
        generation += 1
        batch = min(λ, n_evals - evals)
        arz = [randn(rng, n) for _ in 1:batch]
        arx = [clamp.(m .+ sigma .* (B * (D .* arz[k])), lower, upper) for k in 1:batch]
        fvals = zeros(Float64, batch)
        @threads for k in 1:batch
            fvals[k] = Float64(obj(arx[k]))
        end
        for k in 1:batch
            !isnothing(logger) && log_eval!(logger, arx[k], fvals[k])
            if fvals[k] < f_best
                f_best = fvals[k]
                x_best = copy(arx[k])
            end
        end
        evals += batch

        order = sortperm(fvals)
        μsel = min(μ, batch)
        wsel = copy(weights[1:μsel])
        wsel ./= sum(wsel)
        x_old = copy(m)
        m .= 0.0
        z_mean = zeros(Float64, n)
        for i in 1:μsel
            idx = order[i]
            m .+= wsel[i] .* arx[idx]
            z_mean .+= wsel[i] .* arz[idx]
        end
        y_mean = (m .- x_old) ./ sigma
        ps = (1 - cs) .* ps .+ sqrt(cs * (2 - cs) * μeff) .* (B * z_mean)
        hsig = norm(ps) / sqrt(1 - (1 - cs)^(2 * generation)) / chi_n < (1.4 + 2 / (n + 1))
        pc = (1 - cc) .* pc .+ (hsig ? 1.0 : 0.0) .* sqrt(cc * (2 - cc) * μeff) .* y_mean

        C .*= (1 - c1 - cmu + (hsig ? 0.0 : c1 * cc * (2 - cc)))
        C .+= c1 .* (pc * pc')
        for i in 1:μsel
            idx = order[i]
            y = (arx[idx] .- x_old) ./ sigma
            C .+= cmu .* wsel[i] .* (y * y')
        end
        C .= 0.5 .* (C .+ C')

        eig = eigen(Symmetric(C))
        vals = clamp.(eig.values, 1e-12, Inf)
        B = Matrix(eig.vectors)
        D = sqrt.(vals)
        invsqrtC = B * Diagonal(1.0 ./ D) * B'
        sigma *= exp((cs / damps) * (norm(ps) / chi_n - 1))
        sigma = clamp(sigma, 1e-6, 1e3)
    end

    return x_best, f_best
end

function threaded_finite_difference_gradient!(
    G::AbstractVector{Float64},
    f::Function,
    x::AbstractVector{Float64},
    lower::AbstractVector{Float64},
    upper::AbstractVector{Float64},
)
    fx = Float64(f(x))
    @threads for j in eachindex(x)
        xj = x[j]
        scale = max(abs(xj), 1.0)
        h = max(sqrt(eps(Float64)) * scale, 1e-6)
        xpj = min(xj + h, upper[j])
        if xpj == xj
            xpj = max(xj - h, lower[j])
        end
        step = xpj - xj
        if step == 0.0
            G[j] = 0.0
            continue
        end
        xp = copy(x)
        xp[j] = xpj
        fp = Float64(f(xp))
        G[j] = (fp - fx) / step
    end
    return G
end

function threaded_forwarddiff_hessian!(H::AbstractMatrix{Float64}, f::Function, x::AbstractVector{Float64})
    Htmp = ForwardDiff.hessian(f, x)
    copyto!(H, Htmp)
    return H
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
        tumor_auc = Float64(m.tumor_auc),
        auc_tdbc = Float64(m.auc_tdbc),
    ))
end

function run_methods(case::MethodBenchmarkCase, out_dir::AbstractString)
    lower = copy(case.lower_mg)
    upper = copy(case.upper_mg)
    x0 = clamp.(copy(case.x0), lower, upper)
    loss_primal = make_objective(case; sensealg = nothing)
    loss_dto = make_objective(case; sensealg = DTO_PASSTHROUGH)

    nm_iters = parse(Int, get(ENV, "SINGLE_BENCH_NM_ITERS", "60"))
    lbfgs_iters = parse(Int, get(ENV, "SINGLE_BENCH_LBFGS_ITERS", "25"))
    rs_evals = parse(Int, get(ENV, "SINGLE_BENCH_RS_EVALS", "180"))
    mc_steps = parse(Int, get(ENV, "SINGLE_BENCH_MC_STEPS", "260"))
    cma_evals = parse(Int, get(ENV, "SINGLE_BENCH_CMA_EVALS", string(rs_evals)))
    filter_env = strip(get(ENV, "SINGLE_BENCH_METHOD_FILTER", ""))
    method_filter = isempty(filter_env) ? nothing : Set(lowercase.(split(filter_env, ",")))

    trace_path = artifact_path(out_dir, "optimization_traces")
    summary_path = artifact_path(out_dir, "optimization_summary")
    progress_path = artifact_path(out_dir, "optimization_progress"; ext = "json")
    resume = get(ENV, "SINGLE_BENCH_RESUME", "0") == "1"
    shard_mode = shard_mode_enabled()
    if !resume
        if shard_mode
            method_slugs = if isnothing(method_filter)
                nothing
            else
                Set(method_slug(strip(m)) for m in method_filter)
            end
            for f in readdir(out_dir; join = true)
                base = basename(f)
                if startswith(base, "optimization_summary__") || startswith(base, "optimization_traces__") || startswith(base, "optimization_progress__")
                    if isnothing(method_slugs)
                        rm(f; force = true)
                    else
                        parts = splitext(base)[1]
                        if occursin("__", parts)
                            slug = split(parts, "__", limit = 2)[2]
                            slug in method_slugs && rm(f; force = true)
                        end
                    end
                end
            end
        else
            rm(trace_path; force = true)
            rm(summary_path; force = true)
            rm(progress_path; force = true)
        end
    end

    summary_rows = NamedTuple[]
    trace_rows = NamedTuple[]
    completed_methods = Set{String}()
    if resume
        if shard_mode
            for f in readdir(out_dir; join = true)
                base = basename(f)
                if startswith(base, "optimization_summary__") && endswith(base, ".csv") && filesize(f) > 0
                    existing_summary = DataFrame(CSV.File(f))
                    append!(summary_rows, NamedTuple.(eachrow(existing_summary)))
                    union!(completed_methods, String.(existing_summary.method))
                elseif startswith(base, "optimization_traces__") && endswith(base, ".csv") && filesize(f) > 0
                    existing_traces = DataFrame(CSV.File(f))
                    append!(trace_rows, NamedTuple.(eachrow(existing_traces)))
                end
            end
        else
            if isfile(summary_path) && filesize(summary_path) > 0
                existing_summary = DataFrame(CSV.File(summary_path))
                append!(summary_rows, NamedTuple.(eachrow(existing_summary)))
                union!(completed_methods, String.(existing_summary.method))
            end
            if isfile(trace_path) && filesize(trace_path) > 0
                existing_traces = DataFrame(CSV.File(trace_path))
                append!(trace_rows, NamedTuple.(eachrow(existing_traces)))
            end
        end
    end

    function finish_method(method::String, logger::TraceLogger, x_best::Vector{Float64}, status::String)
        append!(trace_rows, logger.rows)
        runtime_s = time() - logger.t0
        pre_len = length(summary_rows)
        push_summary_row!(summary_rows, case, method, x_best, runtime_s, status)
        row = summary_rows[pre_len + 1]
        append_row!(artifact_path(out_dir, "optimization_summary"; method = method), row)
        write_progress!(artifact_path(out_dir, "optimization_progress"; method = method, ext = "json"), progress_payload(logger; status = status))
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
        method in completed_methods && return false
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
        x_best = copy(x0)
        status = "ok"
        try
            x_best, _ = do_random_search(
                loss_primal,
                x0,
                lower,
                upper;
                n_evals = rs_evals,
                local_prob = 0.75,
                local_sigma_frac = 0.15,
                rng_seed = case.random_seed + 17,
                logger = logger,
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
            x_best = isempty(logger.best_x) ? clamp.(Optim.minimizer(res), lower, upper) : copy(logger.best_x)
            status = string(Optim.converged(res) ? "ok" : "not_converged")
        catch err
            x_best = isempty(logger.best_x) ? copy(x0) : copy(logger.best_x)
            status = err isa MaxEvalStop ? "capped_eval_$(logger.max_evals)" : "error: " * sprint(showerror, err)
        end
        finish_method("nelder_mead", logger, x_best, status)
    end

    # 4) CMA-ES
    if wants("cma_es")
        println("running method=cma_es objective=$(case.objective_name)")
        flush(stdout)
        logger = make_trace_logger(case.objective_name, "cma_es", out_dir)
        x_best = copy(x0)
        status = "ok"
        try
            x_best, _ = do_cma_es(
                loss_primal,
                x0,
                lower,
                upper;
                n_evals = cma_evals,
                rng_seed = case.random_seed + 23,
                sigma_frac = 0.20,
                logger = logger,
            )
        catch err
            x_best = isempty(logger.best_x) ? copy(x0) : copy(logger.best_x)
            status = "error: " * sprint(showerror, err)
        end
        finish_method("cma_es", logger, x_best, status)
    end

    function run_gradient_method(method::String, optimizer, grad_fun; hess_fun = nothing)
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
            default_f_calls = method == "finite_diff_lbfgs" ? "400" : "0"
            default_f_calls = method == "finite_diff_bfgs" ? "400" : default_f_calls
            default_f_calls = method == "dto_forward_ad_bfgs" ? "400" : default_f_calls
            env_key = if method == "finite_diff_bfgs"
                "SINGLE_BENCH_FD_BFGS_F_CALLS"
            elseif method == "dto_forward_ad_bfgs"
                "SINGLE_BENCH_BFGS_F_CALLS"
            elseif method == "finite_diff_lbfgs"
                "SINGLE_BENCH_FD_LBFGS_F_CALLS"
            else
                "SINGLE_BENCH_LBFGS_F_CALLS"
            end
            lbfgs_f_calls = parse(Int, get(ENV, env_key, default_f_calls))
            opts = Optim.Options(
                iterations = lbfgs_iters,
                f_calls_limit = lbfgs_f_calls > 0 ? lbfgs_f_calls : typemax(Int),
                show_trace = false,
                store_trace = false,
            )
            res = if isnothing(hess_fun)
                optimize(f_logged, g!, lower, upper, copy(x0), Fminbox(optimizer), opts)
            else
                h! = (H, x) -> begin
                    hess_fun(H, x)
                    return H
                end
                optimize(f_logged, g!, h!, lower, upper, copy(x0), optimizer, opts)
            end
            x_best = isempty(logger.best_x) ? clamp.(Optim.minimizer(res), lower, upper) : copy(logger.best_x)
            status = string(Optim.converged(res) ? "ok" : "not_converged")
        catch err
            x_best = isempty(logger.best_x) ? copy(x0) : copy(logger.best_x)
            status = err isa MaxEvalStop ? "capped_eval_$(logger.max_evals)" : "error: " * sprint(showerror, err)
        end
        finish_method(method, logger, x_best, status)
    end

    run_gradient_method("finite_diff_bfgs", BFGS(), (G, x) -> threaded_finite_difference_gradient!(G, loss_primal, x, lower, upper))
    run_gradient_method("dto_forward_ad_lbfgs", LBFGS(), (G, x) -> ForwardDiff.gradient!(G, loss_dto, x))
    run_gradient_method("dto_forward_ad_bfgs", BFGS(), (G, x) -> ForwardDiff.gradient!(G, loss_dto, x))
    run_gradient_method(
        "dto_forward_ad_ipnewton",
        IPNewton(),
        (G, x) -> ForwardDiff.gradient!(G, loss_dto, x);
        hess_fun = (H, x) -> threaded_forwarddiff_hessian!(H, loss_dto, x),
    )
    run_gradient_method("dto_reverse_ad_lbfgs", LBFGS(), (G, x) -> copyto!(G, reverse_adjoint_dose_gradient(case, x, InterpolatingAdjoint(autojacvec = EnzymeVJP()))[2]))
    run_gradient_method("otd_reverse_ad_lbfgs", LBFGS(), (G, x) -> copyto!(G, reverse_adjoint_dose_gradient(case, x, QuadratureAdjoint(autojacvec = EnzymeVJP()))[2]))

    return DataFrame(summary_rows), DataFrame(trace_rows)
end

function write_meta(case::MethodBenchmarkCase, out_dir::AbstractString)
    meta = Dict(
        "objective_name" => case.objective_name,
        "objective_label" => case.objective_label,
        "sample_id" => case.sample_id,
        "sample_selection_mode" => case.sample_selection_mode,
        "cohort_source" => case.cohort_source,
        "cohort_dir" => case.cohort_dir,
        "decision_labels" => case.decision_labels,
        "x0_mg" => case.x0,
        "x0_mode" => case.x0_mode,
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
        "remission_refs" => case.remission_refs,
        "paper_stepup_refs" => case.paper_stepup_refs,
        "tracking_target_name" => get(case.target_trajectory, :name, nothing),
        "tracking_target_label" => get(case.target_trajectory, :label, nothing),
    )
    open(joinpath(out_dir, "optimization_meta.json"), "w") do io
        JSON3.pretty(io, meta)
    end
    return nothing
end

function make_case(objective_name::String)
    scenario = load_benchmark_scenario()
    samples, cohort_source, cohort_dir = load_benchmark_samples(scenario.horizon_days)
    seed = parse(Int, get(ENV, "SINGLE_BENCH_RANDOM_SEED", "20260311"))
    sample_selection_mode = isempty(strip(get(ENV, "SINGLE_BENCH_SAMPLE_ID", ""))) ? "random" : "explicit"
    sample = pick_sample(samples, seed)

    summary = JSON3.read(read(joinpath(REPO_ROOT, "generated", "figures", "optimization", "clinical_lhs_100_cohort_dose_lbfgs_20260310", "optimization_summary.json"), String))
    x0_default = Float64.(summary["initial_random_decision_doses_mg"])
    decision_labels, decision_groups, clinical_decision = build_benchmark_decision_scheme(scenario.dose_times_days, scenario.dose_mg, scenario.scenario_name)
    lower = isnothing(scenario.lower_mg) ? Float64.(summary["lower_mg"]) : Float64.(scenario.lower_mg)
    upper = isnothing(scenario.upper_mg) ? Float64.(summary["upper_mg"]) : Float64.(scenario.upper_mg)
    x0, x0_mode = method_specific_x0(x0_default, lower, upper, clinical_decision, seed, scenario.scenario_name)

    alg, solver_label = parse_solver(get(ENV, "SINGLE_BENCH_SOLVER", "rodas4p"))
    saveat = build_saveat(scenario.horizon_days, parse(Float64, get(ENV, "SINGLE_BENCH_SAVEAT_DT", "1.0")))
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
        dose_times_days = scenario.dose_times_days,
        decision_labels = decision_labels,
        decision_groups = decision_groups,
        x_clinical = Float64.(clinical_decision),
        x0 = x0,
        lower_mg = lower,
        upper_mg = upper,
        bw_kg = bw_kg,
        horizon_days = scenario.horizon_days,
        saveat = saveat,
        tox_window_days = tox_window_days,
        tox_tau = tox_tau,
        alg = alg,
        abstol = MMC.SOLVER_ABSTOL,
        reltol = MMC.SOLVER_RELTOL,
        maxiters = 1_000_000,
        post_event_proposed_dt = post_event_proposed_dt,
        target_trajectory = nothing,
        simple_scales = nothing,
        clinical_scales = nothing,
        clinical_weights = normalize_weights(LossWeights(0.2, 0.25, 0.15, 0.2, 0.2)),
        remission_refs = nothing,
        objective_name = objective_name,
        objective_label = objective_name,
        scenario_label = scenario.scenario_label,
        solver_label = solver_label,
        random_seed = seed,
        sample_selection_mode = sample_selection_mode,
        x0_mode = x0_mode,
        cohort_source = cohort_source,
        cohort_dir = cohort_dir,
    )

    m0 = trajectory_metrics(case_stub, x0; sensealg = nothing)
    tracking_target_name = lowercase(strip(get(ENV, "SINGLE_BENCH_TRACKING_TARGET", scenario.scenario_name == "paper84" ? "paper_recommended" : "clinical")))
    target_label = "Nominal clinical decision"
    target_decision = copy(clinical_decision)
    if tracking_target_name != "clinical"
        catalog = paper_regimen_catalog()
        haskey(catalog, tracking_target_name) || error("Unsupported SINGLE_BENCH_TRACKING_TARGET=$tracking_target_name")
        target_decision = Float64.(catalog[tracking_target_name].decision_doses_mg)
        target_label = catalog[tracking_target_name].label
    end
    mclin, anchor_solver_label = safe_anchor_metrics(case_stub, target_decision)
    paper_anchor_decision = Float64[1.6, 10.0, 10.0, 20.0]
    rp2d_anchor_decision = copy(PAPER84_ACTUAL_RP2D_DOSES_MG)
    paper_stepup_metrics, paper_stepup_solver_label = safe_anchor_metrics(case_stub, paper_anchor_decision)
    rp2d_stepup_metrics, rp2d_stepup_solver_label = safe_anchor_metrics(case_stub, rp2d_anchor_decision)
    il6_ref = cohort_source == "topparam10" ? load_topparam10_il6_ref() : max(maximum(Float64.(mclin.il6)), 1.0)
    remission_lb_paper, remission_li_paper = remission_il6_components(case_stub, mclin.t, mclin.bt, mclin.il6, mclin.bt0; il6_ref = il6_ref)
    remission_refs = (
        il6_ref = Float64(il6_ref),
        lb_paper = Float64(remission_lb_paper),
        li_paper = Float64(remission_li_paper),
    )
    paper_stepup_refs = (
        paper = make_stepup_anchor_ref(
            case_stub,
            paper_stepup_metrics;
            anchor_regimen = "paper_best",
            anchor_label = "Paper double-step (1.6/10/10->20)",
            anchor_solver_label = paper_stepup_solver_label,
            anchor_doses_mg = expand_decision_doses_generic(paper_anchor_decision, decision_groups, length(scenario.dose_times_days)),
        ),
        rp2d = make_stepup_anchor_ref(
            case_stub,
            rp2d_stepup_metrics;
            anchor_regimen = "rp2d",
            anchor_label = "RP2D 1/2/60/60/30/30",
            anchor_solver_label = rp2d_stepup_solver_label,
            anchor_doses_mg = rp2d_anchor_decision,
        ),
    )
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
        name = tracking_target_name,
        label = target_label,
        t = mclin.t,
        bt = Float64.(mclin.bt),
        bt0 = Float64(mclin.bt0),
        il6_log = log1p.(smooth_pos.(Float64.(mclin.il6)) ./ max(maximum(smooth_pos.(Float64.(mclin.il6))), 1.0)),
        tdbc_log = log1p.(Float64.(mclin.tdbc) ./ max(maximum(Float64.(mclin.tdbc)), 1e-6)),
        il6_scale = max(maximum(smooth_pos.(Float64.(mclin.il6))), 1.0),
        tdbc_scale = max(maximum(Float64.(mclin.tdbc)), 1e-6),
    )
    patient_ideal_target = if objective_name == "patient_ideal_tracking"
        make_patient_ideal_target(case_stub, paper_stepup_metrics, rp2d_stepup_metrics)
    elseif objective_name in ACHIEVABLE_IDEAL_OBJECTIVES
        make_achievable_ideal_target(case_stub, objective_name)
    else
        nothing
    end
    labels = Dict(
        "simple" => "Simple toxicity + tumor proxy",
        "auc_combo" => "Normalized BTumor AUC + IL6combo AUC",
        "remission_il6" => "Remission-target + early IL6 loss",
        "clinical" => "Clinical 5-term proxy loss",
        "tracking" => "Quadratic tracking loss to $(target_label) trajectory",
        "stepup_paper" => "Paper-anchored early IL6 minimization with one-sided efficacy penalties",
        "rp2d_stepup" => "RP2D-anchored early IL6 minimization with one-sided efficacy penalties",
        "safety_constrained_ti" => "Safety-constrained efficacy-gain step-up loss",
        "patient_ideal_tracking" => "Patient-specific ideal tumor/IL6 trajectory tracking loss",
        "patient_achievable_ideal_tracking" => "Patient-specific achievable max-tumor/1-mg-IL6 tracking loss",
        "global_best_achievable_ideal_tracking" => "Global-best achievable max-tumor/1-mg-IL6 tracking loss",
        "patient_achievable_ideal_tracking_calibrated" => "Calibrated patient-specific achievable max-tumor/1-mg-IL6 tracking loss",
        "achievable_endpoint_cycle20" => "Patient-specific endpoint global-IL6/day-84 tumor loss with 20 mg/cycle target",
        "achievable_endpoint_cycle60" => "Patient-specific endpoint global-IL6/day-84 tumor loss with 60 mg/cycle target",
        "patient_achievable_ideal_tracking_calibrated_cycle20" => "Calibrated patient-specific achievable tracking loss with 20 mg/cycle target",
        "patient_achievable_ideal_tracking_calibrated_cycle60" => "Calibrated patient-specific achievable tracking loss with 60 mg/cycle target",
    )

    return MethodBenchmarkCase(
        sample_id = sample.sample_id,
        params = sample.params,
        pvec = pvec,
        p_idx = p_idx,
        u_idx = u_idx,
        dose_times_days = scenario.dose_times_days,
        decision_labels = decision_labels,
        decision_groups = decision_groups,
        x_clinical = Float64.(clinical_decision),
        x0 = x0,
        lower_mg = lower,
        upper_mg = upper,
        bw_kg = bw_kg,
        horizon_days = scenario.horizon_days,
        saveat = saveat,
        tox_window_days = tox_window_days,
        tox_tau = tox_tau,
        alg = alg,
        abstol = MMC.SOLVER_ABSTOL,
        reltol = MMC.SOLVER_RELTOL,
        maxiters = 1_000_000,
        post_event_proposed_dt = post_event_proposed_dt,
        target_trajectory = is_ideal_tracking_objective(objective_name) ? patient_ideal_target : target_trajectory,
        simple_scales = simple_scales,
        clinical_scales = clinical_scales,
        clinical_weights = normalize_weights(LossWeights(0.2, 0.25, 0.15, 0.2, 0.2)),
        remission_refs = remission_refs,
        paper_stepup_refs = paper_stepup_refs,
        objective_name = objective_name,
        objective_label = labels[objective_name],
        scenario_label = scenario.scenario_label,
        solver_label = solver_label,
        random_seed = seed,
        sample_selection_mode = sample_selection_mode,
        x0_mode = x0_mode,
        cohort_source = cohort_source,
        cohort_dir = cohort_dir,
    )
end

function main()
    objective_name = lowercase(get(ENV, "SINGLE_BENCH_OBJECTIVE", "simple"))
    objective_name in ("simple", "auc_combo", "remission_il6", "clinical", "tracking", "stepup_paper", "rp2d_stepup", "safety_constrained_ti", "patient_ideal_tracking", "patient_achievable_ideal_tracking", "global_best_achievable_ideal_tracking", "patient_achievable_ideal_tracking_calibrated", "achievable_endpoint_cycle20", "achievable_endpoint_cycle60", "patient_achievable_ideal_tracking_calibrated_cycle20", "patient_achievable_ideal_tracking_calibrated_cycle60") || error("Unsupported SINGLE_BENCH_OBJECTIVE=$objective_name")
    out_root_env = get(ENV, "SINGLE_BENCH_OUT_ROOT", joinpath("generated", "figures", "optimization", "single_sample_dose_method_benchmark_20260311"))
    out_root = isabspath(out_root_env) ? out_root_env : joinpath(REPO_ROOT, out_root_env)
    out_dir = joinpath(out_root, objective_name)
    mkpath(out_dir)

    case = make_case(objective_name)
    println(@sprintf("benchmarking objective=%s sample_id=%d solver=%s", case.objective_name, case.sample_id, case.solver_label))
    write_meta(case, out_dir)
    results, traces = run_methods(case, out_dir)
    if !shard_mode_enabled()
        CSV.write(joinpath(out_dir, "optimization_summary.csv"), results)
        CSV.write(joinpath(out_dir, "optimization_traces.csv"), traces)
        println("Wrote ", joinpath(out_dir, "optimization_summary.csv"))
        println("Wrote ", joinpath(out_dir, "optimization_traces.csv"))
    else
        println("Wrote method-sharded summary/trace files in ", out_dir)
    end
    println("Wrote ", joinpath(out_dir, "optimization_meta.json"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
