using CSV
using DataFrames
using JSON3
using Printf
using Random
using Base.Threads

include(joinpath(@__DIR__, "optimize_clinical_lhs_cohort_dosing.jl"))

const REPO_ROOT = TCellEngagerQSP.REPO_ROOT

function default_shared_out_dir()
    return joinpath(REPO_ROOT, "generated", "figures", "optimization", "clinical_lhs_100_cohort_dose_lbfgs_20260310")
end

function default_indiv_out_dir()
    return joinpath(REPO_ROOT, "generated", "figures", "optimization", "clinical_lhs_100_individual_dosing_compare_20260310")
end

function decision_vector_from_summary(summary, key::String)
    return Float64.(summary[key])
end

function expanded_vector(prob::CohortOptimizationProblem, decision_doses::Vector{Float64})
    return expand_decision_doses(decision_doses, prob.decision_groups, length(prob.dose_times_days))
end

function strategy_rows(df::DataFrame, strategy::String)
    out = DataFrame(
        strategy = fill(strategy, nrow(df)),
        sample_id = Int.(df.sample_id),
        best_spd_pct = Float64.(df.best_spd_pct),
        peak_il6 = Float64.(df.peak_il6),
        day_of_global_peak_il6 = Float64.(df.day_of_global_peak_il6),
        auc_tdbc = Float64.(df.auc_tdbc),
        tox_peak_mean = Float64.(df.tox_peak_mean),
        tox_peak_max = Float64.(df.tox_peak_max),
        tox_auc = Float64.(df.tox_auc),
        il6_auc_abs = Float64.(df.il6_auc_abs),
        tumor_terminal = Float64.(df.tumor_terminal),
        tumor_auc = Float64.(df.tumor_auc),
        tumor_auc_abs = Float64.(df.tumor_auc_abs),
    )
    return out
end

function evaluate_decision_vector(prob::CohortOptimizationProblem, sample::CohortSample, decision_doses::Vector{Float64})
    full_doses = expanded_vector(prob, decision_doses)
    regimen = dose_vector_to_regimen(prob.dose_times_days, full_doses, prob.bw_kg)
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
    sol.retcode == SciMLBase.ReturnCode.Success || error("retcode=$(sol.retcode)")
    cache = MMC.zero_observables_cache()
    t, bt, il6, tdbc = solution_series(sol, sample.params, cache)
    raw = raw_metrics_from_series(prob, t, bt, il6)
    return (
        sample_id = sample.sample_id,
        best_spd_pct = Float64(raw.best_spd_pct),
        peak_il6 = Float64(raw.peak_il6),
        day_of_global_peak_il6 = Float64(raw.day_of_global_peak_il6),
        auc_tdbc = trapz(t, max.(tdbc, 0.0)),
        tox_peak_mean = Float64(raw.tox_peak_mean),
        tox_peak_max = Float64(raw.tox_peak_max),
        tox_auc = Float64(raw.tox_auc),
        il6_auc_abs = Float64(raw.il6_auc_abs),
        tumor_terminal = Float64(raw.tumor_terminal),
        tumor_auc = Float64(raw.tumor_auc),
        tumor_auc_abs = Float64(raw.tumor_auc_abs),
    )
end

function evaluate_individual_random(prob::CohortOptimizationProblem)
    n = length(prob.cohort)
    rows = Vector{NamedTuple}(undef, n)
    @threads for i in eachindex(prob.cohort)
        sample = prob.cohort[i]
        rng = MersenneTwister(prob.seed + 1000 + sample.sample_id)
        decision = [prob.lower_mg[j] + rand(rng) * (prob.upper_mg[j] - prob.lower_mg[j]) for j in eachindex(prob.lower_mg)]
        row = evaluate_decision_vector(prob, sample, decision)
        rows[i] = merge(row, (strategy = "individual_random",))
    end
    return DataFrame(rows)
end

function evaluate_individualized_optimum(prob::CohortOptimizationProblem, indiv_ok::DataFrame)
    dose_map = Dict(Int(r.sample_id) => [Float64(r.dose1_mg), Float64(r.dose2_mg), Float64(r.dose3_mg), Float64(r.dose4_mg)] for r in eachrow(indiv_ok))
    sample_lookup = Dict(sample.sample_id => sample for sample in prob.cohort)
    rows = Vector{NamedTuple}(undef, nrow(indiv_ok))
    @threads for i in 1:nrow(indiv_ok)
        sample_id = Int(indiv_ok.sample_id[i])
        rows[i] = merge(evaluate_decision_vector(prob, sample_lookup[sample_id], dose_map[sample_id]), (strategy = "individualized_optimum",))
    end
    return DataFrame(rows)
end

function main()
    shared_dir_env = get(ENV, "STRAT_EVAL_SHARED_OUT_DIR", default_shared_out_dir())
    indiv_dir_env = get(ENV, "STRAT_EVAL_INDIV_OUT_DIR", default_indiv_out_dir())
    out_dir_env = get(ENV, "STRAT_EVAL_OUT_DIR", joinpath(REPO_ROOT, "generated", "figures", "optimization", "clinical_lhs_100_dose_strategy_overlay_20260310"))
    shared_dir = isabspath(shared_dir_env) ? shared_dir_env : joinpath(REPO_ROOT, shared_dir_env)
    indiv_dir = isabspath(indiv_dir_env) ? indiv_dir_env : joinpath(REPO_ROOT, indiv_dir_env)
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(REPO_ROOT, out_dir_env)
    mkpath(out_dir)

    summary = JSON3.read(read(joinpath(shared_dir, "optimization_summary.json"), String))
    prob, regimen_label = build_problem()

    indiv_cmp = DataFrame(CSV.File(joinpath(indiv_dir, "individual_vs_shared_comparison.csv")))
    indiv_ok = indiv_cmp[indiv_cmp.status .== "ok", :]
    nrow(indiv_ok) == length(prob.cohort) || @warn("Individualized comparison has $(nrow(indiv_ok)) rows, expected $(length(prob.cohort))")

    println(@sprintf("strategy evaluation: cohort_n=%d threads=%d solver=%s", length(prob.cohort), Threads.nthreads(), string(prob.alg)))

    min_decision = copy(prob.lower_mg)
    max_decision = copy(prob.upper_mg)
    initial_decision = decision_vector_from_summary(summary, "initial_random_decision_doses_mg")
    cohort_decision = decision_vector_from_summary(summary, "optimized_decision_doses_mg")
    min_metrics, _ = simulate_cohort_trajectories(prob, min_decision)
    max_metrics, _ = simulate_cohort_trajectories(prob, max_decision)
    initial_metrics, _ = simulate_cohort_trajectories(prob, initial_decision)
    cohort_metrics, _ = simulate_cohort_trajectories(prob, cohort_decision)
    individual_random_df = evaluate_individual_random(prob)
    individualized_df = evaluate_individualized_optimum(prob, indiv_ok)

    overlay_df = vcat(
        strategy_rows(min_metrics, "min_dose"),
        strategy_rows(max_metrics, "max_dose"),
        individual_random_df,
        strategy_rows(initial_metrics, "initial_random"),
        strategy_rows(cohort_metrics, "cohort_optimum"),
        individualized_df,
    )

    CSV.write(joinpath(out_dir, "strategy_endpoint_distributions.csv"), overlay_df)

    schedule_rows = DataFrame(
        strategy = String[],
        decision_label = String[],
        decision_dose_mg = Float64[],
        expanded_schedule = String[],
    )

    function add_schedule_rows!(strategy::String, decision_doses::Vector{Float64})
        expanded = expanded_vector(prob, decision_doses)
        expanded_txt = join([@sprintf("%.3f", x) for x in expanded], ",")
        for (lab, val) in zip(prob.decision_labels, decision_doses)
            push!(schedule_rows, (strategy, lab, Float64(val), expanded_txt))
        end
    end

    add_schedule_rows!("min_dose", min_decision)
    add_schedule_rows!("max_dose", max_decision)
    add_schedule_rows!("initial_random", initial_decision)
    add_schedule_rows!("cohort_optimum", cohort_decision)
    random_mean = zeros(Float64, length(prob.decision_labels))
    for j in eachindex(prob.decision_labels)
        random_mean[j] = mean([begin
            rng = MersenneTwister(prob.seed + 1000 + s.sample_id)
            prob.lower_mg[j] + rand(rng) * (prob.upper_mg[j] - prob.lower_mg[j])
        end for s in prob.cohort])
    end
    add_schedule_rows!("individual_random_mean", Float64.(random_mean))

    dose_summary = DataFrame(CSV.File(joinpath(indiv_dir, "individual_dose_summary.csv")))
    for r in eachrow(dose_summary)
        push!(schedule_rows, (
            "individualized_optimum_mean",
            String(r.decision),
            Float64(r.individualized_mean_mg),
            join([@sprintf("%.3f", x) for x in expanded_vector(prob, Float64.(dose_summary.individualized_mean_mg))], ","),
        ))
        push!(schedule_rows, (
            "individualized_optimum_median",
            String(r.decision),
            Float64(r.individualized_median_mg),
            join([@sprintf("%.3f", x) for x in expanded_vector(prob, Float64.(dose_summary.individualized_median_mg))], ","),
        ))
    end
    CSV.write(joinpath(out_dir, "strategy_dose_schedules.csv"), schedule_rows)

    summary_out = Dict(
        "regimen_label" => regimen_label,
        "decision_labels" => prob.decision_labels,
        "strategies" => [
            "min_dose",
            "max_dose",
            "individual_random",
            "initial_random",
            "cohort_optimum",
            "individualized_optimum",
        ],
        "notes" => Dict(
            "min_dose" => "All 4 decision doses set to lower bound",
            "max_dose" => "All 4 decision doses set to upper bound",
            "individual_random" => "Per-parameterization random decision vector within bounds",
            "initial_random" => "Initial random cohort-wide decision vector",
            "cohort_optimum" => "Shared cohort-wide optimized decision vector",
            "individualized_optimum" => "Per-parameterization individualized solution from local search seeded at cohort optimum",
        ),
        "horizon_days" => prob.horizon_days,
        "solver" => string(prob.alg),
        "objective_dt_days" => prob.objective_dt,
        "sample_dt_days" => prob.sample_dt,
        "cohort_size" => length(prob.cohort),
        "shared_summary_path" => joinpath(shared_dir, "optimization_summary.json"),
        "individual_compare_path" => joinpath(indiv_dir, "individual_vs_shared_comparison.csv"),
    )
    open(joinpath(out_dir, "strategy_overlay_summary.json"), "w") do io
        JSON3.pretty(io, summary_out)
    end

    println(joinpath(out_dir, "strategy_endpoint_distributions.csv"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
