using CSV
using DataFrames
using JSON3
using Printf
using Random
using Statistics
using Base.Threads

include(joinpath(@__DIR__, "optimize_clinical_lhs_cohort_dosing.jl"))

const REPO_ROOT = MMC.REPO_ROOT

function loss_scales_from_summary(summary)
    s = summary["scales"]
    return LossScales(
        Float64(s["tox_peak_mean"]),
        Float64(s["tox_peak_max"]),
        Float64(s["tox_auc"]),
        Float64(s["tumor_terminal"]),
        Float64(s["tumor_auc"]),
    )
end

function loss_weights_from_summary(summary)
    w = summary["weights"]
    return normalize_weights(
        LossWeights(
            Float64(w["tox_peak_mean"]),
            Float64(w["tox_peak_max"]),
            Float64(w["tox_auc"]),
            Float64(w["tumor_terminal"]),
            Float64(w["tumor_auc"]),
        ),
    )
end

function default_shared_out_dir()
    return joinpath(REPO_ROOT, "generated", "figures", "optimization", "clinical_lhs_100_cohort_dose_lbfgs_20260310")
end

function evaluate_single_sample(
    prob::CohortOptimizationProblem,
    sample::CohortSample,
    decision_doses::AbstractVector{<:Real},
    scales::LossScales,
    weights::LossWeights;
    saveat = prob.saveat,
)
    full_doses = expand_decision_doses(decision_doses, prob.decision_groups, length(prob.dose_times_days))
    regimen = dose_vector_to_regimen(prob.dose_times_days, full_doses, prob.bw_kg)
    built = MMC.build_problem(
        regimen,
        sample.params;
        tspan = (0.0, prob.horizon_days),
        saveat = saveat,
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
    scaled = scaled_loss(raw, scales, weights)
    auc_tdbc = trapz(t, max.(tdbc, 0.0))
    return (
        raw = raw,
        scaled = scaled,
        auc_tdbc = auc_tdbc,
    )
end

function single_objective(
    prob::CohortOptimizationProblem,
    sample::CohortSample,
    scales::LossScales,
    weights::LossWeights,
)
    return function (x)
        xv = Float64.(x)
        if any(!isfinite, xv)
            return 1.0e12
        end
        viol = max.(prob.lower_mg .- xv, 0.0) .+ max.(xv .- prob.upper_mg, 0.0)
        if any(viol .> 0.0)
            return 1.0e9 + 1.0e6 * sum(abs2, viol)
        end
        try
            result = evaluate_single_sample(prob, sample, xv, scales, weights; saveat = prob.objective_saveat)
            return Float64(result.scaled.loss)
        catch
            return 1.0e12
        end
    end
end

function optimize_single_sample(
    prob::CohortOptimizationProblem,
    sample::CohortSample,
    x_start::Vector{Float64},
    scales::LossScales,
    weights::LossWeights;
    seed::Int,
    n_random::Int,
    local_sigma_frac::Float64,
    cd_sweeps::Int,
    cd_step0::Float64,
    cd_tol::Float64,
)
    rng = MersenneTwister(seed)
    obj = single_objective(prob, sample, scales, weights)
    x_best = clamp.(copy(x_start), prob.lower_mg, prob.upper_mg)
    f_best = Float64(obj(x_best))

    for _ in 1:n_random
        x = similar(x_best)
        if rand(rng) < 0.8
            for j in eachindex(x)
                span = prob.upper_mg[j] - prob.lower_mg[j]
                x[j] = clamp(x_best[j] + randn(rng) * local_sigma_frac * span, prob.lower_mg[j], prob.upper_mg[j])
            end
        else
            for j in eachindex(x)
                x[j] = prob.lower_mg[j] + rand(rng) * (prob.upper_mg[j] - prob.lower_mg[j])
            end
        end
        f = Float64(obj(x))
        if f < f_best
            x_best = copy(x)
            f_best = f
        end
    end

    step = cd_step0
    for _ in 1:cd_sweeps
        improved = false
        for j in eachindex(x_best)
            x_plus = copy(x_best)
            x_plus[j] = min(prob.upper_mg[j], x_plus[j] + step)
            f_plus = Float64(obj(x_plus))

            x_minus = copy(x_best)
            x_minus[j] = max(prob.lower_mg[j], x_minus[j] - step)
            f_minus = Float64(obj(x_minus))

            if f_plus < f_best && f_plus <= f_minus
                x_best = x_plus
                f_best = f_plus
                improved = true
            elseif f_minus < f_best
                x_best = x_minus
                f_best = f_minus
                improved = true
            end
        end
        step *= 0.5
        if step < cd_tol && !improved
            break
        end
    end

    return x_best, f_best
end

function individual_compare_main()
    shared_out_dir_env = get(ENV, "INDIV_COMPARE_SHARED_OUT_DIR", default_shared_out_dir())
    shared_out_dir = isabspath(shared_out_dir_env) ? shared_out_dir_env : joinpath(REPO_ROOT, shared_out_dir_env)
    shared_summary_path = joinpath(shared_out_dir, "optimization_summary.json")
    isfile(shared_summary_path) || error("Shared optimization summary not found: $shared_summary_path")

    summary = JSON3.read(read(shared_summary_path, String))
    prob, regimen_label = build_problem()
    out_dir_default = joinpath(REPO_ROOT, "generated", "figures", "optimization", "clinical_lhs_100_individual_dosing_compare_20260310")
    out_dir_env = get(ENV, "INDIV_COMPARE_OUT_DIR", out_dir_default)
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(REPO_ROOT, out_dir_env)
    mkpath(out_dir)

    scales = loss_scales_from_summary(summary)
    weights = loss_weights_from_summary(summary)
    x_initial = Float64.(summary["initial_random_decision_doses_mg"])
    x_shared = Float64.(summary["optimized_decision_doses_mg"])

    n_random = parse(Int, get(ENV, "INDIV_COMPARE_RANDOM_EVALS", "24"))
    local_sigma_frac = parse(Float64, get(ENV, "INDIV_COMPARE_LOCAL_SIGMA_FRAC", "0.20"))
    cd_sweeps = parse(Int, get(ENV, "INDIV_COMPARE_CD_SWEEPS", "3"))
    cd_step0 = parse(Float64, get(ENV, "INDIV_COMPARE_CD_STEP0_MG", "10.0"))
    cd_tol = parse(Float64, get(ENV, "INDIV_COMPARE_CD_TOL_MG", "0.5"))

    n = length(prob.cohort)
    rows = Vector{NamedTuple}(undef, n)
    completed = Atomic{Int}(0)

    println(@sprintf("individual compare: cohort_n=%d threads=%d solver=%s", n, nthreads(), string(prob.alg)))
    warm_sample = prob.cohort[1]
    println("warming up objective and local search on sample $(warm_sample.sample_id)")
    _ = evaluate_single_sample(prob, warm_sample, x_initial, scales, weights)
    _ = evaluate_single_sample(prob, warm_sample, x_shared, scales, weights)
    _ = optimize_single_sample(
        prob,
        warm_sample,
        x_shared,
        scales,
        weights;
        seed = prob.seed + warm_sample.sample_id,
        n_random = min(n_random, 4),
        local_sigma_frac = local_sigma_frac,
        cd_sweeps = 1,
        cd_step0 = cd_step0,
        cd_tol = cd_tol,
    )

    @threads for i in eachindex(prob.cohort)
        sample = prob.cohort[i]
        row = try
            initial_eval = evaluate_single_sample(prob, sample, x_initial, scales, weights)
            shared_eval = evaluate_single_sample(prob, sample, x_shared, scales, weights)
            t0 = time()
            x_indiv, obj_indiv = optimize_single_sample(
                prob,
                sample,
                x_shared,
                scales,
                weights;
                seed = prob.seed + sample.sample_id,
                n_random = n_random,
                local_sigma_frac = local_sigma_frac,
                cd_sweeps = cd_sweeps,
                cd_step0 = cd_step0,
                cd_tol = cd_tol,
            )
            runtime_s = time() - t0
            indiv_eval = evaluate_single_sample(prob, sample, x_indiv, scales, weights)
            (
                sample_id = sample.sample_id,
                status = "ok",
                runtime_s = runtime_s,
                initial_loss = Float64(initial_eval.scaled.loss),
                shared_loss = Float64(shared_eval.scaled.loss),
                individual_loss = Float64(indiv_eval.scaled.loss),
                individual_gain_vs_shared = Float64(shared_eval.scaled.loss - indiv_eval.scaled.loss),
                initial_best_spd_pct = Float64(initial_eval.raw.best_spd_pct),
                shared_best_spd_pct = Float64(shared_eval.raw.best_spd_pct),
                individual_best_spd_pct = Float64(indiv_eval.raw.best_spd_pct),
                initial_peak_il6 = Float64(initial_eval.raw.peak_il6),
                shared_peak_il6 = Float64(shared_eval.raw.peak_il6),
                individual_peak_il6 = Float64(indiv_eval.raw.peak_il6),
                initial_tox_peak_mean = Float64(initial_eval.raw.tox_peak_mean),
                shared_tox_peak_mean = Float64(shared_eval.raw.tox_peak_mean),
                individual_tox_peak_mean = Float64(indiv_eval.raw.tox_peak_mean),
                initial_tox_peak_max = Float64(initial_eval.raw.tox_peak_max),
                shared_tox_peak_max = Float64(shared_eval.raw.tox_peak_max),
                individual_tox_peak_max = Float64(indiv_eval.raw.tox_peak_max),
                initial_tox_auc = Float64(initial_eval.raw.tox_auc),
                shared_tox_auc = Float64(shared_eval.raw.tox_auc),
                individual_tox_auc = Float64(indiv_eval.raw.tox_auc),
                initial_il6_auc_abs = Float64(initial_eval.raw.il6_auc_abs),
                shared_il6_auc_abs = Float64(shared_eval.raw.il6_auc_abs),
                individual_il6_auc_abs = Float64(indiv_eval.raw.il6_auc_abs),
                initial_auc_tdbc = Float64(initial_eval.auc_tdbc),
                shared_auc_tdbc = Float64(shared_eval.auc_tdbc),
                individual_auc_tdbc = Float64(indiv_eval.auc_tdbc),
                initial_tumor_terminal = Float64(initial_eval.raw.tumor_terminal),
                shared_tumor_terminal = Float64(shared_eval.raw.tumor_terminal),
                individual_tumor_terminal = Float64(indiv_eval.raw.tumor_terminal),
                initial_tumor_auc = Float64(initial_eval.raw.tumor_auc),
                shared_tumor_auc = Float64(shared_eval.raw.tumor_auc),
                individual_tumor_auc = Float64(indiv_eval.raw.tumor_auc),
                initial_tumor_auc_abs = Float64(initial_eval.raw.tumor_auc_abs),
                shared_tumor_auc_abs = Float64(shared_eval.raw.tumor_auc_abs),
                individual_tumor_auc_abs = Float64(indiv_eval.raw.tumor_auc_abs),
                dose1_mg = Float64(x_indiv[1]),
                dose2_mg = Float64(x_indiv[2]),
                dose3_mg = Float64(x_indiv[3]),
                dose4_mg = Float64(x_indiv[4]),
            )
        catch err
            (
                sample_id = sample.sample_id,
                status = "error: " * sprint(showerror, err),
                runtime_s = NaN,
                initial_loss = NaN,
                shared_loss = NaN,
                individual_loss = NaN,
                individual_gain_vs_shared = NaN,
                initial_best_spd_pct = NaN,
                shared_best_spd_pct = NaN,
                individual_best_spd_pct = NaN,
                initial_peak_il6 = NaN,
                shared_peak_il6 = NaN,
                individual_peak_il6 = NaN,
                initial_tox_peak_mean = NaN,
                shared_tox_peak_mean = NaN,
                individual_tox_peak_mean = NaN,
                initial_tox_peak_max = NaN,
                shared_tox_peak_max = NaN,
                individual_tox_peak_max = NaN,
                initial_tox_auc = NaN,
                shared_tox_auc = NaN,
                individual_tox_auc = NaN,
                initial_il6_auc_abs = NaN,
                shared_il6_auc_abs = NaN,
                individual_il6_auc_abs = NaN,
                initial_auc_tdbc = NaN,
                shared_auc_tdbc = NaN,
                individual_auc_tdbc = NaN,
                initial_tumor_terminal = NaN,
                shared_tumor_terminal = NaN,
                individual_tumor_terminal = NaN,
                initial_tumor_auc = NaN,
                shared_tumor_auc = NaN,
                individual_tumor_auc = NaN,
                initial_tumor_auc_abs = NaN,
                shared_tumor_auc_abs = NaN,
                individual_tumor_auc_abs = NaN,
                dose1_mg = NaN,
                dose2_mg = NaN,
                dose3_mg = NaN,
                dose4_mg = NaN,
            )
        end
        rows[i] = row
        done = atomic_add!(completed, 1) + 1
        if done == 1 || done % 10 == 0 || done == n
            println(@sprintf("completed %d/%d", done, n))
        end
    end

    df = DataFrame(rows)
    CSV.write(joinpath(out_dir, "individual_vs_shared_comparison.csv"), df)

    ok = df[df.status .== "ok", :]
    n_ok = nrow(ok)
    n_ok > 0 || error("No successful individual optimizations")

    dose_summary = DataFrame(
        decision = ["C1D1", "C1D8", "C1D15_C2D1", "C3plus_q3w"],
        shared_optimized_mg = x_shared,
        individualized_mean_mg = [mean(ok[!, Symbol("dose$(i)_mg")]) for i in 1:4],
        individualized_median_mg = [median(ok[!, Symbol("dose$(i)_mg")]) for i in 1:4],
        individualized_min_mg = [minimum(ok[!, Symbol("dose$(i)_mg")]) for i in 1:4],
        individualized_max_mg = [maximum(ok[!, Symbol("dose$(i)_mg")]) for i in 1:4],
    )
    CSV.write(joinpath(out_dir, "individual_dose_summary.csv"), dose_summary)

    summary_out = Dict(
        "regimen_label" => regimen_label,
        "cohort_size" => n,
        "successful_optimizations" => n_ok,
        "threads" => nthreads(),
        "solver" => string(prob.alg),
        "objective_dt_days" => prob.objective_dt,
        "sample_dt_days" => prob.sample_dt,
        "decision_labels" => prob.decision_labels,
        "initial_random_decision_doses_mg" => x_initial,
        "shared_optimized_decision_doses_mg" => x_shared,
        "individualized_mean_decision_doses_mg" => [mean(ok[!, Symbol("dose$(i)_mg")]) for i in 1:4],
        "individualized_median_decision_doses_mg" => [median(ok[!, Symbol("dose$(i)_mg")]) for i in 1:4],
        "mean_initial_loss" => mean(ok.initial_loss),
        "mean_shared_loss" => mean(ok.shared_loss),
        "mean_individual_loss" => mean(ok.individual_loss),
        "median_initial_loss" => median(ok.initial_loss),
        "median_shared_loss" => median(ok.shared_loss),
        "median_individual_loss" => median(ok.individual_loss),
        "mean_gain_vs_shared" => mean(ok.individual_gain_vs_shared),
        "median_gain_vs_shared" => median(ok.individual_gain_vs_shared),
        "mean_runtime_s" => mean(ok.runtime_s),
        "median_runtime_s" => median(ok.runtime_s),
    )
    open(joinpath(out_dir, "comparison_summary.json"), "w") do io
        JSON3.pretty(io, summary_out)
    end

    println(joinpath(out_dir, "comparison_summary.json"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    individual_compare_main()
end
