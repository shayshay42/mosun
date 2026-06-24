# Checkpointed narrow Susilo/Fig. 5 eFAST campaign.
#
# This runner reuses the expanded Susilo/Fig. 5 parameter universe from
# run_susilo_fig5_il6_dummy_efast.jl, but narrows the sensitivity target to a
# single fixed-regimen dose tradeoff scalar:
#
#   dose_problem_scalar =
#       delta_log10_first_peak_il6_0_7d_treated_minus_nodose
#     - lambda * clamp(delta_day84_tumor_pct_nodose_minus_treated / 100, 0, cap)
#
# Larger values mean a worse safety/efficacy tradeoff: a larger treatment-
# induced IL6 peak after baseline correction, with less day-84 tumor benefit.

using CSV
using DataFrames
using Dates
using JSON3
using Base.Threads
using SciMLBase

include(joinpath(@__DIR__, "run_susilo_fig5_il6_dummy_efast.jl"))

env_bool(name::String, default::Bool) = begin
    value = lowercase(strip(get(ENV, name, "")))
    isempty(value) && return default
    value in ("1", "true", "yes", "y", "on")
end

function minimal_metric_row(seed::Int, eval_idx::Int, parameter_block::Int, sample_in_block::Int,
    spec::Fig5RegimenSpec, status::String, t::Vector{Float64}, bt::Vector{Float64}, il6::Vector{Float64})

    if isempty(bt)
        return (
            seed = seed,
            efast_eval_idx = eval_idx,
            parameter_block = parameter_block,
            sample_in_block = sample_in_block,
            regimen = spec.name,
            regimen_label = spec.label,
            status = status,
            first_peak_il6_0_7d = NaN,
            day84_tumor_size_change_pct = NaN,
        )
    end
    bt0 = bt[1]
    bt84 = nearest_value(t, bt, 84.0)
    (
        seed = seed,
        efast_eval_idx = eval_idx,
        parameter_block = parameter_block,
        sample_in_block = sample_in_block,
        regimen = spec.name,
        regimen_label = spec.label,
        status = status,
        first_peak_il6_0_7d = max_window(t, il6, 0.0, 7.0),
        day84_tumor_size_change_pct = percent_change(bt84, bt0),
    )
end

function simulate_minimal_regimen_metrics(seed::Int, eval_idx::Int, parameter_block::Int, sample_in_block::Int,
    p, spec::Fig5RegimenSpec, alg, save_times::Vector{Float64}, horizon_days::Float64,
    post_event_dt::Float64, bw_kg::Float64)

    regimen = regimen_from_spec(spec, bw_kg)
    try
        _, sol = MMC.solve_regimen(
            regimen,
            p,
            alg;
            tspan = (0.0, horizon_days),
            saveat = save_times,
            callback_mode = :callback,
            post_event_proposed_dt = post_event_dt,
            abstol = MMC.SOLVER_ABSTOL,
            reltol = MMC.SOLVER_RELTOL,
            save_everystep = false,
            maxiters = 1_000_000,
        )
        status = SciMLBase.successful_retcode(sol.retcode) ? "success" : String(Symbol(sol.retcode))
        cache = MMC.zero_observables_cache()
        t = collect(Float64, sol.t)
        bt = Vector{Float64}(undef, length(t))
        il6 = similar(bt)
        for i in eachindex(t)
            MMC.update_observables!(cache, sol.u[i], p, t[i])
            bt[i] = Float64(cache.Btumor_perml)
            il6[i] = Float64(cache.IL6combo)
        end
        minimal_metric_row(seed, eval_idx, parameter_block, sample_in_block, spec, status, t, bt, il6)
    catch err
        empty = Float64[]
        minimal_metric_row(seed, eval_idx, parameter_block, sample_in_block, spec, "error:$(typeof(err)):$(sprint(showerror, err))", empty, empty, empty)
    end
end

function atomic_csv_write(path::AbstractString, df::DataFrame)
    mkpath(dirname(path))
    tmp = string(path, ".tmp.", getpid(), ".", threadid())
    CSV.write(tmp, df)
    mv(tmp, path; force = true)
    path
end

function shard_dir(out_dir::AbstractString, seed::Int)
    joinpath(out_dir, "shards", "seed_$(seed)")
end

function shard_paths(out_dir::AbstractString, seed::Int, parameter_block::Int)
    dir = shard_dir(out_dir, seed)
    stem = "block_$(lpad(string(parameter_block), 3, '0'))"
    (
        metrics = joinpath(dir, "$(stem)_metrics.csv"),
        paired = joinpath(dir, "$(stem)_paired.csv"),
        done = joinpath(dir, "$(stem).done.json"),
    )
end

function shard_is_done(paths)
    isfile(paths.metrics) && isfile(paths.paired) && isfile(paths.done)
end

function mark_shard_done(path::AbstractString, payload)
    mkpath(dirname(path))
    tmp = string(path, ".tmp.", getpid(), ".", threadid())
    open(tmp, "w") do io
        JSON3.pretty(io, payload)
        println(io)
    end
    mv(tmp, path; force = true)
end

function finite_or_nan(x)
    try
        Float64(x)
    catch
        NaN
    end
end

function build_minimal_paired_metrics(metrics::DataFrame; lambda::Float64 = 1.0, tumor_benefit_cap::Float64 = 5.0)
    no_dose = metrics[metrics.regimen .== "no_dose", :]
    treated = metrics[metrics.regimen .!= "no_dose", :]
    key = [:seed, :efast_eval_idx, :parameter_block, :sample_in_block]
    joined = innerjoin(treated, no_dose; on = key, makeunique = true, renamecols = "_treated" => "_nodose")
    rows = NamedTuple[]
    for row in eachrow(joined)
        treated_peak = finite_or_nan(row.first_peak_il6_0_7d_treated)
        nodose_peak = finite_or_nan(row.first_peak_il6_0_7d_nodose)
        treated_day84 = finite_or_nan(row.day84_tumor_size_change_pct_treated)
        nodose_day84 = finite_or_nan(row.day84_tumor_size_change_pct_nodose)
        il6_delta_log10 = log10(max(treated_peak, EPS)) - log10(max(nodose_peak, EPS))
        tumor_benefit_pct = nodose_day84 - treated_day84
        tumor_benefit_scaled = clamp(tumor_benefit_pct / 100.0, 0.0, tumor_benefit_cap)
        dose_problem_scalar = il6_delta_log10 - lambda * tumor_benefit_scaled
        push!(rows, (
            seed = row.seed,
            efast_eval_idx = row.efast_eval_idx,
            parameter_block = row.parameter_block,
            sample_in_block = row.sample_in_block,
            treated_regimen = row.regimen_treated,
            treated_regimen_label = row.regimen_label_treated,
            untreated_regimen = row.regimen_nodose,
            untreated_regimen_label = row.regimen_label_nodose,
            treated_status = row.status_treated,
            untreated_status = row.status_nodose,
            treated_first_peak_il6_0_7d = treated_peak,
            nodose_first_peak_il6_0_7d = nodose_peak,
            delta_log10_first_peak_il6_0_7d_treated_minus_nodose = il6_delta_log10,
            treated_day84_tumor_size_change_pct = treated_day84,
            nodose_day84_tumor_size_change_pct = nodose_day84,
            delta_day84_tumor_pct_nodose_minus_treated = tumor_benefit_pct,
            scaled_day84_tumor_benefit = tumor_benefit_scaled,
            dose_problem_lambda = lambda,
            dose_problem_tumor_benefit_cap = tumor_benefit_cap,
            dose_problem_scalar = dose_problem_scalar,
        ))
    end
    DataFrame(rows)
end

function analyze_checkpointed_tradeoff(paired::DataFrame, universe::DataFrame, seeds::Vector{Int}, n::Int, d::Int, m::Int, omega::Int; analyze_components::Bool = false)
    rows = NamedTuple[]
    param_names = String.(universe.parameter)
    families = String.(universe.family)
    endpoint_specs = [
        (metric = "dose_problem_scalar", family = "dose_tradeoff", transform = "identity"),
    ]
    if analyze_components
        append!(endpoint_specs, [
            (metric = "delta_log10_first_peak_il6_0_7d_treated_minus_nodose", family = "IL6", transform = "identity"),
            (metric = "delta_day84_tumor_pct_nodose_minus_treated", family = "tumor_response", transform = "identity"),
            (metric = "treated_first_peak_il6_0_7d", family = "IL6", transform = "log10_positive"),
        ])
    end

    for seed in seeds
        seed_paired = paired[paired.seed .== seed, :]
        for regimen in unique(String.(seed_paired.treated_regimen))
            subset = seed_paired[seed_paired.treated_regimen .== regimen, :]
            sort!(subset, [:parameter_block, :sample_in_block])
            nrow(subset) == n * d || error("Expected $(n*d) rows for seed=$(seed), regimen=$(regimen), got $(nrow(subset)).")
            for spec in endpoint_specs
                metric = spec.metric
                vals = Vector{Float64}(undef, nrow(subset))
                col = subset[!, Symbol(metric)]
                for (i, x) in enumerate(col)
                    xf = finite_or_nan(x)
                    vals[i] = spec.transform == "log10_positive" && isfinite(xf) ? log10(max(xf, EPS)) : xf
                end
                s1, st = salib_fast_analyze(vals, n, d, m, omega)
                ranks = rank_desc(st)
                endpoint = "$(regimen)__$(metric)"
                for j in 1:d
                    push!(rows, (
                        seed = seed,
                        endpoint = endpoint,
                        endpoint_family = spec.family,
                        regimen = regimen,
                        metric = metric,
                        endpoint_source = "paired_tradeoff",
                        transform = spec.transform,
                        parameter = param_names[j],
                        parameter_family = families[j],
                        S1 = s1[j],
                        ST = st[j],
                        rank_ST = ranks[j],
                    ))
                end
            end
        end
    end
    DataFrame(rows)
end

function write_checkpoint_manifest(out_dir::AbstractString, seeds::Vector{Int}, d::Int)
    rows = NamedTuple[]
    for seed in seeds
        for block in 1:d
            paths = shard_paths(out_dir, seed, block)
            push!(rows, (
                seed = seed,
                parameter_block = block,
                metrics_path = paths.metrics,
                paired_path = paths.paired,
                done_path = paths.done,
                metrics_exists = isfile(paths.metrics),
                paired_exists = isfile(paths.paired),
                done_exists = isfile(paths.done),
            ))
        end
    end
    manifest = DataFrame(rows)
    CSV.write(joinpath(out_dir, "checkpoint_manifest.csv"), manifest)
    manifest
end

function collect_shards(out_dir::AbstractString, seeds::Vector{Int}, d::Int)
    metric_dfs = DataFrame[]
    paired_dfs = DataFrame[]
    for seed in seeds
        for block in 1:d
            paths = shard_paths(out_dir, seed, block)
            shard_is_done(paths) || error("Missing completed shard for seed=$(seed), parameter_block=$(block).")
            push!(metric_dfs, CSV.read(paths.metrics, DataFrame))
            push!(paired_dfs, CSV.read(paths.paired, DataFrame))
        end
    end
    (
        metrics = vcat(metric_dfs...; cols = :union),
        paired = vcat(paired_dfs...; cols = :union),
    )
end

function run_checkpointed_shards!(out_dir::AbstractString, universe::DataFrame, regimens::Vector{Fig5RegimenSpec},
    sample_values_by_seed, parameter_block_by_seed, sample_in_block_by_seed,
    seeds::Vector{Int}, base_params, alg, save_times::Vector{Float64}, horizon_days::Float64,
    post_event_dt::Float64, bw_kg::Float64, lambda::Float64, tumor_benefit_cap::Float64)

    shard_jobs = NamedTuple[]
    d = nrow(universe)
    for seed in seeds
        for block in 1:d
            push!(shard_jobs, (seed = seed, parameter_block = block))
        end
    end
    progress_lock = ReentrantLock()
    completed = Ref(0)
    skipped = Ref(0)
    total = length(shard_jobs)

    @threads for job_idx in eachindex(shard_jobs)
        job = shard_jobs[job_idx]
        paths = shard_paths(out_dir, job.seed, job.parameter_block)
        if shard_is_done(paths)
            lock(progress_lock)
            try
                skipped[] += 1
                println("Skipping completed shard $(skipped[] + completed[]) / $(total): seed=$(job.seed) block=$(job.parameter_block)")
                flush(stdout)
            finally
                unlock(progress_lock)
            end
            continue
        end

        sample_values = sample_values_by_seed[job.seed]
        parameter_block = parameter_block_by_seed[job.seed]
        sample_in_block = sample_in_block_by_seed[job.seed]
        eval_indices = findall(==(job.parameter_block), parameter_block)
        sort!(eval_indices; by = idx -> sample_in_block[idx])

        rows = NamedTuple[]
        for eval_idx in eval_indices
            p = apply_sample_to_params(base_params, universe, collect(sample_values[eval_idx, :]))
            for spec in regimens
                push!(rows, simulate_minimal_regimen_metrics(
                    job.seed,
                    eval_idx,
                    parameter_block[eval_idx],
                    sample_in_block[eval_idx],
                    p,
                    spec,
                    alg,
                    save_times,
                    horizon_days,
                    post_event_dt,
                    bw_kg,
                ))
            end
        end
        metrics = DataFrame(rows)
        paired = build_minimal_paired_metrics(metrics; lambda = lambda, tumor_benefit_cap = tumor_benefit_cap)
        atomic_csv_write(paths.metrics, metrics)
        atomic_csv_write(paths.paired, paired)
        mark_shard_done(paths.done, Dict(
            "completed_at" => string(Dates.now()),
            "seed" => job.seed,
            "parameter_block" => job.parameter_block,
            "n_metric_rows" => nrow(metrics),
            "n_paired_rows" => nrow(paired),
            "lambda" => lambda,
            "tumor_benefit_cap" => tumor_benefit_cap,
        ))

        lock(progress_lock)
        try
            completed[] += 1
            println("Completed shard $(skipped[] + completed[]) / $(total): seed=$(job.seed) block=$(job.parameter_block) rows=$(nrow(metrics)) paired=$(nrow(paired)) at $(Dates.now())")
            flush(stdout)
        finally
            unlock(progress_lock)
        end
    end
end

function main_checkpointed()
    out_dir = abspath(env_string("SUSILO_FIG5_CKPT_EFAST_OUT_DIR", joinpath(MosunModelCoreSupport.REPO_ROOT, "generated", "figures", "sensitivity", "susilo_fig5_il6_day84_checkpointed_efast_20260504")))
    mkpath(out_dir)
    bounds_path = abspath(env_string("SUSILO_FIG5_CKPT_EFAST_BOUNDS_CSV", joinpath(MosunModelCoreSupport.REPO_ROOT, "generated", "figures", "reference", "dlbcl_stack_parameter_ranges_all", "dlbcl_stack_parameter_ranges_all.csv")))
    n = env_int("SUSILO_FIG5_CKPT_EFAST_N", 1024)
    m = env_int("SUSILO_FIG5_CKPT_EFAST_M", 4)
    seeds = env_seed_list("SUSILO_FIG5_CKPT_EFAST_SEEDS", [20260504, 20260505])
    saveat_dt = env_float("SUSILO_FIG5_CKPT_EFAST_SAVEAT_DT", 0.125)
    horizon_days = env_float("SUSILO_FIG5_CKPT_EFAST_HORIZON_DAYS", 84.0)
    post_event_dt = env_float("SUSILO_FIG5_CKPT_EFAST_POST_EVENT_DT", 0.01)
    bw_kg = env_float("SUSILO_FIG5_CKPT_EFAST_BW_KG", 70.0)
    solver = env_string("SUSILO_FIG5_CKPT_EFAST_SOLVER", "rodas4p")
    param_limit = env_int("SUSILO_FIG5_CKPT_EFAST_PARAM_LIMIT", 0)
    regimen_limit = env_int("SUSILO_FIG5_CKPT_EFAST_REGIMEN_LIMIT", 0)
    lambda = env_float("SUSILO_FIG5_CKPT_EFAST_LAMBDA", 1.0)
    tumor_benefit_cap = env_float("SUSILO_FIG5_CKPT_EFAST_TUMOR_BENEFIT_CAP", 5.0)
    analyze_components = env_bool("SUSILO_FIG5_CKPT_EFAST_ANALYZE_COMPONENTS", false)

    universe = parameter_universe(bounds_path; param_limit = param_limit)
    d = nrow(universe)
    d > 0 || error("Parameter universe is empty.")
    CSV.write(joinpath(out_dir, "parameter_universe.csv"), universe)

    regimens = fig5_regimens()
    if regimen_limit > 0 && regimen_limit < length(regimens)
        regimens = regimens[1:regimen_limit]
    end
    any(r -> r.name == "no_dose", regimens) || error("no_dose regimen must be included for paired deltas.")

    all_sample_rows = DataFrame[]
    sample_values_by_seed = Dict{Int, Matrix{Float64}}()
    parameter_block_by_seed = Dict{Int, Vector{Int}}()
    sample_in_block_by_seed = Dict{Int, Vector{Int}}()
    omega_by_seed = Dict{Int, Int}()
    for seed in seeds
        u, parameter_block, sample_in_block, omega = salib_fast_sample(n, d, m, seed)
        sample_values = transform_unit_sample(u, universe)
        sample_values_by_seed[seed] = sample_values
        parameter_block_by_seed[seed] = parameter_block
        sample_in_block_by_seed[seed] = sample_in_block
        omega_by_seed[seed] = omega
        df = DataFrame(sample_values, Symbol.(universe.parameter))
        df.seed .= seed
        df.efast_eval_idx = collect(1:size(sample_values, 1))
        df.parameter_block = parameter_block
        df.sample_in_block = sample_in_block
        select!(df, [:seed, :efast_eval_idx, :parameter_block, :sample_in_block, Symbol.(universe.parameter)...])
        push!(all_sample_rows, df)
    end
    sample_matrix = vcat(all_sample_rows...; cols = :union)
    CSV.write(joinpath(out_dir, "efast_sample_matrix.csv"), sample_matrix)

    base_params = build_shipped_dlbcl_base_params()
    MMC.has_parameter("end_time") && MMC.set_param!(base_params, "end_time", horizon_days)
    alg = MosunModelCoreSupport.make_solver_alg(solver)
    save_times = collect(0.0:saveat_dt:horizon_days)

    println("Starting checkpointed IL6/day84 eFAST: D=$(d), N=$(n), seeds=$(seeds), blocks=$(d * length(seeds)), regimens=$(length(regimens)), threads=$(nthreads()), out=$(out_dir)")
    println("Primary scalar: delta_log10_first_peak_il6_0_7d_treated_minus_nodose - $(lambda) * clamp(delta_day84_tumor_pct_nodose_minus_treated / 100, 0, $(tumor_benefit_cap))")
    flush(stdout)

    run_checkpointed_shards!(
        out_dir,
        universe,
        regimens,
        sample_values_by_seed,
        parameter_block_by_seed,
        sample_in_block_by_seed,
        seeds,
        base_params,
        alg,
        save_times,
        horizon_days,
        post_event_dt,
        bw_kg,
        lambda,
        tumor_benefit_cap,
    )

    manifest = write_checkpoint_manifest(out_dir, seeds, d)
    all(manifest.done_exists) || error("Not all checkpoint shards completed.")
    collected = collect_shards(out_dir, seeds, d)
    metrics = collected.metrics
    paired = collected.paired
    CSV.write(joinpath(out_dir, "simulation_metrics_minimal_long.csv"), metrics)
    CSV.write(joinpath(out_dir, "paired_delta_metrics_minimal.csv"), paired)

    omega_values = unique(collect(Base.values(omega_by_seed)))
    length(omega_values) == 1 || error("Unexpected different omega values by seed: $(omega_values)")
    indices = analyze_checkpointed_tradeoff(paired, universe, seeds, n, d, m, omega_values[1]; analyze_components = analyze_components)
    CSV.write(joinpath(out_dir, "efast_indices_by_endpoint.csv"), indices)
    thresholds = build_thresholds(indices)
    CSV.write(joinpath(out_dir, "dummy_null_thresholds.csv"), thresholds)
    screen = build_significant_screen(indices, thresholds, seeds)
    CSV.write(joinpath(out_dir, "significant_parameter_screen.csv"), screen)
    stability = build_seed_stability(indices, thresholds)
    CSV.write(joinpath(out_dir, "seed_stability_summary.csv"), stability)

    meta = Dict(
        "created_at" => string(Dates.now()),
        "created_by" => "julia/run_susilo_fig5_il6_day84_checkpointed_efast.jl",
        "output_root" => out_dir,
        "bounds_csv" => bounds_path,
        "n" => n,
        "m" => m,
        "seeds" => seeds,
        "n_parameters" => d,
        "n_regimens" => length(regimens),
        "regimens" => [Dict("name" => r.name, "label" => r.label, "dose_mg" => r.dose_mg) for r in regimens],
        "dummy_parameters" => DUMMY_PARAMS,
        "forced_susilo_response_variables" => ["tumor_burden_factor", "BT_ratio_tumor_init", "kBtumorprolif", "kIL6prod", "thalfIL6", "IL6_tiss_contribution"],
        "primary_scalar" => Dict(
            "name" => "dose_problem_scalar",
            "formula" => "delta_log10_first_peak_il6_0_7d_treated_minus_nodose - lambda * clamp(delta_day84_tumor_pct_nodose_minus_treated / 100, 0, tumor_benefit_cap)",
            "lambda" => lambda,
            "tumor_benefit_cap" => tumor_benefit_cap,
            "interpretation" => "Larger values mean worse treatment-induced IL6 safety cost after untreated baseline correction and less day-84 tumor benefit.",
        ),
        "component_analysis_enabled" => analyze_components,
        "checkpointing" => Dict(
            "shard_unit" => "seed x parameter_block",
            "skip_completed_shards" => true,
            "manifest" => joinpath(out_dir, "checkpoint_manifest.csv"),
        ),
        "simulation_settings" => Dict(
            "cohort_variant_ids" => DLBCL_VARIANT_IDS,
            "solver" => solver,
            "horizon_days" => horizon_days,
            "saveat_dt" => saveat_dt,
            "post_event_proposed_dt" => post_event_dt,
            "bw_kg" => bw_kg,
            "julia_threads" => nthreads(),
        ),
        "significance_rule" => Dict(
            "median_ST_across_seeds_above" => "max dummy-null median ST for endpoint",
            "seed_stability_required" => length(seeds) >= 3 ? "2/3 seeds above dummy threshold" : "1 seed above dummy threshold",
        ),
        "row_counts" => Dict(
            "simulation_metrics_minimal_long" => nrow(metrics),
            "paired_delta_metrics_minimal" => nrow(paired),
            "efast_indices_by_endpoint" => nrow(indices),
            "significant_parameter_screen" => nrow(screen),
        ),
        "outputs" => Dict(
            "parameter_universe" => joinpath(out_dir, "parameter_universe.csv"),
            "efast_sample_matrix" => joinpath(out_dir, "efast_sample_matrix.csv"),
            "simulation_metrics_minimal_long" => joinpath(out_dir, "simulation_metrics_minimal_long.csv"),
            "paired_delta_metrics_minimal" => joinpath(out_dir, "paired_delta_metrics_minimal.csv"),
            "efast_indices_by_endpoint" => joinpath(out_dir, "efast_indices_by_endpoint.csv"),
            "dummy_null_thresholds" => joinpath(out_dir, "dummy_null_thresholds.csv"),
            "significant_parameter_screen" => joinpath(out_dir, "significant_parameter_screen.csv"),
            "seed_stability_summary" => joinpath(out_dir, "seed_stability_summary.csv"),
        ),
    )
    write_meta(joinpath(out_dir, "susilo_fig5_il6_day84_checkpointed_efast_meta.json"), meta)
    println("Completed checkpointed IL6/day84 eFAST outputs in $(out_dir)")
    flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_checkpointed()
end
