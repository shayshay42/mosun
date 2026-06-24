using Base.Threads
using CSV
using DataFrames
using Dates
using JSON3
using SciMLBase

include(joinpath(@__DIR__, "run_susilo_fig5_il6_dummy_efast.jl"))

const VPOP_FIG5_REGIMENS = [r for r in fig5_regimens() if r.name != "no_dose"]

env_bool2(name::String, default::Bool) = begin
    value = lowercase(strip(get(ENV, name, "")))
    isempty(value) && return default
    value in ("1", "true", "yes", "y", "on")
end

function vpop_out_dir()
    raw = get(
        ENV,
        "SUSILO_EFAST_HOSSEINI_VPOP_OUT_DIR",
        joinpath(MosunModelCoreSupport.REPO_ROOT, "generated", "figures", "vpop_pruning", "susilo_efast_hosseini2020_fig5_vpop250_20260505"),
    )
    isabspath(raw) ? raw : joinpath(MosunModelCoreSupport.REPO_ROOT, raw)
end

function resim_candidate_csv(out_dir::AbstractString)
    raw = get(ENV, "SUSILO_EFAST_HOSSEINI_VPOP_CANDIDATES_CSV", joinpath(out_dir, "trajectory_resim_candidate_parameters.csv"))
    isabspath(raw) ? raw : joinpath(MosunModelCoreSupport.REPO_ROOT, raw)
end

function resim_shard_paths(out_dir::AbstractString, shard_id::Int)
    dir = joinpath(out_dir, "trajectory_resim_shards")
    stem = "shard_$(lpad(string(shard_id), 4, '0'))"
    (
        traj = joinpath(dir, "$(stem)_trajectories.csv"),
        waterfall = joinpath(dir, "$(stem)_waterfall.csv"),
        done = joinpath(dir, "$(stem).done.json"),
    )
end

function resim_shard_done(paths)
    isfile(paths.traj) && isfile(paths.waterfall) && isfile(paths.done)
end

function atomic_write_csv(path::AbstractString, df::DataFrame)
    mkpath(dirname(path))
    tmp = string(path, ".tmp.", getpid(), ".", threadid())
    CSV.write(tmp, df)
    mv(tmp, path; force = true)
end

function atomic_write_json(path::AbstractString, payload)
    mkpath(dirname(path))
    tmp = string(path, ".tmp.", getpid(), ".", threadid())
    open(tmp, "w") do io
        JSON3.pretty(io, payload)
        println(io)
    end
    mv(tmp, path; force = true)
end

function concat_csv_files(paths::Vector{String}, out_path::AbstractString)
    isempty(paths) && error("No shard files to concatenate for $(out_path)")
    mkpath(dirname(out_path))
    tmp = string(out_path, ".tmp")
    open(tmp, "w") do out
        first_file = true
        for path in paths
            open(path, "r") do input
                header = readline(input)
                if first_file
                    println(out, header)
                    first_file = false
                end
                for line in eachline(input)
                    println(out, line)
                end
            end
        end
    end
    mv(tmp, out_path; force = true)
end

function response_category_fig5(pct_change::Float64)
    if pct_change <= -99.9
        "cr"
    elseif pct_change <= -30.0
        "pr"
    elseif pct_change < 20.0
        "sd"
    else
        "pd"
    end
end

function candidate_params(base_params, universe::DataFrame, row)
    values = Vector{Float64}(undef, nrow(universe))
    for (i, pname) in enumerate(String.(universe.parameter))
        values[i] = Float64(row[Symbol(pname)])
    end
    p = apply_sample_to_params(base_params, universe, values)
    for pname in ("KBptumor", "KTrptumor")
        sym = Symbol(pname)
        if hasproperty(row, sym) && MMC.has_parameter(pname)
            value = Float64(getproperty(row, sym))
            isfinite(value) && value > 0.0 && MMC.set_param!(p, pname, value)
        end
    end
    p
end

function simulate_candidate_regimen(row, p, spec::Fig5RegimenSpec, alg, saveat::Vector{Float64}, horizon_days::Float64, profile_horizon_days::Float64, post_event_dt::Float64, bw_kg::Float64)
    regimen = regimen_from_spec(spec, bw_kg)
    try
        built, sol = MMC.solve_regimen(
            regimen,
            p,
            alg;
            tspan = (0.0, horizon_days),
            saveat = saveat,
            callback_mode = :callback,
            post_event_proposed_dt = post_event_dt,
            abstol = MMC.SOLVER_ABSTOL,
            reltol = MMC.SOLVER_RELTOL,
            save_everystep = false,
            maxiters = 1_000_000,
        )
        status = SciMLBase.successful_retcode(sol.retcode) ? "success" : String(Symbol(sol.retcode))
        cache = MMC.zero_observables_cache()
        traj_rows = NamedTuple[]
        for j in eachindex(sol.t)
            t = Float64(sol.t[j])
            t < profile_horizon_days || continue
            MMC.update_observables!(cache, sol.u[j], p, t)
            push!(traj_rows, (
                candidate_id = String(row.candidate_id),
                resim_candidate_rank = Int(row.resim_candidate_rank),
                seed = Int(row.seed),
                efast_eval_idx = Int(row.efast_eval_idx),
                parameter_block = Int(row.parameter_block),
                sample_in_block = Int(row.sample_in_block),
                regimen = spec.name,
                regimen_label = spec.label,
                status = status,
                time_day = t,
                il6combo = Float64(cache.IL6combo),
                tafraction_pb_pct = 100.0 * Float64(cache.Tafraction_pb),
                btumor_perml = Float64(cache.Btumor_perml),
            ))
        end
        MMC.update_observables!(cache, built.initial_u, p, 0.0)
        bt0 = Float64(cache.Btumor_perml)
        MMC.update_observables!(cache, sol.u[end], p, Float64(sol.t[end]))
        bt84 = Float64(cache.Btumor_perml)
        pct_change = 100.0 * (bt84 / max(bt0, 1e-12) - 1.0)
        waterfall = (
            candidate_id = String(row.candidate_id),
            resim_candidate_rank = Int(row.resim_candidate_rank),
            seed = Int(row.seed),
            efast_eval_idx = Int(row.efast_eval_idx),
            parameter_block = Int(row.parameter_block),
            sample_in_block = Int(row.sample_in_block),
            regimen = spec.name,
            regimen_label = spec.label,
            status = status,
            btumor_initial_perml = bt0,
            btumor_day84_perml = bt84,
            tumor_size_change_pct = pct_change,
            response_category = response_category_fig5(pct_change),
            responder_gt50 = pct_change <= -50.0,
        )
        traj_rows, waterfall
    catch err
        status = "error:$(typeof(err)):$(sprint(showerror, err))"
        traj_rows = NamedTuple[]
        waterfall = (
            candidate_id = String(row.candidate_id),
            resim_candidate_rank = Int(row.resim_candidate_rank),
            seed = Int(row.seed),
            efast_eval_idx = Int(row.efast_eval_idx),
            parameter_block = Int(row.parameter_block),
            sample_in_block = Int(row.sample_in_block),
            regimen = spec.name,
            regimen_label = spec.label,
            status = status,
            btumor_initial_perml = NaN,
            btumor_day84_perml = NaN,
            tumor_size_change_pct = NaN,
            response_category = "error",
            responder_gt50 = false,
        )
        traj_rows, waterfall
    end
end

function run_shards!(candidates::DataFrame, universe::DataFrame, out_dir::AbstractString, shard_size::Int, alg, saveat::Vector{Float64}, horizon_days::Float64, profile_horizon_days::Float64, post_event_dt::Float64, bw_kg::Float64)
    n = nrow(candidates)
    shard_count = cld(n, shard_size)
    base_params = build_shipped_dlbcl_base_params()
    MMC.has_parameter("end_time") && MMC.set_param!(base_params, "end_time", horizon_days)
    progress_lock = ReentrantLock()
    completed = Ref(0)
    skipped = Ref(0)

    @threads for shard_id in 1:shard_count
        paths = resim_shard_paths(out_dir, shard_id)
        if resim_shard_done(paths)
            lock(progress_lock)
            try
                skipped[] += 1
                println("Skipping completed trajectory shard $(skipped[] + completed[]) / $(shard_count): shard=$(shard_id)")
                flush(stdout)
            finally
                unlock(progress_lock)
            end
            continue
        end
        lo = (shard_id - 1) * shard_size + 1
        hi = min(shard_id * shard_size, n)
        traj_rows = NamedTuple[]
        waterfall_rows = NamedTuple[]
        for row in eachrow(candidates[lo:hi, :])
            p = candidate_params(base_params, universe, row)
            for spec in VPOP_FIG5_REGIMENS
                trs, wf = simulate_candidate_regimen(row, p, spec, alg, saveat, horizon_days, profile_horizon_days, post_event_dt, bw_kg)
                append!(traj_rows, trs)
                push!(waterfall_rows, wf)
            end
        end
        atomic_write_csv(paths.traj, DataFrame(traj_rows))
        atomic_write_csv(paths.waterfall, DataFrame(waterfall_rows))
        atomic_write_json(paths.done, Dict(
            "completed_at" => string(Dates.now()),
            "shard_id" => shard_id,
            "row_start" => lo,
            "row_end" => hi,
            "n_candidates" => hi - lo + 1,
            "n_trajectory_rows" => length(traj_rows),
            "n_waterfall_rows" => length(waterfall_rows),
        ))
        lock(progress_lock)
        try
            completed[] += 1
            println("Completed trajectory shard $(skipped[] + completed[]) / $(shard_count): shard=$(shard_id) candidates=$(hi - lo + 1) traj_rows=$(length(traj_rows)) waterfall_rows=$(length(waterfall_rows)) at $(Dates.now())")
            flush(stdout)
        finally
            unlock(progress_lock)
        end
    end
end

function write_resim_manifest(out_dir::AbstractString, shard_count::Int)
    rows = NamedTuple[]
    for shard_id in 1:shard_count
        paths = resim_shard_paths(out_dir, shard_id)
        push!(rows, (
            shard_id = shard_id,
            trajectories_path = paths.traj,
            waterfall_path = paths.waterfall,
            done_path = paths.done,
            trajectories_exists = isfile(paths.traj),
            waterfall_exists = isfile(paths.waterfall),
            done_exists = isfile(paths.done),
        ))
    end
    manifest = DataFrame(rows)
    CSV.write(joinpath(out_dir, "trajectory_resim_manifest.csv"), manifest)
    manifest
end

function main()
    out_dir = vpop_out_dir()
    mkpath(out_dir)
    candidates_path = resim_candidate_csv(out_dir)
    isfile(candidates_path) || error("Missing candidates CSV: $(candidates_path)")
    candidates = CSV.read(candidates_path, DataFrame)
    sort!(candidates, :resim_candidate_rank)

    bounds_path = abspath(get(ENV, "SUSILO_EFAST_HOSSEINI_VPOP_BOUNDS_CSV", joinpath(MosunModelCoreSupport.REPO_ROOT, "generated", "figures", "reference", "dlbcl_stack_parameter_ranges_all", "dlbcl_stack_parameter_ranges_all.csv")))
    universe = parameter_universe(bounds_path)
    for pname in String.(universe.parameter)
        hasproperty(candidates, Symbol(pname)) || error("Candidate parameter CSV missing column $(pname)")
    end

    n_limit = parse(Int, get(ENV, "SUSILO_EFAST_HOSSEINI_VPOP_RESIM_LIMIT", "0"))
    if n_limit > 0 && n_limit < nrow(candidates)
        candidates = candidates[1:n_limit, :]
    end
    shard_size = parse(Int, get(ENV, "SUSILO_EFAST_HOSSEINI_VPOP_SHARD_SIZE", "50"))
    horizon_days = parse(Float64, get(ENV, "SUSILO_EFAST_HOSSEINI_VPOP_HORIZON_DAYS", "84.0"))
    profile_horizon_days = parse(Float64, get(ENV, "SUSILO_EFAST_HOSSEINI_VPOP_PROFILE_HORIZON_DAYS", "42.0"))
    saveat_dt = parse(Float64, get(ENV, "SUSILO_EFAST_HOSSEINI_VPOP_SAVEAT_DT", "0.125"))
    post_event_dt = parse(Float64, get(ENV, "SUSILO_EFAST_HOSSEINI_VPOP_POST_EVENT_DT", "0.01"))
    bw_kg = parse(Float64, get(ENV, "SUSILO_EFAST_HOSSEINI_VPOP_BW_KG", "70.0"))
    solver = get(ENV, "SUSILO_EFAST_HOSSEINI_VPOP_SOLVER", "rodas4p")
    aggregate = env_bool2("SUSILO_EFAST_HOSSEINI_VPOP_AGGREGATE", true)
    alg = MosunModelCoreSupport.make_solver_alg(solver)
    saveat = collect(0.0:saveat_dt:horizon_days)
    if isempty(saveat) || saveat[end] < horizon_days
        push!(saveat, horizon_days)
    end
    shard_count = cld(nrow(candidates), shard_size)
    println("Starting eFAST-derived Hosseini Fig.5 candidate resimulation: candidates=$(nrow(candidates)), shards=$(shard_count), shard_size=$(shard_size), threads=$(nthreads()), out=$(out_dir)")
    flush(stdout)

    run_shards!(candidates, universe, out_dir, shard_size, alg, saveat, horizon_days, profile_horizon_days, post_event_dt, bw_kg)
    manifest = write_resim_manifest(out_dir, shard_count)
    all(manifest.done_exists) || error("Not all trajectory resimulation shards completed.")

    if aggregate
        traj_paths = String.(manifest.trajectories_path)
        wf_paths = String.(manifest.waterfall_path)
        concat_csv_files(traj_paths, joinpath(out_dir, "candidate_trajectory_summary.csv"))
        concat_csv_files(wf_paths, joinpath(out_dir, "candidate_waterfall.csv"))
    end
    meta = Dict(
        "created_by" => "julia/export_susilo_efast_hosseini2020_fig5_candidates.jl",
        "created_at" => string(Dates.now()),
        "output_root" => out_dir,
        "candidate_parameters_csv" => candidates_path,
        "n_candidates" => nrow(candidates),
        "n_shards" => shard_count,
        "shard_size" => shard_size,
        "aggregate_outputs" => aggregate,
        "simulation_settings" => Dict(
            "solver" => solver,
            "horizon_days" => horizon_days,
            "profile_horizon_days" => profile_horizon_days,
            "saveat_dt" => saveat_dt,
            "post_event_proposed_dt" => post_event_dt,
            "bw_kg" => bw_kg,
            "julia_threads" => nthreads(),
        ),
        "outputs" => Dict(
            "trajectory_resim_manifest" => joinpath(out_dir, "trajectory_resim_manifest.csv"),
            "candidate_trajectory_summary" => joinpath(out_dir, "candidate_trajectory_summary.csv"),
            "candidate_waterfall" => joinpath(out_dir, "candidate_waterfall.csv"),
        ),
    )
    atomic_write_json(joinpath(out_dir, "candidate_resim_meta.json"), meta)
    println("Completed eFAST-derived Hosseini Fig.5 candidate resimulation in $(out_dir)")
    flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
