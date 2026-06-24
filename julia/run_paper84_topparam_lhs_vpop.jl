using CSV
using DataFrames
using JSON3
using Random
using Statistics
using Base.Threads
using SciMLBase

include(joinpath(@__DIR__, "src", "MosunModelCoreSupport.jl"))
using .MosunModelCoreSupport

const REPO_ROOT = MosunModelCoreSupport.REPO_ROOT
const MMC = MosunModelCoreSupport.MMC
const BW_KG = parse(Float64, get(ENV, "TOPPARAM_VPOP_BW_KG", "70.0"))
const DLBCL_VARIANT_IDS = [5, 9, 14, 20, 24, 25, 27, 28]
const TOP_PARAM_NAMES = [
    "KBptumor",
    "kBkill",
    "Vtumor",
    "KTrptumor",
    "fKmTB_kill",
    "kBprolif",
    "fBkill",
    "Vpb",
    "Vtissue3",
    "Cl_tdb",
]

function trapz(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    n == length(y) || throw(ArgumentError("x and y lengths differ"))
    n < 2 && return 0.0
    acc = 0.0
    @inbounds for i in 1:(n - 1)
        dx = Float64(x[i + 1] - x[i])
        acc += 0.5 * dx * (Float64(y[i]) + Float64(y[i + 1]))
    end
    return acc
end

function struct_to_param_dict(p::MMC.MosunParams)
    out = Dict{String, Float64}()
    for sym in MMC.PARAMETER_NAMES
        out[String(sym)] = Float64(getproperty(p, sym))
    end
    return out
end

function lhs_unit(n::Int, d::Int, rng::AbstractRNG)
    U = Matrix{Float64}(undef, n, d)
    for j in 1:d
        perm = randperm(rng, n)
        for i in 1:n
            U[i, j] = (perm[i] - rand(rng)) / n
        end
    end
    return U
end

function map_lhs_log10_bounds(U::Matrix{Float64}, bounds::Vector{Tuple{Float64, Float64}})
    n, d = size(U)
    X = Matrix{Float64}(undef, n, d)
    for j in 1:d
        lb, ub = bounds[j]
        lb > 0.0 || error("Log-space LHS requires positive lower bound for $(TOP_PARAM_NAMES[j]); got $lb")
        ub > lb || error("Invalid bounds for $(TOP_PARAM_NAMES[j]): ($lb, $ub)")
        llb = log10(lb)
        lub = log10(ub)
        @inbounds for i in 1:n
            X[i, j] = 10.0^(llb + (lub - llb) * U[i, j])
        end
    end
    return X
end

function selected_bounds()
    csv_path = joinpath(
        REPO_ROOT,
        "generated",
        "figures",
        "sensitivity",
        "paper84_loginputs_log10out_auc_btumor_n4000_seed20260313",
        "efast_auc_btumor_indices.csv",
    )
    df = DataFrame(CSV.File(csv_path))
    rows = Dict(String(r.parameter) => (Float64(r.lower_bound), Float64(r.upper_bound)) for r in eachrow(df))
    return [rows[name] for name in TOP_PARAM_NAMES]
end

function build_base_params(; horizon_days::Float64)
    p = deepcopy(MMC.default_params())
    MosunModelCoreSupport.apply_variant_overrides!(
        p,
        DLBCL_VARIANT_IDS,
        MosunModelCoreSupport.load_variant_overrides(),
    )
    if MMC.has_parameter("PKflag")
        MMC.set_param!(p, "PKflag", 1.0)
    end
    if MMC.has_parameter("VPid")
        MMC.set_param!(p, "VPid", 1.0)
    end
    if MMC.has_parameter("fvalidation")
        MMC.set_param!(p, "fvalidation", 0.0)
    end
    if MMC.has_parameter("end_time")
        MMC.set_param!(p, "end_time", horizon_days)
    end
    return struct_to_param_dict(p)
end

function make_regimen_specs(rng::AbstractRNG)
    random_doses = 60.0 .* rand(rng, 6)
    return [
        (name = "no_dose", label = "No dose", dose_days = [0.0, 7.0, 14.0, 21.0, 42.0, 63.0], dose_mg = zeros(Float64, 6)),
        (name = "random_shared", label = "Random shared dose", dose_days = [0.0, 7.0, 14.0, 21.0, 42.0, 63.0], dose_mg = random_doses),
        (name = "max_dose", label = "Max dose 60/60/60/60/60/60", dose_days = [0.0, 7.0, 14.0, 21.0, 42.0, 63.0], dose_mg = fill(60.0, 6)),
        (name = "paper_best", label = "Paper best 1.6/10/10/20/20/20", dose_days = [0.0, 7.0, 14.0, 21.0, 42.0, 63.0], dose_mg = [1.6, 10.0, 10.0, 20.0, 20.0, 20.0]),
        (name = "clinical_stepup", label = "1/2/60/60/30/30", dose_days = [0.0, 7.0, 14.0, 21.0, 42.0, 63.0], dose_mg = [1.0, 2.0, 60.0, 60.0, 30.0, 30.0]),
    ]
end

function regimen_from_mg(dose_days::Vector{Float64}, dose_mg::Vector{Float64}; bw_kg::Float64 = BW_KG, target::Symbol = :TDBc_ugperkg)
    length(dose_days) == length(dose_mg) || throw(ArgumentError("dose_days and dose_mg lengths differ"))
    events = MMC.MosunRegimenEvent{Float64}[]
    for (t, d) in zip(dose_days, dose_mg)
        amt = d * 1000.0 / bw_kg
        iszero(amt) && continue
        push!(events, MMC.MosunRegimenEvent{Float64}(target, Float64(t), Float64(amt), 0.0))
    end
    return MMC.MosunRegimen{Float64}(events)
end

function simulate_metrics(
    sample_id::Int,
    regimen_name::String,
    regimen_label::String,
    dose_days::Vector{Float64},
    dose_mg::Vector{Float64},
    pmap::Dict{String, Float64},
    alg,
    saveat::Vector{Float64},
    horizon_days::Float64,
)
    params = MMC.params_from_dict(pmap; ignore_unknown = true)
    regimen = regimen_from_mg(dose_days, dose_mg; bw_kg = BW_KG, target = :TDBc_ugperkg)
    _, sol = MMC.solve_regimen(
        regimen,
        params,
        alg;
        tspan = (0.0, horizon_days),
        saveat = saveat,
        callback_mode = :callback,
        abstol = SOLVER_ABSTOL,
        reltol = SOLVER_RELTOL,
    )
    if sol.retcode != SciMLBase.ReturnCode.Success
    return (
        sample_id = sample_id,
        regimen = regimen_name,
        regimen_label = regimen_label,
        status = String(sol.retcode),
        btumor_auc = NaN,
        il6combo_auc = NaN,
        tdbc_auc = NaN,
        peak_il6 = NaN,
        best_spd_pct = NaN,
        final_tumor = NaN,
    )
end

    bt_idx = MMC.dynamic_state_index(:Btumor)
    cache = MMC.zero_observables_cache()
    n = length(sol.t)
    bt = Vector{Float64}(undef, n)
    il6 = Vector{Float64}(undef, n)
    tdbc = Vector{Float64}(undef, n)
    @inbounds for i in eachindex(sol.t)
        t = Float64(sol.t[i])
        u = sol.u[i]
        bt[i] = Float64(u[bt_idx])
        il6[i] = MMC.value_at(u, params, t, :IL6combo, cache)
        tdbc[i] = MMC.value_at(u, params, t, :TDBc_ugperml, cache)
    end
    bt0 = max(bt[1], 1e-12)
    spd = @. 100.0 * (bt / bt0 - 1.0)
    spd_eval = length(spd) > 1 ? spd[2:end] : spd
    return (
        sample_id = sample_id,
        regimen = regimen_name,
        regimen_label = regimen_label,
        status = "ok",
        btumor_auc = trapz(sol.t, max.(bt, 0.0)),
        il6combo_auc = trapz(sol.t, max.(il6, 0.0)),
        tdbc_auc = trapz(sol.t, max.(tdbc, 0.0)),
        peak_il6 = maximum(il6),
        best_spd_pct = minimum(spd_eval),
        final_tumor = bt[end],
    )
end

function main()
    out_dir_env = get(ENV, "TOPPARAM_VPOP_OUT_DIR", joinpath("generated", "figures", "vpop", "paper84_topparam_loglhs_vpop200_20260316"))
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(REPO_ROOT, out_dir_env)
    mkpath(out_dir)

    n_samples = parse(Int, get(ENV, "TOPPARAM_VPOP_N", "200"))
    seed = parse(Int, get(ENV, "TOPPARAM_VPOP_SEED", "20260316"))
    save_dt = parse(Float64, get(ENV, "TOPPARAM_VPOP_SAVE_DT_DAYS", "0.125"))
    horizon_days = parse(Float64, get(ENV, "TOPPARAM_VPOP_HORIZON_DAYS", "84.0"))
    solver_name = lowercase(strip(get(ENV, "TOPPARAM_VPOP_SOLVER", "cvode_bdf")))
    alg = make_solver_alg(solver_name)

    rng = MersenneTwister(seed)
    base_map = build_base_params(; horizon_days = horizon_days)
    bounds = selected_bounds()
    regimen_specs = make_regimen_specs(rng)

    U = lhs_unit(n_samples, length(TOP_PARAM_NAMES), rng)
    X = map_lhs_log10_bounds(U, bounds)
    saveat = collect(0.0:save_dt:horizon_days)
    if saveat[end] != horizon_days
        push!(saveat, horizon_days)
    end

    sample_rows = Vector{NamedTuple}(undef, n_samples)
    @threads for i in 1:n_samples
        row = Dict{Symbol, Any}()
        row[:sample_id] = i
        for (j, name) in enumerate(TOP_PARAM_NAMES)
            row[Symbol(name)] = X[i, j]
        end
        sample_rows[i] = (; row...)
    end
    samples_df = DataFrame(sample_rows)

    n_jobs = n_samples * length(regimen_specs)
    results = Vector{NamedTuple}(undef, n_jobs)
    @threads for idx in 1:n_jobs
        sample_idx = 1 + ((idx - 1) % n_samples)
        regimen_idx = 1 + ((idx - 1) ÷ n_samples)
        pmap = copy(base_map)
        for (j, name) in enumerate(TOP_PARAM_NAMES)
            pmap[name] = X[sample_idx, j]
        end
        spec = regimen_specs[regimen_idx]
        results[idx] = simulate_metrics(
            sample_idx,
            spec.name,
            spec.label,
            spec.dose_days,
            spec.dose_mg,
            pmap,
            alg,
            saveat,
            horizon_days,
        )
    end

    metrics_df = DataFrame(results)
    CSV.write(joinpath(out_dir, "topparam_vpop_samples.csv"), samples_df)
    CSV.write(joinpath(out_dir, "topparam_vpop_metrics.csv"), metrics_df)

    regimen_df = DataFrame(
        regimen = String[],
        regimen_label = String[],
        dose_index = Int[],
        dose_day = Float64[],
        dose_mg = Float64[],
    )
    for spec in regimen_specs
        for (k, (t, d)) in enumerate(zip(spec.dose_days, spec.dose_mg))
            push!(regimen_df, (spec.name, spec.label, k, t, d))
        end
    end
    CSV.write(joinpath(out_dir, "topparam_vpop_regimens.csv"), regimen_df)

    meta = Dict(
        "out_dir" => out_dir,
        "n_samples" => n_samples,
        "seed" => seed,
        "horizon_days" => horizon_days,
        "save_dt_days" => save_dt,
        "solver" => solver_name,
        "threads" => Threads.nthreads(),
        "bw_kg" => BW_KG,
        "sampled_parameters" => TOP_PARAM_NAMES,
        "sampling_mode" => "lhs_log10_space_then_backtransform",
        "bounds" => Dict(TOP_PARAM_NAMES[i] => Dict("lower" => bounds[i][1], "upper" => bounds[i][2]) for i in eachindex(TOP_PARAM_NAMES)),
        "regimens" => [
            Dict("name" => spec.name, "label" => spec.label, "dose_days" => spec.dose_days, "dose_mg" => spec.dose_mg)
            for spec in regimen_specs
        ],
    )
    open(joinpath(out_dir, "topparam_vpop_meta.json"), "w") do io
        JSON3.pretty(io, meta)
    end

    n_ok = count(==("ok"), metrics_df.status)
    println("Wrote $(joinpath(out_dir, "topparam_vpop_samples.csv"))")
    println("Wrote $(joinpath(out_dir, "topparam_vpop_metrics.csv"))")
    println("Wrote $(joinpath(out_dir, "topparam_vpop_regimens.csv"))")
    println("n_ok=$n_ok / $(nrow(metrics_df))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
