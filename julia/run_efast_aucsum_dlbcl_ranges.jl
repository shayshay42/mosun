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
const DLBCL_VARIANT_IDS = [5, 9, 14, 20, 24, 25, 27, 28]
const SWITCH_PARAMS = Set(["tumor_on", "tissue2on", "tissue3on"])

function salib_fast_sample(
        bounds::Vector{Tuple{Float64, Float64}},
        names::Vector{String},
        N::Int,
        M::Int,
        rng::AbstractRNG;
        input_scale_mode::String = "linear")
    D = length(names)
    if N <= 4 * M^2
        error("Sample size N must satisfy N > 4*M^2 for FAST/eFAST (got N=$N, M=$M).")
    end
    omega = zeros(Int, D)
    omega[1] = floor(Int, (N - 1) / (2 * M))
    m = floor(Int, omega[1] / (2 * M))
    if m < 1
        error("Computed m=$m; increase N or reduce M.")
    end
    if m >= (D - 1)
        omega[2:end] = floor.(Int, range(1, m; length = D - 1))
    else
        omega[2:end] = [mod(i - 1, m) + 1 for i in 1:(D - 1)]
    end

    s = (2 * π / N) .* collect(0:(N - 1))
    X = zeros(Float64, N * D, D)
    omega2 = zeros(Int, D)

    for i in 1:D
        fill!(omega2, 0)
        omega2[i] = omega[1]
        rem_idx = [j for j in 1:D if j != i]
        omega2[rem_idx] = omega[2:end]
        z0 = (i - 1) * N
        ϕ = 2 * π * rand(rng)
        for j in 1:D
            g = @. 0.5 + (1 / π) * asin(sin(omega2[j] * s + ϕ))
            X[(z0 + 1):(z0 + N), j] = g
        end
    end

    param_input_scales = fill("linear", D)
    for j in 1:D
        lb, ub = bounds[j]
        mode = lowercase(input_scale_mode)
        if mode in ("log_positive", "auto_log_positive") && lb > 0.0 && ub > 0.0
            X[:, j] .= 10.0 .^ (log10(lb) .+ X[:, j] .* (log10(ub) - log10(lb)))
            param_input_scales[j] = "log10"
        else
            X[:, j] .= lb .+ X[:, j] .* (ub - lb)
            param_input_scales[j] = "linear"
        end
    end
    return X, omega[1], param_input_scales
end

function salib_fast_analyze(Y::Vector{Float64}, D::Int, N::Int, M::Int)
    if length(Y) != N * D
        error("Expected length(Y)=N*D=$(N * D), got $(length(Y))")
    end
    omega0 = floor(Int, (N - 1) / (2 * M))
    S1 = fill(NaN, D)
    ST = fill(NaN, D)
    for i in 1:D
        y = Y[((i - 1) * N + 1):(i * N)]
        sp = fast_power_spectrum(y)
        V = 2.0 * sum(sp)
        if !(isfinite(V) && V > 0)
            continue
        end
        idx_d1 = (collect(1:M) .* omega0)
        idx_d1 = idx_d1[idx_d1 .<= length(sp)]
        d1 = 2.0 * sum(sp[idx_d1])
        idx_dt = 1:floor(Int, omega0 / 2)
        idx_dt = idx_dt[idx_dt .<= length(sp)]
        dt = 2.0 * sum(sp[idx_dt])
        S1[i] = d1 / V
        ST[i] = 1.0 - dt / V
    end
    return S1, ST
end

function fast_power_spectrum(y::AbstractVector{<:Real})
    n = length(y)
    kmax = ceil(Int, n / 2) - 1
    sp = Vector{Float64}(undef, kmax)
    @inbounds for k in 1:kmax
        re = 0.0
        im = 0.0
        for j in 0:(n - 1)
            θ = -2.0 * π * k * j / n
            yj = Float64(y[j + 1])
            re += yj * cos(θ)
            im += yj * sin(θ)
        end
        sp[k] = (hypot(re, im) / n)^2
    end
    return sp
end

function trapz(t::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(t)
    n == length(y) || throw(ArgumentError("trapz input lengths differ"))
    n < 2 && return 0.0
    acc = 0.0
    @inbounds for i in 1:(n - 1)
        dt = Float64(t[i + 1] - t[i])
        acc += 0.5 * dt * (Float64(y[i]) + Float64(y[i + 1]))
    end
    return acc
end

function build_regimen_8cycle_stepup()
    dose_days = Float64[0.0, 7.0, 14.0, 21.0, 42.0, 63.0, 84.0, 105.0, 126.0, 147.0]
    dose_mg = Float64[1.0, 2.0, 60.0, 60.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0]
    return dose_days, dose_mg
end

function resolve_regimen()
    days_env = strip(get(ENV, "EFAST_DOSE_DAYS", ""))
    mg_env = strip(get(ENV, "EFAST_DOSE_MG", ""))
    if !isempty(days_env) || !isempty(mg_env)
        dose_days = parse_number_vector(days_env)
        dose_mg = parse_number_vector(mg_env)
        isempty(dose_days) && error("EFAST_DOSE_DAYS provided but no valid numeric times were parsed")
        length(dose_days) == length(dose_mg) || error("EFAST_DOSE_DAYS and EFAST_DOSE_MG must have equal lengths")
        return dose_days, dose_mg, "custom_env"
    end
    dose_days, dose_mg = build_regimen_8cycle_stepup()
    return dose_days, dose_mg, "default_8cycle_stepup"
end

function build_shipped_dlbcl_base_params(; horizon_days::Float64)
    variant_overrides = load_variant_overrides()
    p = deepcopy(MMC.default_params())
    apply_variant_overrides!(p, DLBCL_VARIANT_IDS, variant_overrides)
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
    return p
end

function load_dlbcl_bounds(bounds_csv::String)
    df = DataFrame(CSV.File(bounds_csv))
    names = String[]
    bounds = Tuple{Float64, Float64}[]
    finals = Dict{String, Float64}()
    for row in eachrow(df)
        name = String(row.parameter)
        haskey(finals, name) && error("Duplicate parameter row for $name in $bounds_csv")
        lb = try
            Float64(row.dlbcl_min_value)
        catch
            NaN
        end
        ub = try
            Float64(row.dlbcl_max_value)
        catch
            NaN
        end
        fv = try
            Float64(row.dlbcl_final_value)
        catch
            NaN
        end
        if isfinite(fv)
            finals[name] = fv
        end
        if !(isfinite(lb) && isfinite(ub)) || !(ub > lb)
            continue
        end
        if name in SWITCH_PARAMS
            continue
        end
        push!(names, name)
        push!(bounds, (lb, ub))
    end
    return names, bounds, finals
end

function build_params_for_sample(base::MMC.MosunParams, param_names::Vector{String}, rowvals::AbstractVector{<:Real}; horizon_days::Float64)
    p = deepcopy(base)
    @inbounds for (j, nm) in enumerate(param_names)
        MMC.set_param!(p, nm, Float64(rowvals[j]))
    end
    if MMC.has_parameter("end_time")
        MMC.set_param!(p, "end_time", horizon_days)
    end
    return p
end

function eval_aucsum_for_sample(
        rowvals::AbstractVector{<:Real},
        param_names::Vector{String},
        base::MMC.MosunParams,
        regimen::MMC.MosunRegimen,
        horizon_days::Float64,
        save_times::Vector{Float64},
        alg)
    p = build_params_for_sample(base, param_names, rowvals; horizon_days = horizon_days)
    built, sol = MMC.solve_regimen(
        regimen,
        p,
        alg;
        tspan = (0.0, horizon_days),
        saveat = save_times,
        callback_mode = :callback,
        abstol = SOLVER_ABSTOL,
        reltol = SOLVER_RELTOL,
    )
    if sol.retcode != SciMLBase.ReturnCode.Success
        return (NaN, NaN, NaN, String(sol.retcode))
    end
    bt_idx = MMC.dynamic_state_index(:Btumor)
    bt = Vector{Float64}(undef, length(sol.t))
    il6 = solution_observable(sol, p, :IL6combo)
    @inbounds for i in eachindex(sol.t)
        bt[i] = Float64(sol.u[i][bt_idx])
    end
    auc_bt = trapz(sol.t, bt)
    auc_il6 = trapz(sol.t, il6)
    return (auc_bt + auc_il6, auc_bt, auc_il6, "ok")
end

function apply_output_transform(y::Float64, transform::AbstractString)
    mode = lowercase(String(transform))
    if !isfinite(y)
        return NaN
    elseif mode == "identity"
        return y
    elseif mode == "log10"
        return log10(max(y, eps(Float64)))
    elseif mode in ("ln", "log")
        return log(max(y, eps(Float64)))
    elseif mode == "log1p"
        return log1p(max(y, 0.0))
    else
        error("Unsupported EFAST_OUTPUT_TRANSFORM=$transform. Use identity, log10, ln, or log1p.")
    end
end

function objective_spec(kind::AbstractString)
    k = lowercase(String(kind))
    if k == "auc_total"
        return (
            kind = k,
            label = "AUC(Btumor) + AUC(IL6combo)",
            stem = "efast_aucsum",
            raw_col = "auc_total_raw",
            used_col = "auc_total_used",
            fail_env = "EFAST_FAIL_FILL_AUCSUM",
        )
    elseif k == "auc_btumor"
        return (
            kind = k,
            label = "AUC(Btumor)",
            stem = "efast_auc_btumor",
            raw_col = "objective_raw",
            used_col = "objective_used",
            fail_env = "EFAST_FAIL_FILL_AUC_BTUMOR",
        )
    elseif k == "auc_il6combo"
        return (
            kind = k,
            label = "AUC(IL6combo)",
            stem = "efast_auc_il6combo",
            raw_col = "objective_raw",
            used_col = "objective_used",
            fail_env = "EFAST_FAIL_FILL_AUC_IL6COMBO",
        )
    else
        error("Unsupported EFAST_OUTPUT_KIND=$kind. Use auc_total, auc_btumor, or auc_il6combo.")
    end
end

function main()
    bounds_csv = get(ENV, "EFAST_BOUNDS_CSV", joinpath(REPO_ROOT, "generated", "figures", "reference", "dlbcl_stack_parameter_ranges_all", "dlbcl_stack_parameter_ranges_all.csv"))
    N = parse(Int, get(ENV, "EFAST_N", "65"))
    M = parse(Int, get(ENV, "EFAST_M", "4"))
    seed = parse(Int, get(ENV, "EFAST_SEED", "20260313"))
    bw_kg = parse(Float64, get(ENV, "EFAST_BW_KG", "70.0"))
    horizon_days = parse(Float64, get(ENV, "EFAST_HORIZON_DAYS", "168.0"))
    save_dt = parse(Float64, get(ENV, "EFAST_SAVE_DT", "0.25"))
    max_evals = parse(Int, get(ENV, "EFAST_MAX_EVALS", "0"))
    param_limit = parse(Int, get(ENV, "EFAST_PARAM_LIMIT", "0"))
    input_scale_mode = get(ENV, "EFAST_INPUT_SCALE_MODE", "linear")
    output_transform = get(ENV, "EFAST_OUTPUT_TRANSFORM", "identity")
    solver_name = get(ENV, "TCE_SOLVER", "qndf")
    alg = make_solver_alg(solver_name)
    obj = objective_spec(get(ENV, "EFAST_OUTPUT_KIND", "auc_total"))
    out_dir_default = joinpath(REPO_ROOT, "generated", "figures", "sensitivity", "$(obj.kind)_dlbcl_ranges")
    out_dir_env = get(ENV, "EFAST_OUT_DIR", out_dir_default)
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(REPO_ROOT, out_dir_env)
    mkpath(out_dir)
    fail_fill_env = get(ENV, obj.fail_env, "")

    param_names, bounds, finals = load_dlbcl_bounds(bounds_csv)
    if param_limit > 0
        param_names = param_names[1:min(param_limit, length(param_names))]
        bounds = bounds[1:length(param_names)]
    end
    D = length(param_names)
    D > 0 || error("No varying DLBCL bounds loaded from $bounds_csv")

    rng = MersenneTwister(seed)
    X, omega0, param_input_scales = salib_fast_sample(bounds, param_names, N, M, rng; input_scale_mode = input_scale_mode)
    n_total = size(X, 1)
    n_eval = max_evals > 0 ? min(max_evals, n_total) : n_total

    dose_days, dose_mg, regimen_label = resolve_regimen()
    regimen = bolus_regimen_from_mg(dose_days, dose_mg; bw_kg = bw_kg, target = :TDBc_ugperkg)
    save_times = dense_save_times(horizon_days, dose_days, save_dt)
    base = build_shipped_dlbcl_base_params(; horizon_days = horizon_days)
    for sw in SWITCH_PARAMS
        if haskey(finals, sw)
            MMC.set_param!(base, sw, finals[sw])
        end
    end

    Y_raw = fill(NaN, n_total)
    auc_bt = fill(NaN, n_total)
    auc_il6 = fill(NaN, n_total)
    status = fill("not_run", n_total)
    done = Atomic{Int}(0)

    println("eFAST run:")
    println("  D=$D, N=$N, total_evals=$n_total, evals_to_run=$n_eval, omega0=$omega0")
    println("  solver=$(solver_name), threads=$(nthreads())")
    println("  output=$(obj.label)")
    println("  input_scale_mode=$(input_scale_mode), output_transform=$(output_transform)")
    println("  bounds_csv=$(bounds_csv)")

    @threads for i in 1:n_eval
        y, bt_auc, il6_auc, st = try
            eval_aucsum_for_sample(X[i, :], param_names, base, regimen, horizon_days, save_times, alg)
        catch err
            @warn "Sample failed" idx = i err
            (NaN, NaN, NaN, "exception")
        end
        Y_raw[i] = y
        auc_bt[i] = bt_auc
        auc_il6[i] = il6_auc
        status[i] = st
        k = atomic_add!(done, 1) + 1
        if k % 25 == 0 || k == n_eval
            println("  progress $k / $n_eval")
        end
    end

    if obj.kind == "auc_btumor"
        Y_raw .= auc_bt
    elseif obj.kind == "auc_il6combo"
        Y_raw .= auc_il6
    end

    Y_used = [apply_output_transform(Float64(y), output_transform) for y in Y_raw]
    finite_y = Y_used[1:n_eval][isfinite.(Y_used[1:n_eval])]
    fail_fill = if !isempty(fail_fill_env)
        parse(Float64, fail_fill_env)
    elseif isempty(finite_y)
        0.0
    else
        maximum(finite_y)
    end

    Y = copy(Y_used)
    failures = 0
    for i in 1:n_eval
        if !isfinite(Y[i])
            Y[i] = fail_fill
            status[i] = status[i] == "not_run" ? "failed_fill" : string(status[i], "_fill")
            failures += 1
        end
    end
    if n_eval < n_total
        Y[(n_eval + 1):end] .= fail_fill
        status[(n_eval + 1):end] .= "not_run_fill"
    end

    S1, ST = salib_fast_analyze(Y, D, N, M)
    ord = sortperm(ST; rev = true)
    rank = zeros(Int, D)
    for (r, idx) in enumerate(ord)
        rank[idx] = r
    end

    idx_rows = DataFrame(
        parameter = param_names,
        lower_bound = [b[1] for b in bounds],
        upper_bound = [b[2] for b in bounds],
        S1 = S1,
        ST = ST,
        rank_ST = rank,
    )
    sort!(idx_rows, :rank_ST)

    sample_rows = DataFrame(sample_idx = collect(1:n_total))
    for (j, nm) in enumerate(param_names)
        sample_rows[!, nm] = X[:, j]
    end
    sample_rows[!, obj.raw_col] = Y_raw
    sample_rows[!, obj.used_col] = Y
    sample_rows[!, "auc_btumor"] = auc_bt
    sample_rows[!, "auc_il6combo"] = auc_il6
    sample_rows[!, "status"] = status

    idx_path = joinpath(out_dir, "$(obj.stem)_indices.csv")
    samples_path = joinpath(out_dir, "$(obj.stem)_samples.csv")
    meta_path = joinpath(out_dir, "$(obj.stem)_meta.json")

    CSV.write(idx_path, idx_rows)
    CSV.write(samples_path, sample_rows)

    finite_bt = auc_bt[1:n_eval][isfinite.(auc_bt[1:n_eval])]
    finite_il6 = auc_il6[1:n_eval][isfinite.(auc_il6[1:n_eval])]

    meta = Dict(
        "method" => "eFAST (SALib-compatible sampling/analysis implementation)",
        "output_kind" => obj.kind,
        "output" => obj.label,
        "output_transform" => String(output_transform),
        "output_components" => ["AUC(Btumor)", "AUC(IL6combo)"],
        "horizon_days" => horizon_days,
        "dose_days" => dose_days,
        "dose_mg" => dose_mg,
        "regimen_label" => regimen_label,
        "n_params" => D,
        "n_base_samples" => N,
        "n_total_evals" => n_total,
        "n_eval_executed" => n_eval,
        "M" => M,
        "omega0" => omega0,
        "solver" => String(solver_name),
        "threads" => nthreads(),
        "failures" => failures,
        "fail_fill_value" => fail_fill,
        "bw_kg" => bw_kg,
        "save_dt_days" => save_dt,
        "save_times_count" => length(save_times),
        "bounds_csv" => String(bounds_csv),
        "variant_ids" => DLBCL_VARIANT_IDS,
        "fixed_non_sampled_values_source" => "Exact shipped MATLAB DLBCL stack final values from phase1_human_variantset_no_dose",
        "sampled_parameter_count" => D,
        "sampled_parameter_names" => param_names,
        "bounds" => [Dict("name" => nm, "lb" => b[1], "ub" => b[2]) for (nm, b) in zip(param_names, bounds)],
        "input_scale_mode" => String(input_scale_mode),
        "sampled_parameter_input_scales" => [Dict("name" => nm, "scale" => sc) for (nm, sc) in zip(param_names, param_input_scales)],
        "auc_btumor_summary" => isempty(finite_bt) ? Dict() : Dict("min" => minimum(finite_bt), "median" => median(finite_bt), "max" => maximum(finite_bt)),
        "auc_il6combo_summary" => isempty(finite_il6) ? Dict() : Dict("min" => minimum(finite_il6), "median" => median(finite_il6), "max" => maximum(finite_il6)),
        "note" => obj.kind == "auc_total" ? "Raw auc_total is expected to be dominated by AUC(Btumor) unless IL6 exposure reaches comparable scale." : "Component-specific eFAST run with all other settings matched to the DLBCL-range regimen run.",
    )
    open(meta_path, "w") do io
        JSON3.pretty(io, meta)
    end

    println("Wrote:")
    println("  " * idx_path)
    println("  " * samples_path)
    println("  " * meta_path)
    println("Top ST parameters:")
    for r in eachrow(idx_rows[1:min(10, nrow(idx_rows)), :])
        println("  ", r.parameter, " : ST=", round(r.ST, digits = 4), " (S1=", round(r.S1, digits = 4), ")")
    end
end

main()
