using CSV
using DataFrames
using JSON3
using Random
using SciMLBase
using DifferentialEquations
using Sundials
using Optim
using Statistics
using Base.Threads

include(joinpath(@__DIR__, "src", "MosunModelCore.jl"))
using .MosunModelCore

const MMC = MosunModelCore
const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const SOLVER_ABSTOL = parse(Float64, get(ENV, "TCE_ABSTOL", "1e-8"))
const SOLVER_RELTOL = parse(Float64, get(ENV, "TCE_RELTOL", "1e-5"))
const VARIANT_IDS = [5, 9, 14, 20, 24, 25, 27, 28]
const FIXED_SWITCH_PARAMS = Dict("tissue2on" => 1.0, "tissue3on" => 1.0, "tumor_on" => 1.0)

struct CohortSample
    sample_idx::Int
    params::MMC.MosunParams
end

struct RefLossConfig
    eps_tumor::Float64
    eps_il6::Float64
    beta_tumor::Float64
    beta_il6::Float64
    rho::Float64
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

function make_solver_alg(solver_name::AbstractString)
    solver = lowercase(String(solver_name))
    if solver == "cvode_bdf"
        return CVODE_BDF()
    elseif solver == "qndf"
        return QNDF(autodiff = false)
    elseif solver == "rodas4p"
        return Rodas4P(autodiff = false)
    elseif solver == "kencarp4"
        return KenCarp4(autodiff = false)
    elseif solver == "vcabm"
        return VCABM()
    elseif solver == "tsit5"
        return Tsit5()
    else
        error("Unsupported solver=$solver")
    end
end

function load_param_defaults()
    path = joinpath(REPO_ROOT, "generated", "model_reconstruction_bundle", "parameters.tsv")
    df = DataFrame(CSV.File(path; delim = '\t'))
    return String.(df.name), Float64.(df.value)
end

function load_variant_overrides()
    path = joinpath(REPO_ROOT, "generated", "model_reconstruction_bundle", "variants.tsv")
    df = DataFrame(CSV.File(path; delim = '\t'))
    out = Dict{Int, Vector{Pair{String, Float64}}}()
    for row in eachrow(df)
        action = String(row.action)
        is_numeric = Int(row.is_numeric)
        if action != "parameter" || is_numeric != 1
            continue
        end
        vid = Int(row.variant_idx)
        pname = String(row.class)
        pval = Float64(row.value_numeric)
        push!(get!(out, vid, Pair{String, Float64}[]), pname => pval)
    end
    return out
end

function trapz_auc(t::Vector{Float64}, y::Vector{Float64})
    n = length(t)
    n == length(y) || error("trapz input lengths differ")
    n < 2 && return 0.0
    acc = 0.0
    @inbounds for i in 1:(n - 1)
        dt = t[i + 1] - t[i]
        acc += 0.5 * dt * (y[i] + y[i + 1])
    end
    return acc
end

function dense_save_times(horizon_days::Float64, dose_days::Vector{Float64}, save_dt::Float64)
    base = collect(0.0:save_dt:horizon_days)
    push!(base, horizon_days)
    append!(base, dose_days)
    sort!(base)
    return unique(base)
end

function regimen_from_mg(dose_days::Vector{Float64}, dose_mg::Vector{Float64}; bw_kg::Float64)
    length(dose_days) == length(dose_mg) || error("dose_days and dose_mg must have equal lengths")
    dose_map = Dict{Float64, Float64}()
    for (t, mg) in zip(dose_days, dose_mg)
        amt = Float64(mg) * 1000.0 / bw_kg
        if !iszero(amt)
            dose_map[Float64(t)] = amt
        end
    end
    return MMC.bolus_regimen(:TDBc_ugperkg, dose_map)
end

function build_base_params(; horizon_days::Float64)
    param_names, param_defaults = load_param_defaults()
    p = MMC.params_from_named_values(param_names, param_defaults)
    variant_overrides = load_variant_overrides()
    for vid in VARIANT_IDS
        for (nm, val) in get(variant_overrides, vid, Pair{String, Float64}[])
            if MMC.has_parameter(nm)
                MMC.set_param!(p, nm, val)
            end
        end
    end
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
    for (nm, val) in FIXED_SWITCH_PARAMS
        if MMC.has_parameter(nm)
            MMC.set_param!(p, nm, val)
        end
    end
    return p
end

function apply_sample_to_params(base::MMC.MosunParams, sample_row::NamedTuple, sampled_param_names::Vector{String}; horizon_days::Float64)
    p = deepcopy(base)
    for nm in sampled_param_names
        MMC.set_param!(p, nm, Float64(getproperty(sample_row, Symbol(nm))))
    end
    if MMC.has_parameter("end_time")
        MMC.set_param!(p, "end_time", horizon_days)
    end
    return p
end

function load_lhs_cohort(sample_csv::AbstractString, horizon_days::Float64; max_samples::Int)
    sample_df = DataFrame(CSV.File(sample_csv))
    sampled_param_names = [String(nm) for nm in names(sample_df) if String(nm) != "sample_idx" && MMC.has_parameter(String(nm))]
    selected_df = sample_df[1:min(max_samples, nrow(sample_df)), :]
    rows = [NamedTuple(r) for r in eachrow(selected_df)]
    base = build_base_params(; horizon_days = horizon_days)
    cohort = CohortSample[]
    for row in rows
        sample_idx = Int(round(Float64(getproperty(row, :sample_idx))))
        params = apply_sample_to_params(base, row, sampled_param_names; horizon_days = horizon_days)
        push!(cohort, CohortSample(sample_idx, params))
    end
    return cohort, sampled_param_names
end

function solve_regimen(params::MMC.MosunParams, dose_days::Vector{Float64}, dose_mg::Vector{Float64}, horizon_days::Float64, save_dt::Float64, alg; bw_kg::Float64)
    regimen = regimen_from_mg(dose_days, dose_mg; bw_kg = bw_kg)
    saveat = dense_save_times(horizon_days, dose_days, save_dt)
    built = MMC.build_problem(regimen, params; tspan = (0.0, horizon_days), saveat = saveat, callback_mode = :callback)
    sol = MMC.solve_problem(built, alg; abstol = SOLVER_ABSTOL, reltol = SOLVER_RELTOL)
    sol.retcode == SciMLBase.ReturnCode.Success || error("Solve failed with retcode=$(sol.retcode)")
    return sol
end

function extract_metrics(sol, params::MMC.MosunParams)
    bt_idx = MMC.dynamic_state_index(:Btumor)
    cache = MMC.zero_observables_cache()
    n = length(sol.t)
    t = Vector{Float64}(undef, n)
    bt = Vector{Float64}(undef, n)
    il6 = Vector{Float64}(undef, n)
    @inbounds for i in eachindex(sol.t)
        t[i] = Float64(sol.t[i])
        u = sol.u[i]
        bt[i] = Float64(u[bt_idx])
        MMC.update_observables!(cache, u, params, t[i])
        il6[i] = Float64(cache.IL6combo)
    end
    return (
        btumor_auc = trapz_auc(t, bt),
        il6combo_auc = trapz_auc(t, il6),
        peak_il6combo = maximum(il6),
        final_btumor = bt[end],
        best_spd_pct = 100.0 * (minimum(bt) / max(bt[1], 1e-30) - 1.0),
    )
end

function evaluate_cohort_metrics(
    cohort::Vector{CohortSample},
    dose_days::Vector{Float64},
    dose_mg::Vector{Float64},
    horizon_days::Float64,
    save_dt::Float64,
    alg;
    bw_kg::Float64,
)
    n = length(cohort)
    bt_auc = fill(NaN, n)
    il6_auc = fill(NaN, n)
    peak_il6 = fill(NaN, n)
    final_bt = fill(NaN, n)
    best_spd = fill(NaN, n)
    ok = falses(n)
    @threads for i in 1:n
        sample = cohort[i]
        try
            sol = solve_regimen(sample.params, dose_days, dose_mg, horizon_days, save_dt, alg; bw_kg = bw_kg)
            m = extract_metrics(sol, sample.params)
            bt_auc[i] = m.btumor_auc
            il6_auc[i] = m.il6combo_auc
            peak_il6[i] = m.peak_il6combo
            final_bt[i] = m.final_btumor
            best_spd[i] = m.best_spd_pct
            ok[i] = true
        catch
            ok[i] = false
        end
    end
    return (ok = ok, btumor_auc = bt_auc, il6combo_auc = il6_auc, peak_il6 = peak_il6, final_btumor = final_bt, best_spd_pct = best_spd)
end

function stable_softplus(x::Float64)
    if x > 30.0
        return x
    elseif x < -30.0
        return exp(x)
    else
        return log1p(exp(x))
    end
end

softplus_vec(x::AbstractVector{<:Real}) = [stable_softplus(Float64(v)) for v in x]

function inv_softplus(y::Float64)
    y <= 1e-8 && return log(max(y, 1e-12))
    if y > 30.0
        return y
    end
    return log(expm1(y))
end

inv_softplus_vec(y::AbstractVector{<:Real}) = [inv_softplus(Float64(v)) for v in y]

function entropic_risk(x::Vector{Float64}, beta::Float64)
    if isempty(x)
        return Inf
    elseif beta <= 1e-12
        return mean(x)
    end
    m = maximum(x)
    return m + log(mean(exp.(beta .* (x .- m)))) / beta
end

function objective_pair(candidate_metrics, ref_metrics, cfg::RefLossConfig)
    tumor_logratio = log.((candidate_metrics.btumor_auc .+ cfg.eps_tumor) ./ (ref_metrics.btumor_auc .+ cfg.eps_tumor))
    il6_logratio = log.((candidate_metrics.il6combo_auc .+ cfg.eps_il6) ./ (ref_metrics.il6combo_auc .+ cfg.eps_il6))
    ft = entropic_risk(tumor_logratio, cfg.beta_tumor)
    fc = entropic_risk(il6_logratio, cfg.beta_il6)
    return (
        tumor_loss = ft,
        il6_loss = fc,
        tumor_mean_logratio = mean(tumor_logratio),
        il6_mean_logratio = mean(il6_logratio),
    )
end

function scalarize(loss_pair, wt::Float64, cfg::RefLossConfig)
    wc = 1.0 - wt
    return max(wt * loss_pair.tumor_loss, wc * loss_pair.il6_loss) +
           cfg.rho * (wt * loss_pair.tumor_loss + wc * loss_pair.il6_loss)
end

function jittered_start(rng::AbstractRNG, ref_doses::Vector{Float64}, scale::Float64)
    return ref_doses .* exp.(scale .* randn(rng, length(ref_doses)))
end

function nondominated_mask(ft::Vector{Float64}, fc::Vector{Float64})
    n = length(ft)
    keep = trues(n)
    for i in 1:n
        for j in 1:n
            i == j && continue
            if ft[j] <= ft[i] && fc[j] <= fc[i] && (ft[j] < ft[i] || fc[j] < fc[i])
                keep[i] = false
                break
            end
        end
    end
    return keep
end

function main()
    sample_csv = get(ENV, "COHORT_PARETO_SAMPLE_CSV", joinpath(REPO_ROOT, "generated", "figures", "reference", "untreated_dlbcl_variant_range_lhs200", "untreated_dlbcl_lhs_samples.csv"))
    out_dir = get(ENV, "COHORT_PARETO_OUT_DIR", joinpath(REPO_ROOT, "generated", "figures", "optimization", "lhs200_reference_pareto_blackbox"))
    horizon_days = parse(Float64, get(ENV, "COHORT_PARETO_HORIZON_DAYS", "84.0"))
    dose_days = parse_float_list(get(ENV, "COHORT_PARETO_DOSE_DAYS", "0,7,14,21,42,63"))
    ref_doses = parse_float_list(get(ENV, "COHORT_PARETO_REFERENCE_DOSES_MG", "1.6,10,10,20,20,20"))
    length(dose_days) == length(ref_doses) || error("dose_days/reference doses length mismatch")
    save_dt = parse(Float64, get(ENV, "COHORT_PARETO_SAVE_DT", "0.5"))
    bw_kg = parse(Float64, get(ENV, "COHORT_PARETO_BW_KG", "70.0"))
    solver_name = get(ENV, "COHORT_PARETO_SOLVER", "qndf")
    alg = make_solver_alg(solver_name)
    max_samples = parse(Int, get(ENV, "COHORT_PARETO_MAX_SAMPLES", "200"))
    seed = parse(Int, get(ENV, "COHORT_PARETO_SEED", "20260316"))
    weight_grid = parse_float_list(get(ENV, "COHORT_PARETO_WEIGHTS", "0.05,0.2,0.35,0.5,0.65,0.8,0.95"))
    nm_iters = parse(Int, get(ENV, "COHORT_PARETO_NM_ITERS", "80"))
    n_starts = parse(Int, get(ENV, "COHORT_PARETO_NSTARTS", "5"))
    start_jitter = parse(Float64, get(ENV, "COHORT_PARETO_START_JITTER", "0.45"))
    beta_tumor = parse(Float64, get(ENV, "COHORT_PARETO_BETA_TUMOR", "4.0"))
    beta_il6 = parse(Float64, get(ENV, "COHORT_PARETO_BETA_IL6", "6.0"))
    rho = parse(Float64, get(ENV, "COHORT_PARETO_RHO", "1e-3"))
    mkpath(out_dir)

    println("Loading cohort...")
    cohort, sampled_param_names = load_lhs_cohort(sample_csv, horizon_days; max_samples = max_samples)
    println("cohort_n=$(length(cohort)) solver=$(solver_name) threads=$(Threads.nthreads())")

    println("Evaluating reference regimen...")
    ref_metrics = evaluate_cohort_metrics(cohort, dose_days, ref_doses, horizon_days, save_dt, alg; bw_kg = bw_kg)
    n_failed_ref = count(!, ref_metrics.ok)
    n_failed_ref == 0 || error("Reference regimen failed for $(n_failed_ref) cohort members")

    eps_tumor = begin
        pos = filter(>(0.0), ref_metrics.btumor_auc)
        isempty(pos) ? 1e-12 : max(minimum(pos) / 10, 1e-12)
    end
    eps_il6 = begin
        pos = filter(>(0.0), ref_metrics.il6combo_auc)
        isempty(pos) ? 1e-12 : max(minimum(pos) / 10, 1e-12)
    end
    loss_cfg = RefLossConfig(eps_tumor, eps_il6, beta_tumor, beta_il6, rho)

    ref_df = DataFrame(
        sample_idx = getfield.(cohort, :sample_idx),
        btumor_auc = ref_metrics.btumor_auc,
        il6combo_auc = ref_metrics.il6combo_auc,
        peak_il6combo = ref_metrics.peak_il6,
        final_btumor = ref_metrics.final_btumor,
        best_spd_pct = ref_metrics.best_spd_pct,
    )
    CSV.write(joinpath(out_dir, "reference_metrics.csv"), ref_df)

    rng = MersenneTwister(seed)
    ref_x = inv_softplus_vec(ref_doses)
    result_rows = NamedTuple[]
    trace_rows = NamedTuple[]

    for wt in weight_grid
        0.0 <= wt <= 1.0 || error("weights must lie in [0,1]")
        best_scalar = Inf
        best_row = nothing
        for start_idx in 1:n_starts
            start_doses = if start_idx == 1
                copy(ref_doses)
            elseif start_idx == 2
                fill(mean(ref_doses), length(ref_doses))
            else
                jittered_start(rng, ref_doses, start_jitter)
            end
            x0 = inv_softplus_vec(start_doses)
            eval_counter = Ref(0)
            function objective_x(x)
                doses = softplus_vec(x)
                mets = evaluate_cohort_metrics(cohort, dose_days, doses, horizon_days, save_dt, alg; bw_kg = bw_kg)
                n_failed = count(!, mets.ok)
                if n_failed > 0 || any(!isfinite, mets.btumor_auc) || any(!isfinite, mets.il6combo_auc)
                    val = 1e9 + 1e7 * n_failed
                    push!(trace_rows, (weight_tumor = wt, start_idx = start_idx, eval = (eval_counter[] += 1), scalar_loss = val, tumor_loss = NaN, il6_loss = NaN, n_failed = n_failed, doses_mg = join(round.(doses; digits = 6), ",")))
                    return val
                end
                lp = objective_pair(mets, ref_metrics, loss_cfg)
                sval = scalarize(lp, wt, loss_cfg)
                push!(trace_rows, (weight_tumor = wt, start_idx = start_idx, eval = (eval_counter[] += 1), scalar_loss = sval, tumor_loss = lp.tumor_loss, il6_loss = lp.il6_loss, n_failed = 0, doses_mg = join(round.(doses; digits = 6), ",")))
                return sval
            end
            opts = Optim.Options(iterations = nm_iters, show_trace = false, store_trace = false)
            res = optimize(objective_x, x0, NelderMead(), opts)
            xbest = Optim.minimizer(res)
            doses_best = softplus_vec(xbest)
            mets_best = evaluate_cohort_metrics(cohort, dose_days, doses_best, horizon_days, save_dt, alg; bw_kg = bw_kg)
            n_failed_best = count(!, mets_best.ok)
            if n_failed_best > 0
                row = (
                    weight_tumor = wt,
                    start_idx = start_idx,
                    status = "failed",
                    scalar_loss = 1e9 + 1e7 * n_failed_best,
                    tumor_loss = NaN,
                    il6_loss = NaN,
                    tumor_mean_logratio = NaN,
                    il6_mean_logratio = NaN,
                    n_failed = n_failed_best,
                    doses_mg = doses_best,
                )
            else
                lp = objective_pair(mets_best, ref_metrics, loss_cfg)
                sval = scalarize(lp, wt, loss_cfg)
                row = (
                    weight_tumor = wt,
                    start_idx = start_idx,
                    status = string(Optim.converged(res) ? "ok" : "not_converged"),
                    scalar_loss = sval,
                    tumor_loss = lp.tumor_loss,
                    il6_loss = lp.il6_loss,
                    tumor_mean_logratio = lp.tumor_mean_logratio,
                    il6_mean_logratio = lp.il6_mean_logratio,
                    n_failed = 0,
                    doses_mg = doses_best,
                )
                if sval < best_scalar
                    best_scalar = sval
                    best_row = row
                end
            end
            push!(result_rows, row)
        end
        if best_row !== nothing
            println("weight_tumor=$(round(wt; digits=3)) best scalar=$(round(best_row.scalar_loss; digits=5)) tumor_loss=$(round(best_row.tumor_loss; digits=5)) il6_loss=$(round(best_row.il6_loss; digits=5)) doses=$(round.(best_row.doses_mg; digits=3))")
        else
            println("weight_tumor=$(round(wt; digits=3)) no successful starts")
        end
    end

    starts_df = DataFrame(
        weight_tumor = [r.weight_tumor for r in result_rows],
        start_idx = [r.start_idx for r in result_rows],
        status = [r.status for r in result_rows],
        scalar_loss = [r.scalar_loss for r in result_rows],
        tumor_loss = [r.tumor_loss for r in result_rows],
        il6_loss = [r.il6_loss for r in result_rows],
        tumor_mean_logratio = [r.tumor_mean_logratio for r in result_rows],
        il6_mean_logratio = [r.il6_mean_logratio for r in result_rows],
        n_failed = [r.n_failed for r in result_rows],
    )
    for j in eachindex(ref_doses)
        starts_df[!, Symbol("dose$(j)_mg")] = [r.doses_mg[j] for r in result_rows]
    end
    CSV.write(joinpath(out_dir, "optimization_multistart_results.csv"), starts_df)

    best_rows = NamedTuple[]
    for wt in weight_grid
        sub = filter(r -> r.weight_tumor == wt && isfinite(r.scalar_loss), result_rows)
        isempty(sub) && continue
        best = sub[argmin(getfield.(sub, :scalar_loss))]
        push!(best_rows, best)
    end
    ft = [r.tumor_loss for r in best_rows]
    fc = [r.il6_loss for r in best_rows]
    nd = nondominated_mask(ft, fc)
    pareto_df = DataFrame(
        weight_tumor = [r.weight_tumor for r in best_rows],
        scalar_loss = [r.scalar_loss for r in best_rows],
        tumor_loss = ft,
        il6_loss = fc,
        tumor_mean_logratio = [r.tumor_mean_logratio for r in best_rows],
        il6_mean_logratio = [r.il6_mean_logratio for r in best_rows],
        nondominated = nd,
    )
    for j in eachindex(ref_doses)
        pareto_df[!, Symbol("dose$(j)_mg")] = [r.doses_mg[j] for r in best_rows]
    end
    CSV.write(joinpath(out_dir, "pareto_candidates.csv"), pareto_df)

    trace_df = DataFrame(trace_rows)
    CSV.write(joinpath(out_dir, "optimization_trace.csv"), trace_df)

    meta = Dict(
        "sample_csv" => sample_csv,
        "cohort_size" => length(cohort),
        "sampled_parameter_names" => sampled_param_names,
        "reference_regimen_dose_days" => dose_days,
        "reference_regimen_doses_mg" => ref_doses,
        "horizon_days" => horizon_days,
        "solver" => solver_name,
        "save_dt" => save_dt,
        "bw_kg" => bw_kg,
        "solver_abstol" => SOLVER_ABSTOL,
        "solver_reltol" => SOLVER_RELTOL,
        "blackbox_optimizer" => "NelderMead multistart on unconstrained softplus dose parametrization",
        "nm_iterations" => nm_iters,
        "n_starts" => n_starts,
        "weight_grid" => weight_grid,
        "reference_relative_objectives" => Dict(
            "tumor_loss" => "Entropic risk of log((Btumor_AUC + eps_tumor)/(Btumor_AUC_ref + eps_tumor)) across cohort",
            "il6_loss" => "Entropic risk of log((IL6combo_AUC + eps_il6)/(IL6combo_AUC_ref + eps_il6)) across cohort",
            "scalarization" => "Weighted Tchebycheff: max(w*tumor_loss, (1-w)*il6_loss) + rho*(w*tumor_loss + (1-w)*il6_loss)",
        ),
        "loss_parameters" => Dict(
            "eps_tumor" => eps_tumor,
            "eps_il6" => eps_il6,
            "beta_tumor" => beta_tumor,
            "beta_il6" => beta_il6,
            "rho" => rho,
        ),
        "seed" => seed,
    )
    open(joinpath(out_dir, "optimization_meta.json"), "w") do io
        JSON3.pretty(io, meta)
    end

    println("Wrote:")
    println("  " * joinpath(out_dir, "reference_metrics.csv"))
    println("  " * joinpath(out_dir, "optimization_multistart_results.csv"))
    println("  " * joinpath(out_dir, "pareto_candidates.csv"))
    println("  " * joinpath(out_dir, "optimization_trace.csv"))
    println("  " * joinpath(out_dir, "optimization_meta.json"))
end

main()
