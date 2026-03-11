using CSV
using DataFrames
using DifferentialEquations
using DiffEqCallbacks
using ModelingToolkit
using RuntimeGeneratedFunctions
using SciMLBase
using SparseArrays
using Statistics
using Sundials

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

const MMC = TCellEngagerQSP.MosunModelCore
const REPO_ROOT = TCellEngagerQSP.REPO_ROOT

RuntimeGeneratedFunctions.init(@__MODULE__)

elapsed_s(f::Function) = (t0 = time_ns(); val = f(); ((time_ns() - t0) / 1e9, val))

function parse_float_list(txt::AbstractString)
    out = Float64[]
    for tok in split(String(txt), ",")
        s = strip(tok)
        isempty(s) && continue
        push!(out, parse(Float64, s))
    end
    return out
end

function build_base_param_map()
    p = joinpath(REPO_ROOT, "generated", "phase1_design_standard", "dlbcl_param_overrides.csv")
    df = DataFrame(CSV.File(p))
    out = Dict{String, Float64}()
    for r in eachrow(df)
        out[String(r.name)] = Float64(r.value)
    end
    return out
end

function row_to_param_overrides(row, cols)
    out = Dict{String, Float64}()
    for cname in cols
        cname == "patient_id" && continue
        raw = row[Symbol(cname)]
        if !ismissing(raw)
            v = try
                Float64(raw)
            catch
                NaN
            end
            if isfinite(v)
                out[cname] = v
            end
        end
    end
    return out
end

function merged_param_map(base::Dict{String,Float64}, overrides::Dict{String,Float64}, end_time::Float64)
    p = copy(base)
    merge!(p, overrides)
    p["PKflag"] = 1.0
    p["fvalidation"] = 0.0
    p["VPid"] = 1.0
    p["end_time"] = end_time
    return p
end

function realdata_regimen_spec()
    ref = DataFrame(CSV.File(joinpath(REPO_ROOT, "generated", "reference", "musun_2022_88patients_SPD_from_vpop_generation_3a50cd1.csv")))
    row = ref[1, :]
    dose_days = parse_float_list(String(row.dose_days))
    dose_mg = parse_float_list(String(row.dose_mg))
    bw_kg = 70.0
    dose_map = Dict{Float64, Float64}()
    for (t, mg) in zip(dose_days, dose_mg)
        dose_map[t] = get(dose_map, t, 0.0) + mg * 1000.0 / bw_kg
    end
    return (
        scenario = String(row.scenario),
        regimen_label = String(row.regimen_label),
        dose_days = dose_days,
        dose_mg = dose_mg,
        dose_map = dose_map,
        horizon_days = Float64(row.horizon_days),
    )
end

function build_parameter_cases(horizon_days::Float64)
    base = build_base_param_map()
    phase1_patients = DataFrame(CSV.File(joinpath(REPO_ROOT, "generated", "phase1_design", "patients.csv")))
    selected = DataFrame(CSV.File(joinpath(REPO_ROOT, "generated", "vpop_pruning", "vpop_wasserstein_spd88_postopt_20260306", "selected_patients.csv")))
    mid_idx = cld(nrow(selected), 2)

    cases = NamedTuple[]
    push!(cases, (
        case_key = "dlbcl_baseline",
        source = "generated/phase1_design_standard/dlbcl_param_overrides.csv",
        pmap = merged_param_map(base, Dict{String, Float64}(), horizon_days),
    ))

    prow = phase1_patients[1, :]
    push!(cases, (
        case_key = "phase1_patient1",
        source = "generated/phase1_design/patients.csv:1",
        pmap = merged_param_map(base, row_to_param_overrides(prow, names(phase1_patients)), horizon_days),
    ))

    srow = selected[mid_idx, :]
    push!(cases, (
        case_key = "selected_patient_mid",
        source = "generated/vpop_pruning/vpop_wasserstein_spd88_postopt_20260306/selected_patients.csv:$(mid_idx + 1)",
        pmap = merged_param_map(base, row_to_param_overrides(srow, names(selected)), horizon_days),
    ))
    return cases
end

function max_norm_abs_err(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real})
    @assert size(X) == size(Y)
    mask = isfinite.(X) .& isfinite.(Y)
    any(mask) || return Inf
    amp = max(maximum(abs.(Y[mask])), 1e-12)
    return maximum(abs.(X[mask] .- Y[mask])) / amp
end

function best_spd_from_bt(bt::AbstractVector{<:Real})
    bt0 = max(Float64(bt[1]), 1e-12)
    spd_pct = @. 100.0 * (Float64(bt) / bt0 - 1.0)
    spd_eval = length(spd_pct) > 1 ? spd_pct[2:end] : spd_pct
    return minimum(spd_eval)
end

function sample_outputs(sol, params::MMC.MosunParams)
    cache = MMC.zero_observables_cache()
    n = length(sol.t)
    bt = Vector{Float64}(undef, n)
    il6 = Vector{Float64}(undef, n)
    drug = Vector{Float64}(undef, n)
    auc = Vector{Float64}(undef, n)
    for i in 1:n
        u = sol.u[i]
        t = sol.t[i]
        bt[i] = MMC.value_at(u, params, t, :Btumor, cache)
        il6[i] = MMC.value_at(u, params, t, :IL6combo, cache)
        drug[i] = MMC.value_at(u, params, t, :TDBc_ugperml, cache)
        auc[i] = MMC.value_at(u, params, t, :TDBc_ugperml_AUC, cache)
    end
    return bt, il6, drug, auc
end

function make_daily_saveat(horizon_days::Float64)
    n_days = Int(round(horizon_days))
    return collect(0.0:1.0:Float64(n_days))
end

function solver_configurations()
    return [
        (label = "production_ctx_cvode_bdf", kind = :production_ctx, alg = CVODE_BDF(), jacobian = :none, post_event_dt = nothing),
        (label = "vector_cvode_bdf", kind = :vector, alg = CVODE_BDF(), jacobian = :none, post_event_dt = nothing),
        (label = "vector_cvode_bdf_postdt_1e-2", kind = :vector, alg = CVODE_BDF(), jacobian = :none, post_event_dt = 1e-2),
        (label = "vector_cvode_bdf_postdt_1e-3", kind = :vector, alg = CVODE_BDF(), jacobian = :none, post_event_dt = 1e-3),
        (label = "vector_cvode_bdf_postdt_1e-4", kind = :vector, alg = CVODE_BDF(), jacobian = :none, post_event_dt = 1e-4),
        (label = "vector_qndf", kind = :vector, alg = QNDF(autodiff = false), jacobian = :none, post_event_dt = nothing),
        (label = "vector_qndf_postdt_1e-2", kind = :vector, alg = QNDF(autodiff = false), jacobian = :none, post_event_dt = 1e-2),
        (label = "vector_qndf_postdt_1e-3", kind = :vector, alg = QNDF(autodiff = false), jacobian = :none, post_event_dt = 1e-3),
        (label = "vector_qndf_postdt_1e-4", kind = :vector, alg = QNDF(autodiff = false), jacobian = :none, post_event_dt = 1e-4),
        (label = "vector_qndf_mtk_dense", kind = :vector, alg = QNDF(autodiff = false), jacobian = :dense, post_event_dt = nothing),
        (label = "vector_qndf_mtk_sparse", kind = :vector, alg = QNDF(autodiff = false), jacobian = :sparse, post_event_dt = nothing),
        (label = "vector_trbdf2", kind = :vector, alg = TRBDF2(autodiff = false), jacobian = :none, post_event_dt = nothing),
        (label = "vector_rodas4p", kind = :vector, alg = Rodas4P(autodiff = false), jacobian = :none, post_event_dt = nothing),
        (label = "vector_rodas4p_postdt_1e-3", kind = :vector, alg = Rodas4P(autodiff = false), jacobian = :none, post_event_dt = 1e-3),
        (label = "vector_rodas4p_mtk_dense", kind = :vector, alg = Rodas4P(autodiff = false), jacobian = :dense, post_event_dt = nothing),
        (label = "vector_rodas4p_mtk_sparse", kind = :vector, alg = Rodas4P(autodiff = false), jacobian = :sparse, post_event_dt = nothing),
        (label = "vector_kencarp4", kind = :vector, alg = KenCarp4(autodiff = false), jacobian = :none, post_event_dt = nothing),
        (label = "vector_kencarp4_postdt_1e-3", kind = :vector, alg = KenCarp4(autodiff = false), jacobian = :none, post_event_dt = 1e-3),
        (label = "vector_kencarp4_mtk_sparse", kind = :vector, alg = KenCarp4(autodiff = false), jacobian = :sparse, post_event_dt = nothing),
        (label = "vector_tsit5", kind = :vector, alg = Tsit5(), jacobian = :none, post_event_dt = nothing),
    ]
end

function build_vector_callback(target_idx::Int, dose_map::Dict{Float64,Float64}; post_event_proposed_dt = nothing)
    local_dose_map = copy(dose_map)
    t0_dose = get(local_dose_map, 0.0, 0.0)
    if t0_dose != 0.0
        delete!(local_dose_map, 0.0)
    end
    cb_times = sort(collect(keys(local_dose_map)))
    callback = nothing
    if !isempty(cb_times)
        function affect!(integrator)
            integrator.u[target_idx] += local_dose_map[Float64(integrator.t)]
            if !isnothing(post_event_proposed_dt)
                SciMLBase.set_proposed_dt!(integrator, post_event_proposed_dt)
            end
        end
        callback = PresetTimeCallback(cb_times, affect!; save_positions = (false, false))
    end
    return t0_dose, callback, cb_times
end

function build_vector_problem(params::MMC.MosunParams, dose_map::Dict{Float64,Float64}, horizon_days::Float64, saveat; jac_bundle = nothing, post_event_proposed_dt = nothing)
    u0 = MMC.pack_state(MMC.initial_state(params))
    target_idx = MMC.dynamic_state_index(:TDBc_ugperkg)
    t0_dose, callback, cb_times = build_vector_callback(target_idx, dose_map; post_event_proposed_dt = post_event_proposed_dt)
    t0_dose != 0.0 && (u0[target_idx] += t0_dose)
    pvec = MMC.pack_params(params)
    rhs =
        if isnothing(jac_bundle)
            MMC.mosun_rhs_vector!
        elseif isnothing(jac_bundle.proto)
            ODEFunction(MMC.mosun_rhs_vector!; jac = jac_bundle.jac!)
        else
            ODEFunction(MMC.mosun_rhs_vector!; jac = jac_bundle.jac!, jac_prototype = copy(jac_bundle.proto))
        end
    prob = ODEProblem(rhs, u0, (0.0, horizon_days), pvec)
    return (
        prob = prob,
        callback = callback,
        tstops = cb_times,
        d_discontinuities = copy(cb_times),
        saveat = saveat,
    )
end

function solve_vector_problem(built, alg; abstol, reltol)
    return solve(
        built.prob,
        alg;
        abstol = abstol,
        reltol = reltol,
        callback = built.callback,
        tstops = built.tstops,
        d_discontinuities = built.d_discontinuities,
        saveat = built.saveat,
        save_everystep = false,
        save_start = true,
    )
end

function build_production_problem(regimen, params::MMC.MosunParams, horizon_days::Float64, saveat; post_event_proposed_dt = nothing)
    return MMC.build_problem(
        regimen,
        params;
        tspan = (0.0, horizon_days),
        saveat = saveat,
        callback_mode = :callback,
        post_event_proposed_dt = post_event_proposed_dt,
    )
end

function solve_production_problem(built, alg; abstol, reltol)
    return MMC.solve_problem(
        built,
        alg;
        abstol = abstol,
        reltol = reltol,
        save_everystep = false,
        save_start = true,
    )
end

function make_dense_jacobian(u0_template::Vector{Float64}, p_template::Vector{Float64})
    probe_prob = ODEProblem(MMC.mosun_rhs_vector!, copy(u0_template), (0.0, 1.0), copy(p_template))
    sys = ModelingToolkit.modelingtoolkitize(probe_prob)
    (_, jac_expr_inplace) = ModelingToolkit.generate_jacobian(sys; sparse = false, simplify = false)
    jac_vec! = RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, TCellEngagerQSP.normalize_function_expr(jac_expr_inplace))
    n = length(u0_template)
    jac_work = zeros(Float64, n * n)

    function jac!(J, u, p, t)
        jac_vec!(jac_work, u, p, t)
        if J isa AbstractMatrix
            @inbounds for j in 1:n, i in 1:n
                v = jac_work[(j - 1) * n + i]
                J[i, j] = isfinite(v) ? v : 0.0
            end
        else
            @inbounds for i in eachindex(jac_work)
                v = jac_work[i]
                J[i] = isfinite(v) ? v : 0.0
            end
        end
        return nothing
    end

    return (jac! = jac!, proto = nothing, n = n, nnz = n * n)
end

function make_sparse_jacobian(u0_template::Vector{Float64}, p_template::Vector{Float64})
    probe_prob = ODEProblem(MMC.mosun_rhs_vector!, copy(u0_template), (0.0, 1.0), copy(p_template))
    sys = ModelingToolkit.modelingtoolkitize(probe_prob)
    (jac_expr_outofplace, jac_expr_inplace) = ModelingToolkit.generate_jacobian(sys; sparse = true, simplify = false)
    jac_sparse! = RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, TCellEngagerQSP.normalize_function_expr(jac_expr_inplace))
    jac_sparse_outofplace = RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, TCellEngagerQSP.normalize_function_expr(jac_expr_outofplace))
    proto = jac_sparse_outofplace(copy(u0_template), copy(p_template), 0.0)
    proto isa SparseMatrixCSC || (proto = sparse(proto))
    fill!(proto.nzval, 0.0)

    function jac!(J, u, p, t)
        if J isa SparseMatrixCSC
            fill!(J.nzval, 0.0)
            jac_sparse!(J, u, p, t)
            @inbounds for i in eachindex(J.nzval)
                v = J.nzval[i]
                J.nzval[i] = isfinite(v) ? v : 0.0
            end
        elseif J isa AbstractMatrix
            Jsp = copy(proto)
            fill!(Jsp.nzval, 0.0)
            jac_sparse!(Jsp, u, p, t)
            @inbounds for i in eachindex(Jsp.nzval)
                v = Jsp.nzval[i]
                Jsp.nzval[i] = isfinite(v) ? v : 0.0
            end
            copyto!(J, Matrix(Jsp))
        else
            Jsp = copy(proto)
            fill!(Jsp.nzval, 0.0)
            jac_sparse!(Jsp, u, p, t)
            @inbounds for i in eachindex(Jsp.nzval)
                v = Jsp.nzval[i]
                Jsp.nzval[i] = isfinite(v) ? v : 0.0
            end
            copyto!(J, vec(Matrix(Jsp)))
        end
        return nothing
    end

    return (jac! = jac!, proto = proto, n = size(proto, 1), nnz = nnz(proto))
end

function build_symbolic_bundles(u0_template::Vector{Float64}, p_template::Vector{Float64})
    setup = DataFrame(
        jacobian = String[],
        status = String[],
        setup_time_s = Float64[],
        n = Int[],
        nnz = Int[],
        note = String[],
    )
    dense_bundle = nothing
    sparse_bundle = nothing

    dense_time, dense_result = elapsed_s() do
        make_dense_jacobian(u0_template, p_template)
    end
    dense_bundle = dense_result
    push!(setup, ("dense", "ok", dense_time, dense_bundle.n, dense_bundle.nnz, "ModelingToolkit dense symbolic Jacobian for mosun_rhs_vector!"))

    sparse_time, sparse_result = elapsed_s() do
        make_sparse_jacobian(u0_template, p_template)
    end
    sparse_bundle = sparse_result
    push!(setup, ("sparse", "ok", sparse_time, sparse_bundle.n, sparse_bundle.nnz, "ModelingToolkit sparse symbolic Jacobian for mosun_rhs_vector!"))
    return setup, Dict(:dense => dense_bundle, :sparse => sparse_bundle)
end

function trajectory_matrix(bt::Vector{Float64}, il6::Vector{Float64}, drug::Vector{Float64})
    return hcat(bt, il6, drug)
end

function main()
    out_dir_env = get(ENV, "EXPLICIT_CORE_BENCH_OUT_DIR", joinpath("generated", "benchmarks", "explicit_core_solver_sweep_20260310"))
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(REPO_ROOT, out_dir_env)
    mkpath(out_dir)

    repeats = parse(Int, get(ENV, "EXPLICIT_CORE_BENCH_REPEATS", "3"))
    abstol = parse(Float64, get(ENV, "EXPLICIT_CORE_BENCH_ABSTOL", string(TCellEngagerQSP.SOLVER_ABSTOL)))
    reltol = parse(Float64, get(ENV, "EXPLICIT_CORE_BENCH_RELTOL", string(TCellEngagerQSP.SOLVER_RELTOL)))

    regimen_spec = realdata_regimen_spec()
    regimen = MMC.bolus_regimen(:TDBc_ugperkg, regimen_spec.dose_map)
    saveat = make_daily_saveat(regimen_spec.horizon_days)
    parameter_cases = build_parameter_cases(regimen_spec.horizon_days)

    template_params = MMC.params_from_dict(parameter_cases[1].pmap; ignore_unknown = true)
    template_u0 = MMC.pack_state(MMC.initial_state(template_params))
    template_p = MMC.pack_params(template_params)

    setup_df, jacobian_bundles = build_symbolic_bundles(template_u0, template_p)
    CSV.write(joinpath(out_dir, "setup_summary.csv"), setup_df)

    runs = DataFrame(
        config = String[],
        kind = String[],
        jacobian = String[],
        parameter_case = String[],
        parameter_source = String[],
        rep = Int[],
        is_warmup = Bool[],
        solve_time_s = Float64[],
        status = String[],
        retcode = String[],
        best_spd_pct = Float64[],
        peak_il6 = Float64[],
        final_auc_tdbc = Float64[],
        max_rel_err_ref = Float64[],
        best_spd_abs_diff_ref = Float64[],
        peak_il6_abs_diff_ref = Float64[],
        auc_abs_diff_ref = Float64[],
    )

    summary = DataFrame(
        config = String[],
        kind = String[],
        jacobian = String[],
        parameter_case = String[],
        parameter_source = String[],
        n_success = Int[],
        warmup_time_s = Float64[],
        mean_solve_time_s = Float64[],
        median_solve_time_s = Float64[],
        min_solve_time_s = Float64[],
        max_solve_time_s = Float64[],
        mean_best_spd_pct = Float64[],
        mean_peak_il6 = Float64[],
        mean_final_auc_tdbc = Float64[],
        mean_max_rel_err_ref = Float64[],
        mean_best_spd_abs_diff_ref = Float64[],
        mean_peak_il6_abs_diff_ref = Float64[],
        mean_auc_abs_diff_ref = Float64[],
    )

    reference_outputs = Dict{String, NamedTuple}()
    configs = solver_configurations()

    for case in parameter_cases
        params = MMC.params_from_dict(case.pmap; ignore_unknown = true)
        ref_built = build_production_problem(regimen, params, regimen_spec.horizon_days, saveat)
        ref_sol = solve_production_problem(ref_built, CVODE_BDF(); abstol = abstol, reltol = reltol)
        ref_sol.retcode == SciMLBase.ReturnCode.Success || error("Reference production solve failed for $(case.case_key) with retcode=$(ref_sol.retcode)")
        ref_bt, ref_il6, ref_drug, ref_auc = sample_outputs(ref_sol, params)
        reference_outputs[case.case_key] = (
            matrix = trajectory_matrix(ref_bt, ref_il6, ref_drug),
            best_spd = best_spd_from_bt(ref_bt),
            peak_il6 = maximum(ref_il6),
            final_auc = ref_auc[end],
        )
    end

    for cfg in configs
        for case in parameter_cases
            params = MMC.params_from_dict(case.pmap; ignore_unknown = true)
            jac_bundle = get(jacobian_bundles, cfg.jacobian, nothing)
            build_result =
                if cfg.kind == :production_ctx
                    build_production_problem(regimen, params, regimen_spec.horizon_days, saveat; post_event_proposed_dt = cfg.post_event_dt)
                elseif cfg.kind == :vector
                    build_vector_problem(params, regimen_spec.dose_map, regimen_spec.horizon_days, saveat; jac_bundle = jac_bundle, post_event_proposed_dt = cfg.post_event_dt)
                else
                    error("Unsupported benchmark kind $(cfg.kind)")
                end

            solve_fn =
                if cfg.kind == :production_ctx
                    () -> solve_production_problem(build_result, cfg.alg; abstol = abstol, reltol = reltol)
                else
                    () -> solve_vector_problem(build_result, cfg.alg; abstol = abstol, reltol = reltol)
                end

            warmup_time, warmup_sol = elapsed_s(solve_fn)
            ref = reference_outputs[case.case_key]

            function metrics_from_sol(sol)
                bt, il6, drug, auc = sample_outputs(sol, params)
                traj = trajectory_matrix(bt, il6, drug)
                return (
                    best_spd = best_spd_from_bt(bt),
                    peak_il6 = maximum(il6),
                    final_auc = auc[end],
                    max_rel_err = max_norm_abs_err(traj, ref.matrix),
                )
            end

            if warmup_sol.retcode == SciMLBase.ReturnCode.Success
                m = metrics_from_sol(warmup_sol)
                push!(runs, (
                    cfg.label,
                    String(cfg.kind),
                    String(cfg.jacobian),
                    case.case_key,
                    case.source,
                    0,
                    true,
                    warmup_time,
                    "ok",
                    string(warmup_sol.retcode),
                    m.best_spd,
                    m.peak_il6,
                    m.final_auc,
                    m.max_rel_err,
                    abs(m.best_spd - ref.best_spd),
                    abs(m.peak_il6 - ref.peak_il6),
                    abs(m.final_auc - ref.final_auc),
                ))
            else
                push!(runs, (
                    cfg.label,
                    String(cfg.kind),
                    String(cfg.jacobian),
                    case.case_key,
                    case.source,
                    0,
                    true,
                    warmup_time,
                    "solve_error",
                    string(warmup_sol.retcode),
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                ))
            end

            for rep in 1:repeats
                solve_time, sol = elapsed_s(solve_fn)
                if sol.retcode == SciMLBase.ReturnCode.Success
                    m = metrics_from_sol(sol)
                    push!(runs, (
                        cfg.label,
                        String(cfg.kind),
                        String(cfg.jacobian),
                        case.case_key,
                        case.source,
                        rep,
                        false,
                        solve_time,
                        "ok",
                        string(sol.retcode),
                        m.best_spd,
                        m.peak_il6,
                        m.final_auc,
                        m.max_rel_err,
                        abs(m.best_spd - ref.best_spd),
                        abs(m.peak_il6 - ref.peak_il6),
                        abs(m.final_auc - ref.final_auc),
                    ))
                else
                    push!(runs, (
                        cfg.label,
                        String(cfg.kind),
                        String(cfg.jacobian),
                        case.case_key,
                        case.source,
                        rep,
                        false,
                        solve_time,
                        "solve_error",
                        string(sol.retcode),
                        NaN,
                        NaN,
                        NaN,
                        NaN,
                        NaN,
                        NaN,
                        NaN,
                    ))
                end
            end
        end
    end

    grouped = groupby(runs, [:config, :kind, :jacobian, :parameter_case, :parameter_source])
    for sub in grouped
        warm = sub[sub.is_warmup .== true, :]
        timed = sub[(sub.is_warmup .== false) .& (sub.status .== "ok"), :]
        push!(summary, (
            String(sub.config[1]),
            String(sub.kind[1]),
            String(sub.jacobian[1]),
            String(sub.parameter_case[1]),
            String(sub.parameter_source[1]),
            nrow(timed),
            isempty(warm) ? NaN : Float64(warm.solve_time_s[1]),
            isempty(timed) ? NaN : mean(timed.solve_time_s),
            isempty(timed) ? NaN : median(timed.solve_time_s),
            isempty(timed) ? NaN : minimum(timed.solve_time_s),
            isempty(timed) ? NaN : maximum(timed.solve_time_s),
            isempty(timed) ? NaN : mean(timed.best_spd_pct),
            isempty(timed) ? NaN : mean(timed.peak_il6),
            isempty(timed) ? NaN : mean(timed.final_auc_tdbc),
            isempty(timed) ? NaN : mean(timed.max_rel_err_ref),
            isempty(timed) ? NaN : mean(timed.best_spd_abs_diff_ref),
            isempty(timed) ? NaN : mean(timed.peak_il6_abs_diff_ref),
            isempty(timed) ? NaN : mean(timed.auc_abs_diff_ref),
        ))
    end

    aggregate = combine(
        groupby(summary, [:config, :kind, :jacobian]),
        :n_success => sum => :n_success_total,
        :mean_solve_time_s => mean => :mean_solve_time_s,
        :median_solve_time_s => mean => :mean_median_solve_time_s,
        :min_solve_time_s => minimum => :best_case_solve_time_s,
        :max_solve_time_s => maximum => :worst_case_solve_time_s,
        :warmup_time_s => mean => :mean_warmup_time_s,
        :mean_max_rel_err_ref => mean => :mean_max_rel_err_ref,
        :mean_best_spd_abs_diff_ref => mean => :mean_best_spd_abs_diff_ref,
        :mean_peak_il6_abs_diff_ref => mean => :mean_peak_il6_abs_diff_ref,
        :mean_auc_abs_diff_ref => mean => :mean_auc_abs_diff_ref,
    )
    sort!(aggregate, :mean_solve_time_s)

    CSV.write(joinpath(out_dir, "runs.csv"), runs)
    CSV.write(joinpath(out_dir, "summary_by_case.csv"), summary)
    CSV.write(joinpath(out_dir, "summary_aggregate.csv"), aggregate)

    meta = DataFrame(
        key = [
            "scenario",
            "regimen_label",
            "dose_days",
            "dose_mg",
            "horizon_days",
            "saveat_days",
            "repeats",
            "abstol",
            "reltol",
        ],
        value = [
            regimen_spec.scenario,
            regimen_spec.regimen_label,
            join(regimen_spec.dose_days, ","),
            join(regimen_spec.dose_mg, ","),
            string(regimen_spec.horizon_days),
            join(saveat, ","),
            string(repeats),
            string(abstol),
            string(reltol),
        ],
    )
    CSV.write(joinpath(out_dir, "meta.csv"), meta)

    println("Wrote benchmark results to $(out_dir)")
end

main()
