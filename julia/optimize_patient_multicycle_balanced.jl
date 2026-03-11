using CSV
using DataFrames
using DifferentialEquations
using FiniteDiff
using ForwardDiff
using JSON3
using Optim
using Printf
using Random
using RuntimeGeneratedFunctions
using SciMLBase
using SciMLSensitivity
using Statistics

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

RuntimeGeneratedFunctions.init(@__MODULE__)

function smoothmax(v::AbstractVector, tau::Real)
    m = maximum(v)
    return m + tau * log(sum(exp.((v .- m) ./ tau)))
end

function parse_float_list(s::AbstractString)
    vals = Float64[]
    for tok in split(String(s), ",")
        st = strip(tok)
        isempty(st) && continue
        push!(vals, parse(Float64, st))
    end
    return vals
end

function robust_scale(x::Vector{Float64}, floor_val::Float64, fallback::Float64)
    vals = filter(v -> isfinite(v) && v > 0.0, x)
    isempty(vals) && return max(floor_val, fallback)
    s = median(vals)
    return max(floor_val, s)
end

function trapz(x::Vector{Float64}, y::AbstractVector{T}) where {T}
    n = length(x)
    n == length(y) || error("trapz x/y length mismatch: $(length(x)) vs $(length(y))")
    n <= 1 && return zero(T)
    s = zero(T)
    for i in 1:(n - 1)
        dt = x[i + 1] - x[i]
        s += T(0.5 * dt) * (y[i + 1] + y[i])
    end
    return s
end

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

struct PatientMultiCycleProblem
    patient_id::Int
    mdl
    rhs_fun::Function
    il6_obs_fun::Function
    alg
    bt_idx::Int
    e_tdbc::Vector{Float64}
    mg_to_ugkg::Float64
    dose_times_days::Vector{Float64}
    fixed_doses_mg::Vector{Float64}
    lower_mg::Vector{Float64}
    upper_mg::Vector{Float64}
    horizon_days::Float64
    tox_window_days::Float64
    tox_tau::Float64
    sample_dt::Float64
    abstol::Float64
    reltol::Float64
    maxiters::Int
end

mutable struct TraceLogger
    method::String
    t0::Float64
    eval::Int
    rows::Vector{NamedTuple}
end

function log_eval!(logger::TraceLogger, x::AbstractVector, loss::Real)
    logger.eval += 1
    payload = (
        method = logger.method,
        eval = logger.eval,
        elapsed_s = time() - logger.t0,
        loss = Float64(loss),
    )
    d = Dict{Symbol, Any}(pairs(payload))
    for i in eachindex(x)
        d[Symbol("dose$(i)_mg")] = Float64(x[i])
    end
    push!(logger.rows, (; d...))
end

function make_logged_objective(base_obj::Function, logger::TraceLogger)
    return function (x::AbstractVector)
        val = base_obj(x)
        if eltype(x) <: Real && !(x[1] isa ForwardDiff.Dual)
            log_eval!(logger, x, val)
        end
        return val
    end
end

function build_ad_rhs_no_cast(mdl)
    name_to_idx = mdl.name_to_idx
    repeated_rule_defs = mdl.repeated_rule_defs
    rate_exprs = mdl.rate_exprs
    stoich = mdl.stoich

    repeated_assign_exprs = Any[]
    for (lhs_idx, rhs0) in repeated_rule_defs
        rhs =
            if occursin("PK_v26(", rhs0)
                "ifelse(TDBc_ugperkg / Vc_tdb > 1e-5, TDBc_ugperkg / Vc_tdb, 0.0)"
            else
                rhs0
            end
        rhs_expr = TCellEngagerQSP.compile_formula_expr(rhs, name_to_idx)
        push!(repeated_assign_exprs, :(z[$lhs_idx] = real($rhs_expr)))
    end

    reaction_blocks = Any[]
    for j in eachindex(rate_exprs)
        rvar = gensym(:rate)
        rexpr = TCellEngagerQSP.compile_formula_expr(rate_exprs[j], name_to_idx)
        block = Expr(:block, :($rvar = real($rexpr)))
        for (i, coeff) in stoich[j]
            push!(block.args, :(du[$i] += $(Float64(coeff)) * $rvar))
        end
        push!(reaction_blocks, block)
    end

    np = length(mdl.pvals)
    zlen = length(mdl.name_to_idx)
    fexpr = quote
        (du, u, pvals, t) -> begin
            ns = length(u)
            T = eltype(u)
            z = Vector{T}(undef, $zlen)
            @inbounds begin
                z[1:ns] .= u
                z[ns+1:ns+$np] .= pvals
                $(repeated_assign_exprs...)

                fill!(du, zero(T))
                $(reaction_blocks...)
            end
            return nothing
        end
    end

    rhs_rgf = RuntimeGeneratedFunction(
        TCellEngagerQSP,
        TCellEngagerQSP,
        TCellEngagerQSP.normalize_function_expr(fexpr),
    )
    pvals0 = copy(mdl.pvals)
    return (du, u, p, t) -> rhs_rgf(du, u, pvals0, t)
end

function build_symbol_observer_no_cast(mdl, symbol_name::String)
    name_to_idx = mdl.name_to_idx
    repeated_rule_defs = mdl.repeated_rule_defs
    zlen = length(name_to_idx)
    target_idx = name_to_idx[symbol_name]
    np = length(mdl.pvals)
    repeated_assign_exprs = Any[]
    for (lhs_idx, rhs0) in repeated_rule_defs
        rhs =
            if occursin("PK_v26(", rhs0)
                "ifelse(TDBc_ugperkg / Vc_tdb > 1e-5, TDBc_ugperkg / Vc_tdb, 0.0)"
            else
                rhs0
            end
        rhs_expr = TCellEngagerQSP.compile_formula_expr(rhs, name_to_idx)
        push!(repeated_assign_exprs, :(z[$lhs_idx] = real($rhs_expr)))
    end

    fexpr = quote
        (u, t, pvals) -> begin
            ns = length(u)
            T = eltype(u)
            z = Vector{T}(undef, $zlen)
            @inbounds begin
                z[1:ns] .= u
                z[ns+1:ns+$np] .= pvals
                $(repeated_assign_exprs...)
            end
            return z[$target_idx]
        end
    end
    obs_rgf = RuntimeGeneratedFunction(
        TCellEngagerQSP,
        TCellEngagerQSP,
        TCellEngagerQSP.normalize_function_expr(fexpr),
    )
    pvals0 = copy(mdl.pvals)
    return (u, t) -> obs_rgf(u, t, pvals0)
end

function normalize_weights(w::LossWeights)
    s = w.tox_peak_mean + w.tox_peak_max + w.tox_auc + w.tumor_terminal + w.tumor_auc
    s <= 0.0 && error("Loss weights must sum to > 0.")
    return LossWeights(
        w.tox_peak_mean / s,
        w.tox_peak_max / s,
        w.tox_auc / s,
        w.tumor_terminal / s,
        w.tumor_auc / s,
    )
end

function build_dose_schedule(
    mode::String,
    n_cycles::Int;
    step_d1_mg::Float64,
    step_d8_mg::Float64,
    target_mg::Float64,
)
    n_cycles >= 2 || error("n_cycles must be >=2 for clinically meaningful multi-cycle simulation.")
    mode_l = lowercase(mode)
    times = Float64[]
    doses = Float64[]

    if mode_l == "stepup_then_q3w"
        append!(times, [0.0, 7.0, 14.0])
        append!(doses, [step_d1_mg, step_d8_mg, target_mg])
        for c in 2:n_cycles
            push!(times, 21.0 * (c - 1))
            push!(doses, target_mg)
        end
    elseif mode_l == "three_per_cycle"
        for c in 1:n_cycles
            t0 = 21.0 * (c - 1)
            append!(times, [t0, t0 + 7.0, t0 + 14.0])
            append!(doses, [step_d1_mg, step_d8_mg, target_mg])
        end
    else
        error("Unsupported OPTMC_SCHEDULE_MODE=$mode. Use stepup_then_q3w or three_per_cycle.")
    end
    return times, doses
end

function build_patient_problem()
    patients_csv_env = get(ENV, "OPTMC_PATIENTS_CSV", joinpath("assets", "generated_vpop", "selected_patients.csv"))
    patients_csv = isabspath(patients_csv_env) ? patients_csv_env : joinpath(TCellEngagerQSP.REPO_ROOT, patients_csv_env)
    overrides_csv_env = get(ENV, "OPTMC_OVERRIDES_CSV", joinpath("generated", "phase1_design", "dlbcl_param_overrides.csv"))
    overrides_csv = isabspath(overrides_csv_env) ? overrides_csv_env : joinpath(TCellEngagerQSP.REPO_ROOT, overrides_csv_env)

    patients = DataFrame(CSV.File(patients_csv))
    nrow(patients) > 0 || error("No rows found in patients table: $patients_csv")
    patient_id_default = string(Int(patients.patient_id[1]))
    patient_id = parse(Int, get(ENV, "OPTMC_PATIENT_ID", patient_id_default))
    n_cycles = parse(Int, get(ENV, "OPTMC_N_CYCLES", "8"))
    schedule_mode = get(ENV, "OPTMC_SCHEDULE_MODE", "stepup_then_q3w")
    step_d1_mg = parse(Float64, get(ENV, "OPTMC_STEP_D1_MG", "1.0"))
    step_d8_mg = parse(Float64, get(ENV, "OPTMC_STEP_D8_MG", "2.0"))
    target_mg = parse(Float64, get(ENV, "OPTMC_TARGET_MG", "6.0"))

    dose_times_days, fixed_doses_mg = build_dose_schedule(
        schedule_mode,
        n_cycles;
        step_d1_mg = step_d1_mg,
        step_d8_mg = step_d8_mg,
        target_mg = target_mg,
    )
    fixed_override = strip(get(ENV, "OPTMC_FIXED_DOSES_MG", ""))
    if !isempty(fixed_override)
        fixed_doses_mg = parse_float_list(fixed_override)
        length(fixed_doses_mg) == length(dose_times_days) || error(
            "OPTMC_FIXED_DOSES_MG length $(length(fixed_doses_mg)) must match number of dose events $(length(dose_times_days)).",
        )
    end

    lower_scalar = parse(Float64, get(ENV, "OPTMC_DOSE_LOWER_MG", "0.0"))
    upper_scalar = parse(Float64, get(ENV, "OPTMC_DOSE_UPPER_MG", "30.0"))
    lower_mg = fill(lower_scalar, length(dose_times_days))
    upper_mg = fill(upper_scalar, length(dose_times_days))

    horizon_days = parse(Float64, get(ENV, "OPTMC_HORIZON_DAYS", string(21.0 * n_cycles)))
    tox_window_days = parse(Float64, get(ENV, "OPTMC_TOX_WINDOW_DAYS", "2.0"))
    tox_tau = parse(Float64, get(ENV, "OPTMC_TOX_SOFTMAX_TAU", "50.0"))
    sample_dt = parse(Float64, get(ENV, "OPTMC_SAMPLE_DT", "0.25"))
    abstol = parse(Float64, get(ENV, "OPTMC_ABSTOL", "1e-8"))
    reltol = parse(Float64, get(ENV, "OPTMC_RELTOL", "1e-6"))
    maxiters = parse(Int, get(ENV, "OPTMC_MAXITERS_SOLVE", "1000000"))

    overrides = DataFrame(CSV.File(overrides_csv))

    prow_idx = findfirst(patients.patient_id .== patient_id)
    prow_idx === nothing && error("patient_id=$patient_id not found in $patients_csv")
    prow = patients[prow_idx, :]

    pmap = Dict{String, Float64}()
    for rr in eachrow(overrides)
        pmap[String(rr.name)] = Float64(rr.value)
    end

    # Apply patient-specific overrides from the selected VPop row if the parameter exists in the model table.
    for cname in names(patients)
        cname == "patient_id" && continue
        v = try
            Float64(prow[Symbol(cname)])
        catch
            NaN
        end
        isfinite(v) || continue
        if haskey(pmap, cname)
            pmap[cname] = v
        elseif endswith(cname, "_init")
            base = replace(cname, "_init" => "")
            if haskey(pmap, base)
                pmap[base] = v
            end
        end
    end

    # Ensure core per-patient fields are set if present in the row.
    for nm in ("Bpbo_perml", "Bpbref_perml", "Trpbo_perml", "Trpbref_perml", "KBptumor", "KTrptumor", "kBtumorprolif")
        if Symbol(nm) in names(patients)
            pmap[nm] = Float64(prow[Symbol(nm)])
        end
    end
    pmap["PKflag"] = 1.0
    pmap["fvalidation"] = 0.0
    pmap["VPid"] = 1.0
    pmap["end_time"] = horizon_days
    pnames = sort(collect(keys(pmap)))
    pvals = [pmap[n] for n in pnames]

    mdl = TCellEngagerQSP.build_model_with_variant_ids(Int[], pnames, pvals)
    rhs_fun = build_ad_rhs_no_cast(mdl)
    il6_obs_fun = build_symbol_observer_no_cast(mdl, "IL6combo")

    bt_idx = mdl.state_to_idx["Btumor"]
    tdbc_idx = mdl.state_to_idx["TDBc_ugperkg"]
    e_tdbc = zeros(Float64, length(mdl.u0))
    e_tdbc[tdbc_idx] = 1.0
    alg = TCellEngagerQSP.make_solver_alg(lowercase(get(ENV, "OPTMC_SOLVER", "tsit5")))

    prob = PatientMultiCycleProblem(
        patient_id,
        mdl,
        rhs_fun,
        il6_obs_fun,
        alg,
        bt_idx,
        e_tdbc,
        1000.0 / 70.0,
        dose_times_days,
        fixed_doses_mg,
        lower_mg,
        upper_mg,
        horizon_days,
        tox_window_days,
        tox_tau,
        sample_dt,
        abstol,
        reltol,
        maxiters,
    )
    return prob, schedule_mode, n_cycles
end

function build_saveat(t0::Float64, t1::Float64, dt::Float64)
    if t1 <= t0 + 1e-12
        return [t1]
    end
    vals = collect(t0:dt:t1)
    if isempty(vals) || vals[end] < t1 - 1e-12
        push!(vals, t1)
    else
        vals[end] = t1
    end
    return vals
end

function solve_segment(
    prob::PatientMultiCycleProblem,
    u0,
    tspan::Tuple{Float64,Float64},
    saveat::Vector{Float64};
    sensealg = nothing,
)
    ode_prob = ODEProblem(prob.rhs_fun, u0, tspan, [1.0, 1.0, 1.0])
    sol =
        if isnothing(sensealg)
            solve(
                ode_prob,
                prob.alg;
                abstol = prob.abstol,
                reltol = prob.reltol,
                saveat = saveat,
                tstops = [tspan[2]],
                maxiters = prob.maxiters,
            )
        else
            solve(
                ode_prob,
                prob.alg;
                sensealg = sensealg,
                abstol = prob.abstol,
                reltol = prob.reltol,
                saveat = saveat,
                tstops = [tspan[2]],
                maxiters = prob.maxiters,
            )
        end
    if sol.retcode != SciMLBase.ReturnCode.Success
        error("segment solve failed retcode=$(sol.retcode) tspan=$tspan")
    end
    return sol
end

function simulate_trajectory(
    prob::PatientMultiCycleProblem,
    doses_mg::AbstractVector{T};
    sensealg = nothing,
) where {T<:Real}
    nd = length(prob.dose_times_days)
    length(doses_mg) == nd || error("Expected $nd doses, got $(length(doses_mg))")

    d_ugkg = doses_mg .* T(prob.mg_to_ugkg)
    e_tdbc_t = T.(prob.e_tdbc)

    u_curr = T.(prob.mdl.u0)
    bt0 = u_curr[prob.bt_idx]
    t_curr = 0.0

    t_all = Float64[0.0]
    bt_all = T[bt0]
    il6_all = T[prob.il6_obs_fun(u_curr, 0.0)]

    for i in 1:nd
        tdose = prob.dose_times_days[i]
        if tdose > t_curr + 1e-12
            saveat = build_saveat(t_curr, tdose, prob.sample_dt)
            sol = solve_segment(prob, u_curr, (t_curr, tdose), saveat; sensealg = sensealg)
            for j in eachindex(sol.t)
                tj = Float64(sol.t[j])
                if !isempty(t_all) && isapprox(tj, t_all[end]; atol = 1e-12, rtol = 0.0)
                    continue
                end
                uj = sol.u[j]
                push!(t_all, tj)
                push!(bt_all, uj[prob.bt_idx])
                push!(il6_all, prob.il6_obs_fun(uj, sol.t[j]))
            end
            u_curr = sol.u[end]
            t_curr = tdose
        end

        u_curr = u_curr .+ d_ugkg[i] .* e_tdbc_t
        il6_post = prob.il6_obs_fun(u_curr, t_curr)
        if isempty(t_all) || !isapprox(t_curr, t_all[end]; atol = 1e-12, rtol = 0.0)
            push!(t_all, t_curr)
            push!(bt_all, u_curr[prob.bt_idx])
            push!(il6_all, il6_post)
        else
            bt_all[end] = u_curr[prob.bt_idx]
            il6_all[end] = il6_post
        end
    end

    if t_curr < prob.horizon_days - 1e-12
        saveat = build_saveat(t_curr, prob.horizon_days, prob.sample_dt)
        sol = solve_segment(prob, u_curr, (t_curr, prob.horizon_days), saveat; sensealg = sensealg)
        for j in eachindex(sol.t)
            tj = Float64(sol.t[j])
            if !isempty(t_all) && isapprox(tj, t_all[end]; atol = 1e-12, rtol = 0.0)
                continue
            end
            uj = sol.u[j]
            push!(t_all, tj)
            push!(bt_all, uj[prob.bt_idx])
            push!(il6_all, prob.il6_obs_fun(uj, sol.t[j]))
        end
    end

    return (; t = t_all, btumor = bt_all, il6 = il6_all, bt0 = bt0)
end

function collect_window_vals(t::Vector{Float64}, y::AbstractVector{T}, t0::Float64, t1::Float64) where {T}
    vals = T[]
    for i in eachindex(t)
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

function raw_metrics(prob::PatientMultiCycleProblem, traj)
    T = eltype(traj.btumor)
    tox_peaks = T[]
    for tdose in prob.dose_times_days
        vals = collect_window_vals(traj.t, traj.il6, tdose, tdose + prob.tox_window_days)
        push!(tox_peaks, smoothmax(vals, prob.tox_tau))
    end
    tox_peak_mean = sum(tox_peaks) / T(length(tox_peaks))
    tox_peak_max = maximum(tox_peaks)
    tox_auc = trapz(traj.t, max.(traj.il6, zero(T))) / T(prob.horizon_days)
    tumor_terminal = traj.btumor[end] / (traj.bt0 + T(1e-12))
    tumor_auc = trapz(traj.t, max.(traj.btumor, zero(T))) / (T(prob.horizon_days) * (traj.bt0 + T(1e-12)))
    return (; tox_peak_mean, tox_peak_max, tox_auc, tumor_terminal, tumor_auc)
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

function calibrate_scales(
    prob::PatientMultiCycleProblem;
    n_samples::Int,
    rng_seed::Int,
)
    nd = length(prob.fixed_doses_mg)
    rng = MersenneTwister(rng_seed)
    m_tox_mean = Float64[]
    m_tox_max = Float64[]
    m_tox_auc = Float64[]
    m_tum_term = Float64[]
    m_tum_auc = Float64[]

    candidates = Vector{Vector{Float64}}()
    push!(candidates, copy(prob.fixed_doses_mg))
    for _ in 1:n_samples
        x = similar(prob.fixed_doses_mg)
        for j in 1:nd
            x[j] = prob.lower_mg[j] + rand(rng) * (prob.upper_mg[j] - prob.lower_mg[j])
        end
        push!(candidates, x)
    end

    for x in candidates
        traj = simulate_trajectory(prob, x)
        raw = raw_metrics(prob, traj)
        push!(m_tox_mean, Float64(raw.tox_peak_mean))
        push!(m_tox_max, Float64(raw.tox_peak_max))
        push!(m_tox_auc, Float64(raw.tox_auc))
        push!(m_tum_term, Float64(raw.tumor_terminal))
        push!(m_tum_auc, Float64(raw.tumor_auc))
    end

    return LossScales(
        robust_scale(m_tox_mean, 1e-6, m_tox_mean[1]),
        robust_scale(m_tox_max, 1e-6, m_tox_max[1]),
        robust_scale(m_tox_auc, 1e-6, m_tox_auc[1]),
        robust_scale(m_tum_term, 1e-10, m_tum_term[1]),
        robust_scale(m_tum_auc, 1e-10, m_tum_auc[1]),
    )
end

function make_loss_objective(
    prob::PatientMultiCycleProblem,
    scales::LossScales,
    weights::LossWeights;
    sensealg = nothing,
)
    big_penalty = 1.0e9
    return function (x::AbstractVector)
        if !(x[1] isa ForwardDiff.Dual)
            xv = Float64.(x)
            if any(!isfinite, xv)
                return big_penalty
            end
            viol = max.(prob.lower_mg .- xv, 0.0) .+ max.(xv .- prob.upper_mg, 0.0)
            if any(viol .> 0.0)
                return big_penalty + 1.0e6 * sum(abs2, viol)
            end
        end
        try
            traj = simulate_trajectory(prob, x; sensealg = sensealg)
            raw = raw_metrics(prob, traj)
            return scaled_loss(raw, scales, weights).loss
        catch
            if x[1] isa ForwardDiff.Dual
                return one(x[1]) * big_penalty
            end
            return big_penalty
        end
    end
end

function do_coordinate_descent(
    obj::Function,
    x0::Vector{Float64},
    lower::Vector{Float64},
    upper::Vector{Float64};
    sweeps::Int,
    step0::Float64,
    tol::Float64,
)
    x = copy(x0)
    fx = Float64(obj(x))
    step = step0
    for _ in 1:sweeps
        improved = false
        for j in eachindex(x)
            x_plus = copy(x)
            x_plus[j] = min(upper[j], x[j] + step)
            f_plus = Float64(obj(x_plus))

            x_minus = copy(x)
            x_minus[j] = max(lower[j], x[j] - step)
            f_minus = Float64(obj(x_minus))

            if f_plus < fx && f_plus <= f_minus
                x = x_plus
                fx = f_plus
                improved = true
            elseif f_minus < fx
                x = x_minus
                fx = f_minus
                improved = true
            end
        end
        step *= 0.5
        if step < tol && !improved
            break
        end
    end
    return x, fx
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
)
    rng = MersenneTwister(rng_seed)
    nd = length(x0)
    x_best = clamp.(copy(x0), lower, upper)
    f_best = Float64(obj(x_best))
    n_evals <= 1 && return x_best, f_best

    for _ in 2:n_evals
        x = similar(x_best)
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

        f = Float64(obj(x))
        if f < f_best
            x_best = copy(x)
            f_best = f
        end
    end

    return x_best, f_best
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

    n_steps <= 0 && return x_best, f_best, 0.0
    temp0p = max(temp0, 1e-12)
    tempfp = max(tempf, 1e-12)
    denom = max(n_steps - 1, 1)

    for k in 1:n_steps
        α = (k - 1) / denom
        T = exp(log(temp0p) * (1.0 - α) + log(tempfp) * α)

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

function run_all_methods(
    prob::PatientMultiCycleProblem,
    scales::LossScales,
    weights::LossWeights;
    loss_dto = nothing,
    loss_otd = nothing,
)
    lower = copy(prob.lower_mg)
    upper = copy(prob.upper_mg)
    x0 = clamp.(copy(prob.fixed_doses_mg), lower, upper)
    nd = length(x0)

    if isnothing(loss_dto)
        loss_dto = make_loss_objective(prob, scales, weights; sensealg = nothing)
    end
    if isnothing(loss_otd)
        loss_otd = make_loss_objective(prob, scales, weights; sensealg = ForwardSensitivity())
    end

    sa_iters = parse(Int, get(ENV, "OPTMC_SA_ITERS", "120"))
    nm_iters = parse(Int, get(ENV, "OPTMC_NM_ITERS", "80"))
    cd_iters = parse(Int, get(ENV, "OPTMC_CD_ITERS", "60"))
    cd_step0 = parse(Float64, get(ENV, "OPTMC_CD_STEP0_MG", "2.0"))
    cd_tol = parse(Float64, get(ENV, "OPTMC_CD_TOL_MG", "0.05"))
    lbfgs_iters = parse(Int, get(ENV, "OPTMC_LBFGS_ITERS", "80"))
    sa_f_calls = parse(Int, get(ENV, "OPTMC_SA_F_CALLS", "0"))
    nm_f_calls = parse(Int, get(ENV, "OPTMC_NM_F_CALLS", "0"))
    lbfgs_f_calls = parse(Int, get(ENV, "OPTMC_LBFGS_F_CALLS", "0"))

    method_rows = NamedTuple[]
    all_rows = NamedTuple[]

    function push_result!(method::String, best_loss::Float64, x::Vector{Float64}, runtime_s::Float64, status::String)
        d = Dict{Symbol, Any}(
            :method => method,
            :best_loss => best_loss,
            :runtime_s => runtime_s,
            :status => status,
        )
        for i in 1:nd
            d[Symbol("dose$(i)_mg")] = Float64(x[i])
        end
        push!(method_rows, (; d...))
    end

    # 1) Simulated annealing (gradient free)
    Random.seed!(parse(Int, get(ENV, "OPTMC_RANDOM_SEED", "20260303")))
    logger_sa = TraceLogger("simulated_annealing", time(), 0, NamedTuple[])
    f_sa = make_logged_objective(loss_dto, logger_sa)
    status_sa = "ok"
    x_sa = copy(x0)
    loss_sa = Float64(loss_dto(x_sa))
    try
        sa_opts = Optim.Options(
            iterations = sa_iters,
            f_calls_limit = sa_f_calls > 0 ? sa_f_calls : typemax(Int),
            show_trace = false,
            store_trace = false,
        )
        res_sa = optimize(f_sa, lower, upper, copy(x0), SAMIN(), sa_opts)
        x_sa = clamp.(Optim.minimizer(res_sa), lower, upper)
        loss_sa = Float64(loss_dto(x_sa))
    catch err
        status_sa = "error: " * sprint(showerror, err)
    end
    runtime_sa = time() - logger_sa.t0
    append!(all_rows, logger_sa.rows)
    push_result!("simulated_annealing", loss_sa, x_sa, runtime_sa, status_sa)

    # 2) Nelder-Mead + coordinate descent refinement
    logger_nm = TraceLogger("nelder_mead_coordinate_descent", time(), 0, NamedTuple[])
    f_nm = make_logged_objective(loss_dto, logger_nm)
    status_nm = "ok"
    x_nm = copy(x0)
    loss_nm = Float64(loss_dto(x_nm))
    try
        nm_opts = Optim.Options(
            iterations = nm_iters,
            f_calls_limit = nm_f_calls > 0 ? nm_f_calls : typemax(Int),
            show_trace = false,
            store_trace = false,
        )
        res_nm = optimize(f_nm, lower, upper, copy(x0), Fminbox(NelderMead()), nm_opts)
        x_curr = clamp.(Optim.minimizer(res_nm), lower, upper)
        loss_curr = Float64(f_nm(x_curr))
        steps = fill(cd_step0, nd)
        for _ in 1:cd_iters
            improved = false
            for j in 1:nd
                x_plus = copy(x_curr)
                x_plus[j] = clamp(x_plus[j] + steps[j], lower[j], upper[j])
                l_plus = Float64(f_nm(x_plus))

                x_minus = copy(x_curr)
                x_minus[j] = clamp(x_minus[j] - steps[j], lower[j], upper[j])
                l_minus = Float64(f_nm(x_minus))

                if l_plus < loss_curr && l_plus <= l_minus
                    x_curr = x_plus
                    loss_curr = l_plus
                    improved = true
                elseif l_minus < loss_curr
                    x_curr = x_minus
                    loss_curr = l_minus
                    improved = true
                end
            end
            if !improved
                steps .*= 0.5
                if maximum(steps) < cd_tol
                    break
                end
            end
        end
        x_nm = x_curr
        loss_nm = loss_curr
    catch err
        status_nm = "error: " * sprint(showerror, err)
    end
    runtime_nm = time() - logger_nm.t0
    append!(all_rows, logger_nm.rows)
    push_result!("nelder_mead_coordinate_descent", loss_nm, x_nm, runtime_nm, status_nm)

    # 3) Finite-difference + LBFGS
    logger_fd = TraceLogger("finite_diff_lbfgs", time(), 0, NamedTuple[])
    f_fd = make_logged_objective(loss_dto, logger_fd)
    g_fd! = (G, x) -> begin
        FiniteDiff.finite_difference_gradient!(G, loss_dto, x)
        return G
    end
    status_fd = "ok"
    x_fd = copy(x0)
    loss_fd = Float64(loss_dto(x_fd))
    try
        lbfgs_opts = Optim.Options(
            iterations = lbfgs_iters,
            f_calls_limit = lbfgs_f_calls > 0 ? lbfgs_f_calls : typemax(Int),
            show_trace = false,
            store_trace = false,
        )
        res_fd = optimize(f_fd, g_fd!, lower, upper, copy(x0), Fminbox(LBFGS()), lbfgs_opts)
        x_fd = clamp.(Optim.minimizer(res_fd), lower, upper)
        loss_fd = Float64(loss_dto(x_fd))
    catch err
        status_fd = "error: " * sprint(showerror, err)
    end
    runtime_fd = time() - logger_fd.t0
    append!(all_rows, logger_fd.rows)
    push_result!("finite_diff_lbfgs", loss_fd, x_fd, runtime_fd, status_fd)

    # 4) Forward AD DTO + LBFGS
    logger_ad_dto = TraceLogger("forward_ad_dto_lbfgs", time(), 0, NamedTuple[])
    f_ad_dto = make_logged_objective(loss_dto, logger_ad_dto)
    g_ad_dto! = (G, x) -> begin
        ForwardDiff.gradient!(G, loss_dto, x)
        return G
    end
    status_ad_dto = "ok"
    x_ad_dto = copy(x0)
    loss_ad_dto = Float64(loss_dto(x_ad_dto))
    try
        lbfgs_opts = Optim.Options(
            iterations = lbfgs_iters,
            f_calls_limit = lbfgs_f_calls > 0 ? lbfgs_f_calls : typemax(Int),
            show_trace = false,
            store_trace = false,
        )
        res_ad_dto = optimize(f_ad_dto, g_ad_dto!, lower, upper, copy(x0), Fminbox(LBFGS()), lbfgs_opts)
        x_ad_dto = clamp.(Optim.minimizer(res_ad_dto), lower, upper)
        loss_ad_dto = Float64(loss_dto(x_ad_dto))
    catch err
        status_ad_dto = "error: " * sprint(showerror, err)
    end
    runtime_ad_dto = time() - logger_ad_dto.t0
    append!(all_rows, logger_ad_dto.rows)
    push_result!("forward_ad_dto_lbfgs", loss_ad_dto, x_ad_dto, runtime_ad_dto, status_ad_dto)

    # 5) Forward AD OTD + LBFGS
    logger_ad_otd = TraceLogger("forward_ad_otd_lbfgs", time(), 0, NamedTuple[])
    f_ad_otd = make_logged_objective(loss_otd, logger_ad_otd)
    g_ad_otd! = (G, x) -> begin
        ForwardDiff.gradient!(G, loss_otd, x)
        return G
    end
    status_ad_otd = "ok"
    x_ad_otd = copy(x0)
    loss_ad_otd = Float64(loss_otd(x_ad_otd))
    try
        lbfgs_opts = Optim.Options(
            iterations = lbfgs_iters,
            f_calls_limit = lbfgs_f_calls > 0 ? lbfgs_f_calls : typemax(Int),
            show_trace = false,
            store_trace = false,
        )
        res_ad_otd = optimize(f_ad_otd, g_ad_otd!, lower, upper, copy(x0), Fminbox(LBFGS()), lbfgs_opts)
        x_ad_otd = clamp.(Optim.minimizer(res_ad_otd), lower, upper)
        loss_ad_otd = Float64(loss_otd(x_ad_otd))
    catch err
        status_ad_otd = "error: " * sprint(showerror, err)
    end
    runtime_ad_otd = time() - logger_ad_otd.t0
    append!(all_rows, logger_ad_otd.rows)
    push_result!("forward_ad_otd_lbfgs", loss_ad_otd, x_ad_otd, runtime_ad_otd, status_ad_otd)

    return DataFrame(method_rows), DataFrame(all_rows)
end

function to_float_dict(nt)
    out = Dict{String, Float64}()
    for (k, v) in pairs(nt)
        out[String(k)] = Float64(v)
    end
    return out
end

function main()
    prob, schedule_mode, n_cycles = build_patient_problem()
    nd = length(prob.dose_times_days)
    println(@sprintf("patient=%d, schedule=%s, cycles=%d, decision_dim=%d", prob.patient_id, schedule_mode, n_cycles, nd))
    println("fixed doses mg = $(prob.fixed_doses_mg)")

    weights = normalize_weights(
        LossWeights(
            parse(Float64, get(ENV, "OPTMC_W_TOX_PEAK_MEAN", "0.20")),
            parse(Float64, get(ENV, "OPTMC_W_TOX_PEAK_MAX", "0.25")),
            parse(Float64, get(ENV, "OPTMC_W_TOX_AUC", "0.15")),
            parse(Float64, get(ENV, "OPTMC_W_TUMOR_TERMINAL", "0.20")),
            parse(Float64, get(ENV, "OPTMC_W_TUMOR_AUC", "0.20")),
        ),
    )

    n_scale_samples = parse(Int, get(ENV, "OPTMC_SCALE_SAMPLES", "40"))
    scale_seed = parse(Int, get(ENV, "OPTMC_SCALE_SEED", "20260303"))
    scales = calibrate_scales(prob; n_samples = n_scale_samples, rng_seed = scale_seed)
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

    obj_base = make_loss_objective(prob, scales, weights; sensealg = nothing)
    obj_otd = make_loss_objective(prob, scales, weights; sensealg = ForwardSensitivity())
    method = lowercase(get(ENV, "OPTMC_OPT_METHOD", "forward_ad_lbfgs"))
    logger = TraceLogger(method, time(), 0, NamedTuple[])
    obj = make_logged_objective(obj_base, logger)

    x0 = clamp.(copy(prob.fixed_doses_mg), prob.lower_mg, prob.upper_mg)
    out_dir_default = joinpath(TCellEngagerQSP.REPO_ROOT, "generated", "figures", "optimization", "multicycle_patient_$(prob.patient_id)_$(n_cycles)cycles")
    out_dir_env = get(ENV, "OPTMC_OUT_DIR", out_dir_default)
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(TCellEngagerQSP.REPO_ROOT, out_dir_env)
    mkpath(out_dir)

    if method == "all_methods"
        results, traces = run_all_methods(prob, scales, weights)
        ok = results[occursin.("ok", results.status), :]
        nrow(ok) > 0 || error("All methods failed.")
        bidx = argmin(ok.best_loss)
        best = ok[bidx, :]
        x_best = [Float64(best[Symbol("dose$(i)_mg")]) for i in 1:nd]

        traj_fixed = simulate_trajectory(prob, x0)
        traj_best = simulate_trajectory(prob, x_best)
        raw_fixed = raw_metrics(prob, traj_fixed)
        raw_best = raw_metrics(prob, traj_best)
        scaled_fixed = scaled_loss(raw_fixed, scales, weights)
        scaled_best = scaled_loss(raw_best, scales, weights)

        summary_path = joinpath(out_dir, "optimization_summary.csv")
        trace_path = joinpath(out_dir, "optimization_traces.csv")
        traj_path = joinpath(out_dir, "best_dose_trajectory.csv")
        dose_path = joinpath(out_dir, "dose_schedule.csv")
        meta_path = joinpath(out_dir, "optimization_meta.json")

        CSV.write(summary_path, results)
        CSV.write(trace_path, traces)
        CSV.write(traj_path, DataFrame(time_day = traj_best.t, Btumor = Float64.(traj_best.btumor), IL6combo = Float64.(traj_best.il6)))
        CSV.write(
            dose_path,
            DataFrame(
                dose_idx = collect(1:nd),
                time_day = prob.dose_times_days,
                fixed_dose_mg = x0,
                best_dose_mg = x_best,
                lower_mg = prob.lower_mg,
                upper_mg = prob.upper_mg,
            ),
        )

        meta = Dict(
            "patient_id" => prob.patient_id,
            "schedule_mode" => schedule_mode,
            "n_cycles" => n_cycles,
            "horizon_days" => prob.horizon_days,
            "dose_times_days" => prob.dose_times_days,
            "fixed_doses_mg" => x0,
            "best_doses_mg" => x_best,
            "best_method" => String(best.method),
            "best_loss" => Float64(best.best_loss),
            "solver" => string(prob.alg),
            "loss_formula" => "L=w1*(mean_dose_window_IL6_peak/s1)+w2*(max_dose_window_IL6_peak/s2)+w3*(IL6_AUC/s3)+w4*(Btumor_end/Btumor0/s4)+w5*(Btumor_AUC_ratio/s5)",
            "loss_weights" => Dict(
                "tox_peak_mean" => weights.tox_peak_mean,
                "tox_peak_max" => weights.tox_peak_max,
                "tox_auc" => weights.tox_auc,
                "tumor_terminal" => weights.tumor_terminal,
                "tumor_auc" => weights.tumor_auc,
            ),
            "loss_scales" => Dict(
                "tox_peak_mean" => scales.tox_peak_mean,
                "tox_peak_max" => scales.tox_peak_max,
                "tox_auc" => scales.tox_auc,
                "tumor_terminal" => scales.tumor_terminal,
                "tumor_auc" => scales.tumor_auc,
            ),
            "fixed_raw_terms" => to_float_dict(raw_fixed),
            "fixed_scaled_terms" => to_float_dict(scaled_fixed),
            "best_raw_terms" => to_float_dict(raw_best),
            "best_scaled_terms" => to_float_dict(scaled_best),
            "objective_fixed" => Float64(scaled_fixed.loss),
            "objective_best" => Float64(scaled_best.loss),
            "scale_calibration_samples" => n_scale_samples,
            "scale_calibration_seed" => scale_seed,
        )
        open(meta_path, "w") do io
            JSON3.pretty(io, meta)
        end

        println("Best method: ", best.method)
        println(@sprintf("Best loss: %.6g", best.best_loss))
        println("Wrote $summary_path")
        println("Wrote $trace_path")
        println("Wrote $traj_path")
        println("Wrote $dose_path")
        println("Wrote $meta_path")
        return
    end

    x_opt = copy(x0)
    status = "ok"
    t_opt0 = time()
    try
        if method == "forward_ad_lbfgs" || method == "forward_ad_dto_lbfgs"
            g! = function (G, x)
                ForwardDiff.gradient!(G, obj_base, x)
            end
            opt_opts = Optim.Options(
                iterations = parse(Int, get(ENV, "OPTMC_LBFGS_ITERS", "160")),
                show_trace = false,
                store_trace = false,
                g_tol = parse(Float64, get(ENV, "OPTMC_LBFGS_GTOL", "1e-8")),
                f_tol = parse(Float64, get(ENV, "OPTMC_LBFGS_FTOL", "1e-10")),
            )
            res = optimize(obj, g!, prob.lower_mg, prob.upper_mg, copy(x0), Fminbox(LBFGS()), opt_opts)
            x_opt = clamp.(Optim.minimizer(res), prob.lower_mg, prob.upper_mg)
            status = string(Optim.converged(res) ? "ok" : "not_converged")
        elseif method == "forward_ad_otd_lbfgs"
            obj2 = make_logged_objective(obj_otd, logger)
            g! = function (G, x)
                ForwardDiff.gradient!(G, obj_otd, x)
            end
            opt_opts = Optim.Options(
                iterations = parse(Int, get(ENV, "OPTMC_LBFGS_ITERS", "160")),
                show_trace = false,
                store_trace = false,
                g_tol = parse(Float64, get(ENV, "OPTMC_LBFGS_GTOL", "1e-8")),
                f_tol = parse(Float64, get(ENV, "OPTMC_LBFGS_FTOL", "1e-10")),
            )
            res = optimize(obj2, g!, prob.lower_mg, prob.upper_mg, copy(x0), Fminbox(LBFGS()), opt_opts)
            x_opt = clamp.(Optim.minimizer(res), prob.lower_mg, prob.upper_mg)
            status = string(Optim.converged(res) ? "ok" : "not_converged")
        elseif method == "finite_diff_lbfgs"
            g! = function (G, x)
                FiniteDiff.finite_difference_gradient!(G, obj_base, x)
            end
            opt_opts = Optim.Options(
                iterations = parse(Int, get(ENV, "OPTMC_LBFGS_ITERS", "160")),
                show_trace = false,
                store_trace = false,
                g_tol = parse(Float64, get(ENV, "OPTMC_LBFGS_GTOL", "1e-8")),
                f_tol = parse(Float64, get(ENV, "OPTMC_LBFGS_FTOL", "1e-10")),
            )
            res = optimize(obj, g!, prob.lower_mg, prob.upper_mg, copy(x0), Fminbox(LBFGS()), opt_opts)
            x_opt = clamp.(Optim.minimizer(res), prob.lower_mg, prob.upper_mg)
            status = string(Optim.converged(res) ? "ok" : "not_converged")
        elseif method == "simulated_annealing"
            sa_iters = parse(Int, get(ENV, "OPTMC_SA_ITERS", "120"))
            sa_f_calls = parse(Int, get(ENV, "OPTMC_SA_F_CALLS", "0"))
            opt_opts = Optim.Options(
                iterations = sa_iters,
                f_calls_limit = sa_f_calls > 0 ? sa_f_calls : typemax(Int),
                show_trace = false,
                store_trace = false,
            )
            res = optimize(obj, prob.lower_mg, prob.upper_mg, copy(x0), SAMIN(), opt_opts)
            x_opt = clamp.(Optim.minimizer(res), prob.lower_mg, prob.upper_mg)
            status = string(Optim.converged(res) ? "ok" : "not_converged")
        elseif method == "random_search"
            rs_evals = parse(Int, get(ENV, "OPTMC_RS_EVALS", "2000"))
            rs_local_prob = parse(Float64, get(ENV, "OPTMC_RS_LOCAL_PROB", "0.70"))
            rs_sigma_frac = parse(Float64, get(ENV, "OPTMC_RS_SIGMA_FRAC", "0.12"))
            rs_seed = parse(Int, get(ENV, "OPTMC_RANDOM_SEED", "20260303"))
            x_opt, _ = do_random_search(
                obj,
                x0,
                prob.lower_mg,
                prob.upper_mg;
                n_evals = rs_evals,
                local_prob = rs_local_prob,
                local_sigma_frac = rs_sigma_frac,
                rng_seed = rs_seed,
            )
            status = "ok"
        elseif method == "metropolis_anneal" || method == "mcmc_anneal"
            mc_steps = parse(Int, get(ENV, "OPTMC_MC_STEPS", "2500"))
            mc_temp0 = parse(Float64, get(ENV, "OPTMC_MC_TEMP0", "0.03"))
            mc_tempf = parse(Float64, get(ENV, "OPTMC_MC_TEMPF", "2e-4"))
            mc_sigma_frac = parse(Float64, get(ENV, "OPTMC_MC_SIGMA_FRAC", "0.08"))
            mc_seed = parse(Int, get(ENV, "OPTMC_RANDOM_SEED", "20260303"))
            x_opt, _, acc_rate = do_metropolis_anneal(
                obj,
                x0,
                prob.lower_mg,
                prob.upper_mg;
                n_steps = mc_steps,
                temp0 = mc_temp0,
                tempf = mc_tempf,
                proposal_sigma_frac = mc_sigma_frac,
                rng_seed = mc_seed,
            )
            status = @sprintf("ok_accept=%.3f", acc_rate)
        else
            x_opt, _ = do_coordinate_descent(
                obj,
                x0,
                prob.lower_mg,
                prob.upper_mg;
                sweeps = parse(Int, get(ENV, "OPTMC_CD_SWEEPS", "6")),
                step0 = parse(Float64, get(ENV, "OPTMC_CD_STEP0_MG", "2.0")),
                tol = parse(Float64, get(ENV, "OPTMC_CD_TOL_MG", "0.05")),
            )
            status = "ok"
        end
    catch err
        status = "fail: $(typeof(err))"
    end
    runtime_s = time() - t_opt0

    traj_fixed = simulate_trajectory(prob, x0)
    raw_fixed = raw_metrics(prob, traj_fixed)
    scaled_fixed = scaled_loss(raw_fixed, scales, weights)

    traj_opt = simulate_trajectory(prob, x_opt)
    raw_opt = raw_metrics(prob, traj_opt)
    scaled_opt = scaled_loss(raw_opt, scales, weights)

    dose_rows = [
        (
            dose_idx = i,
            time_day = prob.dose_times_days[i],
            fixed_dose_mg = x0[i],
            optimized_dose_mg = x_opt[i],
            lower_mg = prob.lower_mg[i],
            upper_mg = prob.upper_mg[i],
        ) for i in 1:nd
    ]
    dose_path = joinpath(out_dir, "dose_schedule.csv")
    CSV.write(dose_path, DataFrame(dose_rows))

    trace_path = joinpath(out_dir, "loss_trace.csv")
    if !isempty(logger.rows)
        CSV.write(trace_path, DataFrame(logger.rows))
    else
        CSV.write(trace_path, DataFrame(method = String[], eval = Int[], elapsed_s = Float64[], loss = Float64[]))
    end

    traj_fixed_path = joinpath(out_dir, "trajectory_fixed.csv")
    traj_opt_path = joinpath(out_dir, "trajectory_optimized.csv")
    CSV.write(
        traj_fixed_path,
        DataFrame(t_day = traj_fixed.t, btumor = Float64.(traj_fixed.btumor), il6combo = Float64.(traj_fixed.il6)),
    )
    CSV.write(
        traj_opt_path,
        DataFrame(t_day = traj_opt.t, btumor = Float64.(traj_opt.btumor), il6combo = Float64.(traj_opt.il6)),
    )

    summary = Dict(
        "patient_id" => prob.patient_id,
        "schedule_mode" => schedule_mode,
        "n_cycles" => n_cycles,
        "horizon_days" => prob.horizon_days,
        "dose_times_days" => prob.dose_times_days,
        "fixed_doses_mg" => x0,
        "optimized_doses_mg" => x_opt,
        "status" => status,
        "runtime_s" => runtime_s,
        "solver" => string(prob.alg),
        "optimizer_method" => method,
        "loss_weights" => Dict(
            "tox_peak_mean" => weights.tox_peak_mean,
            "tox_peak_max" => weights.tox_peak_max,
            "tox_auc" => weights.tox_auc,
            "tumor_terminal" => weights.tumor_terminal,
            "tumor_auc" => weights.tumor_auc,
        ),
        "loss_scales" => Dict(
            "tox_peak_mean" => scales.tox_peak_mean,
            "tox_peak_max" => scales.tox_peak_max,
            "tox_auc" => scales.tox_auc,
            "tumor_terminal" => scales.tumor_terminal,
            "tumor_auc" => scales.tumor_auc,
        ),
        "loss_formula" => "L=w1*(mean_dose_window_IL6_peak/s1)+w2*(max_dose_window_IL6_peak/s2)+w3*(IL6_AUC/s3)+w4*(Btumor_end/Btumor0/s4)+w5*(Btumor_AUC_ratio/s5)",
        "fixed_raw_terms" => to_float_dict(raw_fixed),
        "fixed_scaled_terms" => to_float_dict(scaled_fixed),
        "optimized_raw_terms" => to_float_dict(raw_opt),
        "optimized_scaled_terms" => to_float_dict(scaled_opt),
        "objective_fixed" => Float64(scaled_fixed.loss),
        "objective_optimized" => Float64(scaled_opt.loss),
        "scale_calibration_samples" => n_scale_samples,
        "scale_calibration_seed" => scale_seed,
    )
    summary_path = joinpath(out_dir, "optimization_summary.json")
    open(summary_path, "w") do io
        JSON3.pretty(io, summary)
    end

    println("Wrote $dose_path")
    println("Wrote $trace_path")
    println("Wrote $traj_fixed_path")
    println("Wrote $traj_opt_path")
    println("Wrote $summary_path")
    println(@sprintf("objective fixed=%.6g optimized=%.6g", summary["objective_fixed"], summary["objective_optimized"]))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
