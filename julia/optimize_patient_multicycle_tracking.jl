using CSV
using DataFrames
using JSON3
using Printf
using SciMLSensitivity

include(joinpath(@__DIR__, "optimize_patient_multicycle_balanced.jl"))

struct TrackingSpec
    t_ref::Vector{Float64}
    bt_target::Vector{Float64}
    il6_target::Vector{Float64}
    il6_scale::Float64
    w_tumor::Float64
    w_tox::Float64
    w_dose::Float64
end

function norm3(a::Float64, b::Float64, c::Float64)
    s = a + b + c
    s <= 0.0 && error("Tracking weights must sum to > 0.")
    return a / s, b / s, c / s
end

function build_ref_grid(horizon::Float64, dt::Float64)
    vals = collect(0.0:dt:horizon)
    if isempty(vals) || vals[end] < horizon - 1e-12
        push!(vals, horizon)
    else
        vals[end] = horizon
    end
    return vals
end

function build_tracking_spec(prob::PatientMultiCycleProblem)
    tumor_end_frac = parse(Float64, get(ENV, "OPTMC_TRACK_TUMOR_END_FRAC", "0.05"))
    il6_base = parse(Float64, get(ENV, "OPTMC_TRACK_IL6_BASE", "0.0"))
    il6_pulse = parse(Float64, get(ENV, "OPTMC_TRACK_IL6_PULSE", "120.0"))
    il6_sigma = parse(Float64, get(ENV, "OPTMC_TRACK_IL6_SIGMA_DAYS", "0.40"))
    il6_scale = parse(Float64, get(ENV, "OPTMC_TRACK_IL6_SCALE", "300.0"))

    w_tumor_raw = parse(Float64, get(ENV, "OPTMC_TRACK_W_TUMOR", "0.70"))
    w_tox_raw = parse(Float64, get(ENV, "OPTMC_TRACK_W_TOX", "0.25"))
    w_dose_raw = parse(Float64, get(ENV, "OPTMC_TRACK_W_DOSE", "0.05"))
    w_tumor, w_tox, w_dose = norm3(w_tumor_raw, w_tox_raw, w_dose_raw)

    t_ref = build_ref_grid(prob.horizon_days, prob.sample_dt)
    bt0 = Float64(prob.mdl.u0[prob.bt_idx])
    k = log(max(tumor_end_frac, 1e-6)) / prob.horizon_days
    bt_target = [bt0 * exp(k * t) for t in t_ref]

    il6_target = fill(il6_base, length(t_ref))
    σ = max(il6_sigma, 1e-6)
    for tdose in prob.dose_times_days
        center = tdose + 0.5
        il6_target .+= il6_pulse .* exp.(-0.5 .* ((t_ref .- center) ./ σ) .^ 2)
    end

    spec = TrackingSpec(t_ref, bt_target, il6_target, il6_scale, w_tumor, w_tox, w_dose)
    params = Dict(
        "tumor_end_frac" => tumor_end_frac,
        "il6_base" => il6_base,
        "il6_pulse" => il6_pulse,
        "il6_sigma_days" => il6_sigma,
        "il6_scale" => il6_scale,
        "w_tumor" => w_tumor,
        "w_tox" => w_tox,
        "w_dose" => w_dose,
    )
    return spec, params
end

function interp_linear(tsrc::Vector{Float64}, ysrc::AbstractVector{T}, tq::Vector{Float64}) where {T}
    n = length(tsrc)
    n == length(ysrc) || error("interp_linear source length mismatch")
    out = Vector{T}(undef, length(tq))
    j = 1
    for i in eachindex(tq)
        t = tq[i]
        if t <= tsrc[1]
            out[i] = ysrc[1]
            continue
        end
        if t >= tsrc[end]
            out[i] = ysrc[end]
            continue
        end
        while j < n - 1 && tsrc[j + 1] < t
            j += 1
        end
        t0 = tsrc[j]
        t1 = tsrc[j + 1]
        y0 = ysrc[j]
        y1 = ysrc[j + 1]
        α = (t - t0) / (t1 - t0 + 1e-12)
        out[i] = y0 + (y1 - y0) * α
    end
    return out
end

function tracking_terms(prob::PatientMultiCycleProblem, traj, doses::AbstractVector, spec::TrackingSpec)
    T = eltype(traj.btumor)
    bt = interp_linear(traj.t, traj.btumor, spec.t_ref)
    il6 = interp_linear(traj.t, traj.il6, spec.t_ref)

    bt0 = traj.bt0 + T(1e-12)
    bt_rel_log = log.(bt ./ bt0 .+ T(1e-12))
    bt_target_rel_log = log.(T.(spec.bt_target) ./ bt0 .+ T(1e-12))
    tumor_mse = trapz(spec.t_ref, (bt_rel_log .- bt_target_rel_log) .^ 2) / T(prob.horizon_days)

    il6_err = (il6 .- T.(spec.il6_target)) ./ T(spec.il6_scale)
    tox_mse = trapz(spec.t_ref, il6_err .^ 2) / T(prob.horizon_days)

    dose_reg = sum((doses ./ prob.upper_mg) .^ 2) / T(length(doses))
    loss = spec.w_tumor * tumor_mse + spec.w_tox * tox_mse + spec.w_dose * dose_reg

    return (; tumor_mse, tox_mse, dose_reg, loss)
end

function make_tracking_objective(prob::PatientMultiCycleProblem, spec::TrackingSpec; sensealg = nothing)
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
            return tracking_terms(prob, traj, x, spec).loss
        catch
            if x[1] isa ForwardDiff.Dual
                return one(x[1]) * big_penalty
            end
            return big_penalty
        end
    end
end

function main()
    prob, schedule_mode, n_cycles = build_patient_problem()
    nd = length(prob.dose_times_days)
    println(@sprintf("tracking objective run: patient=%d, schedule=%s, cycles=%d, decision_dim=%d", prob.patient_id, schedule_mode, n_cycles, nd))

    spec, spec_params = build_tracking_spec(prob)

    # Keep placeholders for run_all_methods schema compatibility.
    scales_stub = LossScales(1.0, 1.0, 1.0, 1.0, 1.0)
    weights_stub = LossWeights(0.2, 0.2, 0.2, 0.2, 0.2)

    loss_dto = make_tracking_objective(prob, spec; sensealg = nothing)
    loss_otd = make_tracking_objective(prob, spec; sensealg = ForwardSensitivity())

    out_dir_default = joinpath(TCellEngagerQSP.REPO_ROOT, "generated", "figures", "optimization", "multicycle_tracking_patient_$(prob.patient_id)_$(n_cycles)cycles")
    out_dir_env = get(ENV, "OPTMC_OUT_DIR", out_dir_default)
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(TCellEngagerQSP.REPO_ROOT, out_dir_env)
    mkpath(out_dir)

    x_fixed = clamp.(copy(prob.fixed_doses_mg), prob.lower_mg, prob.upper_mg)
    method = lowercase(get(ENV, "OPTMC_OPT_METHOD", "all_methods"))

    results = DataFrame()
    traces = DataFrame()
    best_method = method
    best_loss = Inf
    x_best = copy(x_fixed)
    best_status = "ok"

    if method == "all_methods"
        results, traces = run_all_methods(prob, scales_stub, weights_stub; loss_dto = loss_dto, loss_otd = loss_otd)
        ok = results[occursin.("ok", results.status), :]
        nrow(ok) > 0 || error("All methods failed.")
        bidx = argmin(ok.best_loss)
        best = ok[bidx, :]
        x_best = [Float64(best[Symbol("dose$(i)_mg")]) for i in 1:nd]
        best_method = String(best.method)
        best_loss = Float64(best.best_loss)
        best_status = String(best.status)
    else
        logger = TraceLogger(method, time(), 0, NamedTuple[])
        obj = make_logged_objective(loss_dto, logger)
        x_opt = copy(x_fixed)
        status = "ok"
        t_opt0 = time()
        try
            if method == "forward_ad_lbfgs" || method == "forward_ad_dto_lbfgs"
                g! = (G, x) -> ForwardDiff.gradient!(G, loss_dto, x)
                opt_opts = Optim.Options(
                    iterations = parse(Int, get(ENV, "OPTMC_LBFGS_ITERS", "80")),
                    show_trace = false,
                    store_trace = false,
                )
                res = optimize(obj, g!, prob.lower_mg, prob.upper_mg, copy(x_fixed), Fminbox(LBFGS()), opt_opts)
                x_opt = clamp.(Optim.minimizer(res), prob.lower_mg, prob.upper_mg)
                status = string(Optim.converged(res) ? "ok" : "not_converged")
            elseif method == "forward_ad_otd_lbfgs"
                obj2 = make_logged_objective(loss_otd, logger)
                g! = (G, x) -> ForwardDiff.gradient!(G, loss_otd, x)
                opt_opts = Optim.Options(
                    iterations = parse(Int, get(ENV, "OPTMC_LBFGS_ITERS", "80")),
                    show_trace = false,
                    store_trace = false,
                )
                res = optimize(obj2, g!, prob.lower_mg, prob.upper_mg, copy(x_fixed), Fminbox(LBFGS()), opt_opts)
                x_opt = clamp.(Optim.minimizer(res), prob.lower_mg, prob.upper_mg)
                status = string(Optim.converged(res) ? "ok" : "not_converged")
            elseif method == "finite_diff_lbfgs"
                g! = (G, x) -> FiniteDiff.finite_difference_gradient!(G, loss_dto, x)
                opt_opts = Optim.Options(
                    iterations = parse(Int, get(ENV, "OPTMC_LBFGS_ITERS", "80")),
                    show_trace = false,
                    store_trace = false,
                )
                res = optimize(obj, g!, prob.lower_mg, prob.upper_mg, copy(x_fixed), Fminbox(LBFGS()), opt_opts)
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
                res = optimize(obj, prob.lower_mg, prob.upper_mg, copy(x_fixed), SAMIN(), opt_opts)
                x_opt = clamp.(Optim.minimizer(res), prob.lower_mg, prob.upper_mg)
                status = string(Optim.converged(res) ? "ok" : "not_converged")
            elseif method == "random_search"
                x_opt, _ = do_random_search(
                    obj,
                    x_fixed,
                    prob.lower_mg,
                    prob.upper_mg;
                    n_evals = parse(Int, get(ENV, "OPTMC_RS_EVALS", "1000")),
                    local_prob = parse(Float64, get(ENV, "OPTMC_RS_LOCAL_PROB", "0.70")),
                    local_sigma_frac = parse(Float64, get(ENV, "OPTMC_RS_SIGMA_FRAC", "0.12")),
                    rng_seed = parse(Int, get(ENV, "OPTMC_RANDOM_SEED", "20260303")),
                )
                status = "ok"
            elseif method == "metropolis_anneal" || method == "mcmc_anneal"
                x_opt, _, acc_rate = do_metropolis_anneal(
                    obj,
                    x_fixed,
                    prob.lower_mg,
                    prob.upper_mg;
                    n_steps = parse(Int, get(ENV, "OPTMC_MC_STEPS", "1200")),
                    temp0 = parse(Float64, get(ENV, "OPTMC_MC_TEMP0", "0.03")),
                    tempf = parse(Float64, get(ENV, "OPTMC_MC_TEMPF", "2e-4")),
                    proposal_sigma_frac = parse(Float64, get(ENV, "OPTMC_MC_SIGMA_FRAC", "0.08")),
                    rng_seed = parse(Int, get(ENV, "OPTMC_RANDOM_SEED", "20260303")),
                )
                status = @sprintf("ok_accept=%.3f", acc_rate)
            elseif method == "nelder_mead_coordinate_descent"
                nm_iters = parse(Int, get(ENV, "OPTMC_NM_ITERS", "60"))
                cd_iters = parse(Int, get(ENV, "OPTMC_CD_ITERS", "30"))
                cd_step0 = parse(Float64, get(ENV, "OPTMC_CD_STEP0_MG", "2.0"))
                cd_tol = parse(Float64, get(ENV, "OPTMC_CD_TOL_MG", "0.05"))
                nm_f_calls = parse(Int, get(ENV, "OPTMC_NM_F_CALLS", "0"))
                nm_opts = Optim.Options(
                    iterations = nm_iters,
                    f_calls_limit = nm_f_calls > 0 ? nm_f_calls : typemax(Int),
                    show_trace = false,
                    store_trace = false,
                )
                res_nm = optimize(obj, prob.lower_mg, prob.upper_mg, copy(x_fixed), Fminbox(NelderMead()), nm_opts)
                x_curr = clamp.(Optim.minimizer(res_nm), prob.lower_mg, prob.upper_mg)
                loss_curr = Float64(obj(x_curr))
                steps = fill(cd_step0, nd)
                for _ in 1:cd_iters
                    improved = false
                    for j in 1:nd
                        x_plus = copy(x_curr)
                        x_plus[j] = clamp(x_plus[j] + steps[j], prob.lower_mg[j], prob.upper_mg[j])
                        l_plus = Float64(obj(x_plus))
                        x_minus = copy(x_curr)
                        x_minus[j] = clamp(x_minus[j] - steps[j], prob.lower_mg[j], prob.upper_mg[j])
                        l_minus = Float64(obj(x_minus))
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
                x_opt = x_curr
                status = "ok"
            else
                x_opt, _ = do_coordinate_descent(
                    obj,
                    x_fixed,
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
        best_loss = Float64(loss_dto(x_opt))
        x_best = x_opt
        best_method = method
        best_status = status

        row = Dict{Symbol, Any}(
            :method => method,
            :best_loss => best_loss,
            :runtime_s => runtime_s,
            :status => status,
        )
        for i in 1:nd
            row[Symbol("dose$(i)_mg")] = Float64(x_opt[i])
        end
        results = DataFrame([(; row...)])
        traces = isempty(logger.rows) ? DataFrame(method = String[], eval = Int[], elapsed_s = Float64[], loss = Float64[]) : DataFrame(logger.rows)
    end

    traj_fixed = simulate_trajectory(prob, x_fixed)
    traj_best = simulate_trajectory(prob, x_best)
    term_fixed = tracking_terms(prob, traj_fixed, x_fixed, spec)
    term_best = tracking_terms(prob, traj_best, x_best, spec)

    summary_path = joinpath(out_dir, "optimization_summary.csv")
    trace_path = joinpath(out_dir, "optimization_traces.csv")
    dose_path = joinpath(out_dir, "dose_schedule.csv")
    meta_path = joinpath(out_dir, "optimization_meta.json")
    traj_fixed_path = joinpath(out_dir, "trajectory_fixed.csv")
    traj_best_path = joinpath(out_dir, "trajectory_best.csv")
    target_path = joinpath(out_dir, "trajectory_target.csv")

    CSV.write(summary_path, results)
    CSV.write(trace_path, traces)
    CSV.write(
        dose_path,
        DataFrame(
            dose_idx = collect(1:nd),
            time_day = prob.dose_times_days,
            fixed_dose_mg = x_fixed,
            best_dose_mg = x_best,
            lower_mg = prob.lower_mg,
            upper_mg = prob.upper_mg,
        ),
    )
    CSV.write(traj_fixed_path, DataFrame(time_day = traj_fixed.t, Btumor = Float64.(traj_fixed.btumor), IL6combo = Float64.(traj_fixed.il6)))
    CSV.write(traj_best_path, DataFrame(time_day = traj_best.t, Btumor = Float64.(traj_best.btumor), IL6combo = Float64.(traj_best.il6)))
    CSV.write(target_path, DataFrame(time_day = spec.t_ref, Btumor_target = spec.bt_target, IL6_target = spec.il6_target))

    meta = Dict(
        "patient_id" => prob.patient_id,
        "schedule_mode" => schedule_mode,
        "n_cycles" => n_cycles,
        "horizon_days" => prob.horizon_days,
        "dose_times_days" => prob.dose_times_days,
        "fixed_doses_mg" => x_fixed,
        "best_doses_mg" => x_best,
        "best_method" => best_method,
        "best_loss" => best_loss,
        "best_status" => best_status,
        "solver" => string(prob.alg),
        "loss_mode" => "trajectory_tracking_quadratic",
        "loss_formula" => "L=w_tumor*mean((log(B/B0)-log(Btarget/B0))^2)+w_tox*mean(((IL6-IL6target)/il6_scale)^2)+w_dose*mean((dose/ub)^2)",
        "tracking_spec" => spec_params,
        "fixed_terms" => Dict("tumor_mse" => Float64(term_fixed.tumor_mse), "tox_mse" => Float64(term_fixed.tox_mse), "dose_reg" => Float64(term_fixed.dose_reg), "loss" => Float64(term_fixed.loss)),
        "best_terms" => Dict("tumor_mse" => Float64(term_best.tumor_mse), "tox_mse" => Float64(term_best.tox_mse), "dose_reg" => Float64(term_best.dose_reg), "loss" => Float64(term_best.loss)),
    )
    open(meta_path, "w") do io
        JSON3.pretty(io, meta)
    end

    println("Best method: ", best_method)
    println(@sprintf("Best loss: %.6g", best_loss))
    println("Wrote $summary_path")
    println("Wrote $trace_path")
    println("Wrote $traj_fixed_path")
    println("Wrote $traj_best_path")
    println("Wrote $target_path")
    println("Wrote $dose_path")
    println("Wrote $meta_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
