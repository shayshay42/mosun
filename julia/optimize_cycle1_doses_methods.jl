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

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

RuntimeGeneratedFunctions.init(@__MODULE__)

function smoothmax(v::AbstractVector, tau::Real)
    m = maximum(v)
    return m + tau * log(sum(exp.((v .- m) ./ tau)))
end

struct CycleOptProblem
    mdl
    rhs_fun::Function
    il6_obs_fun::Function
    alg
    lower_mg::Vector{Float64}
    upper_mg::Vector{Float64}
    x0_mg::Vector{Float64}
    mg_to_ugkg::Float64
    e_tdbc::Vector{Float64}
    bt_idx::Int
    horizon_days::Float64
    tox_window_days::Float64
    tox_grid::Vector{Float64}
    tox_tau::Float64
    tox_scale::Float64
    tumor_scale::Float64
    w_tox::Float64
    w_tumor::Float64
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
    push!(
        logger.rows,
        (
            method = logger.method,
            eval = logger.eval,
            elapsed_s = time() - logger.t0,
            loss = Float64(loss),
            dose1_mg = Float64(x[1]),
            dose2_mg = Float64(x[2]),
            dose3_mg = Float64(x[3]),
        ),
    )
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
                # Phase-1 setup fixes PKflag=1, so this simplification is exact on the optimization path.
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

function build_problem(; patient_id::Int = 1)
    repo_root = TCellEngagerQSP.REPO_ROOT
    design_dir = joinpath(repo_root, "generated", "phase1_design")
    patients = DataFrame(CSV.File(joinpath(design_dir, "patients.csv")))
    overrides = DataFrame(CSV.File(joinpath(design_dir, "dlbcl_param_overrides.csv")))

    prow_idx = findfirst(patients.patient_id .== patient_id)
    prow_idx === nothing && error("patient_id=$patient_id not found in $(joinpath(design_dir, "patients.csv"))")
    prow = patients[prow_idx, :]

    pmap = Dict{String, Float64}()
    for rr in eachrow(overrides)
        pmap[String(rr.name)] = Float64(rr.value)
    end
    pmap["Bpbo_perml"] = Float64(prow.Bpbo_perml)
    pmap["Bpbref_perml"] = Float64(prow.Bpbref_perml)
    pmap["Trpbo_perml"] = Float64(prow.Trpbo_perml)
    pmap["Trpbref_perml"] = Float64(prow.Trpbref_perml)
    pmap["KBptumor"] = Float64(prow.KBptumor)
    pmap["KTrptumor"] = Float64(prow.KTrptumor)
    pmap["kBtumorprolif"] = Float64(prow.kBtumorprolif)
    pmap["PKflag"] = 1.0
    pmap["fvalidation"] = 0.0
    pmap["VPid"] = 1.0
    pmap["end_time"] = 21.0
    pnames = sort(collect(keys(pmap)))
    pvals = [pmap[n] for n in pnames]

    mdl = TCellEngagerQSP.build_model_with_variant_ids(Int[], pnames, pvals)
    rhs_fun = build_ad_rhs_no_cast(mdl)
    il6_obs_fun = build_symbol_observer_no_cast(mdl, "IL6combo")

    bt_idx = mdl.state_to_idx["Btumor"]
    tdbc_idx = mdl.state_to_idx["TDBc_ugperkg"]

    x0_mg = [
        parse(Float64, get(ENV, "OPT_DOSE0_MG", "0.8")),
        parse(Float64, get(ENV, "OPT_DOSE7_MG", "2.0")),
        parse(Float64, get(ENV, "OPT_DOSE14_MG", "6.0")),
    ]
    lower_mg = fill(parse(Float64, get(ENV, "OPT_DOSE_LOWER_MG", "0.0")), 3)
    upper_mg = fill(parse(Float64, get(ENV, "OPT_DOSE_UPPER_MG", "30.0")), 3)

    mg_to_ugkg = 1000.0 / 70.0
    e_tdbc = zeros(Float64, length(mdl.u0))
    e_tdbc[tdbc_idx] = 1.0
    alg = TCellEngagerQSP.make_solver_alg(lowercase(get(ENV, "OPT_SOLVER", "tsit5")))

    tox_window_days = 2.0
    tox_grid = collect(0.0:0.1:tox_window_days)
    horizon_days = 21.0
    tox_tau = parse(Float64, get(ENV, "OPT_TOX_SOFTMAX_TAU", "50.0"))
    w_tox = parse(Float64, get(ENV, "OPT_LOSS_W_TOX", "0.5"))
    w_tumor = parse(Float64, get(ENV, "OPT_LOSS_W_TUMOR", "0.5"))
    abstol = parse(Float64, get(ENV, "OPT_ABSTOL", "1e-8"))
    reltol = parse(Float64, get(ENV, "OPT_RELTOL", "1e-6"))
    maxiters = parse(Int, get(ENV, "OPT_MAXITERS_SOLVE", "1000000"))

    # Reference scales from the gradient-check setup.
    tox_scale = parse(Float64, get(ENV, "OPT_TOX_SCALE", "152.22612188617114"))
    tumor_scale = parse(Float64, get(ENV, "OPT_TUMOR_SCALE", "0.7050236467345111"))

    return CycleOptProblem(
        mdl,
        rhs_fun,
        il6_obs_fun,
        alg,
        lower_mg,
        upper_mg,
        x0_mg,
        mg_to_ugkg,
        e_tdbc,
        bt_idx,
        horizon_days,
        tox_window_days,
        tox_grid,
        tox_tau,
        tox_scale,
        tumor_scale,
        w_tox,
        w_tumor,
        abstol,
        reltol,
        maxiters,
    )
end

function solve_segment(prob::CycleOptProblem, u0, tspan::Tuple{Float64,Float64}, saveat::Vector{Float64}; sensealg = nothing)
    ode_prob = ODEProblem(prob.rhs_fun, u0, tspan, [1.0, 1.0, 1.0])
    if isnothing(sensealg)
        sol = solve(
            ode_prob,
            prob.alg;
            abstol = prob.abstol,
            reltol = prob.reltol,
            saveat = saveat,
            tstops = [tspan[2]],
            maxiters = prob.maxiters,
        )
    else
        sol = solve(
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
        error("segment solve failed with retcode=$(sol.retcode) tspan=$tspan")
    end
    return sol
end

function simulate_cycle(prob::CycleOptProblem, doses_mg::AbstractVector{T}; save_trajectory::Bool = false, sensealg = nothing) where {T<:Real}
    length(doses_mg) == 3 || error("Expected 3 doses")
    d_ugkg = doses_mg .* T(prob.mg_to_ugkg)
    e_tdbc_t = T.(prob.e_tdbc)

    u1 = T.(prob.mdl.u0) .+ d_ugkg[1] .* e_tdbc_t
    bt0 = u1[prob.bt_idx]

    seg1_save = save_trajectory ? collect(0.0:0.1:7.0) : vcat(prob.tox_grid, [7.0])
    sol1 = solve_segment(prob, u1, (0.0, 7.0), seg1_save; sensealg = sensealg)
    il6_vals = [
        prob.il6_obs_fun(sol1.u[i], sol1.t[i]) for i in eachindex(sol1.t)
        if Float64(sol1.t[i]) <= prob.tox_window_days + 1e-12
    ]
    tox_proxy = smoothmax(il6_vals, prob.tox_tau)

    u2 = sol1.u[end] .+ d_ugkg[2] .* e_tdbc_t
    seg2_save = save_trajectory ? collect(7.0:0.1:14.0) : [14.0]
    sol2 = solve_segment(prob, u2, (7.0, 14.0), seg2_save; sensealg = sensealg)

    u3 = sol2.u[end] .+ d_ugkg[3] .* e_tdbc_t
    seg3_save = save_trajectory ? collect(14.0:0.1:21.0) : [21.0]
    sol3 = solve_segment(prob, u3, (14.0, 21.0), seg3_save; sensealg = sensealg)

    bt_end = sol3.u[end][prob.bt_idx]
    tumor_proxy = bt_end / (bt0 + T(1e-12))
    loss = T(prob.w_tox) * (tox_proxy / T(max(prob.tox_scale, 1e-12))) +
           T(prob.w_tumor) * (tumor_proxy / T(max(prob.tumor_scale, 1e-12)))

    if !save_trajectory
        return (; tox_proxy = tox_proxy, tumor_proxy = tumor_proxy, loss = loss)
    end

    t_all = Float64[]
    bt_all = Float64[]
    il6_all = Float64[]
    for sol in (sol1, sol2, sol3)
        for i in eachindex(sol.t)
            t = Float64(sol.t[i])
            if !isempty(t_all) && isapprox(t, t_all[end]; atol = 1e-12, rtol = 0.0)
                continue
            end
            u = sol.u[i]
            push!(t_all, t)
            push!(bt_all, Float64(u[prob.bt_idx]))
            push!(il6_all, Float64(prob.il6_obs_fun(u, sol.t[i])))
        end
    end
    return (; tox_proxy = tox_proxy, tumor_proxy = tumor_proxy, loss = loss, t = t_all, btumor = bt_all, il6combo = il6_all)
end

function make_loss_objective(prob::CycleOptProblem; sensealg = nothing)
    big_penalty = 1.0e9
    return function (x::AbstractVector)
        T = eltype(x)
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
            return simulate_cycle(prob, x; save_trajectory = false, sensealg = sensealg).loss
        catch
            if x[1] isa ForwardDiff.Dual
                return one(x[1]) * big_penalty
            end
            return big_penalty
        end
    end
end

function run_all_methods(prob::CycleOptProblem)
    loss_dto = make_loss_objective(prob; sensealg = nothing)
    loss_otd = make_loss_objective(prob; sensealg = ForwardSensitivity())
    lower = copy(prob.lower_mg)
    upper = copy(prob.upper_mg)
    x0 = clamp.(copy(prob.x0_mg), lower, upper)

    sa_iters = parse(Int, get(ENV, "OPT_SA_ITERS", "220"))
    nm_iters = parse(Int, get(ENV, "OPT_NM_ITERS", "180"))
    cd_iters = parse(Int, get(ENV, "OPT_CD_ITERS", "80"))
    cd_step0 = parse(Float64, get(ENV, "OPT_CD_STEP0_MG", "2.0"))
    cd_tol = parse(Float64, get(ENV, "OPT_CD_TOL_MG", "0.02"))
    lbfgs_iters = parse(Int, get(ENV, "OPT_LBFGS_ITERS", "120"))

    results = DataFrame(
        method = String[],
        best_loss = Float64[],
        dose1_mg = Float64[],
        dose2_mg = Float64[],
        dose3_mg = Float64[],
        runtime_s = Float64[],
        status = String[],
    )
    all_rows = NamedTuple[]

    # 1) Gradient-free: Simulated Annealing.
    Random.seed!(parse(Int, get(ENV, "OPT_RANDOM_SEED", "20260302")))
    logger_sa = TraceLogger("simulated_annealing", time(), 0, NamedTuple[])
    f_sa = make_logged_objective(loss_dto, logger_sa)
    sa_opt = Optim.Options(iterations = sa_iters, show_trace = false, store_trace = false)
    status_sa = "ok"
    x_sa = copy(x0)
    loss_sa = loss_dto(x_sa)
    try
        res_sa = optimize(f_sa, lower, upper, copy(x0), SAMIN(), sa_opt)
        x_sa = clamp.(Optim.minimizer(res_sa), lower, upper)
        loss_sa = loss_dto(x_sa)
    catch err
        status_sa = "error: " * sprint(showerror, err)
    end
    runtime_sa = time() - logger_sa.t0
    append!(all_rows, logger_sa.rows)
    push!(results, ("simulated_annealing", loss_sa, x_sa[1], x_sa[2], x_sa[3], runtime_sa, status_sa))

    # 2) Gradient-free: Nelder-Mead + coordinate descent refinement.
    logger_nm = TraceLogger("nelder_mead_coordinate_descent", time(), 0, NamedTuple[])
    f_nm = make_logged_objective(loss_dto, logger_nm)
    nm_opt = Optim.Options(iterations = nm_iters, show_trace = false, store_trace = false)
    status_nm = "ok"
    x_nm = copy(x0)
    loss_nm = loss_dto(x_nm)
    try
        res_nm = optimize(f_nm, lower, upper, copy(x0), Fminbox(NelderMead()), nm_opt)
        x_curr = clamp.(Optim.minimizer(res_nm), lower, upper)
        loss_curr = f_nm(x_curr)
        steps = fill(cd_step0, 3)
        for _ in 1:cd_iters
            improved = false
            for j in 1:3
                x_plus = copy(x_curr)
                x_plus[j] = clamp(x_plus[j] + steps[j], lower[j], upper[j])
                l_plus = f_nm(x_plus)

                x_minus = copy(x_curr)
                x_minus[j] = clamp(x_minus[j] - steps[j], lower[j], upper[j])
                l_minus = f_nm(x_minus)

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
    push!(results, ("nelder_mead_coordinate_descent", loss_nm, x_nm[1], x_nm[2], x_nm[3], runtime_nm, status_nm))

    # 3) Finite-difference + L-BFGS.
    logger_fd = TraceLogger("finite_diff_lbfgs", time(), 0, NamedTuple[])
    f_fd = make_logged_objective(loss_dto, logger_fd)
    g_fd! = (G, x) -> begin
        copyto!(G, FiniteDiff.finite_difference_gradient(loss_dto, x))
        return G
    end
    lbfgs_opt = Optim.Options(iterations = lbfgs_iters, show_trace = false, store_trace = false)
    status_fd = "ok"
    x_fd = copy(x0)
    loss_fd = loss_dto(x_fd)
    try
        res_fd = optimize(f_fd, g_fd!, lower, upper, copy(x0), Fminbox(LBFGS()), lbfgs_opt)
        x_fd = clamp.(Optim.minimizer(res_fd), lower, upper)
        loss_fd = loss_dto(x_fd)
    catch err
        status_fd = "error: " * sprint(showerror, err)
    end
    runtime_fd = time() - logger_fd.t0
    append!(all_rows, logger_fd.rows)
    push!(results, ("finite_diff_lbfgs", loss_fd, x_fd[1], x_fd[2], x_fd[3], runtime_fd, status_fd))

    # 4) Forward-mode AD DTO + L-BFGS.
    logger_ad_dto = TraceLogger("forward_ad_dto_lbfgs", time(), 0, NamedTuple[])
    f_ad_dto = make_logged_objective(loss_dto, logger_ad_dto)
    g_ad_dto! = (G, x) -> begin
        ForwardDiff.gradient!(G, loss_dto, x)
        return G
    end
    status_ad_dto = "ok"
    x_ad_dto = copy(x0)
    loss_ad_dto = loss_dto(x_ad_dto)
    try
        res_ad_dto = optimize(f_ad_dto, g_ad_dto!, lower, upper, copy(x0), Fminbox(LBFGS()), lbfgs_opt)
        x_ad_dto = clamp.(Optim.minimizer(res_ad_dto), lower, upper)
        loss_ad_dto = loss_dto(x_ad_dto)
    catch err
        status_ad_dto = "error: " * sprint(showerror, err)
    end
    runtime_ad_dto = time() - logger_ad_dto.t0
    append!(all_rows, logger_ad_dto.rows)
    push!(results, ("forward_ad_dto_lbfgs", loss_ad_dto, x_ad_dto[1], x_ad_dto[2], x_ad_dto[3], runtime_ad_dto, status_ad_dto))

    # 5) Forward-mode AD OTD + L-BFGS.
    logger_ad_otd = TraceLogger("forward_ad_otd_lbfgs", time(), 0, NamedTuple[])
    f_ad_otd = make_logged_objective(loss_otd, logger_ad_otd)
    g_ad_otd! = (G, x) -> begin
        ForwardDiff.gradient!(G, loss_otd, x)
        return G
    end
    status_ad_otd = "ok"
    x_ad_otd = copy(x0)
    loss_ad_otd = loss_otd(x_ad_otd)
    try
        res_ad_otd = optimize(f_ad_otd, g_ad_otd!, lower, upper, copy(x0), Fminbox(LBFGS()), lbfgs_opt)
        x_ad_otd = clamp.(Optim.minimizer(res_ad_otd), lower, upper)
        loss_ad_otd = loss_otd(x_ad_otd)
    catch err
        status_ad_otd = "error: " * sprint(showerror, err)
    end
    runtime_ad_otd = time() - logger_ad_otd.t0
    append!(all_rows, logger_ad_otd.rows)
    push!(results, ("forward_ad_otd_lbfgs", loss_ad_otd, x_ad_otd[1], x_ad_otd[2], x_ad_otd[3], runtime_ad_otd, status_ad_otd))

    traces = DataFrame(all_rows)
    return results, traces
end

function main()
    patient_id = parse(Int, get(ENV, "OPT_PATIENT_ID", "1"))
    prob = build_problem(patient_id = patient_id)
    results, traces = run_all_methods(prob)

    ok_rows = results[occursin.("ok", results.status), :]
    if nrow(ok_rows) == 0
        error("All optimization methods failed.")
    end
    best_idx = argmin(ok_rows.best_loss)
    best_row = ok_rows[best_idx, :]
    best_doses = [best_row.dose1_mg, best_row.dose2_mg, best_row.dose3_mg]
    sim_best = simulate_cycle(prob, best_doses; save_trajectory = true)

    out_dir_default = joinpath(TCellEngagerQSP.REPO_ROOT, "generated", "figures", "optimization", "patient_$(patient_id)")
    out_dir_env = get(ENV, "OPT_OUT_DIR", out_dir_default)
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(TCellEngagerQSP.REPO_ROOT, out_dir_env)
    mkpath(out_dir)

    trace_path = joinpath(out_dir, "optimization_traces.csv")
    summary_path = joinpath(out_dir, "optimization_summary.csv")
    traj_path = joinpath(out_dir, "best_dose_trajectory.csv")
    meta_path = joinpath(out_dir, "optimization_meta.json")

    CSV.write(trace_path, traces)
    CSV.write(summary_path, results)
    traj_df = DataFrame(
        time_day = sim_best.t,
        Btumor = sim_best.btumor,
        IL6combo = sim_best.il6combo,
    )
    CSV.write(traj_path, traj_df)

    loss_formula = @sprintf(
        "L=%.2f*(smoothmax(IL6combo[0-2],tau=%.1f)/%.6f) + %.2f*((Btumor21/Btumor0)/%.6f)",
        prob.w_tox,
        prob.tox_tau,
        prob.tox_scale,
        prob.w_tumor,
        prob.tumor_scale,
    )

    meta = Dict(
        "patient_id" => patient_id,
        "dose_times_days" => [0.0, 7.0, 14.0],
        "horizon_days" => prob.horizon_days,
        "loss_formula" => loss_formula,
        "tox_proxy" => "smoothmax(IL6combo day 0-2)",
        "tumor_proxy" => "Btumor(day21)/Btumor(day0)",
        "best_method" => String(best_row.method),
        "best_loss" => Float64(best_row.best_loss),
        "best_doses_mg" => best_doses,
        "solver" => string(prob.alg),
    )
    open(meta_path, "w") do io
        JSON3.pretty(io, meta)
    end

    println("Best method: ", best_row.method)
    println(@sprintf("Best doses (mg): [%.4f, %.4f, %.4f]", best_doses[1], best_doses[2], best_doses[3]))
    println(@sprintf("Best loss: %.6f", best_row.best_loss))
    println("Wrote ", trace_path)
    println("Wrote ", summary_path)
    println("Wrote ", traj_path)
    println("Wrote ", meta_path)
end

main()
