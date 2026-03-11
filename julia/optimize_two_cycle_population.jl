using CSV
using DataFrames
using DifferentialEquations
using ForwardDiff
using JSON3
using Optim
using Printf
using RuntimeGeneratedFunctions
using SciMLBase
using Statistics

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

RuntimeGeneratedFunctions.init(@__MODULE__)

function smoothmax(v::AbstractVector, tau::Real)
    m = maximum(v)
    return m + tau * log(sum(exp.((v .- m) ./ tau)))
end

function robust_scale(x::AbstractVector{Float64}, floor_val::Float64, fallback::Float64)
    vals = filter(v -> isfinite(v) && v > 0.0, x)
    isempty(vals) && return max(floor_val, fallback)
    s = median(vals)
    return max(floor_val, s)
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

function parse_fixed_doses()
    raw = get(ENV, "OPT2C_FIXED_DOSES_MG", "0.8,2.0,6.0,6.0,6.0,6.0")
    x = parse_float_list(raw)
    if length(x) == 3
        return vcat(x, x)
    elseif length(x) == 6
        return x
    end
    error("OPT2C_FIXED_DOSES_MG must contain 3 or 6 comma-separated values.")
end

function parse_int_list(s::AbstractString)
    vals = Int[]
    for tok in split(String(s), ",")
        st = strip(tok)
        isempty(st) && continue
        push!(vals, parse(Int, st))
    end
    return vals
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
                # Optimization uses PKflag=1 so this exact simplification is valid on-path.
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

struct PatientProblem
    patient_id::Int
    mdl
    rhs_fun::Function
    il6_obs_fun::Function
    alg
    e_tdbc::Vector{Float64}
    bt_idx::Int
    dose_times_days::Vector{Float64}
    fixed_doses_mg::Vector{Float64}
    lower_mg::Vector{Float64}
    upper_mg::Vector{Float64}
    mg_to_ugkg::Float64
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

function solve_segment(prob::PatientProblem, u0, tspan::Tuple{Float64,Float64}, saveat::Vector{Float64})
    ode_prob = ODEProblem(prob.rhs_fun, u0, tspan, [1.0, 1.0, 1.0])
    sol = solve(
        ode_prob,
        prob.alg;
        abstol = prob.abstol,
        reltol = prob.reltol,
        saveat = saveat,
        tstops = [tspan[2]],
        maxiters = prob.maxiters,
    )
    if sol.retcode != SciMLBase.ReturnCode.Success
        error("segment solve failed with retcode=$(sol.retcode) tspan=$tspan")
    end
    return sol
end

function simulate_patient(prob::PatientProblem, doses_mg::AbstractVector{T}; save_trajectory::Bool = false) where {T<:Real}
    nd = length(prob.dose_times_days)
    length(doses_mg) == nd || error("Expected $nd doses, got $(length(doses_mg))")

    d_ugkg = doses_mg .* T(prob.mg_to_ugkg)
    e_tdbc_t = T.(prob.e_tdbc)
    dose_times = prob.dose_times_days

    u_curr = T.(prob.mdl.u0)
    t_curr = 0.0
    tox_vals = T[]
    bt0 = zero(T)

    t_all = Float64[]
    bt_all = Float64[]
    il6_all = Float64[]

    for i in 1:nd
        tdose = dose_times[i]
        if tdose > t_curr + 1e-12
            saveat =
                if save_trajectory
                    collect(t_curr:0.1:tdose)
                else
                    tox_pts = [tt for tt in prob.tox_grid if tt > t_curr + 1e-12 && tt <= min(tdose, prob.tox_window_days) + 1e-12]
                    sort(unique(vcat([tdose], tox_pts)))
                end

            sol_pre = solve_segment(prob, u_curr, (t_curr, tdose), saveat)
            for j in eachindex(sol_pre.t)
                tj = Float64(sol_pre.t[j])
                uj = sol_pre.u[j]
                il6v = prob.il6_obs_fun(uj, sol_pre.t[j])
                if tj <= prob.tox_window_days + 1e-12
                    push!(tox_vals, il6v)
                end
                if save_trajectory
                    if !isempty(t_all) && isapprox(tj, t_all[end]; atol = 1e-12, rtol = 0.0)
                        continue
                    end
                    push!(t_all, tj)
                    push!(bt_all, Float64(uj[prob.bt_idx]))
                    push!(il6_all, Float64(il6v))
                end
            end
            u_curr = sol_pre.u[end]
            t_curr = tdose
        end

        u_curr = u_curr .+ d_ugkg[i] .* e_tdbc_t
        if i == 1
            bt0 = u_curr[prob.bt_idx]
        end
        if t_curr <= prob.tox_window_days + 1e-12
            push!(tox_vals, prob.il6_obs_fun(u_curr, t_curr))
        end

        if save_trajectory
            if isempty(t_all) || !isapprox(t_curr, t_all[end]; atol = 1e-12, rtol = 0.0)
                push!(t_all, t_curr)
                push!(bt_all, Float64(u_curr[prob.bt_idx]))
                push!(il6_all, Float64(prob.il6_obs_fun(u_curr, t_curr)))
            end
        end
    end

    if t_curr < prob.horizon_days - 1e-12
        saveat =
            if save_trajectory
                collect(t_curr:0.1:prob.horizon_days)
            else
                tox_pts = [tt for tt in prob.tox_grid if tt > t_curr + 1e-12 && tt <= prob.tox_window_days + 1e-12]
                sort(unique(vcat([prob.horizon_days], tox_pts)))
            end
        sol_post = solve_segment(prob, u_curr, (t_curr, prob.horizon_days), saveat)
        for j in eachindex(sol_post.t)
            tj = Float64(sol_post.t[j])
            uj = sol_post.u[j]
            il6v = prob.il6_obs_fun(uj, sol_post.t[j])
            if tj <= prob.tox_window_days + 1e-12
                push!(tox_vals, il6v)
            end
            if save_trajectory
                if !isempty(t_all) && isapprox(tj, t_all[end]; atol = 1e-12, rtol = 0.0)
                    continue
                end
                push!(t_all, tj)
                push!(bt_all, Float64(uj[prob.bt_idx]))
                push!(il6_all, Float64(il6v))
            end
        end
        u_curr = sol_post.u[end]
    end

    isempty(tox_vals) && error("No IL6 samples captured in toxicity window.")
    tox_proxy = smoothmax(tox_vals, prob.tox_tau)
    tumor_proxy = u_curr[prob.bt_idx] / (bt0 + T(1e-12))
    loss = T(prob.w_tox) * (tox_proxy / T(max(prob.tox_scale, 1e-12))) +
           T(prob.w_tumor) * (tumor_proxy / T(max(prob.tumor_scale, 1e-12)))

    if !save_trajectory
        return (; tox_proxy = tox_proxy, tumor_proxy = tumor_proxy, loss = loss)
    end
    return (; tox_proxy = tox_proxy, tumor_proxy = tumor_proxy, loss = loss, t = t_all, btumor = bt_all, il6combo = il6_all)
end

function make_loss_objective(prob::PatientProblem)
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
            return simulate_patient(prob, x; save_trajectory = false).loss
        catch
            if x[1] isa ForwardDiff.Dual
                return one(x[1]) * big_penalty
            end
            return big_penalty
        end
    end
end

function build_patient_model_rows()
    repo_root = TCellEngagerQSP.REPO_ROOT
    design_dir = joinpath(repo_root, "generated", "phase1_design")
    patients = DataFrame(CSV.File(joinpath(design_dir, "patients.csv")))
    overrides = DataFrame(CSV.File(joinpath(design_dir, "dlbcl_param_overrides.csv")))

    patient_ids_env = strip(get(ENV, "OPT2C_PATIENT_IDS", ""))
    if !isempty(patient_ids_env)
        keep_ids = Set(parse_int_list(patient_ids_env))
        sel = [Int(pid) in keep_ids for pid in patients.patient_id]
        patients = patients[sel, :]
        sort!(patients, :patient_id)
    else
        max_patients = parse(Int, get(ENV, "OPT2C_MAX_PATIENTS", string(nrow(patients))))
        if max_patients < nrow(patients)
            patients = patients[1:max_patients, :]
        end
    end

    rows = NamedTuple[]
    for prow in eachrow(patients)
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
        pmap["end_time"] = 42.0

        pnames = sort(collect(keys(pmap)))
        pvals = [pmap[n] for n in pnames]

        mdl = TCellEngagerQSP.build_model_with_variant_ids(Int[], pnames, pvals)
        rhs_fun = build_ad_rhs_no_cast(mdl)

        bt_idx = mdl.state_to_idx["Btumor"]
        il6_obs_fun = build_symbol_observer_no_cast(mdl, "IL6combo")
        tdbc_idx = mdl.state_to_idx["TDBc_ugperkg"]
        e_tdbc = zeros(Float64, length(mdl.u0))
        e_tdbc[tdbc_idx] = 1.0

        push!(rows, (patient_id = Int(prow.patient_id), mdl = mdl, rhs_fun = rhs_fun, il6_obs_fun = il6_obs_fun, bt_idx = bt_idx, e_tdbc = e_tdbc))
    end
    return rows
end

function instantiate_problem(
    row;
    dose_times_days::Vector{Float64},
    fixed_doses_mg::Vector{Float64},
    lower_mg::Vector{Float64},
    upper_mg::Vector{Float64},
    tox_scale::Float64,
    tumor_scale::Float64,
    alg,
    horizon_days::Float64,
    tox_window_days::Float64,
    tox_grid::Vector{Float64},
    tox_tau::Float64,
    w_tox::Float64,
    w_tumor::Float64,
    abstol::Float64,
    reltol::Float64,
    maxiters::Int,
)
    return PatientProblem(
        row.patient_id,
        row.mdl,
        row.rhs_fun,
        row.il6_obs_fun,
        alg,
        row.e_tdbc,
        row.bt_idx,
        dose_times_days,
        fixed_doses_mg,
        lower_mg,
        upper_mg,
        1000.0 / 70.0,
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

function main()
    rows = build_patient_model_rows()
    n_patients = length(rows)
    n_patients == 0 && error("No patients loaded.")

    dose_times_days = parse_float_list(get(ENV, "OPT2C_DOSE_TIMES_DAYS", "0,7,14,21,28,35"))
    fixed_doses_mg = parse_fixed_doses()
    length(dose_times_days) == length(fixed_doses_mg) || error("Dose times and fixed dose vectors must have equal length.")

    horizon_days = parse(Float64, get(ENV, "OPT2C_HORIZON_DAYS", "42.0"))
    tox_window_days = parse(Float64, get(ENV, "OPT2C_TOX_WINDOW_DAYS", "2.0"))
    tox_grid = collect(0.0:0.1:tox_window_days)
    tox_tau = parse(Float64, get(ENV, "OPT2C_TOX_SOFTMAX_TAU", "50.0"))
    w_tox = parse(Float64, get(ENV, "OPT2C_LOSS_W_TOX", "0.5"))
    w_tumor = parse(Float64, get(ENV, "OPT2C_LOSS_W_TUMOR", "0.5"))
    abstol = parse(Float64, get(ENV, "OPT2C_ABSTOL", "1e-8"))
    reltol = parse(Float64, get(ENV, "OPT2C_RELTOL", "1e-6"))
    maxiters = parse(Int, get(ENV, "OPT2C_MAXITERS_SOLVE", "1000000"))
    opt_method = lowercase(get(ENV, "OPT2C_OPT_METHOD", "coordinate_descent"))
    lbfgs_iters = parse(Int, get(ENV, "OPT2C_LBFGS_ITERS", "80"))
    cd_sweeps = parse(Int, get(ENV, "OPT2C_CD_SWEEPS", "4"))
    cd_step0 = parse(Float64, get(ENV, "OPT2C_CD_STEP0_MG", "2.0"))
    cd_tol = parse(Float64, get(ENV, "OPT2C_CD_TOL_MG", "0.05"))
    alg = TCellEngagerQSP.make_solver_alg(lowercase(get(ENV, "OPT2C_SOLVER", "tsit5")))

    lower_scalar = parse(Float64, get(ENV, "OPT2C_DOSE_LOWER_MG", "0.0"))
    upper_scalar = parse(Float64, get(ENV, "OPT2C_DOSE_UPPER_MG", "30.0"))
    lower_mg = fill(lower_scalar, length(dose_times_days))
    upper_mg = fill(upper_scalar, length(dose_times_days))

    println(@sprintf("patients=%d, doses_per_patient=%d", n_patients, length(dose_times_days)))
    println("fixed_doses_mg=$(fixed_doses_mg)")

    # Pass 1: fixed-regimen proxies to calibrate loss scales.
    fixed_tox = Float64[]
    fixed_tumor = Float64[]
    for row in rows
        ptmp = instantiate_problem(
            row;
            dose_times_days = dose_times_days,
            fixed_doses_mg = fixed_doses_mg,
            lower_mg = lower_mg,
            upper_mg = upper_mg,
            tox_scale = 1.0,
            tumor_scale = 1.0,
            alg = alg,
            horizon_days = horizon_days,
            tox_window_days = tox_window_days,
            tox_grid = tox_grid,
            tox_tau = tox_tau,
            w_tox = w_tox,
            w_tumor = w_tumor,
            abstol = abstol,
            reltol = reltol,
            maxiters = maxiters,
        )
        sim = simulate_patient(ptmp, fixed_doses_mg; save_trajectory = false)
        push!(fixed_tox, Float64(sim.tox_proxy))
        push!(fixed_tumor, Float64(sim.tumor_proxy))
    end

    tox_scale_floor = parse(Float64, get(ENV, "OPT2C_TOX_SCALE_FLOOR", "1e-3"))
    tumor_scale_floor = parse(Float64, get(ENV, "OPT2C_TUMOR_SCALE_FLOOR", "1e-12"))
    tox_scale_default = robust_scale(fixed_tox, tox_scale_floor, 100.0)
    tumor_scale_default = robust_scale(fixed_tumor, tumor_scale_floor, 0.5)
    tox_scale = parse(Float64, get(ENV, "OPT2C_TOX_SCALE", string(tox_scale_default)))
    tumor_scale = parse(Float64, get(ENV, "OPT2C_TUMOR_SCALE", string(tumor_scale_default)))
    println(@sprintf("loss_scales: tox_scale=%.6g tumor_scale=%.6g", tox_scale, tumor_scale))

    out_dir_default = joinpath(TCellEngagerQSP.REPO_ROOT, "generated", "figures", "optimization", "two_cycle_population")
    out_dir_env = get(ENV, "OPT2C_OUT_DIR", out_dir_default)
    out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(TCellEngagerQSP.REPO_ROOT, out_dir_env)
    mkpath(out_dir)

    summary_rows = NamedTuple[]
    nd = length(dose_times_days)
    dcols = Symbol[]
    for i in 1:nd
        push!(dcols, Symbol("dose_day$(Int(round(dose_times_days[i])))_mg"))
    end

    for (k, row) in enumerate(rows)
        prob = instantiate_problem(
            row;
            dose_times_days = dose_times_days,
            fixed_doses_mg = fixed_doses_mg,
            lower_mg = lower_mg,
            upper_mg = upper_mg,
            tox_scale = tox_scale,
            tumor_scale = tumor_scale,
            alg = alg,
            horizon_days = horizon_days,
            tox_window_days = tox_window_days,
            tox_grid = tox_grid,
            tox_tau = tox_tau,
            w_tox = w_tox,
            w_tumor = w_tumor,
            abstol = abstol,
            reltol = reltol,
            maxiters = maxiters,
        )

        fixed = simulate_patient(prob, fixed_doses_mg; save_trajectory = false)

        obj = make_loss_objective(prob)

        x0 = clamp.(copy(fixed_doses_mg), lower_mg, upper_mg)
        status = "ok"
        t0 = time()
        x_opt = copy(x0)
        loss_opt = Float64(fixed.loss)
        tox_opt = Float64(fixed.tox_proxy)
        tumor_opt = Float64(fixed.tumor_proxy)
        try
            if opt_method == "lbfgs"
                g! = function (G, x)
                    ForwardDiff.gradient!(G, obj, x)
                end
                opt_opts = Optim.Options(iterations = lbfgs_iters, show_trace = false, store_trace = false)
                res = optimize(obj, g!, lower_mg, upper_mg, copy(x0), Fminbox(LBFGS()), opt_opts)
                x_opt = clamp.(Optim.minimizer(res), lower_mg, upper_mg)
            else
                x_curr = copy(x0)
                loss_curr = Float64(obj(x_curr))
                step = cd_step0
                for _ in 1:cd_sweeps
                    improved = false
                    for j in eachindex(x_curr)
                        x_plus = copy(x_curr)
                        x_plus[j] = min(upper_mg[j], x_curr[j] + step)
                        l_plus = Float64(obj(x_plus))

                        x_minus = copy(x_curr)
                        x_minus[j] = max(lower_mg[j], x_curr[j] - step)
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
                    step *= 0.5
                    if step < cd_tol && !improved
                        break
                    end
                end
                x_opt = x_curr
            end
            sim_opt = simulate_patient(prob, x_opt; save_trajectory = false)
            loss_opt = Float64(sim_opt.loss)
            tox_opt = Float64(sim_opt.tox_proxy)
            tumor_opt = Float64(sim_opt.tumor_proxy)
        catch err
            status = "fail: $(typeof(err))"
        end
        runtime_s = time() - t0

        dose_nt = NamedTuple{Tuple(dcols)}(Tuple(Float64.(x_opt)))
        push!(
            summary_rows,
            merge(
                (
                    patient_id = prob.patient_id,
                    status = status,
                    runtime_s = runtime_s,
                    loss_fixed = Float64(fixed.loss),
                    loss_opt = loss_opt,
                    tox_fixed = Float64(fixed.tox_proxy),
                    tox_opt = tox_opt,
                    tumor_fixed = Float64(fixed.tumor_proxy),
                    tumor_opt = tumor_opt,
                ),
                dose_nt,
            ),
        )

        if k % 10 == 0 || k == n_patients
            println(@sprintf("processed %d/%d patients", k, n_patients))
        end
    end

    summary_df = DataFrame(summary_rows)
    summary_path = joinpath(out_dir, "two_cycle_population_optimization_summary.csv")
    CSV.write(summary_path, summary_df)

    meta = Dict(
        "n_patients" => n_patients,
        "patient_ids" => [r.patient_id for r in rows],
        "dose_times_days" => dose_times_days,
        "fixed_doses_mg" => fixed_doses_mg,
        "label_regimen_name" => get(ENV, "OPT2C_LABEL_REGIMEN_NAME", "step_0.8_2_6mg"),
        "horizon_days" => horizon_days,
        "tox_window_days" => tox_window_days,
        "loss_formula" => @sprintf(
            "L=%.2f*(smoothmax(IL6combo[0-%.1f],tau=%.1f)/%.6g) + %.2f*((Btumor_day%.0f/Btumor_day0)/%.6g)",
            w_tox,
            tox_window_days,
            tox_tau,
            tox_scale,
            w_tumor,
            horizon_days,
            tumor_scale,
        ),
        "solver" => string(alg),
        "optimization_method" => opt_method,
        "coordinate_descent" => Dict(
            "sweeps" => cd_sweeps,
            "step0_mg" => cd_step0,
            "tol_mg" => cd_tol,
        ),
        "lbfgs_iterations" => lbfgs_iters,
        "dose_bounds_mg" => Dict("lower" => lower_scalar, "upper" => upper_scalar),
        "n_failures" => sum(.!occursin.("ok", summary_df.status)),
    )
    meta_path = joinpath(out_dir, "two_cycle_population_optimization_meta.json")
    open(meta_path, "w") do io
        JSON3.pretty(io, meta)
    end

    println("Wrote $summary_path")
    println("Wrote $meta_path")
end

main()
