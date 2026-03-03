using CSV
using DataFrames
using DifferentialEquations
using DiffEqCallbacks
using Sundials
using SciMLBase

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

repo_root = TCellEngagerQSP.REPO_ROOT
design_dir_default = joinpath(repo_root, "generated", "phase1_design")
design_dir_env = get(ENV, "PHASE1_DESIGN_DIR", design_dir_default)
design_dir = isabspath(design_dir_env) ? design_dir_env : joinpath(repo_root, design_dir_env)
patients_path = joinpath(design_dir, "patients.csv")
reg_path = joinpath(design_dir, "regimen_events.csv")

if !isfile(patients_path) || !isfile(reg_path)
    error("Design files missing. Run scripts/generate_phase1_design.py first.")
end

patients = DataFrame(CSV.File(patients_path))
reg_events = DataFrame(CSV.File(reg_path))

out_dir_default = joinpath(repo_root, "generated", "phase1_julia")
out_dir_env = get(ENV, "PHASE1_JULIA_OUT_DIR", out_dir_default)
out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(repo_root, out_dir_env)
mkpath(out_dir)

phase1_active_variant_ids = [5, 9, 14, 20, 24, 25, 27, 28]
variants_mode = lowercase(get(ENV, "PHASE1_VARIANTS_MODE", "matlab_empty"))
variant_ids =
    if variants_mode == "matlab_empty"
        # Match MATLAB script behavior in run_phase1_mosu_matlab_reference.m:
        # sbiosimulate(model, cs, [], d) -> no variants applied.
        Int[]
    elseif variants_mode == "active_phase1"
        phase1_active_variant_ids
    else
        error("Unsupported PHASE1_VARIANTS_MODE=$variants_mode (use matlab_empty or active_phase1)")
    end
bw_kg_for_dose_conversion = 70.0
dlbcl_overrides_path = joinpath(design_dir, "dlbcl_param_overrides.csv")
dlbcl_overrides = isfile(dlbcl_overrides_path) ? DataFrame(CSV.File(dlbcl_overrides_path)) : DataFrame(name = String[], value = Float64[])

solver_name_full = lowercase(get(ENV, "TCE_SOLVER", "cvode_bdf"))
use_mtk_jac = endswith(solver_name_full, "_mtkjac") || get(ENV, "TCE_MTK_JAC", "0") == "1"
solver_name = replace(solver_name_full, "_mtkjac" => "")
hybrid_mode = startswith(solver_name, "hybrid")
alg = hybrid_mode ? nothing : TCellEngagerQSP.make_solver_alg(solver_name)
hybrid_stiff_alg = hybrid_mode ? TCellEngagerQSP.make_solver_alg(get(ENV, "TCE_HYBRID_STIFF_SOLVER", "cvode_bdf")) : nothing
hybrid_nonstiff_alg = hybrid_mode ? TCellEngagerQSP.make_solver_alg(get(ENV, "TCE_HYBRID_NONSTIFF_SOLVER", "tsit5")) : nothing
hybrid_stiff_window_days = parse(Float64, get(ENV, "TCE_HYBRID_STIFF_WINDOW_DAYS", "2.0"))
engine = Symbol(lowercase(get(ENV, "TCE_ENGINE", "canonical")))
if !(engine in (:legacy, :canonical))
    error("Unsupported TCE_ENGINE=$engine. Use legacy or canonical.")
end

function eval_symbol(mdl, u::Vector{Float64}, t::Float64, name::String)
    ns = length(mdl.state_names)
    np = length(mdl.pvals)
    z = copy(mdl.z)
    z[1:ns] .= u
    z[ns+1:ns+np] .= mdl.pvals
    for r in mdl.repeated_rule_exprs
        z[r.lhs_idx] = Base.invokelatest(r.fn, z, t)
    end
    return z[mdl.name_to_idx[name]]
end

function interp1_linear(x::Vector{Float64}, y::Vector{Float64}, tq::Float64)
    if tq <= x[1]
        return y[1]
    elseif tq >= x[end]
        return y[end]
    end
    i = searchsortedlast(x, tq)
    x1, x2 = x[i], x[i + 1]
    y1, y2 = y[i], y[i + 1]
    if x2 == x1
        return y1
    end
    return y1 + (y2 - y1) * (tq - x1) / (x2 - x1)
end

function window_peak(t::Vector{Float64}, y::Vector{Float64}, t0::Float64, t1::Float64)
    if isempty(t) || isempty(y)
        return NaN
    end
    if t1 < t0
        t0, t1 = t1, t0
    end
    vals = Float64[]
    push!(vals, interp1_linear(t, y, t0))
    push!(vals, interp1_linear(t, y, t1))
    for i in eachindex(t)
        ti = t[i]
        if t0 <= ti <= t1
            push!(vals, y[i])
        end
    end
    isempty(vals) && return NaN
    return maximum(vals)
end

function window_auc(t::Vector{Float64}, y::Vector{Float64}, t0::Float64, t1::Float64)
    if isempty(t) || isempty(y)
        return NaN
    end
    if t1 < t0
        t0, t1 = t1, t0
    end
    if t1 == t0
        return 0.0
    end

    tw = Float64[t0]
    yw = Float64[interp1_linear(t, y, t0)]
    for i in eachindex(t)
        ti = t[i]
        if t0 < ti < t1
            push!(tw, ti)
            push!(yw, y[i])
        end
    end
    push!(tw, t1)
    push!(yw, interp1_linear(t, y, t1))

    ord = sortperm(tw)
    tw = tw[ord]
    yw = yw[ord]

    auc = 0.0
    for i in 1:(length(tw) - 1)
        dt = tw[i + 1] - tw[i]
        if dt > 0.0
            auc += 0.5 * dt * (yw[i] + yw[i + 1])
        end
    end
    return auc
end

function append_solution!(sol_t::Vector{Float64}, sol_u::Vector{Vector{Float64}}, ts::Vector{Float64}, us::Vector)
    isempty(ts) && return
    if isempty(sol_t)
        append!(sol_t, Float64.(ts))
        for u in us
            push!(sol_u, Float64.(u))
        end
        return
    end
    start_idx = isapprox(ts[1], sol_t[end]; atol = 1e-12, rtol = 0.0) ? 2 : 1
    for i in start_idx:length(ts)
        push!(sol_t, Float64(ts[i]))
        push!(sol_u, Float64.(us[i]))
    end
end

function run_segment!(
        sol_t::Vector{Float64},
        sol_u::Vector{Vector{Float64}},
        ode_rhs,
        ctx,
        u_start::Vector{Float64},
        t0::Float64,
        t1::Float64,
        seg_alg)
    if t1 <= t0
        return u_start
    end
    prob = ODEProblem(ode_rhs, u_start, (t0, t1), ctx)
    sol = solve(
        prob,
        seg_alg;
        abstol = TCellEngagerQSP.SOLVER_ABSTOL,
        reltol = TCellEngagerQSP.SOLVER_RELTOL,
        tstops = [t1],
    )
    if sol.retcode != SciMLBase.ReturnCode.Success
        error("segment solve failed with retcode=$(sol.retcode)")
    end
    append_solution!(sol_t, sol_u, sol.t, sol.u)
    return Float64.(sol.u[end])
end

function run_hybrid_segmented(
        ode_rhs,
        ctx,
        u0::Vector{Float64},
        dose_map::Dict{Float64,Float64},
        target_idx::Int,
        tf::Float64,
        stiff_alg,
        nonstiff_alg,
        stiff_window_days::Float64)
    local_dose_map = copy(dose_map)
    dose_at_t0 = get(local_dose_map, 0.0, 0.0)
    if dose_at_t0 != 0.0
        local_dose_map[0.0] = 0.0
    end

    u_curr = copy(u0)
    if dose_at_t0 != 0.0
        u_curr[target_idx] += dose_at_t0
    end

    dose_times = sort([t for t in keys(local_dose_map) if t > 0.0])
    segment_ends = sort(unique(vcat(dose_times, [tf])))
    stiff_starts = Set(vcat([0.0], dose_times))

    sol_t = Float64[]
    sol_u = Vector{Vector{Float64}}()
    t_curr = 0.0

    for te in segment_ends
        if te > t_curr
            stiff_here = any(isapprox(t_curr, ts; atol = 1e-10, rtol = 0.0) for ts in stiff_starts)
            if stiff_here && stiff_window_days > 0.0
                t_mid = min(te, t_curr + stiff_window_days)
                u_curr = run_segment!(sol_t, sol_u, ode_rhs, ctx, u_curr, t_curr, t_mid, stiff_alg)
                if te > t_mid
                    u_curr = run_segment!(sol_t, sol_u, ode_rhs, ctx, u_curr, t_mid, te, nonstiff_alg)
                end
            else
                u_curr = run_segment!(sol_t, sol_u, ode_rhs, ctx, u_curr, t_curr, te, nonstiff_alg)
            end
        end
        if haskey(local_dose_map, te)
            u_curr[target_idx] += local_dose_map[te]
        end
        t_curr = te
    end
    return sol_t, sol_u
end

rows = DataFrame(
    regimen = String[],
    regimen_type = String[],
    patient_id = Int[],
    il6_peak_0_2 = Float64[],
    il6_peak_0_21 = Float64[],
    il6_auc_0_2 = Float64[],
    il6_auc_0_21 = Float64[],
    btumor_t0 = Float64[],
    btumor_day42 = Float64[],
    tumor_resid_day42 = Float64[],
    tumor_cfbl_day42 = Float64[],
    loss_raw = Float64[],
    bt_ratio_tumor_init = Float64[],
    status = String[],
)

regimens = unique(reg_events.regimen)
max_regimens = parse(Int, get(ENV, "PHASE1_MAX_REGIMENS", string(length(regimens))))
max_patients = parse(Int, get(ENV, "PHASE1_MAX_PATIENTS", string(nrow(patients))))
regimens = regimens[1:min(end, max_regimens)]
for reg in regimens
    sub = reg_events[reg_events.regimen .== reg, :]
    sort!(sub, :event_idx)

    dose_map = Dict{Float64, Float64}()
    for rr in eachrow(sub)
        t = Float64(rr.time_day)
        amt = Float64(rr.dose_mg) * 1000.0 / bw_kg_for_dose_conversion
        dose_map[t] = get(dose_map, t, 0.0) + amt
    end
    dose_times = sort(collect(keys(dose_map)))
    reg_type = String(sub.regimen_type[1])

    for prow in eachrow(patients[1:min(end, max_patients), :])
        pmap = Dict{String, Float64}()
        for rr in eachrow(dlbcl_overrides)
            pmap[String(rr.name)] = Float64(rr.value)
        end
        for cname in names(patients)
            if cname == "patient_id"
                continue
            end
            raw = prow[cname]
            if !ismissing(raw)
                v = try
                    Float64(raw)
                catch
                    NaN
                end
                if isfinite(v)
                    pmap[cname] = v
                end
            end
        end
        pmap["PKflag"] = 1.0
        pmap["fvalidation"] = 0.0
        pmap["VPid"] = 1.0
        pmap["end_time"] = 84.0
        pnames = sort(collect(keys(pmap)))
        pvals = [pmap[n] for n in pnames]

        try
            mdl = TCellEngagerQSP.build_model_with_variant_ids(variant_ids, pnames, pvals)
            rhs_fun = TCellEngagerQSP.rhs!
            ctx = nothing
            if engine == :legacy
                ctx = TCellEngagerQSP.SimContext(
                    copy(mdl.z),
                    mdl.pvals,
                    mdl.repeated_rule_exprs,
                    mdl.rate_fns,
                    mdl.stoich,
                    TCellEngagerQSP.InfusionEvent[],
                )
                rhs_fun = TCellEngagerQSP.rhs!
            else
                ctx = TCellEngagerQSP.CanonicalSimContext(
                    copy(mdl.z),
                    mdl.pvals,
                    TCellEngagerQSP.InfusionEvent[],
                )
                rhs_fun = mdl.canonical_rhs
            end
            jac_fun = use_mtk_jac ? TCellEngagerQSP.make_mtk_dense_jacobian(mdl) : nothing
            ode_rhs = isnothing(jac_fun) ? rhs_fun : ODEFunction(rhs_fun; jac = jac_fun)

            target_idx = mdl.state_to_idx["TDBc_ugperkg"]

            local_dose_map = copy(dose_map)
            dose_at_t0 = get(local_dose_map, 0.0, 0.0)
            if dose_at_t0 != 0.0
                local_dose_map[0.0] = 0.0
            end
            cb_times = filter(t -> t > 0.0, dose_times)

            function affect!(integrator)
                t = integrator.t
                amt = get(local_dose_map, t, NaN)
                if isnan(amt)
                    for (tt, vv) in local_dose_map
                        if isapprox(t, tt; atol = 1e-8, rtol = 0.0)
                            amt = vv
                            break
                        end
                    end
                end
                if !isnan(amt)
                    integrator.u[target_idx] += amt
                end
            end

            u0 = copy(mdl.u0)
            if dose_at_t0 != 0.0
                u0[target_idx] += dose_at_t0
            end
            il6_series = Float64[]
            bt_series = Float64[]
            t_series = Float64[]
            if hybrid_mode
                sol_t, sol_u = run_hybrid_segmented(
                    ode_rhs,
                    ctx,
                    u0,
                    local_dose_map,
                    target_idx,
                    84.0,
                    hybrid_stiff_alg,
                    hybrid_nonstiff_alg,
                    hybrid_stiff_window_days,
                )
                for i in eachindex(sol_t)
                    t = Float64(sol_t[i])
                    u = Float64.(sol_u[i])
                    push!(t_series, t)
                    push!(il6_series, eval_symbol(mdl, u, t, "IL6combo"))
                    push!(bt_series, eval_symbol(mdl, u, t, "Btumor"))
                end
            else
                cb = PresetTimeCallback(cb_times, affect!; save_positions = (true, true))
                prob = ODEProblem(ode_rhs, u0, (0.0, 84.0), ctx)
                sol = solve(
                    prob,
                    alg;
                    abstol = TCellEngagerQSP.SOLVER_ABSTOL,
                    reltol = TCellEngagerQSP.SOLVER_RELTOL,
                    callback = cb,
                    tstops = cb_times,
                    d_discontinuities = cb_times,
                )
                if sol.retcode != SciMLBase.ReturnCode.Success
                    error("solve failed with retcode=$(sol.retcode)")
                end
                for i in eachindex(sol.t)
                    t = Float64(sol.t[i])
                    u = Float64.(sol.u[i])
                    push!(t_series, t)
                    push!(il6_series, eval_symbol(mdl, u, t, "IL6combo"))
                    push!(bt_series, eval_symbol(mdl, u, t, "Btumor"))
                end
            end

            bt0 = bt_series[1]
            tuniq = Float64[]
            btuniq = Float64[]
            il6uniq = Float64[]
            seen = Set{Float64}()
            for (t, b, il6v) in zip(t_series, bt_series, il6_series)
                if !(t in seen)
                    push!(seen, t)
                    push!(tuniq, t)
                    push!(btuniq, b)
                    push!(il6uniq, il6v)
                end
            end
            il6_peak_0_2 = window_peak(tuniq, il6uniq, 0.0, 2.0)
            il6_peak_0_21 = window_peak(tuniq, il6uniq, 0.0, 21.0)
            il6_auc_0_2 = window_auc(tuniq, il6uniq, 0.0, 2.0)
            il6_auc_0_21 = window_auc(tuniq, il6uniq, 0.0, 21.0)
            bt42 = interp1_linear(tuniq, btuniq, 42.0)
            tumor_resid_day42 = bt42 / max(bt0, 1e-12)
            tumor_cfbl_day42 = (bt42 - bt0) / max(bt0, 1e-12)
            loss_raw = 0.5 * log10(1 + max(il6_peak_0_2, 0.0)) + 0.5 * tumor_resid_day42

            push!(
                rows,
                (
                    String(reg),
                    reg_type,
                    Int(prow.patient_id),
                    il6_peak_0_2,
                    il6_peak_0_21,
                    il6_auc_0_2,
                    il6_auc_0_21,
                    bt0,
                    bt42,
                    tumor_resid_day42,
                    tumor_cfbl_day42,
                    loss_raw,
                    Float64(prow.BT_ratio_tumor_init),
                    "ok",
                ),
            )
        catch err
            push!(
                rows,
                (
                    String(reg),
                    reg_type,
                    Int(prow.patient_id),
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    Float64(prow.BT_ratio_tumor_init),
                    "error",
                ),
            )
            @warn "Simulation failed" regimen = reg patient = Int(prow.patient_id) err
        end
    end
end

out_path = joinpath(out_dir, "phase1_metrics_julia.csv")
CSV.write(out_path, rows)
println("Wrote " * out_path)
