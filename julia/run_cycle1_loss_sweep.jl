using CSV
using DataFrames
using DifferentialEquations
using DiffEqCallbacks
using SciMLBase

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

function smoothmax(v::Vector{Float64}, tau::Float64)
    m = maximum(v)
    return m + tau * log(sum(exp.((v .- m) ./ tau)))
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

repo_root = TCellEngagerQSP.REPO_ROOT
design_dir_default = joinpath(repo_root, "generated", "phase1_design")
design_dir_env = get(ENV, "PHASE1_DESIGN_DIR", design_dir_default)
design_dir = isabspath(design_dir_env) ? design_dir_env : joinpath(repo_root, design_dir_env)
patients_path = joinpath(design_dir, "patients.csv")
reg_path = joinpath(design_dir, "regimen_events.csv")
dlbcl_overrides_path = joinpath(design_dir, "dlbcl_param_overrides.csv")

if !isfile(patients_path) || !isfile(reg_path) || !isfile(dlbcl_overrides_path)
    error("Design files missing. Run scripts/generate_phase1_design.py first.")
end

patients = DataFrame(CSV.File(patients_path))
reg_events = DataFrame(CSV.File(reg_path))
dlbcl_overrides = DataFrame(CSV.File(dlbcl_overrides_path))

out_dir_default = joinpath(repo_root, "generated", "cycle1_loss")
out_dir_env = get(ENV, "CYCLE1_LOSS_OUT_DIR", out_dir_default)
out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(repo_root, out_dir_env)
mkpath(out_dir)

bw_kg_for_dose_conversion = 70.0
horizon_days = parse(Float64, get(ENV, "CYCLE1_HORIZON_DAYS", "21.0"))
tox_window_days = parse(Float64, get(ENV, "CYCLE1_TOX_WINDOW_DAYS", "2.0"))
tox_tau = parse(Float64, get(ENV, "CYCLE1_TOX_SOFTMAX_TAU", "50.0"))
tox_scale = parse(Float64, get(ENV, "CYCLE1_TOX_SCALE", "152.22612188617114"))
tumor_scale = parse(Float64, get(ENV, "CYCLE1_TUMOR_SCALE", "0.7050236467345111"))
loss_weight_tox = parse(Float64, get(ENV, "CYCLE1_LOSS_W_TOX", "0.5"))
loss_weight_tumor = parse(Float64, get(ENV, "CYCLE1_LOSS_W_TUMOR", "0.5"))

solver_name = lowercase(get(ENV, "TCE_SOLVER", "cvode_bdf"))
alg = TCellEngagerQSP.make_solver_alg(solver_name)
engine = Symbol(lowercase(get(ENV, "TCE_ENGINE", "canonical")))
if !(engine in (:legacy, :canonical))
    error("Unsupported TCE_ENGINE=$engine. Use legacy or canonical.")
end

rows = DataFrame(
    regimen = String[],
    regimen_type = String[],
    patient_id = Int[],
    tox_proxy_smax_day0_2 = Float64[],
    tumor_proxy_day21_over_day0 = Float64[],
    loss_cycle1 = Float64[],
    status = String[],
)

regimens = unique(reg_events.regimen)
max_regimens = parse(Int, get(ENV, "PHASE1_MAX_REGIMENS", string(length(regimens))))
max_patients = parse(Int, get(ENV, "PHASE1_MAX_PATIENTS", string(nrow(patients))))
regimens = regimens[1:min(end, max_regimens)]
maxiters = parse(Int, get(ENV, "CYCLE1_MAXITERS", "2000000"))

regimen_specs = NamedTuple[]
for reg in regimens
    sub = reg_events[reg_events.regimen .== reg, :]
    sort!(sub, :event_idx)
    reg_type = String(sub.regimen_type[1])
    dose_map = Dict{Float64, Float64}()
    for rr in eachrow(sub)
        t = Float64(rr.time_day)
        if t < horizon_days
            amt = Float64(rr.dose_mg) * 1000.0 / bw_kg_for_dose_conversion
            dose_map[t] = get(dose_map, t, 0.0) + amt
        end
    end
    dose_times = sort(collect(keys(dose_map)))
    push!(regimen_specs, (regimen = String(reg), regimen_type = reg_type, dose_map = dose_map, dose_times = dose_times))
end

tox_grid = collect(0.0:0.1:tox_window_days)
saveat = sort(unique(vcat(tox_grid, [horizon_days])))

for prow in eachrow(patients[1:min(end, max_patients), :])
    pmap = Dict{String, Float64}()
    for rr in eachrow(dlbcl_overrides)
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
    pmap["end_time"] = horizon_days
    pnames = sort(collect(keys(pmap)))
    pvals = [pmap[n] for n in pnames]

    try
        mdl = TCellEngagerQSP.build_model_with_variant_ids(Int[], pnames, pvals)
        target_idx = mdl.state_to_idx["TDBc_ugperkg"]
        btumor_idx = mdl.state_to_idx["Btumor"]

        for spec in regimen_specs
            reg = spec.regimen
            reg_type = spec.regimen_type
            dose_map = spec.dose_map
            dose_times = spec.dose_times

            try
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
                bt0 = u0[btumor_idx]

                cb = PresetTimeCallback(cb_times, affect!; save_positions = (true, true))
                prob = ODEProblem(rhs_fun, u0, (0.0, horizon_days), ctx)
                sol = solve(
                    prob,
                    alg;
                    abstol = TCellEngagerQSP.SOLVER_ABSTOL,
                    reltol = TCellEngagerQSP.SOLVER_RELTOL,
                    callback = cb,
                    tstops = cb_times,
                    d_discontinuities = cb_times,
                    saveat = saveat,
                    maxiters = maxiters,
                )
                if sol.retcode != SciMLBase.ReturnCode.Success
                    error("solve failed with retcode=$(sol.retcode)")
                end

                il6_vals = Float64[]
                bt_end = NaN
                for i in eachindex(sol.t)
                    t = Float64(sol.t[i])
                    u = Float64.(sol.u[i])
                    if t <= tox_window_days + 1e-12
                        push!(il6_vals, eval_symbol(mdl, u, t, "IL6combo"))
                    end
                    if isapprox(t, horizon_days; atol = 1e-12, rtol = 0.0)
                        bt_end = u[btumor_idx]
                    end
                end
                if isnan(bt_end)
                    bt_end = Float64(sol.u[end][btumor_idx])
                end
                tox_proxy = smoothmax(il6_vals, tox_tau)
                tumor_proxy = bt_end / max(bt0, 1e-12)

                loss_cycle1 = loss_weight_tox * (tox_proxy / max(tox_scale, 1e-12)) +
                              loss_weight_tumor * (tumor_proxy / max(tumor_scale, 1e-12))

                push!(
                    rows,
                    (
                        reg,
                        reg_type,
                        Int(prow.patient_id),
                        tox_proxy,
                        tumor_proxy,
                        loss_cycle1,
                        "ok",
                    ),
                )
            catch err
                push!(
                    rows,
                    (
                        reg,
                        reg_type,
                        Int(prow.patient_id),
                        NaN,
                        NaN,
                        NaN,
                        "error",
                    ),
                )
                @warn "Cycle-1 loss simulation failed" regimen = reg patient = Int(prow.patient_id) err
            end
        end
    catch err
        for spec in regimen_specs
            push!(
                rows,
                (
                    spec.regimen,
                    spec.regimen_type,
                    Int(prow.patient_id),
                    NaN,
                    NaN,
                    NaN,
                    "error",
                ),
            )
        end
        @warn "Model build failed for patient; marking all regimens as error" patient = Int(prow.patient_id) err
    end
end

out_path = joinpath(out_dir, "cycle1_loss_metrics_julia.csv")
CSV.write(out_path, rows)
println("Wrote " * out_path)
