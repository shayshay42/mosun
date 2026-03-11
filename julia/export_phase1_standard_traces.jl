using CSV
using DataFrames
using DifferentialEquations
using DiffEqCallbacks
using Sundials
using SciMLBase

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

repo_root = TCellEngagerQSP.REPO_ROOT
design_dir_default = joinpath(repo_root, "generated", "phase1_design_standard")
design_dir_env = get(ENV, "PHASE1_DESIGN_DIR", design_dir_default)
design_dir = isabspath(design_dir_env) ? design_dir_env : joinpath(repo_root, design_dir_env)

out_dir_default = joinpath(repo_root, "generated", "phase1_standard_julia_traces")
out_dir_env = get(ENV, "PHASE1_JULIA_TRACE_OUT_DIR", out_dir_default)
out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(repo_root, out_dir_env)
mkpath(out_dir)

patients_path = joinpath(design_dir, "patients.csv")
reg_path = joinpath(design_dir, "regimen_events.csv")
if !isfile(patients_path) || !isfile(reg_path)
    error("Design files missing in $design_dir")
end

patients = DataFrame(CSV.File(patients_path))
reg_events = DataFrame(CSV.File(reg_path))
dlbcl_overrides_path = joinpath(design_dir, "dlbcl_param_overrides.csv")
dlbcl_overrides = isfile(dlbcl_overrides_path) ? DataFrame(CSV.File(dlbcl_overrides_path)) : DataFrame(name = String[], value = Float64[])

phase1_active_variant_ids = [5, 9, 14, 20, 24, 25, 27, 28]
variants_mode = lowercase(get(ENV, "PHASE1_VARIANTS_MODE", "matlab_empty"))
variant_ids =
    if variants_mode == "matlab_empty"
        # Match MATLAB script behavior in export_phase1_standard_traces_matlab.m:
        # sbiosimulate(model, cs, [], d) -> no variants applied.
        Int[]
    elseif variants_mode == "active_phase1"
        phase1_active_variant_ids
    else
        error("Unsupported PHASE1_VARIANTS_MODE=$variants_mode (use matlab_empty or active_phase1)")
    end
bw_kg_for_dose_conversion = 70.0
tgrid = sort(unique(vcat(collect(0.0:0.01:2.0), collect(2.0:0.1:84.0))))
alg = TCellEngagerQSP.make_solver_alg()
engine = Symbol(lowercase(get(ENV, "TCE_ENGINE", "canonical")))
if !(engine in (:legacy, :canonical))
    error("Unsupported TCE_ENGINE=$engine. Use legacy or canonical.")
end

manifest = DataFrame(
    regimen = String[],
    patient_id = Int[],
    trace_csv = String[],
    status = String[],
)

for reg in unique(reg_events.regimen)
    sub = reg_events[reg_events.regimen .== reg, :]
    sort!(sub, :event_idx)

    dose_map = Dict{Float64, Float64}()
    for rr in eachrow(sub)
        t = Float64(rr.time_day)
        amt = Float64(rr.dose_mg) * 1000.0 / bw_kg_for_dose_conversion
        dose_map[t] = get(dose_map, t, 0.0) + amt
    end
    dose_times = sort(collect(keys(dose_map)))

    for prow in eachrow(patients)
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
        pmap["end_time"] = 84.0
        pnames = sort(collect(keys(pmap)))
        pvals = [pmap[n] for n in pnames]

        patient_id = Int(prow.patient_id)
        safe_reg = replace(String(reg), r"[^A-Za-z0-9_]" => "_")
        out_csv = joinpath(out_dir, "regimen_$(safe_reg)_patient_$(patient_id).csv")

        try
            mdl = TCellEngagerQSP.build_model_with_variant_ids(variant_ids, pnames, pvals)
            rhs_fun = TCellEngagerQSP.rhs!
            ctx = nothing
            if engine == :legacy
                ctx = TCellEngagerQSP.make_legacy_context(mdl, TCellEngagerQSP.InfusionEvent[])
                rhs_fun = TCellEngagerQSP.rhs!
            else
                ctx = TCellEngagerQSP.CanonicalSimContext(
                    copy(mdl.z),
                    mdl.pvals,
                    TCellEngagerQSP.InfusionEvent[],
                )
                rhs_fun = mdl.canonical_rhs
            end
            target_idx = TCellEngagerQSP.canonical_state_to_idx(mdl)["TDBc_ugperkg"]

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

            u0 = copy(TCellEngagerQSP.canonical_u0(mdl))
            if dose_at_t0 != 0.0
                u0[target_idx] += dose_at_t0
            end

            cb = PresetTimeCallback(cb_times, affect!; save_positions = (true, true))
            prob = ODEProblem(rhs_fun, u0, (0.0, 84.0), ctx)
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

            il6 = Float64[]
            bt = Float64[]
            for t in tgrid
                u = Float64.(sol(t))
                push!(il6, TCellEngagerQSP.eval_symbol(mdl, u, t, "IL6combo"))
                push!(bt, TCellEngagerQSP.eval_symbol(mdl, u, t, "Btumor"))
            end

            out_df = DataFrame(
                time = tgrid,
                IL6combo = il6,
                Btumor = bt,
                regimen = fill(String(reg), length(tgrid)),
                patient_id = fill(patient_id, length(tgrid)),
            )
            CSV.write(out_csv, out_df)
            push!(manifest, (String(reg), patient_id, out_csv, "ok"))
        catch
            push!(manifest, (String(reg), patient_id, out_csv, "error"))
        end
    end
end

manifest_path = joinpath(out_dir, "manifest.csv")
CSV.write(manifest_path, manifest)
println("Wrote " * manifest_path)
