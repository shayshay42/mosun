using CSV
using DataFrames
using DifferentialEquations
using Sundials
using SciMLBase

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

repo_root = TCellEngagerQSP.REPO_ROOT
design_dir_default = joinpath(repo_root, "generated", "phase1_design_standard")
design_dir_env = get(ENV, "PHASE1_DESIGN_DIR", design_dir_default)
design_dir = isabspath(design_dir_env) ? design_dir_env : joinpath(repo_root, design_dir_env)

out_dir_default = joinpath(repo_root, "generated", "phase1_standard_julia_traces_canonical_segmented")
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
        Int[]
    elseif variants_mode == "active_phase1"
        phase1_active_variant_ids
    else
        error("Unsupported PHASE1_VARIANTS_MODE=$variants_mode (use matlab_empty or active_phase1)")
    end

bw_kg_for_dose_conversion = 70.0
tgrid = sort(unique(vcat(collect(0.0:0.01:2.0), collect(2.0:0.1:84.0))))
alg = TCellEngagerQSP.make_solver_alg()

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

function append_solution!(sol_t::Vector{Float64}, sol_u::Vector{Vector{Float64}}, ts::Vector{Float64}, us::Vector{Vector{Float64}})
    if isempty(ts)
        return
    end
    if isempty(sol_t)
        append!(sol_t, ts)
        append!(sol_u, us)
        return
    end
    start_idx = isapprox(ts[1], sol_t[end]; atol = 1e-12, rtol = 0.0) ? 2 : 1
    for i in start_idx:length(ts)
        push!(sol_t, ts[i])
        push!(sol_u, copy(us[i]))
    end
end

function run_segmented_trace(rhs_fun, ctx, u0::Vector{Float64}, target_idx::Int, dose_map::Dict{Float64,Float64}, tf::Float64, tgrid::Vector{Float64}, alg)
    u_curr = copy(u0)
    t_curr = 0.0
    sol_t = Float64[]
    sol_u = Vector{Vector{Float64}}()

    dose_at_t0 = get(dose_map, 0.0, 0.0)
    if dose_at_t0 != 0.0
        u_curr[target_idx] += dose_at_t0
    end

    event_times = sort(unique([t for t in keys(dose_map) if t > 0.0]))
    all_event_times = sort(unique(vcat(event_times, [tf])))

    for te in all_event_times
        te < t_curr && continue
        seg_save = [t for t in tgrid if t_curr <= t <= te]
        if te > t_curr
            prob = ODEProblem(rhs_fun, u_curr, (t_curr, te), ctx)
            sol_seg = solve(
                prob,
                alg;
                abstol = TCellEngagerQSP.SOLVER_ABSTOL,
                reltol = TCellEngagerQSP.SOLVER_RELTOL,
                saveat = seg_save,
                tstops = [te],
            )
            append_solution!(sol_t, sol_u, sol_seg.t, sol_seg.u)
            u_curr = copy(sol_seg.u[end])
        elseif !isempty(seg_save) && (isempty(sol_t) || !isapprox(sol_t[end], te; atol = 1e-12, rtol = 0.0))
            push!(sol_t, te)
            push!(sol_u, copy(u_curr))
        end

        if haskey(dose_map, te)
            u_curr[target_idx] += dose_map[te]
        end
        t_curr = te
    end

    lookup = Dict{Float64, Vector{Float64}}()
    for i in eachindex(sol_t)
        lookup[sol_t[i]] = sol_u[i]
    end
    aligned_u = [lookup[t] for t in tgrid]
    return aligned_u
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
            rhs_fun = mdl.canonical_rhs
            ctx = TCellEngagerQSP.CanonicalSimContext(
                copy(mdl.z),
                mdl.pvals,
                TCellEngagerQSP.InfusionEvent[],
            )
            target_idx = mdl.state_to_idx["TDBc_ugperkg"]
            u0 = copy(mdl.u0)

            u_aligned = run_segmented_trace(rhs_fun, ctx, u0, target_idx, dose_map, 84.0, tgrid, alg)

            il6 = Float64[]
            bt = Float64[]
            for (i, t) in enumerate(tgrid)
                u = Float64.(u_aligned[i])
                push!(il6, eval_symbol(mdl, u, t, "IL6combo"))
                push!(bt, eval_symbol(mdl, u, t, "Btumor"))
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

