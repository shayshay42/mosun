using CSV
using DataFrames
using DifferentialEquations
using DiffEqCallbacks
using SciMLBase
using Sundials
using Base.Threads

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

const REPO_ROOT = TCellEngagerQSP.REPO_ROOT
const MMC = TCellEngagerQSP.MosunModelCore

function parse_float_list(txt::String)
    vals = Float64[]
    for tok in split(txt, ",")
        s = strip(tok)
        if isempty(s)
            continue
        end
        push!(vals, parse(Float64, s))
    end
    return vals
end

function load_base_param_map()
    p = joinpath(REPO_ROOT, "generated", "phase1_design_standard", "dlbcl_param_overrides.csv")
    if !isfile(p)
        error("Missing DLBCL overrides at $p")
    end
    df = DataFrame(CSV.File(p))
    pmap = Dict{String, Float64}()
    for r in eachrow(df)
        pmap[String(r.name)] = Float64(r.value)
    end
    pmap["Bpbo_perml"] = get(pmap, "Bpbo_perml", get(pmap, "Bpbref_perml", 500_000.0))
    pmap["Trpbo_perml"] = get(pmap, "Trpbo_perml", get(pmap, "Trpbref_perml", 500_000.0))
    pmap["Bpbref_perml"] = get(pmap, "Bpbref_perml", pmap["Bpbo_perml"])
    pmap["Trpbref_perml"] = get(pmap, "Trpbref_perml", pmap["Trpbo_perml"])
    return pmap
end

function eval_best_spd_for_pmap(
    pmap::Dict{String, Float64},
    variant_ids::Vector{Int},
    dose_days::Vector{Float64},
    dose_mg::Vector{Float64},
    horizon_days::Float64,
    save_dt::Float64,
    bw_kg::Float64,
    alg,
    use_mtk_jac::Bool,
)
    pmap["PKflag"] = 1.0
    pmap["fvalidation"] = 0.0
    pmap["VPid"] = 1.0
    pmap["end_time"] = horizon_days

    dose_map = Dict{Float64, Float64}()
    for (t, d) in zip(dose_days, dose_mg)
        amt = d * 1000.0 / bw_kg
        dose_map[t] = get(dose_map, t, 0.0) + amt
    end
    saveat = collect(0.0:save_dt:horizon_days)

    sol =
        if isempty(variant_ids) && !use_mtk_jac
            params = MMC.params_from_dict(pmap; ignore_unknown = true)
            regimen = MMC.bolus_regimen(:TDBc_ugperkg, dose_map)
            _, direct_sol = MMC.solve_regimen(
                regimen,
                params,
                alg;
                tspan = (0.0, horizon_days),
                saveat = saveat,
                callback_mode = :callback,
                abstol = TCellEngagerQSP.SOLVER_ABSTOL,
                reltol = TCellEngagerQSP.SOLVER_RELTOL,
            )
            direct_sol
        else
            pnames = sort(collect(keys(pmap)))
            pvals = [pmap[n] for n in pnames]
            mdl = TCellEngagerQSP.build_model_with_variant_ids(variant_ids, pnames, pvals)

            ctx = TCellEngagerQSP.CanonicalSimContext(
                copy(mdl.z),
                mdl.pvals,
                TCellEngagerQSP.InfusionEvent[],
            )
            rhs_fun = mdl.canonical_rhs
            jac_fun = use_mtk_jac ? TCellEngagerQSP.make_mtk_dense_jacobian(mdl) : nothing
            ode_rhs = isnothing(jac_fun) ? rhs_fun : ODEFunction(rhs_fun; jac = jac_fun)
            state_to_idx = TCellEngagerQSP.canonical_state_to_idx(mdl)
            target_idx = state_to_idx["TDBc_ugperkg"]

            local_dose_map = copy(dose_map)
            dose_at_t0 = get(local_dose_map, 0.0, 0.0)
            local_dose_map[0.0] = 0.0
            cb_times = sort([t for t in keys(local_dose_map) if t > 0.0])

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
                if !isnan(amt) && amt != 0.0
                    integrator.u[target_idx] += amt
                end
            end

            u0 = copy(TCellEngagerQSP.canonical_u0(mdl))
            if dose_at_t0 != 0.0
                u0[target_idx] += dose_at_t0
            end

            cb = PresetTimeCallback(cb_times, affect!; save_positions = (true, true))
            prob = ODEProblem(ode_rhs, u0, (0.0, horizon_days), ctx)
            solve(
                prob,
                alg;
                abstol = TCellEngagerQSP.SOLVER_ABSTOL,
                reltol = TCellEngagerQSP.SOLVER_RELTOL,
                callback = cb,
                tstops = cb_times,
                d_discontinuities = cb_times,
                saveat = saveat,
            )
        end

    if sol.retcode != SciMLBase.ReturnCode.Success
        return NaN
    end
    bt_idx = MMC.dynamic_state_index(:Btumor)
    bt = [Float64(u[bt_idx]) for u in sol.u]
    bt0 = max(bt[1], 1e-12)
    spd_pct = @. 100.0 * (bt / bt0 - 1.0)
    spd_eval = length(spd_pct) > 1 ? spd_pct[2:end] : spd_pct
    return minimum(spd_eval)
end

function main()
    patients_csv_env = get(
        ENV,
        "BESTSPD_PATIENTS_CSV",
        joinpath("generated", "vpop_pruning", "vpop_sens_top5_20260306_w02", "selected_patients.csv"),
    )
    out_csv_env = get(
        ENV,
        "BESTSPD_OUT_CSV",
        joinpath("generated", "figures", "optimization", "vpop_selected_best_spd.csv"),
    )
    patients_csv = isabspath(patients_csv_env) ? patients_csv_env : joinpath(REPO_ROOT, patients_csv_env)
    out_csv = isabspath(out_csv_env) ? out_csv_env : joinpath(REPO_ROOT, out_csv_env)
    mkpath(dirname(out_csv))

    dose_days = parse_float_list(get(ENV, "BESTSPD_DOSE_DAYS", "0,7,14,21,42,63,84,105,126,147"))
    dose_mg = parse_float_list(get(ENV, "BESTSPD_DOSE_MG", "1,2,60,60,30,30,30,30,30,30"))
    if length(dose_days) != length(dose_mg)
        error("BESTSPD_DOSE_DAYS and BESTSPD_DOSE_MG lengths differ.")
    end
    horizon_days = parse(Float64, get(ENV, "BESTSPD_HORIZON_DAYS", "168"))
    save_dt = parse(Float64, get(ENV, "BESTSPD_SAVE_DT", "1.0"))
    bw_kg = parse(Float64, get(ENV, "BESTSPD_BW_KG", "70.0"))
    fail_fill = parse(Float64, get(ENV, "BESTSPD_FAIL_FILL", "100.0"))

    solver_name_full = lowercase(get(ENV, "TCE_SOLVER", "qndf"))
    solver_name = replace(solver_name_full, "_mtkjac" => "")
    use_mtk_jac = endswith(solver_name_full, "_mtkjac") || get(ENV, "TCE_MTK_JAC", "0") == "1"
    alg = TCellEngagerQSP.make_solver_alg(solver_name)

    phase1_active_variant_ids = [5, 9, 14, 20, 24, 25, 27, 28]
    variants_mode = lowercase(get(ENV, "BESTSPD_VARIANTS_MODE", "matlab_empty"))
    variant_ids = if variants_mode == "matlab_empty"
        Int[]
    elseif variants_mode == "active_phase1"
        phase1_active_variant_ids
    else
        error("Unsupported BESTSPD_VARIANTS_MODE=$variants_mode")
    end

    patients = DataFrame(CSV.File(patients_csv))
    if !("patient_id" in names(patients))
        error("Expected patient_id column in $patients_csv")
    end

    base_pmap = load_base_param_map()
    base_keys = Set(keys(base_pmap))

    n = nrow(patients)
    best_spd = fill(NaN, n)
    status = fill("not_run", n)
    patient_ids = Vector{Int}(undef, n)
    progress = Atomic{Int}(0)
    progress_lock = ReentrantLock()

    println("Evaluating best %SPD for $n patients")
    println("  solver=$(solver_name_full), variants_mode=$variants_mode, threads=$(Threads.nthreads())")
    println("  regimen_days=$(join(dose_days, ','))")
    println("  regimen_mg=$(join(dose_mg, ','))")
    println("  horizon_days=$horizon_days")

    Threads.@threads for i in 1:n
        row = patients[i, :]
        pid = Int(round(Float64(row[Symbol("patient_id")])))
        patient_ids[i] = pid

        pmap = copy(base_pmap)
        for c in names(patients)
            if c == "patient_id"
                continue
            end
            cname = String(c)
            if !(cname in base_keys)
                continue
            end
            v = row[c]
            if ismissing(v)
                continue
            end
            fv = try
                Float64(v)
            catch
                NaN
            end
            if isfinite(fv)
                pmap[cname] = fv
            end
        end

        y = try
            eval_best_spd_for_pmap(
                pmap,
                variant_ids,
                dose_days,
                dose_mg,
                horizon_days,
                save_dt,
                bw_kg,
                alg,
                use_mtk_jac,
            )
        catch err
            @warn "Patient evaluation failed" idx = i patient_id = pid err
            NaN
        end

        if isfinite(y)
            best_spd[i] = y
            status[i] = "ok"
        else
            best_spd[i] = fail_fill
            status[i] = "failed_fill"
        end

        done = atomic_add!(progress, 1) + 1
        if done % 20 == 0 || done == n
            lock(progress_lock) do
                println("  progress $done / $n")
            end
        end
    end

    out = DataFrame(
        patient_id = patient_ids,
        best_spd_pct = best_spd,
        status = status,
    )
    sort!(out, :patient_id)
    CSV.write(out_csv, out)
    n_fail = count(x -> x != "ok", out.status)
    println("Wrote $out_csv")
    println("n_fail=$n_fail")
end

main()
