using CSV
using DataFrames
using JSON3
using Optim
using Printf

include(joinpath(@__DIR__, "optimize_patient_multicycle_tracking.jl"))

function resolve_repo_path(rel_or_abs::String)
    return isabspath(rel_or_abs) ? rel_or_abs : joinpath(TCellEngagerQSP.REPO_ROOT, rel_or_abs)
end

function main()
    patients_csv_env = get(ENV, "OPTMC_PATIENTS_CSV", joinpath("assets", "generated_vpop", "selected_patients.csv"))
    patients_csv = resolve_repo_path(patients_csv_env)
    patients = DataFrame(CSV.File(patients_csv))
    nrow(patients) > 0 || error("No rows in $patients_csv")
    sort!(patients, :patient_id)

    max_patients = parse(Int, get(ENV, "OPTVPOP_MAX_PATIENTS", string(nrow(patients))))
    if max_patients < nrow(patients)
        patients = patients[1:max_patients, :]
    end
    patient_ids = Int.(patients.patient_id)
    n_patients = length(patient_ids)

    nm_iters = parse(Int, get(ENV, "OPTVPOP_NM_ITERS", "20"))
    nm_f_calls = parse(Int, get(ENV, "OPTVPOP_NM_F_CALLS", "250"))
    cd_sweeps = parse(Int, get(ENV, "OPTVPOP_CD_SWEEPS", "10"))
    cd_step0 = parse(Float64, get(ENV, "OPTVPOP_CD_STEP0_MG", "2.0"))
    cd_tol = parse(Float64, get(ENV, "OPTVPOP_CD_TOL_MG", "0.05"))
    progress_every = parse(Int, get(ENV, "OPTVPOP_PROGRESS_EVERY", "10"))
    checkpoint_every = parse(Int, get(ENV, "OPTVPOP_CHECKPOINT_EVERY", "5"))

    out_dir_default = joinpath(TCellEngagerQSP.REPO_ROOT, "generated", "figures", "optimization", "vpop_tracking_nmcd")
    out_dir_env = get(ENV, "OPTVPOP_OUT_DIR", out_dir_default)
    out_dir = resolve_repo_path(out_dir_env)
    mkpath(out_dir)
    summary_path = joinpath(out_dir, "vpop_tracking_nmcd_summary.csv")
    meta_path = joinpath(out_dir, "vpop_tracking_nmcd_meta.json")
    resume = lowercase(get(ENV, "OPTVPOP_RESUME", "1")) in ("1", "true", "yes", "y")

    orig_patient_id = get(ENV, "OPTMC_PATIENT_ID", "")
    ENV["OPTMC_PATIENT_ID"] = string(patient_ids[1])
    prob0, schedule_mode, n_cycles = build_patient_problem()
    nd = length(prob0.dose_times_days)

    summary_df =
        if resume && isfile(summary_path)
            DataFrame(CSV.File(summary_path))
        else
            DataFrame()
        end
    done_ids = Set{Int}()
    if nrow(summary_df) > 0 && ("patient_id" in names(summary_df))
        # Keep only the latest row per patient when resuming from an existing summary.
        seen = Set{Int}()
        keep = falses(nrow(summary_df))
        for i in nrow(summary_df):-1:1
            pid = Int(summary_df.patient_id[i])
            if !(pid in seen)
                keep[i] = true
                push!(seen, pid)
            end
        end
        summary_df = summary_df[keep, :]
        for pid in Int.(summary_df.patient_id)
            push!(done_ids, pid)
        end
    end
    dose_cols = Symbol[]
    for i in 1:nd
        push!(dose_cols, Symbol("dose$(i)_mg"))
    end

    println(
        @sprintf(
            "tracking VPop optimize (NM+CD): n_patients=%d, done=%d, cycles=%d, dose_controls=%d, nm_iters=%d, cd_sweeps=%d",
            n_patients,
            length(done_ids),
            n_cycles,
            nd,
            nm_iters,
            cd_sweeps,
        ),
    )

    new_done = 0
    for pid in patient_ids
        if pid in done_ids
            continue
        end
        ENV["OPTMC_PATIENT_ID"] = string(pid)
        prob, _, _ = build_patient_problem()
        length(prob.dose_times_days) == nd || error("Dose event count mismatch for patient $pid")
        spec, _ = build_tracking_spec(prob)
        obj = make_tracking_objective(prob, spec; sensealg = nothing)

        x0 = clamp.(copy(prob.fixed_doses_mg), prob.lower_mg, prob.upper_mg)
        x_opt = copy(x0)
        status = "ok"
        t0 = time()
        try
            nm_opts = Optim.Options(
                iterations = nm_iters,
                f_calls_limit = nm_f_calls > 0 ? nm_f_calls : typemax(Int),
                show_trace = false,
                store_trace = false,
            )
            res_nm = optimize(obj, prob.lower_mg, prob.upper_mg, copy(x0), Fminbox(NelderMead()), nm_opts)
            x_curr = clamp.(Optim.minimizer(res_nm), prob.lower_mg, prob.upper_mg)
            loss_curr = Float64(obj(x_curr))
            steps = fill(cd_step0, nd)
            for _ in 1:cd_sweeps
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
        catch err
            status = "fail: " * sprint(showerror, err)
        end
        runtime_s = time() - t0

        traj_fixed = simulate_trajectory(prob, x0)
        term_fixed = tracking_terms(prob, traj_fixed, x0, spec)
        traj_opt = simulate_trajectory(prob, x_opt)
        term_opt = tracking_terms(prob, traj_opt, x_opt, spec)

        dose_nt = NamedTuple{Tuple(dose_cols)}(Tuple(Float64.(x_opt)))
        row = merge(
            (
                patient_id = pid,
                status = status,
                runtime_s = runtime_s,
                loss_fixed = Float64(term_fixed.loss),
                loss_opt = Float64(term_opt.loss),
                tumor_mse_fixed = Float64(term_fixed.tumor_mse),
                tumor_mse_opt = Float64(term_opt.tumor_mse),
                tox_mse_fixed = Float64(term_fixed.tox_mse),
                tox_mse_opt = Float64(term_opt.tox_mse),
                dose_reg_fixed = Float64(term_fixed.dose_reg),
                dose_reg_opt = Float64(term_opt.dose_reg),
            ),
            dose_nt,
        )
        push!(summary_df, row, cols = :union)
        push!(done_ids, pid)
        new_done += 1

        if checkpoint_every > 0 && (new_done % checkpoint_every == 0)
            CSV.write(summary_path, summary_df)
        end
        done_count = length(done_ids)
        if done_count % progress_every == 0 || done_count == n_patients
            println(@sprintf("processed %d / %d patients", done_count, n_patients))
        end
    end

    if isempty(orig_patient_id)
        delete!(ENV, "OPTMC_PATIENT_ID")
    else
        ENV["OPTMC_PATIENT_ID"] = orig_patient_id
    end

    CSV.write(summary_path, summary_df)

    meta = Dict(
        "n_patients" => n_patients,
        "patient_ids" => patient_ids,
        "schedule_mode" => schedule_mode,
        "n_cycles" => n_cycles,
        "dose_times_days" => prob0.dose_times_days,
        "method" => "nelder_mead_coordinate_descent",
        "nm_iters" => nm_iters,
        "nm_f_calls" => nm_f_calls,
        "cd_sweeps" => cd_sweeps,
        "cd_step0_mg" => cd_step0,
        "cd_tol_mg" => cd_tol,
        "loss_mode" => "trajectory_tracking_quadratic",
        "loss_improvement_mean" => mean(summary_df.loss_fixed .- summary_df.loss_opt),
        "loss_improvement_median" => median(summary_df.loss_fixed .- summary_df.loss_opt),
    )
    open(meta_path, "w") do io
        JSON3.pretty(io, meta)
    end

    println("Wrote $summary_path")
    println("Wrote $meta_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
