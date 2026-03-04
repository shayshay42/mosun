using CSV
using DataFrames
using JSON3
using Printf

include(joinpath(@__DIR__, "optimize_patient_multicycle_tracking.jl"))

function resolve_repo_path(rel_or_abs::String)
    return isabspath(rel_or_abs) ? rel_or_abs : joinpath(TCellEngagerQSP.REPO_ROOT, rel_or_abs)
end

function parse_patient_ids(summary::DataFrame)
    return Int.(summary.patient_id)
end

function extract_dose_vector(row, nd::Int)
    x = zeros(Float64, nd)
    for i in 1:nd
        x[i] = Float64(row[Symbol("dose$(i)_mg")])
    end
    return x
end

function main()
    summary_csv_env = get(ENV, "OPTVPOP_SUMMARY_CSV", joinpath("generated", "figures", "optimization", "vpop_tracking_nmcd", "vpop_tracking_nmcd_summary.csv"))
    summary_csv = resolve_repo_path(summary_csv_env)
    summary = DataFrame(CSV.File(summary_csv))
    nrow(summary) > 0 || error("No rows in $summary_csv")

    orig_patient_id = get(ENV, "OPTMC_PATIENT_ID", "")
    first_pid = Int(summary.patient_id[1])
    ENV["OPTMC_PATIENT_ID"] = string(first_pid)
    prob0, _, _ = build_patient_problem()
    nd = length(prob0.dose_times_days)

    all_rows = NamedTuple[]
    for (k, srow) in enumerate(eachrow(summary))
        pid = Int(srow.patient_id)
        ENV["OPTMC_PATIENT_ID"] = string(pid)
        prob, _, _ = build_patient_problem()
        length(prob.dose_times_days) == nd || error("Dose event count mismatch for patient $pid")

        x_fixed = clamp.(copy(prob.fixed_doses_mg), prob.lower_mg, prob.upper_mg)
        x_opt = extract_dose_vector(srow, nd)

        for (scenario, doses) in (("fixed", x_fixed), ("optimized", x_opt))
            traj = simulate_trajectory(prob, doses)
            for j in eachindex(traj.t)
                push!(
                    all_rows,
                    (
                        patient_id = pid,
                        scenario = scenario,
                        time_day = Float64(traj.t[j]),
                        Btumor = Float64(traj.btumor[j]),
                        IL6combo = Float64(traj.il6[j]),
                    ),
                )
            end
        end
        if k % 10 == 0 || k == nrow(summary)
            println(@sprintf("trajectory export processed %d / %d patients", k, nrow(summary)))
        end
    end

    if isempty(orig_patient_id)
        delete!(ENV, "OPTMC_PATIENT_ID")
    else
        ENV["OPTMC_PATIENT_ID"] = orig_patient_id
    end

    out_dir_default = dirname(summary_csv)
    out_dir_env = get(ENV, "OPTVPOP_OUT_DIR", out_dir_default)
    out_dir = resolve_repo_path(out_dir_env)
    mkpath(out_dir)

    traj_df = DataFrame(all_rows)
    out_csv = joinpath(out_dir, "vpop_tracking_nmcd_trajectories.csv")
    CSV.write(out_csv, traj_df)

    meta = Dict(
        "source_summary_csv" => summary_csv,
        "n_patients" => nrow(summary),
        "dose_times_days" => prob0.dose_times_days,
        "scenarios" => ["fixed", "optimized"],
    )
    meta_path = joinpath(out_dir, "vpop_tracking_nmcd_trajectories_meta.json")
    open(meta_path, "w") do io
        JSON3.pretty(io, meta)
    end

    println("Wrote $out_csv")
    println("Wrote $meta_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
