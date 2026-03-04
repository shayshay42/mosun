using CSV
using DataFrames
using JSON3
using Random
using Statistics

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

    # Use a single random dose vector shared across all patients.
    rand_seed = parse(Int, get(ENV, "OPTVPOP_RANDOM_SEED", "20260304"))
    rng = MersenneTwister(rand_seed)

    orig_patient_id = get(ENV, "OPTMC_PATIENT_ID", "")
    ENV["OPTMC_PATIENT_ID"] = string(patient_ids[1])
    prob0, schedule_mode, n_cycles = build_patient_problem()
    nd = length(prob0.dose_times_days)

    random_override = strip(get(ENV, "OPTVPOP_RANDOM_DOSES_MG", ""))
    random_doses_mg =
        if isempty(random_override)
            [prob0.lower_mg[i] + rand(rng) * (prob0.upper_mg[i] - prob0.lower_mg[i]) for i in 1:nd]
        else
            vals = parse_float_list(random_override)
            length(vals) == nd || error("OPTVPOP_RANDOM_DOSES_MG length $(length(vals)) does not match dose events $nd")
            vals
        end

    rows = NamedTuple[]
    for (k, pid) in enumerate(patient_ids)
        ENV["OPTMC_PATIENT_ID"] = string(pid)
        prob, _, _ = build_patient_problem()
        length(prob.dose_times_days) == nd || error("Dose event count mismatch for patient $pid")
        spec, _ = build_tracking_spec(prob)

        scenarios = (
            ("fixed_min", copy(prob.lower_mg)),
            ("fixed_max", copy(prob.upper_mg)),
            ("random_shared", copy(random_doses_mg)),
        )

        for (scenario, doses) in scenarios
            traj = simulate_trajectory(prob, doses)
            terms = tracking_terms(prob, traj, doses, spec)
            push!(
                rows,
                (
                    patient_id = pid,
                    scenario = scenario,
                    loss = Float64(terms.loss),
                    tumor_mse = Float64(terms.tumor_mse),
                    tox_mse = Float64(terms.tox_mse),
                    dose_reg = Float64(terms.dose_reg),
                ),
            )
        end
        if k % 20 == 0 || k == n_patients
            println("processed $k / $n_patients patients")
        end
    end

    if isempty(orig_patient_id)
        delete!(ENV, "OPTMC_PATIENT_ID")
    else
        ENV["OPTMC_PATIENT_ID"] = orig_patient_id
    end

    out_dir_default = joinpath(TCellEngagerQSP.REPO_ROOT, "generated", "figures", "optimization", "vpop_tracking_loss_scenarios")
    out_dir_env = get(ENV, "OPTVPOP_OUT_DIR", out_dir_default)
    out_dir = resolve_repo_path(out_dir_env)
    mkpath(out_dir)

    df = DataFrame(rows)
    csv_path = joinpath(out_dir, "vpop_tracking_loss_by_scenario.csv")
    CSV.write(csv_path, df)

    summary = combine(groupby(df, :scenario), :loss => mean => :loss_mean, :loss => median => :loss_median, :loss => std => :loss_std)
    summary_path = joinpath(out_dir, "vpop_tracking_loss_summary.csv")
    CSV.write(summary_path, summary)

    meta = Dict(
        "n_patients" => n_patients,
        "patient_ids" => patient_ids,
        "schedule_mode" => schedule_mode,
        "n_cycles" => n_cycles,
        "dose_times_days" => prob0.dose_times_days,
        "random_seed" => rand_seed,
        "random_shared_doses_mg" => random_doses_mg,
        "loss_mode" => "trajectory_tracking_quadratic",
    )
    meta_path = joinpath(out_dir, "vpop_tracking_loss_meta.json")
    open(meta_path, "w") do io
        JSON3.pretty(io, meta)
    end

    println("Wrote $csv_path")
    println("Wrote $summary_path")
    println("Wrote $meta_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
