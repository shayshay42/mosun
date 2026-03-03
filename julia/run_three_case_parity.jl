using CSV
using DataFrames

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

repo_root = TCellEngagerQSP.REPO_ROOT
manifest_path = joinpath(repo_root, "generated", "matlab_reference", "manifest_three_cases.csv")
out_dir = joinpath(repo_root, "generated", "julia_reference")
mkpath(out_dir)
engine = Symbol(lowercase(get(ENV, "TCE_ENGINE", "canonical")))
if !(engine in (:legacy, :canonical))
    error("Unsupported TCE_ENGINE=$engine. Use legacy or canonical.")
end

manifest = DataFrame(CSV.File(manifest_path))
summary_rows = DataFrame(
    sim_key = String[],
    case_no = Int[],
    id = Int[],
    max_abs_err = Float64[],
    max_rel_err = Float64[]
)

for row in eachrow(manifest)
    times, outvec, Xj, ref_df = run_manifest_row(row; engine = engine)

    Xref = hcat([Float64.(ref_df[!, Symbol(name)]) for name in outvec]...)
    abs_err = abs.(Xj .- Xref)
    rel_err = abs_err ./ max.(abs.(Xref), 1e-12)

    max_abs = maximum(abs_err)
    max_rel = maximum(rel_err)

    sim_key = String(row.sim_key)
    case_no = Int(row.case_no)
    id = Int(row.id)

    out_df = DataFrame(time = times)
    for (j, name) in enumerate(outvec)
        out_df[!, Symbol(name)] = Xj[:, j]
    end
    CSV.write(joinpath(out_dir, sim_key * ".csv"), out_df)

    push!(summary_rows, (sim_key, case_no, id, max_abs, max_rel))

    println("$sim_key => max_abs_err=$(max_abs), max_rel_err=$(max_rel)")
end

summary_path = joinpath(out_dir, "parity_summary_" * String(engine) * ".csv")
CSV.write(summary_path, summary_rows)
println("Wrote parity summary: " * summary_path)
