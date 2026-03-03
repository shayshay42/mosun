using CSV
using DataFrames
using Test

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

function max_norm_abs_err(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real})
    @assert size(X) == size(Y)
    mask = isfinite.(X) .& isfinite.(Y)
    if !any(mask)
        return Inf
    end
    d = abs.(X .- Y)
    amp = maximum(abs.(Y[mask]))
    amp = max(amp, 1e-12)
    return maximum(d[mask]) / amp
end

repo_root = TCellEngagerQSP.REPO_ROOT
manifest_path = joinpath(repo_root, "generated", "matlab_reference", "manifest_three_cases.csv")
manifest = DataFrame(CSV.File(manifest_path))

@testset "Canonical RHS Not Worse Than Legacy" begin
    @test nrow(manifest) > 0

    summary = DataFrame(
        sim_key = String[],
        legacy_norm_err = Float64[],
        canonical_norm_err = Float64[],
        canon_minus_legacy = Float64[],
        canon_vs_legacy_norm = Float64[],
    )

    for row in eachrow(manifest)
        _, outvec, X_legacy, ref_df = run_manifest_row(row; engine = :legacy)
        _, _, X_canon, _ = run_manifest_row(row; engine = :canonical)
        X_ref = hcat([Float64.(ref_df[!, Symbol(name)]) for name in outvec]...)

        legacy_norm = max_norm_abs_err(X_legacy, X_ref)
        canon_norm = max_norm_abs_err(X_canon, X_ref)
        canon_vs_legacy = max_norm_abs_err(X_canon, X_legacy)
        margin = canon_norm - legacy_norm

        push!(summary, (String(row.sim_key), legacy_norm, canon_norm, margin, canon_vs_legacy))

        # Canonical path should not degrade reference fit beyond tiny numerical jitter.
        @test canon_norm <= legacy_norm + 1e-6
        # Canonical and legacy implementations should be numerically near-identical.
        @test canon_vs_legacy <= 1e-6
    end

    out_dir = joinpath(repo_root, "generated", "julia_reference")
    mkpath(out_dir)
    out_path = joinpath(out_dir, "canonical_vs_legacy_quality.csv")
    CSV.write(out_path, summary)
    @info "Wrote canonical-vs-legacy summary" out_path
end

