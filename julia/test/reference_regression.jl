using CSV
using DataFrames
using Statistics
using Test

const MATLAB_REFERENCE_REL_TOL = Dict(
    "case10_id2001" => 0.0016,
    "case20_id1001" => 0.0006,
    "case30_id1" => 0.00025,
)

function max_abs_err(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real})
    @assert size(X) == size(Y)
    mask = isfinite.(X) .& isfinite.(Y)
    any(mask) || return Inf
    return maximum(abs.(X[mask] .- Y[mask]))
end

function max_norm_abs_err(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real})
    @assert size(X) == size(Y)
    mask = isfinite.(X) .& isfinite.(Y)
    any(mask) || return Inf
    amp = max(maximum(abs.(Y[mask])), 1e-12)
    return maximum(abs.(X[mask] .- Y[mask])) / amp
end

function rms_norm_abs_err(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real})
    @assert size(X) == size(Y)
    mask = isfinite.(X) .& isfinite.(Y)
    any(mask) || return Inf
    amp = max(maximum(abs.(Y[mask])), 1e-12)
    return sqrt(mean(abs2.(X[mask] .- Y[mask]))) / amp
end

repo_root = TCellEngagerQSP.REPO_ROOT
manifest_path = joinpath(repo_root, "generated", "matlab_reference", "manifest_three_cases.csv")
manifest = DataFrame(CSV.File(manifest_path))

@testset "Manual Core Structure" begin
    artifacts = TCellEngagerQSP.load_artifacts()
    @test length(artifacts.rate_exprs) == 116
    @test length(artifacts.initial_rules) + length(artifacts.repeated_rules) == 101
    @test TCellEngagerQSP.direct_model_layout_matches(artifacts.state_names, artifacts.param_names)
end

@testset "Manual Core Matches MATLAB Reference" begin
    @test nrow(manifest) == length(MATLAB_REFERENCE_REL_TOL)

    summary = DataFrame(
        sim_key = String[],
        max_abs_err = Float64[],
        max_rel_err = Float64[],
        rms_rel_err = Float64[],
        legacy_rel_err = Float64[],
    )

    for row in eachrow(manifest)
        times_legacy, outvec, X_legacy, ref_df = TCellEngagerQSP.run_manifest_row(row; engine = :legacy)
        times_manual, _, X_manual, _ = TCellEngagerQSP.run_manifest_row(row; engine = :canonical, core_mode_override = :legacy_reference)
        X_ref = hcat([Float64.(ref_df[!, Symbol(name)]) for name in outvec]...)

        @test Float64.(times_legacy) == Float64.(ref_df.time)
        @test Float64.(times_manual) == Float64.(ref_df.time)

        max_abs = max_abs_err(X_manual, X_ref)
        max_rel = max_norm_abs_err(X_manual, X_ref)
        rms_rel = rms_norm_abs_err(X_manual, X_ref)
        legacy_rel = max_norm_abs_err(X_manual, X_legacy)
        tol = MATLAB_REFERENCE_REL_TOL[String(row.sim_key)]

        @test legacy_rel <= max(1e-6, tol)
        @test max_rel <= tol
        @test rms_rel <= tol / 2

        push!(summary, (String(row.sim_key), max_abs, max_rel, rms_rel, legacy_rel))
    end

    out_dir = joinpath(repo_root, "generated", "julia_reference")
    mkpath(out_dir)
    out_path = joinpath(out_dir, "manual_core_reference_quality.csv")
    CSV.write(out_path, summary)
    @info "Wrote manual-core reference summary" out_path
end
