using CSV
using DataFrames
using DifferentialEquations
using SciMLBase
using SciMLSensitivity
using ForwardDiff
using ReverseDiff
using Zygote
using Printf
using JSON3
using Statistics
using LinearAlgebra

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

repo_root = TCellEngagerQSP.REPO_ROOT
design_dir_default = joinpath(repo_root, "generated", "phase1_design")
design_dir_env = get(ENV, "PHASE1_DESIGN_DIR", design_dir_default)
design_dir = isabspath(design_dir_env) ? design_dir_env : joinpath(repo_root, design_dir_env)
patients_path = joinpath(design_dir, "patients.csv")
overrides_path = joinpath(design_dir, "dlbcl_param_overrides.csv")

if !isfile(patients_path) || !isfile(overrides_path)
    error("Design files missing. Run scripts/generate_phase1_design.py first.")
end

patients = DataFrame(CSV.File(patients_path))
overrides = DataFrame(CSV.File(overrides_path))
patient_id = parse(Int, get(ENV, "GRAD_PATIENT_ID", "1"))
prow_idx = findfirst(patients.patient_id .== patient_id)
if prow_idx === nothing
    error("GRAD_PATIENT_ID=$patient_id not found in $(patients_path)")
end
prow = patients[prow_idx, :]

function parse_dose_vector(s::String)
    toks = split(strip(s), ",")
    out = Float64[]
    for tok in toks
        st = strip(tok)
        isempty(st) && continue
        push!(out, parse(Float64, st))
    end
    return out
end

function build_patient_param_map(prow_row, overrides_df::DataFrame)
    pmap = Dict{String, Float64}()
    for rr in eachrow(overrides_df)
        pmap[String(rr.name)] = Float64(rr.value)
    end

    pmap["Bpbo_perml"] = Float64(prow_row.Bpbo_perml)
    pmap["Bpbref_perml"] = Float64(prow_row.Bpbref_perml)
    pmap["Trpbo_perml"] = Float64(prow_row.Trpbo_perml)
    pmap["Trpbref_perml"] = Float64(prow_row.Trpbref_perml)
    pmap["KBptumor"] = Float64(prow_row.KBptumor)
    pmap["KTrptumor"] = Float64(prow_row.KTrptumor)
    pmap["kBtumorprolif"] = Float64(prow_row.kBtumorprolif)

    # Phase-1 setup used throughout this repo.
    pmap["PKflag"] = 1.0
    pmap["fvalidation"] = 0.0
    pmap["VPid"] = 1.0
    pmap["end_time"] = 21.0
    return pmap
end

function build_ad_rhs_no_cast(mdl)
    name_to_idx = mdl.name_to_idx
    repeated_rule_defs = mdl.repeated_rule_defs
    rate_exprs = mdl.rate_exprs
    stoich = mdl.stoich

    repeated_assign_exprs = Any[]
    for (lhs_idx, rhs0) in repeated_rule_defs
        rhs =
            if occursin("PK_v26(", rhs0)
                # In this phase-1 setup PKflag is fixed to 1, so PK_v26 simplifies.
                "ifelse(TDBc_ugperkg / Vc_tdb > 1e-5, TDBc_ugperkg / Vc_tdb, 0.0)"
            else
                rhs0
            end
        rhs_expr = TCellEngagerQSP.compile_formula_expr(rhs, name_to_idx)
        push!(repeated_assign_exprs, :(z[$lhs_idx] = real($rhs_expr)))
    end

    reaction_blocks = Any[]
    for j in eachindex(rate_exprs)
        rvar = gensym(:rate)
        rexpr = TCellEngagerQSP.compile_formula_expr(rate_exprs[j], name_to_idx)
        block = Expr(:block, :($rvar = real($rexpr)))
        for (i, coeff) in stoich[j]
            push!(block.args, :(du[$i] += $(Float64(coeff)) * $rvar))
        end
        push!(reaction_blocks, block)
    end

    np = length(mdl.pvals)
    zlen = length(mdl.name_to_idx)
    pvals = copy(mdl.pvals)

    fexpr = quote
        let pvals0 = $pvals
            (du, u, p, t) -> begin
                ns = length(u)
                T = eltype(u)
                z = Vector{T}(undef, $zlen)
                @inbounds begin
                    z[1:ns] .= u
                    z[ns+1:ns+$np] .= T.(pvals0)
                    $(repeated_assign_exprs...)

                    fill!(du, zero(T))
                    $(reaction_blocks...)
                end
                return nothing
            end
        end
    end

    return Base.invokelatest(Core.eval, TCellEngagerQSP, fexpr)
end

function build_symbol_observer_no_cast(mdl, symbol_name::String)
    name_to_idx = mdl.name_to_idx
    repeated_rule_defs = mdl.repeated_rule_defs
    zlen = length(name_to_idx)
    target_idx = name_to_idx[symbol_name]
    np = length(mdl.pvals)
    pvals = copy(mdl.pvals)

    repeated_assign_exprs = Any[]
    for (lhs_idx, rhs0) in repeated_rule_defs
        rhs =
            if occursin("PK_v26(", rhs0)
                "ifelse(TDBc_ugperkg / Vc_tdb > 1e-5, TDBc_ugperkg / Vc_tdb, 0.0)"
            else
                rhs0
            end
        rhs_expr = TCellEngagerQSP.compile_formula_expr(rhs, name_to_idx)
        push!(repeated_assign_exprs, :(z[$lhs_idx] = real($rhs_expr)))
    end

    fexpr = quote
        let pvals0 = $pvals
            (u, t) -> begin
                ns = length(u)
                T = eltype(u)
                z = Vector{T}(undef, $zlen)
                @inbounds begin
                    z[1:ns] .= u
                    z[ns+1:ns+$np] .= T.(pvals0)
                    $(repeated_assign_exprs...)
                end
                return z[$target_idx]
            end
        end
    end
    return Base.invokelatest(Core.eval, TCellEngagerQSP, fexpr)
end

function smoothmax(v::AbstractVector, tau::Float64)
    m = maximum(v)
    return m + tau * log(sum(exp.((v .- m) ./ tau)))
end

pmap = build_patient_param_map(prow, overrides)
pnames = sort(collect(keys(pmap)))
pvals = [pmap[n] for n in pnames]
variant_ids = Int[]
mdl = TCellEngagerQSP.build_model_with_variant_ids(variant_ids, pnames, pvals)
rhs_ad = build_ad_rhs_no_cast(mdl)
il6_obs_fun = build_symbol_observer_no_cast(mdl, "IL6combo")

bt_idx = mdl.state_to_idx["Btumor"]
tdbc_idx = mdl.state_to_idx["TDBc_ugperkg"]

dose_times = [0.0, 7.0, 14.0]
horizon_days = 21.0
tox_eval_days = collect(0.0:0.1:2.0)
saveat_seg1 = vcat(tox_eval_days, [7.0])
prob_p_const = [1.0, 1.0, 1.0]

mg_to_ugkg = 1000.0 / 70.0
loss_weight_tox = 0.5
loss_weight_tum = 0.5
tox_tau = parse(Float64, get(ENV, "GRAD_TOX_SOFTMAX_TAU", "50.0"))
abstol = parse(Float64, get(ENV, "GRAD_ABSTOL", "1e-8"))
reltol = parse(Float64, get(ENV, "GRAD_RELTOL", "1e-6"))
maxiters = parse(Int, get(ENV, "GRAD_MAXITERS", "100000000"))

function solve_segment(rhs_fun, u0, p, tspan::Tuple{Float64,Float64}, saveat::Vector{Float64}, alg;
        sensealg = nothing, abstol::Float64 = 1e-8, reltol::Float64 = 1e-6, maxiters::Int = 100000000)
    prob = ODEProblem(rhs_fun, u0, tspan, p)
    if isnothing(sensealg)
        sol = solve(
            prob,
            alg;
            abstol = abstol,
            reltol = reltol,
            saveat = saveat,
            tstops = [tspan[2]],
            maxiters = maxiters,
        )
    else
        sol = solve(
            prob,
            alg;
            sensealg = sensealg,
            abstol = abstol,
            reltol = reltol,
            saveat = saveat,
            tstops = [tspan[2]],
            maxiters = maxiters,
        )
    end
    if sol.retcode != SciMLBase.ReturnCode.Success
        error("segment solve failed with retcode=$(sol.retcode) over tspan=$tspan")
    end
    return sol
end

function cycle_metrics(dose_mg::AbstractVector{T}; sensealg = nothing, alg = Tsit5()) where {T}
    length(dose_mg) == 3 || error("Expected exactly 3 doses (mg), got $(length(dose_mg))")
    d_ugkg = dose_mg .* T(mg_to_ugkg)
    e_tdbc = zeros(T, length(mdl.u0))
    e_tdbc[tdbc_idx] = one(T)

    u = T.(mdl.u0) .+ d_ugkg[1] .* e_tdbc
    bt0 = u[bt_idx]

    sol1 = solve_segment(
        rhs_ad,
        u,
        prob_p_const,
        (dose_times[1], dose_times[2]),
        saveat_seg1,
        alg;
        sensealg = sensealg,
        abstol = abstol,
        reltol = reltol,
        maxiters = maxiters,
    )
    n_tox = length(tox_eval_days)
    il6_vals = [Base.invokelatest(il6_obs_fun, sol1.u[i], sol1.t[i]) for i in 1:n_tox]
    tox_proxy = smoothmax(il6_vals, tox_tau)

    u2 = sol1.u[end] .+ d_ugkg[2] .* e_tdbc
    sol2 = solve_segment(
        rhs_ad,
        u2,
        prob_p_const,
        (dose_times[2], dose_times[3]),
        [dose_times[3]],
        alg;
        sensealg = sensealg,
        abstol = abstol,
        reltol = reltol,
        maxiters = maxiters,
    )

    u3 = sol2.u[end] .+ d_ugkg[3] .* e_tdbc
    sol3 = solve_segment(
        rhs_ad,
        u3,
        prob_p_const,
        (dose_times[3], horizon_days),
        [horizon_days],
        alg;
        sensealg = sensealg,
        abstol = abstol,
        reltol = reltol,
        maxiters = maxiters,
    )

    bt_end = sol3.u[end][bt_idx]
    tumor_proxy = bt_end / (bt0 + T(1e-12))
    return (; tox_proxy = tox_proxy, tumor_proxy = tumor_proxy)
end

ref_dose_mg = parse_dose_vector(get(ENV, "GRAD_REF_DOSES_MG", "0.8,2.0,6.0"))
length(ref_dose_mg) == 3 || error("GRAD_REF_DOSES_MG must contain exactly 3 comma-separated values")
solver_otd_name = lowercase(get(ENV, "GRAD_SOLVER_OTD", "rodas4p"))
solver_dto_name = lowercase(get(ENV, "GRAD_SOLVER_DTO", "tsit5"))
solver_scale_name = lowercase(get(ENV, "GRAD_SOLVER_SCALE", solver_dto_name))
alg_otd = TCellEngagerQSP.make_solver_alg(solver_otd_name)
alg_dto = TCellEngagerQSP.make_solver_alg(solver_dto_name)
alg_scale = TCellEngagerQSP.make_solver_alg(solver_scale_name)

ref_metrics = cycle_metrics(ref_dose_mg; sensealg = nothing, alg = alg_scale)
tox_scale = max(abs(Float64(ref_metrics.tox_proxy)), 1e-12)
tum_scale = max(abs(Float64(ref_metrics.tumor_proxy)), 1e-12)

function balanced_loss(dose_mg::AbstractVector; sensealg = nothing, alg = Tsit5())
    m = cycle_metrics(dose_mg; sensealg = sensealg, alg = alg)
    tox_norm = m.tox_proxy / tox_scale
    tum_norm = m.tumor_proxy / tum_scale
    return loss_weight_tox * tox_norm + loss_weight_tum * tum_norm
end

function finite_difference_gradient(f, x::Vector{Float64}; rel_step::Float64 = 1e-3, abs_step::Float64 = 1e-5)
    g = zeros(Float64, length(x))
    for i in eachindex(x)
        h = max(abs_step, rel_step * max(abs(x[i]), 1.0))
        xp = copy(x)
        xm = copy(x)
        xp[i] += h
        xm[i] -= h
        g[i] = (f(xp) - f(xm)) / (2h)
    end
    return g
end

function grad_report(name::String, grad_fun, g_fd::Vector{Float64})
    t0 = time()
    try
        g = grad_fun()
        gvec = Float64.(g)
        abs_err = abs.(gvec .- g_fd)
        rel_err = abs_err ./ (abs.(g_fd) .+ 1e-12)
        cos_sim = dot(gvec, g_fd) / (norm(gvec) * norm(g_fd) + 1e-12)
        return Dict(
            "method" => name,
            "status" => "ok",
            "runtime_s" => time() - t0,
            "grad" => gvec,
            "max_abs_err_vs_fd" => maximum(abs_err),
            "max_rel_err_vs_fd" => maximum(rel_err),
            "median_rel_err_vs_fd" => median(rel_err),
            "cosine_similarity_vs_fd" => cos_sim,
            "error" => "",
        )
    catch err
        return Dict(
            "method" => name,
            "status" => "error",
            "runtime_s" => time() - t0,
            "grad" => [NaN, NaN, NaN],
            "max_abs_err_vs_fd" => NaN,
            "max_rel_err_vs_fd" => NaN,
            "median_rel_err_vs_fd" => NaN,
            "cosine_similarity_vs_fd" => NaN,
            "error" => sprint(showerror, err),
        )
    end
end

otd_forward = ForwardSensitivity()
otd_reverse = ReverseDiffAdjoint()
dto_passthrough = SensitivityADPassThrough()

reports = Dict{String, Any}[]

function run_fd_row!(reports_arr::Vector{Dict{String, Any}}, label::String, fd_obj)
    t0 = time()
    g_fd = finite_difference_gradient(fd_obj, copy(ref_dose_mg))
    runtime = time() - t0
    push!(
        reports_arr,
        Dict(
            "method" => label,
            "status" => "ok",
            "runtime_s" => runtime,
            "grad" => g_fd,
            "max_abs_err_vs_fd" => 0.0,
            "max_rel_err_vs_fd" => 0.0,
            "median_rel_err_vs_fd" => 0.0,
            "cosine_similarity_vs_fd" => 1.0,
            "error" => "",
        ),
    )
    return g_fd
end

# OtD forward
fd_otd_f = run_fd_row!(reports, "fd_otd_forward", d -> balanced_loss(d; sensealg = nothing, alg = alg_otd))
push!(
    reports,
    grad_report(
        "otd_forward",
        () -> ForwardDiff.gradient(d -> balanced_loss(d; sensealg = otd_forward, alg = alg_otd), ref_dose_mg),
        fd_otd_f,
    ),
)

# OtD reverse
fd_otd_r = run_fd_row!(reports, "fd_otd_reverse", d -> balanced_loss(d; sensealg = nothing, alg = alg_otd))
push!(
    reports,
    grad_report(
        "otd_reverse",
        () -> first(Zygote.gradient(d -> balanced_loss(d; sensealg = otd_reverse, alg = alg_otd), ref_dose_mg)),
        fd_otd_r,
    ),
)

# DtO forward
fd_dto_f = run_fd_row!(reports, "fd_dto_forward", d -> balanced_loss(d; sensealg = nothing, alg = alg_dto))
push!(
    reports,
    grad_report(
        "dto_forward",
        () -> ForwardDiff.gradient(d -> balanced_loss(d; sensealg = dto_passthrough, alg = alg_dto), ref_dose_mg),
        fd_dto_f,
    ),
)

# DtO reverse
fd_dto_r = run_fd_row!(reports, "fd_dto_reverse", d -> balanced_loss(d; sensealg = nothing, alg = alg_dto))
push!(
    reports,
    grad_report(
        "dto_reverse",
        () -> first(Zygote.gradient(d -> balanced_loss(d; sensealg = dto_passthrough, alg = alg_dto), ref_dose_mg)),
        fd_dto_r,
    ),
)

out_dir_default = joinpath(repo_root, "generated", "gradient_checks")
out_dir_env = get(ENV, "GRAD_OUT_DIR", out_dir_default)
out_dir = isabspath(out_dir_env) ? out_dir_env : joinpath(repo_root, out_dir_env)
mkpath(out_dir)

table_rows = DataFrame(
    method = String[],
    status = String[],
    runtime_s = Float64[],
    grad_dose1 = Float64[],
    grad_dose2 = Float64[],
    grad_dose3 = Float64[],
    max_abs_err_vs_fd = Float64[],
    max_rel_err_vs_fd = Float64[],
    median_rel_err_vs_fd = Float64[],
    cosine_similarity_vs_fd = Float64[],
    error = String[],
)

for r in reports
    g = r["grad"]
    push!(
        table_rows,
        (
            String(r["method"]),
            String(r["status"]),
            Float64(r["runtime_s"]),
            Float64(g[1]),
            Float64(g[2]),
            Float64(g[3]),
            Float64(r["max_abs_err_vs_fd"]),
            Float64(r["max_rel_err_vs_fd"]),
            Float64(r["median_rel_err_vs_fd"]),
            Float64(r["cosine_similarity_vs_fd"]),
            String(r["error"]),
        ),
    )
end

csv_path = joinpath(out_dir, "dose3_cycle_gradient_check.csv")
CSV.write(csv_path, table_rows)

meta = Dict(
    "patient_id" => patient_id,
    "solver_otd" => solver_otd_name,
    "solver_dto" => solver_dto_name,
    "solver_scale" => solver_scale_name,
    "doses_mg" => ref_dose_mg,
    "dose_times_days" => dose_times,
    "horizon_days" => horizon_days,
    "tox_proxy" => "smoothmax(IL6combo, day 0-2, tau=$(tox_tau))",
    "tumor_proxy" => "Btumor(day21)/Btumor(day0)",
    "loss" => "$(loss_weight_tox)*tox/tox_scale + $(loss_weight_tum)*tumor/tum_scale",
    "tox_scale" => tox_scale,
    "tumor_scale" => tum_scale,
    "ref_metrics" => Dict(
        "tox_proxy" => Float64(ref_metrics.tox_proxy),
        "tumor_proxy" => Float64(ref_metrics.tumor_proxy),
    ),
    "abstol" => abstol,
    "reltol" => reltol,
)
json_path = joinpath(out_dir, "dose3_cycle_gradient_check_meta.json")

sanitize_json(x::Float64) = isfinite(x) ? x : nothing
sanitize_json(x::AbstractVector) = [sanitize_json(v) for v in x]
sanitize_json(x::Dict) = Dict(String(k) => sanitize_json(v) for (k, v) in x)
sanitize_json(x) = x

open(json_path, "w") do io
    JSON3.pretty(io, sanitize_json(Dict("meta" => meta, "results" => reports)))
end

println("Reference doses (mg): ", ref_dose_mg)
println("Reference tox proxy: ", Float64(ref_metrics.tox_proxy), " (scale=", tox_scale, ")")
println("Reference tumor proxy: ", Float64(ref_metrics.tumor_proxy), " (scale=", tum_scale, ")")
println()
for r in reports
    g = r["grad"]
    @printf(
        "%-22s status=%-5s grad=[%.6e, %.6e, %.6e] max_rel_err_vs_fd=%.3e cos=%.6f runtime=%.2fs\n",
        String(r["method"]),
        String(r["status"]),
        Float64(g[1]),
        Float64(g[2]),
        Float64(g[3]),
        Float64(r["max_rel_err_vs_fd"]),
        Float64(r["cosine_similarity_vs_fd"]),
        Float64(r["runtime_s"]),
    )
    if String(r["status"]) != "ok"
        println("  error: ", String(r["error"]))
    end
end
println()
println("Wrote ", csv_path)
println("Wrote ", json_path)
