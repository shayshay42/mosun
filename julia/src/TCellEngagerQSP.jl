module TCellEngagerQSP

using CSV
using DataFrames
using JSON3
using SciMLBase
using DifferentialEquations
using Sundials
using ModelingToolkit
using RuntimeGeneratedFunctions
using SparseArrays

RuntimeGeneratedFunctions.init(@__MODULE__)

include(joinpath(@__DIR__, "MosunModelCore.jl"))
using .MosunModelCore

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MODEL_TABLES_DIR = joinpath(REPO_ROOT, "generated", "model_tables")
const RECON_BUNDLE_DIR = joinpath(REPO_ROOT, "generated", "model_reconstruction_bundle")

const SAFE_FUNCS = Set([:+, :-, :*, :/, :max, :min, :log, :exp, :sqrt, :abs, :ifelse])

struct RuleExpr
    lhs_idx::Int
    fn::Function
    lhs_name::String
end

struct DoseDef
    name::String
    target::String
    amount::Vector{Float64}
    rate::Vector{Float64}
    times::Vector{Float64}
end

struct InfusionEvent
    target_idx::Int
    t_start::Float64
    t_end::Float64
    rate::Float64
end

mutable struct SimContext
    z::Vector{Float64}
    pvals::Vector{Float64}
    repeated_rules::Vector{RuleExpr}
    rate_fns::Vector{Function}
    stoich::Vector{Vector{Pair{Int,Float64}}}
    infusions::Vector{InfusionEvent}
end

mutable struct CanonicalSimContext
    z::Vector{Float64}
    pvals::Vector{Float64}
    infusions::Vector{InfusionEvent}
    cache::MosunModelCore.MosunObservablesCache
end

CanonicalSimContext(z::Vector{Float64}, pvals::Vector{Float64}, infusions::Vector{InfusionEvent}) =
    CanonicalSimContext(z, pvals, infusions, MosunModelCore.zero_observables_cache())

struct ModelArtifacts
    state_names::Vector{String}
    state_u0::Vector{Float64}
    param_names::Vector{String}
    param_defaults::Vector{Float64}
    initial_rules::Vector{Tuple{String,String}}
    repeated_rules::Vector{Tuple{String,String}}
    rate_exprs::Vector{String}
    reaction_strs::Vector{String}
    variants::Dict{Int, Vector{Pair{String,Float64}}}
    stoich::Vector{Vector{Pair{Int,Float64}}}
    doses::Dict{String,DoseDef}
end

struct PKDataset
    animal::Vector{Int}
    time::Vector{Float64}
    pk::Vector{Float64}
end

const PK_DATA = Ref{PKDataset}()
const PK_VALIDATION_DATA = Ref{PKDataset}()
const PK_VALIDATION2_DATA = Ref{PKDataset}()
const SOLVER_ABSTOL = parse(Float64, get(ENV, "TCE_ABSTOL", "1e-8"))
const SOLVER_RELTOL = parse(Float64, get(ENV, "TCE_RELTOL", "1e-5"))
const MTK_JAC_CACHE = Dict{UInt64, Function}()
const MTK_SPARSE_PROTOTYPE_CACHE = Dict{UInt64, SparseMatrixCSC{Float64,Int}}()
const DEFAULT_CORE_MODE = let core_mode = lowercase(strip(get(ENV, "TCE_CORE_MODE", "production")))
    if core_mode == "production"
        :production
    elseif core_mode == "legacy_reference"
        :legacy_reference
    elseif core_mode == "auto"
        :auto
    else
        error("Unsupported TCE_CORE_MODE=$core_mode. Use production, legacy_reference, or auto.")
    end
end

function make_solver_alg(solver_name::AbstractString)
    solver = lowercase(String(solver_name))
    if solver == "cvode_bdf"
        return CVODE_BDF()
    elseif solver == "qndf"
        return QNDF(autodiff = false)
    elseif solver == "rodas4p"
        return Rodas4P(autodiff = false)
    elseif solver == "kencarp4"
        return KenCarp4(autodiff = false)
    elseif solver == "vcabm"
        return VCABM()
    elseif solver == "tsit5"
        return Tsit5()
    else
        error("Unsupported solver=$solver. Use cvode_bdf, qndf, rodas4p, kencarp4, vcabm, or tsit5.")
    end
end

function make_solver_alg()
    solver = lowercase(get(ENV, "TCE_SOLVER", "cvode_bdf"))
    return make_solver_alg(solver)
end

function sanitize_name(name::AbstractString)
    return replace(String(name), "unnamed." => "")
end

function sanitize_formula(expr::AbstractString)
    return replace(String(expr), "unnamed." => "")
end

function parse_value_or_nan(s::AbstractString)
    st = strip(s)
    if isempty(st)
        return NaN
    elseif lowercase(st) == "nan"
        return NaN
    else
        return parse(Float64, st)
    end
end

function parse_pk_source(path::String)
    txt = read(path, String)

    function extract_vec(name::String)
        m = match(Regex("data\\.$name\\s*=\\s*\\[(.*?)\\]'\\s*(?:/\\s*([0-9.eE+\\-]+))?\\s*;", "s"), txt)
        m === nothing && error("Could not find data.$name in $path")
        toks = split(strip(m.captures[1]))
        vals = [parse_value_or_nan(tok) for tok in toks]
        if length(m.captures) >= 2 && !(m.captures[2] === nothing)
            divv = parse(Float64, m.captures[2])
            vals = vals ./ divv
        end
        return vals
    end

    animal = Int.(round.(extract_vec("AnimalID")))
    time = extract_vec("Time")
    pk = extract_vec("PK")
    return PKDataset(animal, time, pk)
end

function load_pk_datasets!()
    pk_path = joinpath(REPO_ROOT, "assets", "Supp Matlab Code", "PK.m")
    pkv_path = joinpath(REPO_ROOT, "assets", "Supp Matlab Code", "PK_validation.m")
    pkv2_path = joinpath(REPO_ROOT, "assets", "Supp Matlab Code", "PK_validation_2.m")

    PK_DATA[] = parse_pk_source(pk_path)
    PK_VALIDATION_DATA[] = parse_pk_source(pkv_path)
    PK_VALIDATION2_DATA[] = parse_pk_source(pkv2_path)
end

function unique_sorted_with_first(T::Vector{Float64}, X::Vector{Float64})
    Tu = sort(unique(T))
    Xu = similar(Tu)
    for i in eachindex(Tu)
        idx = findfirst(==(Tu[i]), T)
        Xu[i] = X[idx]
    end
    return Tu, Xu
end

function interp1_linear(x::Vector{Float64}, y::Vector{Float64}, t::Float64)
    if t <= x[1]
        return y[1]
    elseif t >= x[end]
        return y[end]
    end
    i = searchsortedlast(x, t)
    if i >= length(x)
        return y[end]
    end
    x1, x2 = x[i], x[i + 1]
    y1, y2 = y[i], y[i + 1]
    if x2 == x1
        return y1
    end
    return y1 + (y2 - y1) * (t - x1) / (x2 - x1)
end

function pk_core(data::PKDataset, vpid::Int, ticker::Float64, end_time::Float64;
        floor_log::Float64 = -6.0, add_missing_20p8::Bool = false, clip_ticker::Bool = false)

    if clip_ticker && ticker > end_time
        return 1e-6
    end

    idx = findall(==(vpid), data.animal)
    Tvec = data.time[idx]
    PKvec = data.pk[idx]

    keep = .!isnan.(PKvec)
    Tvec = Tvec[keep]
    PKvec = PKvec[keep]

    TvecR, PKvecR = unique_sorted_with_first(Tvec, PKvec)

    PKvecR = max.(1e-6, PKvecR)
    TvecR = vcat(-1.0, TvecR)
    PKvecR = vcat(PKvecR[1], PKvecR)
    PKvecR = log10.(PKvecR)

    if add_missing_20p8 && (vpid == 2501 || vpid == 2502 || vpid == 3001)
        ind1 = findfirst(==(8.0), TvecR)
        ind2 = findfirst(==(13.99), TvecR)
        ind3 = findfirst(==(15.0), TvecR)
        if !(ind1 === nothing || ind2 === nothing || ind3 === nothing)
            m_tmp = (PKvecR[ind2] - PKvecR[ind1]) / (13.99 - 8.0)
            PK_tmp = PKvecR[ind3] + m_tmp * (20.8 - 15.0)
            ind4 = findall(x -> x <= 15.0, TvecR)
            ind5 = findall(x -> x > 15.0, TvecR)
            TvecR = vcat(TvecR[ind4], [20.8], TvecR[ind5])
            PKvecR = vcat(PKvecR[ind4], [PK_tmp], PKvecR[ind5])
        end
    end

    m = (PKvecR[end] - PKvecR[end - 1]) / (TvecR[end] - TvecR[end - 1])
    m = min(m, -0.02)
    Tprev = TvecR[end]
    PKprev = PKvecR[end]
    TvecR = vcat(TvecR, [end_time])
    PKvecR = vcat(PKvecR, [max(floor_log, PKprev + m * (end_time - Tprev))])

    return 10.0 ^ interp1_linear(TvecR, PKvecR, ticker)
end

function PK(TDBc_ugperkg::Float64, Vc::Float64, PKflag::Float64, VPid::Float64, ticker::Float64, end_time::Float64)
    if PKflag == 1.0
        val = TDBc_ugperkg / Vc
        return (val > 1e-5 ? 1.0 : 0.0) * val
    end
    return pk_core(PK_DATA[], Int(round(VPid)), ticker, end_time; floor_log = -6.0, add_missing_20p8 = true, clip_ticker = true)
end

function PK_validation(TDBc_ugperkg::Float64, Vc::Float64, PKflag::Float64, VPid::Float64, ticker::Float64, end_time::Float64)
    if PKflag == 1.0
        val = TDBc_ugperkg / Vc
        return (val > 1e-5 ? 1.0 : 0.0) * val
    end
    return pk_core(PK_VALIDATION_DATA[], Int(round(VPid)), ticker, end_time; floor_log = -9.0, add_missing_20p8 = false, clip_ticker = false)
end

function PK_validation_2(TDBc_ugperkg::Float64, Vc::Float64, PKflag::Float64, VPid::Float64, ticker::Float64, end_time::Float64)
    if PKflag == 1.0
        val = TDBc_ugperkg / Vc
        return (val > 1e-5 ? 1.0 : 0.0) * val
    end

    data = PK_VALIDATION2_DATA[]
    idx = findall(==(Int(round(VPid))), data.animal)
    Tvec = data.time[idx]
    PKvec = data.pk[idx]

    ord = sortperm(Tvec)
    Tvec = Tvec[ord]
    PKvec = PKvec[ord]

    nan_idx = findall(isnan, PKvec)
    if length(nan_idx) > 1
        PKvec[nan_idx[2]] = 0.0015
    end

    keep = .!isnan.(PKvec)
    Tvec = Tvec[keep]
    PKvec = PKvec[keep]

    TvecR, PKvecR = unique_sorted_with_first(Tvec, PKvec)

    PKvecR = max.(1e-6, PKvecR)
    TvecR = vcat(-1.0, TvecR)
    PKvecR = vcat(PKvecR[1], PKvecR)
    PKvecR = log10.(PKvecR)

    m = (PKvecR[end] - PKvecR[end - 1]) / (TvecR[end] - TvecR[end - 1])
    m = min(m, -0.02)
    Tprev = TvecR[end]
    PKprev = PKvecR[end]
    TvecR = vcat(TvecR, [end_time])
    PKvecR = vcat(PKvecR, [max(-9.0, PKprev + m * (end_time - Tprev))])

    return 10.0 ^ interp1_linear(TvecR, PKvecR, ticker)
end

function PK_mean(TDBc_ugperkg::Float64, Vc::Float64, PKflag::Float64, case_no::Float64, ticker::Float64, end_time::Float64)
    if PKflag == 1.0
        val = TDBc_ugperkg / Vc
        return (val > 1e-5 ? 1.0 : 0.0) * val
    elseif PKflag == 2.0
        return PK(10.0, 1.0, 0.0, case_no, ticker, end_time)
    elseif PKflag != 0.0
        return 0.0
    end

    c = Int(round(case_no))
    vids = Int[]
    if c == 103 || c == 113
        vids = [2001, 2002, 2501, 2502]
    elseif c == 104 || c == 114
        vids = [3001, 3002, 3501, 3502]
    elseif c == 105 || c == 115
        vids = [4001, 4002, 4003, 4004, 4501, 4502, 4503, 20013, 20023, 20033, 20043]
    elseif c == 106
        vids = [10012, 10022, 10032, 10015, 10025, 10035, 20015, 20025, 20035, 30015, 30025, 30035]
    end
    isempty(vids) && return 0.0

    pk_mul = 1.0
    cnt = 0
    for vid in vids
        pk_mul *= PK(10.0, 1.0, 0.0, vid, ticker, end_time)
        cnt += 1
    end
    return pk_mul ^ (1.0 / cnt)
end

function PK_v26(TDBc_ugperkg::Float64, Vc::Float64, PKflag::Float64, VPid::Float64, ticker::Float64, end_time::Float64, fvalidation::Float64)
    fv = Int(round(fvalidation))
    if fv == 0
        return PK(TDBc_ugperkg, Vc, PKflag, VPid, ticker, end_time)
    elseif fv == 1
        return PK_validation(TDBc_ugperkg, Vc, PKflag, VPid, ticker, end_time)
    elseif fv == 2
        return PK_validation_2(TDBc_ugperkg, Vc, PKflag, VPid, ticker, end_time)
    elseif fv == 5
        return PK_mean(TDBc_ugperkg, Vc, 0.0, VPid, ticker, end_time)
    else
        return 0.0
    end
end

# Generic numeric fallback for AD types (e.g., Dual numbers) used when PKflag is fixed to 1.
# If PKflag is not 1, we evaluate the existing Float64 implementation.
function PK_v26(TDBc_ugperkg::Real, Vc::Real, PKflag::Real, VPid::Real, ticker::Real, end_time::Real, fvalidation::Real)
    if PKflag == one(PKflag)
        val = TDBc_ugperkg / Vc
        return ifelse(val > 1e-5, val, zero(val))
    end
    return PK_v26(
        Float64(TDBc_ugperkg),
        Float64(Vc),
        Float64(PKflag),
        Float64(VPid),
        Float64(ticker),
        Float64(end_time),
        Float64(fvalidation),
    )
end

# Symbolics/MTK tracing path for Jacobian generation (used only for symbolic probe RHS).
# For phase-1 settings PKflag=1 and this branch matches the on-path behavior.
function PK_v26(TDBc_ugperkg::Num, Vc::Num, PKflag::Num, VPid::Num, ticker::Num, end_time::Num, fvalidation::Num)
    val = TDBc_ugperkg / Vc
    return ifelse(val > 1e-5, val, zero(val))
end

pow_safe(a::Float64, b::Float64) = Float64(real((complex(a) ^ b)))
pow_safe(a::Float64, b::Integer) = Float64(real((complex(a) ^ b)))
pow_safe(a::Real, b::Real) = real((complex(a) ^ b))
pow_safe(a::Num, b::Num) = a ^ b
pow_safe(a::Num, b::Real) = a ^ b
pow_safe(a::Num, b) = a ^ b
pow_safe(a, b) = a ^ b

function normalize_function_expr(ex::Expr)
    Base.remove_linenums!(ex)
    while ex.head == :block && length(ex.args) == 1 && ex.args[1] isa Expr
        ex = ex.args[1]
    end
    return ex
end

const MANUAL_CORE_PATH = joinpath(@__DIR__, "TCellEngagerQSPManualCore.jl")
if isfile(MANUAL_CORE_PATH)
    include(MANUAL_CORE_PATH)
else
    error("Missing manual model core at $MANUAL_CORE_PATH.")
end

function build_core_params(param_names::Vector{String}, pvals::Vector{Float64})
    p = MosunModelCore.default_params()
    for (nm, pv) in zip(param_names, pvals)
        sym = Symbol(nm)
        if hasproperty(p, sym)
            setproperty!(p, sym, Float64(pv))
        end
    end
    return p
end

function canonical_core_mode(param_to_idx::Dict{String,Int}, pvals::Vector{Float64})
    pkflag_idx = get(param_to_idx, "PKflag", 0)
    if pkflag_idx == 0
        return :legacy_reference
    end
    return isapprox(pvals[pkflag_idx], 1.0; atol = 0.0, rtol = 0.0) ? :production : :legacy_reference
end

core_state_names() = [String(sym) for sym in MosunModelCore.DYNAMIC_STATE_NAMES]
core_state_to_idx() = Dict(name => i for (i, name) in enumerate(core_state_names()))
core_observable_names() = Set(String(sym) for sym in MosunModelCore.OBSERVABLE_NAMES)
canonical_u0(mdl) = get(mdl, :core_mode, :legacy_reference) == :production ? mdl.core_u0 : mdl.u0
canonical_state_to_idx(mdl) = get(mdl, :core_mode, :legacy_reference) == :production ? mdl.core_state_to_idx : mdl.state_to_idx
uses_production_layout(mdl, u::AbstractVector) =
    get(mdl, :core_mode, :legacy_reference) == :production && haskey(mdl, :core_u0) && length(u) == length(mdl.core_u0)

function rewrite_expr(ex, name_to_idx::Dict{String,Int})
    if ex isa Symbol
        s = String(ex)
        if s == "time"
            return :t
        elseif haskey(name_to_idx, s)
            return :(z[$(name_to_idx[s])])
        else
            return ex
        end
    elseif ex isa Expr
        if ex.head == :call
            f = ex.args[1]
            if f == :^
                a1 = rewrite_expr(ex.args[2], name_to_idx)
                a2 = rewrite_expr(ex.args[3], name_to_idx)
                return Expr(:call, :pow_safe, a1, a2)
            else
                new_args = Any[f]
                for a in ex.args[2:end]
                    push!(new_args, rewrite_expr(a, name_to_idx))
                end
                return Expr(:call, new_args...)
            end
        else
            return Expr(ex.head, map(a -> rewrite_expr(a, name_to_idx), ex.args)...)
        end
    else
        return ex
    end
end

function compile_formula(formula::String, name_to_idx::Dict{String,Int})
    expr = Meta.parse(sanitize_formula(formula))
    rexpr = rewrite_expr(expr, name_to_idx)
    fexpr = quote
        (z, t) -> begin
            v = $rexpr
            Float64(real(v))
        end
    end
    return RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, normalize_function_expr(fexpr))
end

function populate_value_buffer!(z::Vector{Float64}, mdl, u::Vector{Float64}, t::Float64)
    uses_production_layout(mdl, u) &&
        throw(ArgumentError("populate_value_buffer! is only available for the legacy/reference layout"))
    apply_repeated_rules_direct!(z, u, mdl.pvals, t)
    return z
end

function eval_symbol(mdl, u::Vector{Float64}, t::Float64, name::String)
    if uses_production_layout(mdl, u)
        cache = MosunModelCore.zero_observables_cache()
        MosunModelCore.update_observables!(cache, u, mdl.core_params, t)
        return MosunModelCore.state_or_observable(u, cache, Symbol(sanitize_name(name)))
    end
    z = copy(mdl.z)
    populate_value_buffer!(z, mdl, u, t)
    return z[mdl.name_to_idx[name]]
end

function compile_legacy_runtime(mdl)
    repeated_rule_exprs = RuleExpr[]
    for (lhs_idx, rhs) in mdl.repeated_rule_defs
        push!(repeated_rule_exprs, RuleExpr(lhs_idx, compile_formula(rhs, mdl.name_to_idx), ""))
    end
    rate_fns = [compile_formula(expr, mdl.name_to_idx) for expr in mdl.rate_exprs]
    return repeated_rule_exprs, rate_fns
end

function make_legacy_context(mdl, infusions::Vector{InfusionEvent})
    repeated_rule_exprs, rate_fns = compile_legacy_runtime(mdl)
    return SimContext(copy(mdl.z), mdl.pvals, repeated_rule_exprs, rate_fns, mdl.stoich, infusions)
end

function compile_formula_expr(formula::String, name_to_idx::Dict{String,Int})
    expr = Meta.parse(sanitize_formula(formula))
    return rewrite_expr(expr, name_to_idx)
end

function build_symbolic_probe_rhs_function(
        name_to_idx::Dict{String,Int},
        repeated_rule_defs::Vector{Tuple{Int,String}},
        rate_exprs::Vector{String},
        stoich::Vector{Vector{Pair{Int,Float64}}},
        ns::Int,
        np::Int)

    repeated_assign_exprs = Any[]
    for (lhs_idx, rhs) in repeated_rule_defs
        rhs_expr = compile_formula_expr(rhs, name_to_idx)
        push!(repeated_assign_exprs, :(z[$lhs_idx] = real($rhs_expr)))
    end

    reaction_blocks = Any[]
    for j in eachindex(rate_exprs)
        rvar = gensym(:rate)
        rexpr = compile_formula_expr(rate_exprs[j], name_to_idx)
        block = Expr(:block, :($rvar = real($rexpr)))
        for (i, coeff) in stoich[j]
            push!(block.args, :(du[$i] += $(Float64(coeff)) * $rvar))
        end
        push!(reaction_blocks, block)
    end

    naux = length(name_to_idx) - ns - np
    zlen = ns + np + max(naux, 0)
    fexpr = quote
        (du, u, p, t) -> begin
            z = Vector{eltype(u)}(undef, $zlen)
            @inbounds begin
                z[1:$ns] .= u
                z[$(ns + 1):$(ns + np)] .= p
                $(repeated_assign_exprs...)
                fill!(du, zero(eltype(u)))
                $(reaction_blocks...)
            end
            return nothing
        end
    end
    return RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, normalize_function_expr(fexpr))
end

function mtk_jac_cache_key(mdl; sparse::Bool = false)
    return hash((length(mdl.u0), length(mdl.pvals), mdl.repeated_rule_defs, mdl.rate_exprs, mdl.stoich, sparse))
end

function make_mtk_dense_jacobian(mdl)
    get(mdl, :core_mode, :legacy_reference) == :production &&
        throw(ArgumentError("MTK dense Jacobians are only available on the legacy/reference model path"))
    key = mtk_jac_cache_key(mdl; sparse = false)
    if haskey(MTK_JAC_CACHE, key)
        return MTK_JAC_CACHE[key]
    end

    sym_rhs = build_symbolic_probe_rhs_function(
        mdl.name_to_idx,
        mdl.repeated_rule_defs,
        mdl.rate_exprs,
        mdl.stoich,
        length(mdl.u0),
        length(mdl.pvals),
    )
    probe_prob = ODEProblem(sym_rhs, copy(mdl.u0), (0.0, 1.0), copy(mdl.pvals))
    sys = ModelingToolkit.modelingtoolkitize(probe_prob)
    jac_exprs = ModelingToolkit.generate_jacobian(sys; sparse = false, simplify = false)
    (_, jac_expr_inplace) = jac_exprs
    jac_vec! = RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, normalize_function_expr(jac_expr_inplace))

    n = length(mdl.u0)
    jac_work = zeros(Float64, n * n)

    function jac!(J, u, p, t)
        pvals =
            if p isa CanonicalSimContext
                p.pvals
            elseif p isa SimContext
                p.pvals
            elseif p isa AbstractVector
                p
            else
                throw(ArgumentError("Unsupported parameter container type $(typeof(p)) for MTK Jacobian"))
            end
        jac_vec!(jac_work, u, pvals, t)
        if J isa AbstractMatrix
            @inbounds for j in 1:n, i in 1:n
                v = jac_work[(j - 1) * n + i]
                J[i, j] = isfinite(v) ? v : 0.0
            end
        else
            @inbounds for i in eachindex(jac_work)
                v = jac_work[i]
                J[i] = isfinite(v) ? v : 0.0
            end
        end
        return nothing
    end

    MTK_JAC_CACHE[key] = jac!
    return jac!
end

function make_mtk_sparse_jacobian(mdl)
    get(mdl, :core_mode, :legacy_reference) == :production &&
        throw(ArgumentError("MTK sparse Jacobians are only available on the legacy/reference model path"))
    key = mtk_jac_cache_key(mdl; sparse = true)
    if haskey(MTK_JAC_CACHE, key)
        return MTK_JAC_CACHE[key]
    end

    sym_rhs = build_symbolic_probe_rhs_function(
        mdl.name_to_idx,
        mdl.repeated_rule_defs,
        mdl.rate_exprs,
        mdl.stoich,
        length(mdl.u0),
        length(mdl.pvals),
    )
    probe_prob = ODEProblem(sym_rhs, copy(mdl.u0), (0.0, 1.0), copy(mdl.pvals))
    sys = ModelingToolkit.modelingtoolkitize(probe_prob)
    jac_exprs = ModelingToolkit.generate_jacobian(sys; sparse = true, simplify = false)
    (_, jac_expr_inplace) = jac_exprs
    jac_sparse! = RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, normalize_function_expr(jac_expr_inplace))
    n = length(mdl.u0)
    proto = copy(make_mtk_sparse_jacobian_prototype(mdl))

    function sanitize_sparse_jac!(J::SparseMatrixCSC)
        @inbounds for i in eachindex(J.nzval)
            v = J.nzval[i]
            J.nzval[i] = isfinite(v) ? v : 0.0
        end
        return J
    end

    function jac!(J, u, p, t)
        pvals =
            if p isa CanonicalSimContext
                p.pvals
            elseif p isa SimContext
                p.pvals
            elseif p isa AbstractVector
                p
            else
                throw(ArgumentError("Unsupported parameter container type $(typeof(p)) for MTK sparse Jacobian"))
            end

        if J isa SparseMatrixCSC
            fill!(J.nzval, 0.0)
            jac_sparse!(J, u, pvals, t)
            sanitize_sparse_jac!(J)
        elseif J isa AbstractMatrix
            Jsp = copy(proto)
            fill!(Jsp.nzval, 0.0)
            jac_sparse!(Jsp, u, pvals, t)
            sanitize_sparse_jac!(Jsp)
            copyto!(J, Matrix(Jsp))
        else
            Jsp = copy(proto)
            fill!(Jsp.nzval, 0.0)
            jac_sparse!(Jsp, u, pvals, t)
            sanitize_sparse_jac!(Jsp)
            copyto!(J, vec(Matrix(Jsp)))
        end
        return nothing
    end

    MTK_JAC_CACHE[key] = jac!
    return jac!
end

function make_mtk_sparse_jacobian_prototype(mdl)
    get(mdl, :core_mode, :legacy_reference) == :production &&
        throw(ArgumentError("MTK sparse Jacobian prototypes are only available on the legacy/reference model path"))
    key = mtk_jac_cache_key(mdl; sparse = true)
    if haskey(MTK_SPARSE_PROTOTYPE_CACHE, key)
        return MTK_SPARSE_PROTOTYPE_CACHE[key]
    end

    sym_rhs = build_symbolic_probe_rhs_function(
        mdl.name_to_idx,
        mdl.repeated_rule_defs,
        mdl.rate_exprs,
        mdl.stoich,
        length(mdl.u0),
        length(mdl.pvals),
    )
    probe_prob = ODEProblem(sym_rhs, copy(mdl.u0), (0.0, 1.0), copy(mdl.pvals))
    sys = ModelingToolkit.modelingtoolkitize(probe_prob)
    jac_exprs = ModelingToolkit.generate_jacobian(sys; sparse = true, simplify = false)
    (jac_expr_outofplace, _) = jac_exprs
    jac_sparse_outofplace = RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, normalize_function_expr(jac_expr_outofplace))
    proto = jac_sparse_outofplace(copy(mdl.u0), copy(mdl.pvals), 0.0)
    if !(proto isa SparseMatrixCSC)
        proto = sparse(proto)
    end
    fill!(proto.nzval, 0.0)
    MTK_SPARSE_PROTOTYPE_CACHE[key] = proto
    return proto
end

function extract_identifiers(expr::String)
    ids = String[]
    for m in eachmatch(r"[A-Za-z_][A-Za-z0-9_]*", expr)
        push!(ids, String(m.match))
    end
    return ids
end

function repeated_rule_order(rules::Vector{Tuple{String,String}})
    n = length(rules)
    lhs_names = [sanitize_name(r[1]) for r in rules]
    rhs_exprs = [sanitize_formula(r[2]) for r in rules]
    lhs_to_idx = Dict(lhs => i for (i, lhs) in enumerate(lhs_names))

    incoming = [Set{Int}() for _ in 1:n]
    outgoing = [Int[] for _ in 1:n]

    for i in 1:n
        ids = extract_identifiers(rhs_exprs[i])
        for id in ids
            if haskey(lhs_to_idx, id)
                dep = lhs_to_idx[id]
                dep == i && continue
                if !(dep in incoming[i])
                    push!(incoming[i], dep)
                    push!(outgoing[dep], i)
                end
            end
        end
    end

    indeg = [length(incoming[i]) for i in 1:n]
    q = [i for i in 1:n if indeg[i] == 0]
    order = Int[]

    while !isempty(q)
        i = popfirst!(q)
        push!(order, i)
        for j in outgoing[i]
            indeg[j] -= 1
            if indeg[j] == 0
                push!(q, j)
            end
        end
    end

    if length(order) != n
        return collect(1:n)
    end
    return order
end

function parse_rule(rule_text::String)
    parts = split(rule_text, "=", limit = 2)
    length(parts) == 2 || error("Could not parse rule: $rule_text")
    lhs = strip(parts[1])
    rhs = strip(parts[2])
    return sanitize_name(lhs), rhs
end

function parse_reaction_side(side::AbstractString)
    out = Vector{Pair{String,Float64}}()
    for raw_term in split(side, "+")
        term = strip(raw_term)
        isempty(term) && continue
        lowercase(term) == "null" && continue
        m = match(r"^([0-9]+(?:\.[0-9]+)?)\s+(.+)$", term)
        if m === nothing
            push!(out, term => 1.0)
        else
            coeff = parse(Float64, m.captures[1])
            name = strip(m.captures[2])
            push!(out, name => coeff)
        end
    end
    return out
end

function reaction_to_stoich(reaction::String, state_to_idx::Dict{String,Int}, repeated_lhs::Set{String})
    arrow = occursin("<->", reaction) ? "<->" : "->"
    parts = split(reaction, arrow)
    length(parts) == 2 || error("Could not parse reaction string: $reaction")

    left = parse_reaction_side(strip(parts[1]))
    right = parse_reaction_side(strip(parts[2]))

    delta = Dict{Int,Float64}()

    for (name, coeff) in left
        haskey(state_to_idx, name) || continue
        name in repeated_lhs && continue
        i = state_to_idx[name]
        delta[i] = get(delta, i, 0.0) - coeff
    end
    for (name, coeff) in right
        haskey(state_to_idx, name) || continue
        name in repeated_lhs && continue
        i = state_to_idx[name]
        delta[i] = get(delta, i, 0.0) + coeff
    end

    return collect(delta)
end

function parse_variant_value(v)
    if v isa Number
        return Float64(v)
    end
    s = strip(String(v))
    try
        return parse(Float64, s)
    catch
        return NaN
    end
end

function scalar_or_missing(v)
    if v === missing || isnothing(v)
        return ""
    end
    return string(v)
end

function parse_json_number_vector(v)::Vector{Float64}
    s = strip(scalar_or_missing(v))
    if isempty(s)
        return Float64[]
    end

    try
        parsed = JSON3.read(s)
        if parsed isa Number
            return [Float64(parsed)]
        elseif parsed isa AbstractVector
            out = Float64[]
            for x in parsed
                if x isa Number
                    push!(out, Float64(x))
                elseif x === nothing
                    continue
                else
                    xs = strip(String(x))
                    if isempty(xs)
                        continue
                    end
                    push!(out, parse(Float64, xs))
                end
            end
            return out
        end
    catch
        # Fall back to scalars and comma-separated fields.
    end

    if occursin(",", s)
        toks = split(s, ",")
        out = Float64[]
        for tok in toks
            st = strip(tok)
            isempty(st) && continue
            push!(out, parse(Float64, st))
        end
        return out
    end
    return [parse(Float64, s)]
end

function parse_first_number(v; default::Float64 = 0.0)::Float64
    vec = parse_json_number_vector(v)
    isempty(vec) && return default
    return vec[1]
end

function normalize_vector_length(vals::Vector{Float64}, n::Int, default::Float64 = 0.0)::Vector{Float64}
    if n <= 0
        return Float64[]
    end
    if isempty(vals)
        return fill(default, n)
    elseif length(vals) == 1
        return fill(vals[1], n)
    elseif length(vals) == n
        return vals
    elseif length(vals) > n
        return vals[1:n]
    else
        out = copy(vals)
        append!(out, fill(vals[end], n - length(vals)))
        return out
    end
end

function dose_from_row(r)
    dclass = String(r.class)
    dname = String(r.name)
    target =
        if hasproperty(r, :target_name)
            String(r.target_name)
        else
            String(r.target)
        end

    amount_vec =
        if hasproperty(r, :amount_json)
            parse_json_number_vector(r.amount_json)
        else
            parse_json_number_vector(r.amount)
        end
    rate_vec =
        if hasproperty(r, :rate_json)
            parse_json_number_vector(r.rate_json)
        else
            parse_json_number_vector(r.rate)
        end

    times = Float64[]
    if occursin("RepeatDose", dclass)
        startv =
            if hasproperty(r, :start_time_json)
                parse_first_number(r.start_time_json; default = 0.0)
            else
                parse_first_number(r.starttime; default = 0.0)
            end
        intervalv =
            if hasproperty(r, :interval_json)
                parse_first_number(r.interval_json; default = 0.0)
            else
                parse_first_number(r.interval; default = 0.0)
            end
        repeatcountv =
            if hasproperty(r, :repeat_count_json)
                parse_first_number(r.repeat_count_json; default = 0.0)
            else
                parse_first_number(r.repeatcount; default = 0.0)
            end
        n_doses = max(1, Int(round(repeatcountv)) + 1)
        times = [startv + intervalv * i for i in 0:(n_doses - 1)]
    else
        times =
            if hasproperty(r, :time_json)
                parse_json_number_vector(r.time_json)
            else
                parse_json_number_vector(r.time)
            end
        if isempty(times)
            startv =
                if hasproperty(r, :start_time_json)
                    parse_first_number(r.start_time_json; default = 0.0)
                else
                    parse_first_number(r.starttime; default = 0.0)
                end
            times = [startv]
        end
    end

    n = length(times)
    amounts = normalize_vector_length(amount_vec, n, 0.0)
    rates = normalize_vector_length(rate_vec, n, 0.0)
    return DoseDef(dname, target, amounts, rates, times)
end

function build_stoich_from_table(stoich_df::DataFrame, n_rxn::Int)
    stoich = [Pair{Int,Float64}[] for _ in 1:n_rxn]
    rx_col = hasproperty(stoich_df, :reaction_idx) ? :reaction_idx : :rxn_idx
    coeff_col = hasproperty(stoich_df, :stoich_coeff) ? :stoich_coeff : :stoich
    for r in eachrow(stoich_df)
        rx = Int(getproperty(r, rx_col))
        sp = Int(r.species_idx)
        coeff = Float64(getproperty(r, coeff_col))
        push!(stoich[rx], sp => coeff)
    end
    return stoich
end

function load_artifacts(; model_tables_dir::String = "")
    tables_dir = model_tables_dir
    if isempty(tables_dir)
        tables_dir = isdir(RECON_BUNDLE_DIR) ? RECON_BUNDLE_DIR : MODEL_TABLES_DIR
    end

    species_df = DataFrame(CSV.File(joinpath(tables_dir, "species.tsv"); delim = '\t'))
    params_df = DataFrame(CSV.File(joinpath(tables_dir, "parameters.tsv"); delim = '\t'))
    rules_df = DataFrame(CSV.File(joinpath(tables_dir, "rules.tsv"); delim = '\t'))
    rxn_df = DataFrame(CSV.File(joinpath(tables_dir, "reactions.tsv"); delim = '\t'))
    var_df = DataFrame(CSV.File(joinpath(tables_dir, "variants.tsv"); delim = '\t'))
    stoich_df =
        isfile(joinpath(tables_dir, "stoich_sparse.tsv")) ?
        DataFrame(CSV.File(joinpath(tables_dir, "stoich_sparse.tsv"); delim = '\t')) :
        DataFrame(CSV.File(joinpath(tables_dir, "stoich_truth.tsv"); delim = '\t'))
    doses_df = DataFrame(CSV.File(joinpath(tables_dir, "doses.tsv"); delim = '\t'))

    state_names = String.(species_df.name)
    state_u0 = Float64.(species_df.initial_amount)
    param_names = String.(params_df.name)
    param_defaults = Float64.(params_df.value)

    initial_rules = Tuple{String,String}[]
    repeated_rules = Tuple{String,String}[]

    for r in eachrow(rules_df)
        lhs, rhs = parse_rule(String(r.rule))
        rtype = hasproperty(r, :rule_type) ? String(r.rule_type) : String(r.type)
        if rtype == "initialAssignment"
            push!(initial_rules, (lhs, rhs))
        elseif rtype == "repeatedAssignment"
            push!(repeated_rules, (lhs, rhs))
        end
    end

    reaction_strs = String.(rxn_df.reaction)
    rate_exprs =
        if hasproperty(rxn_df, :reaction_rate)
            String.(rxn_df.reaction_rate)
        else
            String.(rxn_df.rate)
        end

    variants = Dict{Int, Vector{Pair{String,Float64}}}()
    for r in eachrow(var_df)
        idx = hasproperty(r, :variant_idx) ? Int(r.variant_idx) : Int(r.idx)
        action = String(r.action)
        param_name = String(r.class)
        field_name = String(r.name)
        value = if hasproperty(r, :value_numeric)
            r.value_numeric === missing ? NaN : Float64(r.value_numeric)
        else
            parse_variant_value(r.value)
        end
        is_numeric = hasproperty(r, :is_numeric) ? Int(r.is_numeric) == 1 : !isnan(value)
        if action == "parameter" && field_name == "Value" && is_numeric
            if !haskey(variants, idx)
                variants[idx] = Pair{String,Float64}[]
            end
            push!(variants[idx], param_name => value)
        end
    end

    doses = Dict{String,DoseDef}()
    for r in eachrow(doses_df)
        ddef = dose_from_row(r)
        doses[ddef.name] = ddef
    end

    stoich = build_stoich_from_table(stoich_df, nrow(rxn_df))

    return ModelArtifacts(
        state_names,
        state_u0,
        param_names,
        param_defaults,
        initial_rules,
        repeated_rules,
        rate_exprs,
        reaction_strs,
        variants,
        stoich,
        doses,
    )
end

function case_variant_ids(case_no::Int)
    if case_no in (10, 11, 12, 13)
        return [5, 9, 11, 13, 14, 24]
    elseif case_no in (20, 21, 22)
        return [5, 9, 11, 13, 14, 15, 24]
    elseif case_no == 30
        return [5, 9, 11, 13, 14, 20, 22, 24]
    else
        error("Unsupported case_no: $case_no")
    end
end

function dose_lookup(artifacts::ModelArtifacts, name::String)
    haskey(artifacts.doses, name) || error("Dose '$name' not found in exported doses table")
    return artifacts.doses[name]
end

function expand_dose_events(dose::DoseDef)
    nts = length(dose.times)
    amounts = normalize_vector_length(dose.amount, nts, 0.0)
    rates = normalize_vector_length(dose.rate, nts, 0.0)
    return dose.times, amounts, rates
end

function build_model_with_variant_ids(variant_ids::Vector{Int}, param_names::Vector{String}, param_values::Vector{Float64}; core_mode_override::Symbol = DEFAULT_CORE_MODE)
    artifacts = load_artifacts()

    load_pk_datasets!()

    state_names = copy(artifacts.state_names)
    param_names_all = copy(artifacts.param_names)

    aux_names = String[]
    for (lhs, _) in artifacts.initial_rules
        lhss = sanitize_name(lhs)
        if !(lhss in state_names) && !(lhss in param_names_all) && !(lhss in aux_names)
            push!(aux_names, lhss)
        end
    end
    for (lhs, _) in artifacts.repeated_rules
        lhss = sanitize_name(lhs)
        if !(lhss in state_names) && !(lhss in param_names_all) && !(lhss in aux_names)
            push!(aux_names, lhss)
        end
    end

    all_names = vcat(state_names, param_names_all, aux_names)
    name_to_idx = Dict(name => i for (i, name) in enumerate(all_names))

    state_to_idx = Dict(name => i for (i, name) in enumerate(state_names))
    param_to_idx = Dict(name => i for (i, name) in enumerate(param_names_all))
    direct_model_layout_matches(state_names, param_names_all) ||
        error("Manual model core layout does not match exported tables.")

    rep_order = repeated_rule_order(artifacts.repeated_rules)
    repeated_rule_defs = Tuple{Int,String}[]
    for idx in rep_order
        lhs, rhs = artifacts.repeated_rules[idx]
        lhs_s = sanitize_name(lhs)
        lhs_idx = name_to_idx[lhs_s]
        push!(repeated_rule_defs, (lhs_idx, rhs))
    end

    stoich = artifacts.stoich

    pvals = copy(artifacts.param_defaults)

    for vid in variant_ids
        entries = get(artifacts.variants, vid, Pair{String,Float64}[])
        for (pname, pval) in entries
            if haskey(param_to_idx, pname)
                pvals[param_to_idx[pname]] = pval
            end
        end
    end

    for (pn, pv) in zip(param_names, param_values)
        if haskey(param_to_idx, pn)
            pvals[param_to_idx[pn]] = pv
        end
    end

    u0 = copy(artifacts.state_u0)

    z = zeros(Float64, length(all_names))
    z[1:length(state_names)] .= u0
    z[length(state_names)+1:length(state_names)+length(param_names_all)] .= pvals

    apply_initial_rules_direct!(z, u0, pvals)

    core_params = build_core_params(param_names_all, pvals)
    core_mode = core_mode_override == :auto ? canonical_core_mode(param_to_idx, pvals) : core_mode_override
    core_mode in (:production, :legacy_reference) || error("Unsupported core_mode_override=$core_mode_override")
    core_u0 = MosunModelCore.pack_state(MosunModelCore.initial_state(core_params))
    core_state_names_vec = core_state_names()
    core_state_to_idx_map = core_state_to_idx()
    core_observable_names_set = core_observable_names()
    canonical_rhs =
        if core_mode == :production
            let params_fixed = deepcopy(core_params)
                (du, u, ctx, t) -> MosunModelCore.mosun_rhs!(du, u, params_fixed, t, ctx.cache; infusions = ctx.infusions)
            end
        else
            rhs_direct!
        end

    return (
        state_names = state_names,
        state_to_idx = state_to_idx,
        name_to_idx = name_to_idx,
        param_names = param_names_all,
        pvals = pvals,
        u0 = u0,
        core_mode = core_mode,
        core_params = core_params,
        core_u0 = core_u0,
        core_state_names = core_state_names_vec,
        core_state_to_idx = core_state_to_idx_map,
        core_observable_names = core_observable_names_set,
        z = z,
        repeated_rule_exprs = RuleExpr[],
        repeated_rule_defs = repeated_rule_defs,
        rate_fns = Function[],
        rate_exprs = artifacts.rate_exprs,
        stoich = stoich,
        canonical_rhs = canonical_rhs,
    )
end

function build_model(case_no::Int, param_names::Vector{String}, param_values::Vector{Float64}; core_mode_override::Symbol = DEFAULT_CORE_MODE)
    return build_model_with_variant_ids(case_variant_ids(case_no), param_names, param_values; core_mode_override = core_mode_override)
end

function build_dose_events(selected_doses::Vector{String}, state_to_idx::Dict{String,Int}, artifacts::ModelArtifacts)
    discrete_times = Float64[]
    discrete_map = Dict{Float64, Vector{Tuple{Int,Float64}}}()
    infusions = InfusionEvent[]

    for dname in selected_doses
        ddef = dose_lookup(artifacts, dname)
        target_idx = state_to_idx[ddef.target]
        times, amounts, rates = expand_dose_events(ddef)

        for i in eachindex(times)
            t0 = times[i]
            amt = amounts[i]
            rate = rates[i]
            if rate == 0.0
                if !haskey(discrete_map, t0)
                    discrete_map[t0] = Tuple{Int,Float64}[]
                    push!(discrete_times, t0)
                end
                push!(discrete_map[t0], (target_idx, amt))
            else
                t1 = t0 + amt / rate
                push!(infusions, InfusionEvent(target_idx, t0, t1, rate))
            end
        end
    end

    sort!(discrete_times)
    unique!(discrete_times)

    return discrete_times, discrete_map, infusions
end

function build_core_regimen(selected_doses::Vector{String}, artifacts::ModelArtifacts)
    events = MosunModelCore.MosunRegimenEvent[]
    for dname in selected_doses
        ddef = dose_lookup(artifacts, dname)
        times, amounts, rates = expand_dose_events(ddef)
        target = Symbol(ddef.target)
        for i in eachindex(times)
            push!(events, MosunModelCore.MosunRegimenEvent(
                target = target,
                time = Float64(times[i]),
                amount = Float64(amounts[i]),
                rate = Float64(rates[i]),
            ))
        end
    end
    return MosunModelCore.MosunRegimen(events = events)
end

function rhs!(du, u, ctx::SimContext, t)
    ns = length(u)
    np = length(ctx.pvals)

    @inbounds begin
        ctx.z[1:ns] .= u
        ctx.z[ns+1:ns+np] .= ctx.pvals

        for r in ctx.repeated_rules
            ctx.z[r.lhs_idx] = r.fn(ctx.z, t)
        end

        fill!(du, 0.0)

        for inf in ctx.infusions
            if t >= inf.t_start && t <= inf.t_end
                du[inf.target_idx] += inf.rate
            end
        end

        for j in eachindex(ctx.rate_fns)
            rate = ctx.rate_fns[j](ctx.z, t)
            for (i, coeff) in ctx.stoich[j]
                du[i] += coeff * rate
            end
        end
    end

    return nothing
end

function run_simulation(
        case_no::Int,
        param_names::Vector{String},
        param_values::Vector{Float64},
        selected_doses::Vector{String},
        save_times::Vector{Float64};
        engine::Symbol = :canonical,
        core_mode_override::Symbol = DEFAULT_CORE_MODE)
    artifacts = load_artifacts()
    mdl = build_model(case_no, param_names, param_values; core_mode_override = core_mode_override)
    dose_state_to_idx = engine == :canonical && get(mdl, :core_mode, :legacy_reference) == :production ? mdl.core_state_to_idx : mdl.state_to_idx
    discrete_times, discrete_map, infusions = build_dose_events(selected_doses, dose_state_to_idx, artifacts)
    t0 = minimum(save_times)
    tf = maximum(save_times)
    save_times_sorted = sort(unique(save_times))

    if engine == :canonical && get(mdl, :core_mode, :legacy_reference) == :production
        regimen = build_core_regimen(selected_doses, artifacts)
        tstart = min(0.0, t0)
        built = MosunModelCore.build_problem(
            regimen,
            mdl.core_params;
            tspan = (tstart, tf),
            saveat = save_times_sorted,
            callback_mode = :callback,
        )
        sol = MosunModelCore.solve_problem(
            built,
            make_solver_alg();
            abstol = SOLVER_ABSTOL,
            reltol = SOLVER_RELTOL,
        )
        sol.retcode == SciMLBase.ReturnCode.Success || error("solve failed with retcode=$(sol.retcode)")
        return mdl, sol
    end

    rhs_fun = rhs!
    ctx = nothing
    if engine == :legacy
        ctx = make_legacy_context(mdl, infusions)
        rhs_fun = rhs!
    elseif engine == :canonical
        ctx = CanonicalSimContext(copy(mdl.z), mdl.pvals, infusions)
        rhs_fun = mdl.canonical_rhs
    else
        error("Unsupported engine=$engine. Use :legacy or :canonical.")
    end

    u_curr = engine == :canonical && get(mdl, :core_mode, :legacy_reference) == :production ? copy(mdl.core_u0) : copy(mdl.u0)
    t_curr = t0
    alg = make_solver_alg()
    sol_t = Float64[]
    sol_u = Vector{Vector{Float64}}()

    function append_solution!(ts::Vector{Float64}, us::Vector{Vector{Float64}})
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

    all_event_times = sort(unique(vcat(discrete_times, [tf])))

    for te in all_event_times
        te < t_curr && continue

        seg_save = [t for t in save_times_sorted if t_curr <= t <= te]
        if te > t_curr
            prob = ODEProblem(rhs_fun, u_curr, (t_curr, te), ctx)
            sol_seg = solve(
                prob,
                alg;
                abstol = SOLVER_ABSTOL,
                reltol = SOLVER_RELTOL,
                saveat = seg_save,
                tstops = [te],
            )
            append_solution!(sol_seg.t, sol_seg.u)
            u_curr = copy(sol_seg.u[end])
        elseif !isempty(seg_save) && (isempty(sol_t) || !isapprox(sol_t[end], te; atol = 1e-12, rtol = 0.0))
            push!(sol_t, te)
            push!(sol_u, copy(u_curr))
        end

        if haskey(discrete_map, te)
            for (idx, amt) in discrete_map[te]
                u_curr[idx] += amt
            end
        end
        t_curr = te
    end

    if isempty(sol_t) || !isapprox(sol_t[end], tf; atol = 1e-12, rtol = 0.0)
        prob = ODEProblem(rhs_fun, u_curr, (t_curr, tf), ctx)
        sol_seg = solve(prob, alg; abstol = SOLVER_ABSTOL, reltol = SOLVER_RELTOL, saveat = [tf], tstops = [tf])
        append_solution!(sol_seg.t, sol_seg.u)
    end

    # Align outputs to requested save times order.
    lookup = Dict{Float64, Vector{Float64}}()
    for i in eachindex(sol_t)
        lookup[sol_t[i]] = sol_u[i]
    end
    aligned_u = [lookup[t] for t in save_times]

    return mdl, (t = save_times, u = aligned_u)
end

function extract_outputs(sol, mdl, outvec::Vector{String})
    n = length(sol.t)
    m = length(outvec)
    X = zeros(Float64, n, m)
    if n > 0 && uses_production_layout(mdl, sol.u[1])
        cache = MosunModelCore.zero_observables_cache()
        for i in 1:n
            u = sol.u[i]
            t = sol.t[i]
            MosunModelCore.update_observables!(cache, u, mdl.core_params, t)
            for (j, name0) in enumerate(outvec)
                X[i, j] = MosunModelCore.state_or_observable(u, cache, Symbol(sanitize_name(name0)))
            end
        end
        return X
    end

    z = copy(mdl.z)
    for i in 1:n
        u = sol.u[i]
        t = sol.t[i]
        populate_value_buffer!(z, mdl, u, t)
        for (j, name0) in enumerate(outvec)
            name = sanitize_name(name0)
            idx = mdl.name_to_idx[name]
            X[i, j] = z[idx]
        end
    end

    return X
end

function run_manifest_row(row; engine::Symbol = :canonical, core_mode_override::Symbol = DEFAULT_CORE_MODE)
    case_no = Int(row.case_no)
    outvec = Vector{String}(JSON3.read(String(row.outvec_json), Vector{String}))
    param_names = Vector{String}(JSON3.read(String(row.param_names_json), Vector{String}))
    param_values = Float64.(JSON3.read(String(row.param_values_json), Vector{Float64}))
    selected_doses = Vector{String}(JSON3.read(String(row.selected_doses_json), Vector{String}))

    ref_path = joinpath(REPO_ROOT, String(row.sim_csv_relpath))
    ref_df = DataFrame(CSV.File(ref_path))
    save_times = Float64.(ref_df.time)

    mdl, sol = run_simulation(case_no, param_names, param_values, selected_doses, save_times; engine = engine, core_mode_override = core_mode_override)
    X = extract_outputs(sol, mdl, outvec)

    return save_times, outvec, X, ref_df
end

export MosunModelCore, DEFAULT_CORE_MODE, run_manifest_row, REPO_ROOT, build_model_with_variant_ids, run_simulation, make_mtk_dense_jacobian, make_mtk_sparse_jacobian, make_mtk_sparse_jacobian_prototype, canonical_u0, canonical_state_to_idx

end # module
