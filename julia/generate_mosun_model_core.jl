using Printf

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

const OUTPUT_PATH = joinpath(@__DIR__, "src", "MosunModelCore.jl")

function emit_expr(expr)
    return sprint(Base.show_unquoted, expr)
end

function write_symbol_vector(io, name::String, values)
    println(io, "const $name = (")
    for v in values
        println(io, "    :", v, ",")
    end
    println(io, ")")
    println(io)
end

function build_dynamic_state_info(artifacts)
    dyn_idx = sort(unique(first(pair) for pairs in artifacts.stoich for pair in pairs))
    dyn_names = artifacts.state_names[dyn_idx]
    dyn_name_to_new_idx = Dict(name => i for (i, name) in enumerate(dyn_names))
    dyn_name_to_old_idx = Dict(name => dyn_idx[i] for (i, name) in enumerate(dyn_names))
    return dyn_idx, dyn_names, dyn_name_to_new_idx, dyn_name_to_old_idx
end

function build_observable_info(artifacts)
    rep_order = TCellEngagerQSP.repeated_rule_order(artifacts.repeated_rules)
    rep_pairs = Tuple{String,String}[]
    for idx in rep_order
        lhs, rhs = artifacts.repeated_rules[idx]
        push!(rep_pairs, (TCellEngagerQSP.sanitize_name(lhs), rhs))
    end
    return rep_pairs
end

function build_initial_rule_map(artifacts)
    out = Dict{String,String}()
    for (lhs, rhs) in artifacts.initial_rules
        out[TCellEngagerQSP.sanitize_name(lhs)] = rhs
    end
    return out
end

function collect_symbol_refs!(refs::Set{String}, ex)
    if ex isa Symbol
        push!(refs, String(ex))
    elseif ex isa Expr
        for arg in ex.args
            collect_symbol_refs!(refs, arg)
        end
    end
    return refs
end

function compute_rhs_observable_names(repeated_pairs, rate_exprs)
    obs_names = Set(lhs for (lhs, _) in repeated_pairs)
    obs_rhs_map = Dict(lhs => rhs for (lhs, rhs) in repeated_pairs)
    needed = Set{String}()
    for rhs in rate_exprs
        refs = collect_symbol_refs!(Set{String}(), Meta.parse(rhs))
        foreach(ref -> (ref in obs_names) && push!(needed, ref), refs)
    end
    closure = Set{String}()
    stack = collect(needed)
    while !isempty(stack)
        lhs = pop!(stack)
        lhs in closure && continue
        push!(closure, lhs)
        refs = collect_symbol_refs!(Set{String}(), Meta.parse(obs_rhs_map[lhs]))
        for ref in refs
            ref in obs_names && push!(stack, ref)
        end
    end
    return [lhs for (lhs, _) in repeated_pairs if lhs in closure]
end

function rewrite_core_expr(
        ex,
        state_names::Set{String},
        param_names::Set{String},
        observable_names::Set{String};
        observable_style::Symbol = :cache,
        param_style::Symbol = :struct,
        pow_function::Symbol = :pow_safe,
        pk_function::Symbol = :tdb_central_concentration)
    if ex isa Symbol
        s = String(ex)
        if s == "time"
            return :t
        elseif s in state_names
            return Symbol(s)
        elseif s in param_names
            if param_style == :struct
                return Expr(:., :p, QuoteNode(Symbol(s)))
            elseif param_style == :local
                return Symbol(s)
            else
                error("Unsupported param_style=$param_style")
            end
        elseif s in observable_names
            if observable_style == :cache
                return Expr(:., :cache, QuoteNode(Symbol(s)))
            elseif observable_style == :local
                return Symbol(s)
            else
                error("Unsupported observable_style=$observable_style")
            end
        else
            return ex
        end
    elseif ex isa Expr
        if ex.head == :call
            f = ex.args[1]
            if f == :^
                a1 = rewrite_core_expr(ex.args[2], state_names, param_names, observable_names; observable_style = observable_style, param_style = param_style, pow_function = pow_function, pk_function = pk_function)
                a2 = rewrite_core_expr(ex.args[3], state_names, param_names, observable_names; observable_style = observable_style, param_style = param_style, pow_function = pow_function, pk_function = pk_function)
                return Expr(:call, pow_function, a1, a2)
            elseif f == :PK_v26
                new_args = map(arg -> rewrite_core_expr(arg, state_names, param_names, observable_names; observable_style = observable_style, param_style = param_style, pow_function = pow_function, pk_function = pk_function), ex.args[2:end])
                return Expr(:call, pk_function, new_args...)
            else
                new_args = Any[rewrite_core_expr(f, state_names, param_names, observable_names; observable_style = observable_style, param_style = param_style, pow_function = pow_function, pk_function = pk_function)]
                append!(new_args, map(arg -> rewrite_core_expr(arg, state_names, param_names, observable_names; observable_style = observable_style, param_style = param_style, pow_function = pow_function, pk_function = pk_function), ex.args[2:end]))
                return Expr(:call, new_args...)
            end
        else
            return Expr(ex.head, map(arg -> rewrite_core_expr(arg, state_names, param_names, observable_names; observable_style = observable_style, param_style = param_style, pow_function = pow_function, pk_function = pk_function), ex.args)...)
        end
    else
        return ex
    end
end

function write_ifelseif_lookup(io, fname::String, arg_decl::String, entries::Vector{String}, value_expr)
    println(io, "function $fname($arg_decl)")
    for (i, nm) in enumerate(entries)
        prefix = i == 1 ? "    if" : "    elseif"
        println(io, "$prefix name === :", nm)
        println(io, "        return ", value_expr(nm))
    end
    println(io, "    end")
    println(io, "    throw(ArgumentError(\"Unknown name: \$(name)\"))")
    println(io, "end")
    println(io)
end

function write_module(io, artifacts)
    dyn_idx, dyn_names, dyn_name_to_new_idx, _ = build_dynamic_state_info(artifacts)
    repeated_pairs = build_observable_info(artifacts)
    observable_names = [lhs for (lhs, _) in repeated_pairs]
    observable_set = Set(observable_names)
    rhs_observable_names = compute_rhs_observable_names(repeated_pairs, artifacts.rate_exprs)
    rhs_observable_set = Set(rhs_observable_names)
    param_names = copy(artifacts.param_names)
    param_set = Set(param_names)
    state_set = Set(dyn_names)
    initial_map = build_initial_rule_map(artifacts)
    species_default = Dict(name => artifacts.state_u0[idx] for (idx, name) in zip(1:length(artifacts.state_names), artifacts.state_names))
    dead_legacy_names = [nm for nm in artifacts.state_names if !(nm in dyn_names) && !(nm in observable_names)]

    println(io, "module MosunModelCore")
    println(io)
    println(io, "using SciMLBase")
    println(io, "using DifferentialEquations")
    println(io, "using DiffEqCallbacks")
    println(io)
    println(io, "pow_safe(a::Float64, b::Float64) = Float64(real((complex(a) ^ b)))")
    println(io, "pow_safe(a::Float64, b::Integer) = Float64(real((complex(a) ^ b)))")
    println(io, "pow_safe(a::Real, b::Real) = real((complex(a) ^ b))")
    println(io, "pow_safe(a, b) = a ^ b")
    println(io, "pow_safe_symbolic(a::Float64, b::Float64) = pow_safe(a, b)")
    println(io, "pow_safe_symbolic(a::Float64, b::Integer) = pow_safe(a, b)")
    println(io, "pow_safe_symbolic(a, b) = a ^ b")
    println(io)
    println(io, "@inline function tdb_central_concentration(TDBc_ugperkg, Vc_tdb, PKflag, VPid, t, end_time, fvalidation)")
    println(io, "    val = TDBc_ugperkg / Vc_tdb")
    println(io, "    return ifelse(val > 1e-5, val, zero(val))")
    println(io, "end")
    println(io, "@inline tdb_central_concentration_symbolic(TDBc_ugperkg, Vc_tdb, PKflag, VPid, t, end_time, fvalidation) = TDBc_ugperkg / Vc_tdb")
    println(io)

    write_symbol_vector(io, "DYNAMIC_STATE_NAMES", dyn_names)
    write_symbol_vector(io, "OBSERVABLE_NAMES", observable_names)
    write_symbol_vector(io, "RHS_OBSERVABLE_NAMES", rhs_observable_names)
    write_symbol_vector(io, "PARAMETER_NAMES", param_names)
    write_symbol_vector(io, "DEAD_LEGACY_NAMES", dead_legacy_names)
    println(io, "const DYNAMIC_STATE_COUNT = ", length(dyn_names))
    println(io, "const OBSERVABLE_COUNT = ", length(observable_names))
    println(io, "const RHS_OBSERVABLE_COUNT = ", length(rhs_observable_names))
    println(io)

    println(io, "Base.@kwdef mutable struct MosunParams")
    for (nm, val) in zip(param_names, artifacts.param_defaults)
        println(io, "    ", nm, "::Float64 = ", repr(val))
    end
    println(io, "end")
    println(io)

    println(io, "Base.@kwdef struct MosunDynamicState")
    for nm in dyn_names
        println(io, "    ", nm, "::Float64 = 0.0")
    end
    println(io, "end")
    println(io)

    println(io, "Base.@kwdef mutable struct MosunObservablesCache")
    for nm in observable_names
        println(io, "    ", nm, "::Float64 = 0.0")
    end
    println(io, "end")
    println(io)

    println(io, "Base.@kwdef struct MosunRegimenEvent")
    println(io, "    target::Symbol")
    println(io, "    time::Float64")
    println(io, "    amount::Float64 = 0.0")
    println(io, "    rate::Float64 = 0.0")
    println(io, "end")
    println(io)

    println(io, "Base.@kwdef struct MosunRegimen")
    println(io, "    events::Vector{MosunRegimenEvent} = MosunRegimenEvent[]")
    println(io, "end")
    println(io)

    println(io, "Base.@kwdef mutable struct MosunProblemContext")
    println(io, "    params::MosunParams")
    println(io, "    cache::MosunObservablesCache = MosunObservablesCache()")
    println(io, "    active_rates::Vector{Float64} = zeros(Float64, DYNAMIC_STATE_COUNT)")
    println(io, "end")
    println(io)

    println(io, "Base.@kwdef struct MosunBuiltProblem")
    println(io, "    prob")
    println(io, "    ctx::MosunProblemContext")
    println(io, "    initial_u::Vector{Float64}")
    println(io, "    callback")
    println(io, "    tstops::Vector{Float64}")
    println(io, "    d_discontinuities::Vector{Float64}")
    println(io, "    event_map::Dict{Float64, Vector{Tuple{Int,Float64,Float64}}}")
    println(io, "    saveat")
    println(io, "    regimen::MosunRegimen")
    println(io, "    callback_mode::Symbol")
    println(io, "end")
    println(io)

    println(io, "default_params() = MosunParams()")
    println(io, "zero_observables_cache() = MosunObservablesCache()")
    println(io, "has_parameter(name) = Symbol(name) in PARAMETER_NAMES")
    println(io)
    println(io, "function set_param!(p::MosunParams, name, value)")
    println(io, "    sym = Symbol(name)")
    println(io, "    sym in PARAMETER_NAMES || throw(ArgumentError(\"Unknown parameter: \$(name)\"))")
    println(io, "    setproperty!(p, sym, Float64(value))")
    println(io, "    return p")
    println(io, "end")
    println(io)
    println(io, "function params_from_named_values(names, values; base::MosunParams = default_params(), ignore_unknown::Bool = false)")
    println(io, "    length(names) == length(values) || throw(ArgumentError(\"names and values must have equal length\"))")
    println(io, "    p = deepcopy(base)")
    println(io, "    for (name, value) in zip(names, values)")
    println(io, "        if has_parameter(name)")
    println(io, "            set_param!(p, name, value)")
    println(io, "        elseif !ignore_unknown")
    println(io, "            throw(ArgumentError(\"Unknown parameter: \$(name)\"))")
    println(io, "        end")
    println(io, "    end")
    println(io, "    return p")
    println(io, "end")
    println(io)
    println(io, "function params_from_dict(values; base::MosunParams = default_params(), ignore_unknown::Bool = false)")
    println(io, "    p = deepcopy(base)")
    println(io, "    for (name, value) in pairs(values)")
    println(io, "        if has_parameter(name)")
    println(io, "            set_param!(p, name, value)")
    println(io, "        elseif !ignore_unknown")
    println(io, "            throw(ArgumentError(\"Unknown parameter: \$(name)\"))")
    println(io, "        end")
    println(io, "    end")
    println(io, "    return p")
    println(io, "end")
    println(io)
    println(io, "function bolus_regimen(target::Symbol, dose_map)")
    println(io, "    events = MosunRegimenEvent[]")
    println(io, "    for t in sort(collect(keys(dose_map)))")
    println(io, "        amt = Float64(dose_map[t])")
    println(io, "        amt == 0.0 && continue")
    println(io, "        push!(events, MosunRegimenEvent(target = target, time = Float64(t), amount = amt, rate = 0.0))")
    println(io, "    end")
    println(io, "    return MosunRegimen(events = events)")
    println(io, "end")
    println(io)

    println(io, "function initial_state(p::MosunParams)")
    for nm in dyn_names
        if haskey(initial_map, nm)
            expr = Meta.parse(initial_map[nm])
            rewritten = rewrite_core_expr(expr, state_set, param_set, observable_set)
            println(io, "    ", nm, " = Float64(real(", emit_expr(rewritten), "))")
        else
            println(io, "    ", nm, " = ", repr(Float64(species_default[nm])))
        end
    end
    println(io, "    return MosunDynamicState(")
    for nm in dyn_names
        println(io, "        ", nm, " = ", nm, ",")
    end
    println(io, "    )")
    println(io, "end")
    println(io)

    println(io, "function pack_state(s::MosunDynamicState)")
    println(io, "    return Float64[")
    for nm in dyn_names
        println(io, "        s.", nm, ",")
    end
    println(io, "    ]")
    println(io, "end")
    println(io)

    println(io, "function pack_params(p::MosunParams)")
    println(io, "    return Float64[")
    for nm in param_names
        println(io, "        p.", nm, ",")
    end
    println(io, "    ]")
    println(io, "end")
    println(io)

    println(io, "function unpack_state(u::AbstractVector{<:Real})")
    println(io, "    return MosunDynamicState(")
    for (i, nm) in enumerate(dyn_names)
        println(io, "        ", nm, " = Float64(u[$i]),")
    end
    println(io, "    )")
    println(io, "end")
    println(io)

    write_ifelseif_lookup(io, "dynamic_state_index", "name::Symbol", dyn_names, nm -> string(dyn_name_to_new_idx[nm]))
    write_ifelseif_lookup(io, "dynamic_state_value", "u::AbstractVector{<:Real}, name::Symbol", dyn_names, nm -> "Float64(u[$(dyn_name_to_new_idx[nm])])")
    write_ifelseif_lookup(io, "observable", "cache::MosunObservablesCache, name::Symbol", observable_names, nm -> "cache.$nm")

    for nm in observable_names
        println(io, "observable(cache::MosunObservablesCache, ::Val{:$nm}) = cache.$nm")
    end
    println(io)
    for nm in dyn_names
        println(io, "dynamic_state_value(u::AbstractVector{<:Real}, ::Val{:$nm}) = Float64(u[$(dyn_name_to_new_idx[nm])])")
    end
    println(io)

    println(io, "function state_or_observable(u::AbstractVector{<:Real}, cache::MosunObservablesCache, name::Symbol)")
    println(io, "    if name in DYNAMIC_STATE_NAMES")
    println(io, "        return dynamic_state_value(u, name)")
    println(io, "    elseif name in OBSERVABLE_NAMES")
    println(io, "        return observable(cache, name)")
    println(io, "    end")
    println(io, "    throw(ArgumentError(\"Unknown state or observable: \$(name)\"))")
    println(io, "end")
    println(io)

    println(io, "function update_observables!(cache::MosunObservablesCache, u::AbstractVector{<:Real}, p::MosunParams, t)")
    for (i, nm) in enumerate(dyn_names)
        println(io, "    ", nm, " = u[$i]")
    end
    for (lhs, rhs) in repeated_pairs
        expr = rewrite_core_expr(Meta.parse(rhs), state_set, param_set, observable_set)
        println(io, "    cache.", lhs, " = Float64(real(", emit_expr(expr), "))")
    end
    println(io, "    return cache")
    println(io, "end")
    println(io)

    println(io, "function value_at(u::AbstractVector{<:Real}, p::MosunParams, t, name::Symbol, cache::MosunObservablesCache = zero_observables_cache())")
    println(io, "    update_observables!(cache, u, p, t)")
    println(io, "    return state_or_observable(u, cache, name)")
    println(io, "end")
    println(io)

    println(io, "function mosun_rhs!(du, u::AbstractVector{<:Real}, p::MosunParams, t, cache::MosunObservablesCache; active_rates = nothing, infusions = nothing)")
    for (i, nm) in enumerate(dyn_names)
        println(io, "    ", nm, " = u[$i]")
    end
    for (lhs, rhs) in repeated_pairs
        lhs in rhs_observable_set || continue
        expr = rewrite_core_expr(Meta.parse(rhs), state_set, param_set, rhs_observable_set; observable_style = :local)
        println(io, "    ", lhs, " = Float64(real(", emit_expr(expr), "))")
    end
    println(io, "    fill!(du, 0.0)")
    println(io, "    if active_rates !== nothing")
    println(io, "        @inbounds for i in eachindex(du, active_rates)")
    println(io, "            du[i] += active_rates[i]")
    println(io, "        end")
    println(io, "    elseif infusions !== nothing")
    println(io, "        for inf in infusions")
    println(io, "            if t >= inf.t_start && t <= inf.t_end")
    println(io, "                du[inf.target_idx] += inf.rate")
    println(io, "            end")
    println(io, "        end")
    println(io, "    end")
    for (j, rhs) in enumerate(artifacts.rate_exprs)
        expr = rewrite_core_expr(Meta.parse(rhs), state_set, param_set, rhs_observable_set; observable_style = :local)
        println(io, "    rate_", j, " = Float64(real(", emit_expr(expr), "))")
        for (old_idx, coeff) in artifacts.stoich[j]
            name = artifacts.state_names[old_idx]
            new_idx = dyn_name_to_new_idx[name]
            println(io, "    du[$new_idx] += ", repr(Float64(coeff)), " * rate_", j)
        end
    end
    println(io, "    return nothing")
    println(io, "end")
    println(io)

    println(io, "function mosun_rhs_vector!(du, u::AbstractVector, p::AbstractVector, t)")
    for (i, nm) in enumerate(dyn_names)
        println(io, "    ", nm, " = u[$i]")
    end
    for (i, nm) in enumerate(param_names)
        println(io, "    ", nm, " = p[$i]")
    end
    for (lhs, rhs) in repeated_pairs
        lhs in rhs_observable_set || continue
        expr = rewrite_core_expr(Meta.parse(rhs), state_set, param_set, rhs_observable_set; observable_style = :local, param_style = :local, pow_function = :pow_safe_symbolic, pk_function = :tdb_central_concentration_symbolic)
        println(io, "    ", lhs, " = real(", emit_expr(expr), ")")
    end
    println(io, "    fill!(du, zero(eltype(u)))")
    for (j, rhs) in enumerate(artifacts.rate_exprs)
        expr = rewrite_core_expr(Meta.parse(rhs), state_set, param_set, rhs_observable_set; observable_style = :local, param_style = :local, pow_function = :pow_safe_symbolic, pk_function = :tdb_central_concentration_symbolic)
        println(io, "    rate_", j, " = real(", emit_expr(expr), ")")
        for (old_idx, coeff) in artifacts.stoich[j]
            name = artifacts.state_names[old_idx]
            new_idx = dyn_name_to_new_idx[name]
            println(io, "    du[$new_idx] += ", repr(Float64(coeff)), " * rate_", j)
        end
    end
    println(io, "    return nothing")
    println(io, "end")
    println(io)

    println(io, "function mosun_rhs!(du, u::AbstractVector{<:Real}, ctx::MosunProblemContext, t)")
    println(io, "    return mosun_rhs!(du, u, ctx.params, t, ctx.cache; active_rates = ctx.active_rates)")
    println(io, "end")
    println(io)

    println(io, "function push_event_delta!(event_map::Dict{Float64, Vector{Tuple{Int,Float64,Float64}}}, t::Float64, entry::Tuple{Int,Float64,Float64})")
    println(io, "    if !haskey(event_map, t)")
    println(io, "        event_map[t] = Tuple{Int,Float64,Float64}[]")
    println(io, "    end")
    println(io, "    push!(event_map[t], entry)")
    println(io, "    return event_map")
    println(io, "end")
    println(io)

    println(io, "function event_items_at(event_map::Dict{Float64, Vector{Tuple{Int,Float64,Float64}}}, t::Float64)")
    println(io, "    if haskey(event_map, t)")
    println(io, "        return event_map[t]")
    println(io, "    end")
    println(io, "    for (tt, items) in event_map")
    println(io, "        if isapprox(t, tt; atol = 1e-8, rtol = 0.0)")
    println(io, "            return items")
    println(io, "        end")
    println(io, "    end")
    println(io, "    return nothing")
    println(io, "end")
    println(io)

    println(io, "function regimen_to_event_map(regimen::MosunRegimen)")
    println(io, "    event_map = Dict{Float64, Vector{Tuple{Int,Float64,Float64}}}()")
    println(io, "    for ev in regimen.events")
    println(io, "        idx = dynamic_state_index(ev.target)")
    println(io, "        if ev.rate == 0.0")
    println(io, "            push_event_delta!(event_map, ev.time, (idx, ev.amount, 0.0))")
    println(io, "        else")
    println(io, "            ev.rate > 0.0 || throw(ArgumentError(\"Infusion rates must be positive\"))")
    println(io, "            duration = ev.amount / ev.rate")
    println(io, "            duration >= 0.0 || throw(ArgumentError(\"Infusion amount/rate produced negative duration\"))")
    println(io, "            t_end = ev.time + duration")
    println(io, "            push_event_delta!(event_map, ev.time, (idx, 0.0, ev.rate))")
    println(io, "            push_event_delta!(event_map, t_end, (idx, 0.0, -ev.rate))")
    println(io, "        end")
    println(io, "    end")
    println(io, "    return event_map")
    println(io, "end")
    println(io)

    println(io, "function apply_event_deltas!(u::Vector{Float64}, active_rates::Vector{Float64}, event_map::Dict{Float64, Vector{Tuple{Int,Float64,Float64}}}, t::Float64)")
    println(io, "    items = event_items_at(event_map, t)")
    println(io, "    if items !== nothing")
    println(io, "        for (idx, amt_delta, rate_delta) in items")
    println(io, "            u[idx] += amt_delta")
    println(io, "            active_rates[idx] += rate_delta")
    println(io, "        end")
    println(io, "    end")
    println(io, "    return u, active_rates")
    println(io, "end")
    println(io)

    println(io, "function build_problem(regimen::MosunRegimen, p::MosunParams; tspan::Tuple{Float64,Float64} = (0.0, 84.0), saveat = nothing, callback_mode::Symbol = :segmented, post_event_proposed_dt = nothing)")
    println(io, "    callback_mode in (:segmented, :callback) || throw(ArgumentError(\"Unsupported callback_mode=\$(callback_mode)\"))")
    println(io, "    event_map = regimen_to_event_map(regimen)")
    println(io, "    u0 = pack_state(initial_state(p))")
    println(io, "    active_rates = zeros(Float64, DYNAMIC_STATE_COUNT)")
    println(io, "    t0 = tspan[1]")
    println(io, "    if event_items_at(event_map, t0) !== nothing")
    println(io, "        apply_event_deltas!(u0, active_rates, event_map, t0)")
    println(io, "        delete!(event_map, t0)")
    println(io, "    end")
    println(io, "    ctx = MosunProblemContext(params = p, cache = zero_observables_cache(), active_rates = active_rates)")
    println(io, "    prob = ODEProblem(mosun_rhs!, u0, tspan, ctx)")
    println(io, "    cb_times = sort(collect(keys(event_map)))")
    println(io, "    callback = nothing")
    println(io, "    tstops = copy(cb_times)")
    println(io, "    d_discontinuities = copy(cb_times)")
    println(io, "    if !isempty(cb_times)")
    println(io, "        function affect!(integrator)")
    println(io, "            apply_event_deltas!(integrator.u, integrator.p.active_rates, event_map, Float64(integrator.t))")
    println(io, "            if !isnothing(post_event_proposed_dt)")
    println(io, "                SciMLBase.set_proposed_dt!(integrator, post_event_proposed_dt)")
    println(io, "            end")
    println(io, "        end")
    println(io, "        callback = PresetTimeCallback(cb_times, affect!; save_positions = (false, false))")
    println(io, "    end")
    println(io, "    return MosunBuiltProblem(prob = prob, ctx = ctx, initial_u = u0, callback = callback, tstops = tstops, d_discontinuities = d_discontinuities, event_map = event_map, saveat = saveat, regimen = regimen, callback_mode = callback_mode == :segmented ? :callback : callback_mode)")
    println(io, "end")
    println(io)

    println(io, "function solve_problem(built::MosunBuiltProblem, alg; abstol = 1e-8, reltol = 1e-5, kwargs...)")
    println(io, "    if isnothing(built.saveat)")
    println(io, "        return solve(built.prob, alg; abstol = abstol, reltol = reltol, callback = built.callback, tstops = built.tstops, d_discontinuities = built.d_discontinuities, kwargs...)")
    println(io, "    end")
    println(io, "    return solve(built.prob, alg; abstol = abstol, reltol = reltol, callback = built.callback, tstops = built.tstops, d_discontinuities = built.d_discontinuities, saveat = built.saveat, kwargs...)")
    println(io, "end")
    println(io)

    println(io, "function solve_regimen(regimen::MosunRegimen, p::MosunParams, alg; tspan::Tuple{Float64,Float64} = (0.0, 84.0), saveat = nothing, callback_mode::Symbol = :callback, post_event_proposed_dt = nothing, abstol = 1e-8, reltol = 1e-5, kwargs...)")
    println(io, "    built = build_problem(regimen, p; tspan = tspan, saveat = saveat, callback_mode = callback_mode, post_event_proposed_dt = post_event_proposed_dt)")
    println(io, "    sol = solve_problem(built, alg; abstol = abstol, reltol = reltol, kwargs...)")
    println(io, "    return built, sol")
    println(io, "end")
    println(io)

    println(io, "export MosunParams, MosunDynamicState, MosunObservablesCache, MosunRegimenEvent, MosunRegimen, MosunProblemContext, MosunBuiltProblem,")
    println(io, "    DYNAMIC_STATE_NAMES, OBSERVABLE_NAMES, RHS_OBSERVABLE_NAMES, PARAMETER_NAMES, DEAD_LEGACY_NAMES, DYNAMIC_STATE_COUNT, OBSERVABLE_COUNT, RHS_OBSERVABLE_COUNT,")
    println(io, "    default_params, zero_observables_cache, has_parameter, set_param!, params_from_named_values, params_from_dict, bolus_regimen,")
    println(io, "    initial_state, pack_state, pack_params, unpack_state, update_observables!, value_at, mosun_rhs!, mosun_rhs_vector!, build_problem, solve_problem, solve_regimen,")
    println(io, "    observable, dynamic_state_value, state_or_observable, dynamic_state_index, regimen_to_event_map, apply_event_deltas!")
    println(io)
    println(io, "end # module")
end

artifacts = TCellEngagerQSP.load_artifacts()
open(OUTPUT_PATH, "w") do io
    write_module(io, artifacts)
end
@printf("Wrote %s\n", OUTPUT_PATH)
