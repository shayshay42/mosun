using Printf

include(joinpath(@__DIR__, "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

function emit_expr(expr)
    return sprint(Base.show_unquoted, expr)
end

function write_string_vector(io, name::String, values::Vector{String})
    println(io, "const $name = [")
    for v in values
        println(io, "    ", repr(v), ",")
    end
    println(io, "]")
    println(io)
end

function build_name_layout(artifacts)
    state_names = copy(artifacts.state_names)
    param_names = copy(artifacts.param_names)

    aux_names = String[]
    for (lhs, _) in artifacts.initial_rules
        lhss = TCellEngagerQSP.sanitize_name(lhs)
        if !(lhss in state_names) && !(lhss in param_names) && !(lhss in aux_names)
            push!(aux_names, lhss)
        end
    end
    for (lhs, _) in artifacts.repeated_rules
        lhss = TCellEngagerQSP.sanitize_name(lhs)
        if !(lhss in state_names) && !(lhss in param_names) && !(lhss in aux_names)
            push!(aux_names, lhss)
        end
    end

    all_names = vcat(state_names, param_names, aux_names)
    name_to_idx = Dict(name => i for (i, name) in enumerate(all_names))
    state_to_idx = Dict(name => i for (i, name) in enumerate(state_names))
    return state_names, param_names, aux_names, all_names, name_to_idx, state_to_idx
end

function write_direct_core(out_path::String)
    artifacts = TCellEngagerQSP.load_artifacts()
    state_names, param_names, aux_names, _, name_to_idx, state_to_idx = build_name_layout(artifacts)

    rep_order = TCellEngagerQSP.repeated_rule_order(artifacts.repeated_rules)

    open(out_path, "w") do io
        println(io, "# Candidate manual core generated for audit/diff against the checked-in manual source.")
        println(io, "# This file is not included by runtime code.")
        println(io)
        println(io, "const DIRECT_STATE_COUNT = ", length(state_names))
        println(io, "const DIRECT_PARAM_COUNT = ", length(param_names))
        println(io, "const DIRECT_VALUE_COUNT = ", length(state_names) + length(param_names) + length(aux_names))
        println(io)
        write_string_vector(io, "DIRECT_STATE_NAMES", state_names)
        write_string_vector(io, "DIRECT_PARAM_NAMES", param_names)

        println(io, "function direct_model_layout_matches(state_names::Vector{String}, param_names::Vector{String})")
        println(io, "    return state_names == DIRECT_STATE_NAMES && param_names == DIRECT_PARAM_NAMES")
        println(io, "end")
        println(io)

        println(io, "function apply_initial_rules_direct!(z, u0, pvals)")
        println(io, "    @inbounds begin")
        println(io, "        z[1:DIRECT_STATE_COUNT] .= u0")
        println(io, "        z[(DIRECT_STATE_COUNT + 1):(DIRECT_STATE_COUNT + DIRECT_PARAM_COUNT)] .= pvals")
        for (lhs, rhs) in artifacts.initial_rules
            lhs_s = TCellEngagerQSP.sanitize_name(lhs)
            lhs_idx = name_to_idx[lhs_s]
            rhs_expr = TCellEngagerQSP.compile_formula_expr(rhs, name_to_idx)
            rhs_str = emit_expr(rhs_expr)
            println(io, "        z[$lhs_idx] = Float64(real($rhs_str))")
            if haskey(state_to_idx, lhs_s)
                println(io, "        u0[$(state_to_idx[lhs_s])] = z[$lhs_idx]")
            end
        end
        println(io, "    end")
        println(io, "    return nothing")
        println(io, "end")
        println(io)

        println(io, "function apply_repeated_rules_direct!(z, u, pvals, t)")
        println(io, "    @inbounds begin")
        println(io, "        z[1:DIRECT_STATE_COUNT] .= u")
        println(io, "        z[(DIRECT_STATE_COUNT + 1):(DIRECT_STATE_COUNT + DIRECT_PARAM_COUNT)] .= pvals")
        for idx in rep_order
            lhs, rhs = artifacts.repeated_rules[idx]
            lhs_s = TCellEngagerQSP.sanitize_name(lhs)
            lhs_idx = name_to_idx[lhs_s]
            rhs_expr = TCellEngagerQSP.compile_formula_expr(rhs, name_to_idx)
            rhs_str = emit_expr(rhs_expr)
            println(io, "        z[$lhs_idx] = Float64(real($rhs_str))")
        end
        println(io, "    end")
        println(io, "    return z")
        println(io, "end")
        println(io)

        println(io, "function rhs_direct!(du, u, ctx, t)")
        println(io, "    z = ctx.z")
        println(io, "    pvals = ctx.pvals")
        println(io, "    apply_repeated_rules_direct!(z, u, pvals, t)")
        println(io, "    @inbounds begin")
        println(io, "        fill!(du, 0.0)")
        println(io, "        for inf in ctx.infusions")
        println(io, "            if t >= inf.t_start && t <= inf.t_end")
        println(io, "                du[inf.target_idx] += inf.rate")
        println(io, "            end")
        println(io, "        end")
        for j in eachindex(artifacts.rate_exprs)
            rexpr = TCellEngagerQSP.compile_formula_expr(artifacts.rate_exprs[j], name_to_idx)
            rstr = emit_expr(rexpr)
            println(io, "        rate_$j = Float64(real($rstr))")
            for (i, coeff) in artifacts.stoich[j]
                coefff = Float64(coeff)
                println(io, "        du[$i] += ", repr(coefff), " * rate_$j")
            end
        end
        println(io, "    end")
        println(io, "    return nothing")
        println(io, "end")
        println(io)
    end
end

out_dir = joinpath(TCellEngagerQSP.REPO_ROOT, "generated", "manual_core_candidate")
mkpath(out_dir)
out_path = joinpath(out_dir, "TCellEngagerQSPManualCore.generated.jl")
write_direct_core(out_path)
@printf("Wrote %s\n", out_path)
