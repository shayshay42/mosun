module MosunModelCoreSupport

using CSV
using DataFrames
using DifferentialEquations
using Sundials

include(joinpath(@__DIR__, "MosunModelCore.jl"))
using .MosunModelCore

const MMC = MosunModelCore
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MODEL_TABLES_DIR = joinpath(REPO_ROOT, "generated", "model_tables")
const RECON_BUNDLE_DIR = joinpath(REPO_ROOT, "generated", "model_reconstruction_bundle")
const SOLVER_ABSTOL = parse(Float64, get(ENV, "TCE_ABSTOL", "1e-8"))
const SOLVER_RELTOL = parse(Float64, get(ENV, "TCE_RELTOL", "1e-5"))

function effective_model_tables_dir()
    if isfile(joinpath(RECON_BUNDLE_DIR, "doses.tsv")) && isfile(joinpath(RECON_BUNDLE_DIR, "variants.tsv"))
        return RECON_BUNDLE_DIR
    end
    return MODEL_TABLES_DIR
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

parse_float_or_nan(x) = try
    parse(Float64, strip(String(x)))
catch
    NaN
end

function parse_number_vector(txt)
    s = strip(String(txt))
    isempty(s) && return Float64[]
    s = replace(s, '[' => ' ', ']' => ' ', ';' => ' ', '\'' => ' ')
    toks = split(s, [',', ' ', '\t', '\n', '\r']; keepempty = false)
    vals = Float64[]
    for tok in toks
        push!(vals, parse(Float64, tok))
    end
    return vals
end

function parse_first_number(txt; default::Float64 = 0.0)
    vals = parse_number_vector(txt)
    return isempty(vals) ? default : vals[1]
end

function normalize_vector_length(v::Vector{Float64}, n::Int, default::Float64)
    if isempty(v)
        return fill(default, n)
    elseif length(v) == n
        return v
    elseif length(v) == 1
        return fill(v[1], n)
    elseif length(v) > n
        return v[1:n]
    else
        out = Vector{Float64}(undef, n)
        for i in 1:n
            out[i] = i <= length(v) ? v[i] : v[end]
        end
        return out
    end
end

function load_variant_overrides(; model_tables_dir::String = effective_model_tables_dir())
    variants_df = DataFrame(CSV.File(joinpath(model_tables_dir, "variants.tsv"); delim = '\t'))
    variants = Dict{Int, Vector{Pair{String, Float64}}}()
    for r in eachrow(variants_df)
        idx = hasproperty(r, :variant_idx) ? Int(r.variant_idx) : Int(r.idx)
        action = String(r.action)
        field_name = String(r.name)
        value = hasproperty(r, :value_numeric) ?
            (r.value_numeric === missing ? NaN : Float64(r.value_numeric)) :
            parse_float_or_nan(r.value)
        if action == "parameter" && field_name == "Value" && isfinite(value)
            push!(get!(variants, idx, Pair{String, Float64}[]), String(r.class) => value)
        end
    end
    return variants
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

function apply_variant_overrides!(p::MMC.MosunParams, variant_ids::Vector{Int}, variant_overrides::Dict{Int, Vector{Pair{String, Float64}}})
    for vid in variant_ids
        for (name, value) in get(variant_overrides, vid, Pair{String, Float64}[])
            if MMC.has_parameter(name)
                MMC.set_param!(p, name, value)
            end
        end
    end
    return p
end

function build_params_with_variants(param_names, param_values;
        variant_ids::Vector{Int} = Int[],
        variant_overrides::Dict{Int, Vector{Pair{String, Float64}}} = load_variant_overrides(),
        base::MMC.MosunParams = MMC.default_params(),
        ignore_unknown::Bool = true)
    p = deepcopy(base)
    apply_variant_overrides!(p, variant_ids, variant_overrides)
    return MMC.params_from_named_values(param_names, param_values; base = p, ignore_unknown = ignore_unknown)
end

function dose_row_to_entry(row)
    dclass = String(row.class)
    target = hasproperty(row, :target_name) ? String(row.target_name) : String(row.target)
    times = Float64[]
    if occursin("RepeatDose", dclass)
        startv = hasproperty(row, :start_time_json) ? parse_first_number(row.start_time_json; default = 0.0) : parse_first_number(row.starttime; default = 0.0)
        intervalv = hasproperty(row, :interval_json) ? parse_first_number(row.interval_json; default = 0.0) : parse_first_number(row.interval; default = 0.0)
        repeatcountv = hasproperty(row, :repeat_count_json) ? parse_first_number(row.repeat_count_json; default = 0.0) : parse_first_number(row.repeatcount; default = 0.0)
        n_doses = max(1, Int(round(repeatcountv)) + 1)
        times = [startv + intervalv * i for i in 0:(n_doses - 1)]
    else
        times = hasproperty(row, :time_json) ? parse_number_vector(row.time_json) : parse_number_vector(row.time)
        if isempty(times)
            startv = hasproperty(row, :start_time_json) ? parse_first_number(row.start_time_json; default = 0.0) : parse_first_number(row.starttime; default = 0.0)
            times = [startv]
        end
    end

    amount_vec = hasproperty(row, :amount_json) ? parse_number_vector(row.amount_json) : parse_number_vector(row.amount)
    rate_vec = hasproperty(row, :rate_json) ? parse_number_vector(row.rate_json) : parse_number_vector(row.rate)
    amounts = normalize_vector_length(amount_vec, length(times), 0.0)
    rates = normalize_vector_length(rate_vec, length(times), 0.0)
    return (target = target, times = times, amount = amounts, rate = rates)
end

function load_dose_catalog(; model_tables_dir::String = effective_model_tables_dir())
    doses_df = DataFrame(CSV.File(joinpath(model_tables_dir, "doses.tsv"); delim = '\t'))
    catalog = Dict{String, NamedTuple{(:target, :times, :amount, :rate), Tuple{String, Vector{Float64}, Vector{Float64}, Vector{Float64}}}}()
    for row in eachrow(doses_df)
        catalog[String(row.name)] = dose_row_to_entry(row)
    end
    return catalog
end

function regimen_from_dose_catalog(selected_doses::Vector{String}, dose_catalog::Dict{String, <:NamedTuple})
    events = MMC.MosunRegimenEvent[]
    for dname in selected_doses
        haskey(dose_catalog, dname) || error("Dose '$dname' not found in dose catalog")
        ddef = dose_catalog[dname]
        target = Symbol(ddef.target)
        for i in eachindex(ddef.times)
            amt = ddef.amount[i]
            rate = ddef.rate[i]
            if !iszero(amt) || !iszero(rate)
                push!(events, MMC.MosunRegimenEvent(target = target, time = ddef.times[i], amount = amt, rate = rate))
            end
        end
    end
    return MMC.MosunRegimen(events)
end

function bolus_regimen_from_mg(dose_days::Vector{Float64}, dose_mg::Vector{Float64}; bw_kg::Float64 = 70.0, target::Symbol = :TDBc_ugperkg)
    length(dose_days) == length(dose_mg) || throw(ArgumentError("dose_days and dose_mg lengths differ"))
    events = MMC.MosunRegimenEvent[]
    for (t, d) in zip(dose_days, dose_mg)
        amt = d * 1000.0 / bw_kg
        iszero(amt) && continue
        push!(events, MMC.MosunRegimenEvent(target = target, time = t, amount = amt, rate = 0.0))
    end
    return MMC.MosunRegimen(events)
end

function solution_observable(sol, params::MMC.MosunParams, name::Symbol)
    cache = MMC.zero_observables_cache()
    vals = Vector{Float64}(undef, length(sol.t))
    for i in eachindex(sol.t)
        t = Float64(sol.t[i])
        u = sol.u[i]
        MMC.update_observables!(cache, u, params, t)
        vals[i] = Float64(MMC.state_or_observable(u, cache, name))
    end
    return vals
end

function dense_save_times(sim_time::Float64, base_times::Vector{Float64}, save_dt::Float64)
    dense = collect(0.0:save_dt:sim_time)
    all_times = unique(vcat(base_times, dense, [sim_time]))
    sort!(all_times)
    return all_times
end

export MMC, REPO_ROOT, MODEL_TABLES_DIR, RECON_BUNDLE_DIR, SOLVER_ABSTOL, SOLVER_RELTOL,
    make_solver_alg, parse_number_vector, parse_first_number, normalize_vector_length,
    load_variant_overrides, case_variant_ids, apply_variant_overrides!, build_params_with_variants,
    load_dose_catalog, regimen_from_dose_catalog, bolus_regimen_from_mg, effective_model_tables_dir,
    solution_observable, dense_save_times

end # module
