# MosunModelCore Refactor And Validation

## Summary

`MosunModelCore` is now the default production model path for the Julia implementation.

The legacy/reference implementations remain in the repo for two specific purposes:

1. reproducing the original shipped MATLAB-compatible reference behavior
2. auditing parity against the older Julia pathways

The public default is controlled in [TCellEngagerQSP.jl](../julia/src/TCellEngagerQSP.jl):

- `DEFAULT_CORE_MODE = :production`
- override with `TCE_CORE_MODE=legacy_reference` or explicit `core_mode_override = :legacy_reference` only when reference reproduction is required

## What Changed

The main production changes were:

1. Replaced runtime formula/code generation on the hot solve path with explicit checked-in Julia code in [MosunModelCore.jl](../julia/src/MosunModelCore.jl)
2. Separated the model into:
   - `36` dynamic ODE states
   - `89` observables/repeated-assignment quantities
   - typed parameter/state/regimen/cache structs
3. Made production PK use the algebraic central concentration path
4. Moved callback/event dosing into the production core with `PresetTimeCallback`
5. Made callback dose amounts compatible with forward AD
6. Kept the legacy/reference path explicit and opt-in inside [TCellEngagerQSP.jl](../julia/src/TCellEngagerQSP.jl)

## Refactor Process

The implementation sequence was:

1. Remove runtime `invokelatest`/`Core.eval` dependence from the canonical path.
2. Replace runtime-generated equations with direct Julia source.
3. Introduce the standalone `MosunModelCore` module.
4. Reduce the runtime state representation from the old mixed `131`-slot layout to the production `36`-state layout.
5. Keep a compatibility adapter in `TCellEngagerQSP` so old scripts could continue to run while the production path moved to `MosunModelCore`.
6. Split production observables from RHS-only algebra so the hot path does not have to materialize the full observable cache for every solve step.
7. Make callback dosing AD-compatible for dose-valued optimization variables.

## Default Runtime Behavior

The repo now defaults to production:

- [build_model_with_variant_ids](../julia/src/TCellEngagerQSP.jl)
- [build_model](../julia/src/TCellEngagerQSP.jl)
- [run_simulation](../julia/src/TCellEngagerQSP.jl)
- [run_manifest_row](../julia/src/TCellEngagerQSP.jl)

all default to `core_mode_override = DEFAULT_CORE_MODE`, and `DEFAULT_CORE_MODE` resolves to `:production` unless `TCE_CORE_MODE` is set.

Use the legacy/reference implementation only when needed:

```julia
mdl = TCellEngagerQSP.build_model_with_variant_ids(variant_ids, pnames, pvals; core_mode_override = :legacy_reference)
```

or

```powershell
$env:TCE_CORE_MODE = "legacy_reference"
```

## Validation Coverage

Validation is split into three layers:

1. MATLAB reference regression
2. production-vs-legacy parity checks
3. solver and AD benchmark checks

### 1. MATLAB Reference Regression

Primary test file:

- [reference_regression.jl](../julia/test/reference_regression.jl)

What it checks:

1. The model still contains `116` reactions.
2. The model still contains `101` total rules (`initial + repeated`).
3. The explicit Julia/manual reference path still matches the saved MATLAB three-case manifest on the original output grids.

Command:

```powershell
julia --project=./julia julia/test/runtests.jl
```

Reference artifact written by the test:

- `generated/julia_reference/manual_core_reference_quality.csv`

This is intentionally run with `core_mode_override = :legacy_reference`, because two of the shipped MATLAB reference cases depend on the old lookup PK branch rather than the production PK simplification.

### 2. Production Core Regression

Primary test file:

- [production_core_regression.jl](../julia/test/production_core_regression.jl)

What it checks:

1. `MosunModelCore` has no runtime dependency on:
   - `CSV`
   - `DataFrames`
   - `JSON3`
   - `ModelingToolkit`
   - `RuntimeGeneratedFunctions`
2. The production layout is:
   - `36` dynamic states
   - `89` observables
   - `58` RHS-needed observables
   - `6` dead legacy names removed from the production API
3. Production initial state matches the legacy/reference state values for shared state names.
4. Representative observables match the legacy/reference repeated-assignment values at fixed probe states.
5. Production case-30 canonical output matches both:
   - legacy/reference canonical output
   - saved MATLAB reference output
6. A phase-1 workflow regression matches legacy/reference trajectories for:
   - `IL6combo`
   - `Btumor`
   - best `%SPD`
7. ForwardDiff works on the production path for:
   - parameter gradients
   - callback-dose gradients

Command:

```powershell
julia --project=./julia julia/test/runtests.jl
```

### 3. Solver Benchmarks

Benchmark scripts:

- [benchmark_explicit_core_solver_sweep.jl](../julia/benchmark_explicit_core_solver_sweep.jl)
- [benchmark_mosun_model_core_gradient_modes.jl](../julia/benchmark_mosun_model_core_gradient_modes.jl)

Summary markdown:

- `generated/benchmarks/solver_benchmark_summary_20260310.md`

Main forward-solve conclusion from the current benchmark:

1. Fastest plain solve: `CVODE_BDF` with post-event proposed `dt = 1e-2`
2. Best pure-Julia forward-solve fallback: `QNDF`
3. `ModelingToolkit` dense/sparse Jacobian paths were not competitive for this workload

Plain solve artifact directory:

- `generated/benchmarks/explicit_core_solver_sweep_20260310`

### 4. Gradient Benchmarks

Gradient benchmark script:

- [benchmark_mosun_model_core_gradient_modes.jl](../julia/benchmark_mosun_model_core_gradient_modes.jl)

Current conclusion:

1. Forward-mode AD:
   - best speed: `Tsit5`
   - production callback dose AD path is validated in the unit tests
2. Reverse-mode AD:
   - best benchmarked adjoint runtime: `Rodas4P` with `Enzyme` adjoint
   - `QNDF` is the main alternative
3. ReverseDiff-based direct reverse optimization paths are still not the stable default for the long objective stack benchmark

Gradient artifact directory:

- `generated/benchmarks/mosun_core_gradient_modes_20260310`

## Known Scope Boundary

`MosunModelCore` is the production model, not the full historical-reference package.

What remains legacy/reference-only:

1. `PK_v26` lookup behavior from the shipped MATLAB supplementary code
2. MTK Jacobian helper paths inside the compatibility adapter
3. old `131`-slot mixed layout utilities

That split is deliberate:

- production solves should be plain Julia, typed, explicit, and AD-friendly
- MATLAB-paper reproduction should stay opt-in and auditable

## Recommended Usage

Production:

```julia
using .TCellEngagerQSP
const MMC = TCellEngagerQSP.MosunModelCore

p = MMC.default_params()
regimen = MMC.bolus_regimen(:TDBc_ugperkg, Dict(0.0 => 1.0))
built = MMC.build_problem(regimen, p; tspan = (0.0, 84.0), callback_mode = :callback)
sol = MMC.solve_problem(built, TCellEngagerQSP.make_solver_alg("cvode_bdf"))
```

Reference reproduction:

```julia
times, outvec, X, ref_df = TCellEngagerQSP.run_manifest_row(row; engine = :canonical, core_mode_override = :legacy_reference)
```

## Push Checklist

Before pushing:

1. run `julia --project=./julia julia/test/runtests.jl`
2. if you changed solver infrastructure, rerun the solver sweep
3. if you changed AD-sensitive code, rerun the gradient benchmark
4. use `legacy_reference` explicitly for any script whose goal is historical MATLAB parity rather than production execution
