using Test

include(joinpath(@__DIR__, "..", "src", "TCellEngagerQSP.jl"))
using .TCellEngagerQSP

include(joinpath(@__DIR__, "reference_regression.jl"))
include(joinpath(@__DIR__, "production_core_regression.jl"))
