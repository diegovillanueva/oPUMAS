using Test
using oPUMAS

include("golden.jl")
include("processes.jl")

@testset "oPUMAS" begin
    include("test_foundations.jl")
    include("test_size_distributions.jl")
    include("test_warm_rain.jl")
    include("test_ice_processes.jl")
    include("test_step.jl")
    include("test_time_loop.jl")
    include("test_physics.jl")
end

println("🔬 worst relative error per routine")
foreach(((name, err),) -> println(rpad(name, 40), err), REPORT)
