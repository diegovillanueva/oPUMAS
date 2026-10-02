using Test
using oPUMAS

include("golden.jl")
include("processes.jl")

@testset "oPUMAS" begin
    include("test_foundations.jl")
    include("test_size_distributions.jl")
    include("test_warm_rain.jl")
end

println("🔬 worst relative error per routine")
foreach(((name, err),) -> println(rpad(name, 40), err), REPORT)
