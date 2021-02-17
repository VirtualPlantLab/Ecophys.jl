using Ecophys
using Test

@testset "Photosynthesis" begin
    include("Photosynthesis/test_components.jl")
    include("Photosynthesis/test_C3.jl")
    include("Photosynthesis/test_C4.jl")
    include("Photosynthesis/test_gb.jl")
    include("Photosynthesis/test_energybalance.jl")
end
