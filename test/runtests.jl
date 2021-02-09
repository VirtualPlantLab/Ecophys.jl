using Ecophys
using Test

@testset "Photosynthesis" begin
    include("Photosynthesis/test_components.jl")
    include("Photosynthesis/test_C3.jl")
end
