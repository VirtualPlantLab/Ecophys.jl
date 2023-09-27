using Ecophys
using Test
using Documenter
import Aqua

@testset "Photosynthesis" begin

    # Test examples on documentation (jldoctest blocks)
    DocMeta.setdocmeta!(Ecophys,
        :DocTestSetup,
        :(using Ecophys);
        recursive = true)
    doctest(Ecophys)

    # Aqua
    @testset "Aqua" begin
        Aqua.test_all(Ecophys, ambiguities = false, project_extras = false)
        Aqua.test_ambiguities([Ecophys])
    end

    # Unit tests
    include("Photosynthesis/test_components.jl")
    include("Photosynthesis/test_C3.jl")
    include("Photosynthesis/test_C4.jl")
    include("Photosynthesis/test_gb.jl")
    include("Photosynthesis/test_energybalance.jl")
end
