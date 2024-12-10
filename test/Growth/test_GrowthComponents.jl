using Test, Ecophys
GR = Ecophys.Growth

let

    # Component values related to growth respiration given the composition of new tissue
    f = [0.53, 0.25, 0.05, 0.05, 0.06, 0.06]
    ASRQ, CO2PFF, CF = GR.compute_ASRQ(f)
    @testset "Accuracy of compute_ASRQ" begin
        @test abs(ASRQ - 1.656) < 1e-4
        @test abs(CO2PFF - 0.71502) < 1e-4
        @test abs(CF - 0.467373) < 1e-4
    end
    
end