using Test, Ecophys, Unitful
import Unitful: °C, g
GR = Ecophys.Growth

let

    # Growth models without units
    par = GR.Organ()
    t = 20.0
    te = 100.0
    tm = 50.0
    tb = 10.0
    wmax = 0.1
    c = GR.compute_potential_GR(t, te, tm, tb, wmax)

    # Growth models with units
    par_q = GR.OrganQ()
    t_q = 20.0u"°C"
    te_q = 100.0u"°C"
    tm_q = 50.0u"°C"
    tb_q = 10.0u"°C"
    wmax_q = 0.1u"g"
    c_q = GR.compute_potential_GR(par_q)

    @testset "Accuracy differences between unit and unit-free Growth models" begin
        @test abs(c - 0.0031) < 1e-4
        @test ustrip(abs(c_q - 0.0031u"g/d")) < 1e-4
    end
    
end
