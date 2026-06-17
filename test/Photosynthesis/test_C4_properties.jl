
using Test
import Ecophys
PH = Ecophys.Photosynthesis

@testset "C4 photosynthesis properties" begin

    @testset "Monotonic A with Ca" begin
        c4 = PH.C4()
        A_vals = [PH.photosynthesis(c4, Ca = ca).A for ca in 10.0:50.0:2000.0]
        @test all(diff(A_vals) .>= 0)
    end

    @testset "Monotonic A with PAR" begin
        c4 = PH.C4()
        A_vals = [PH.photosynthesis(c4, PAR = par).A for par in 0.0:50.0:2000.0]
        @test all(diff(A_vals) .>= 0)
    end

    @testset "gs increases with RH" begin
        c4 = PH.C4()
        gs_vals = [PH.photosynthesis(c4, RH = rh).gs for rh in 0.10:0.05:0.95]
        @test all(diff(gs_vals) .>= 0)
    end

    @testset "Peaked temperature response (10-40 C)" begin
        c4 = PH.C4()
        A_vals = [PH.photosynthesis(c4, Tleaf = 273.15 + t).A for t in 10.0:1.0:40.0]
        peak_idx = argmax(A_vals)
        @test peak_idx > 1 && peak_idx < length(A_vals)
    end

    @testset "Dark respiration at PAR = 0" begin
        @test PH.photosynthesis(PH.C4(), PAR = 0.0).A < 0
    end

    @testset "CO2 compensation point" begin
        c4 = PH.C4()
        @test PH.photosynthesis(c4, Ca = 1.0).A < 0
        @test PH.photosynthesis(c4, Ca = 400.0).A > 0
    end

    @testset "gs >= gso" begin
        c4 = PH.C4()
        results = [PH.photosynthesis(c4, PAR = par, Ca = ca)
                   for par in 100.0:100.0:2000.0, ca in 200.0:100.0:2000.0]
        @test all(r.gs >= c4.gso for r in results)
    end

end
