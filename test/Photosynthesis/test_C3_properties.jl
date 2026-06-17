
using Test
import Ecophys
PH = Ecophys.Photosynthesis

@testset "C3 photosynthesis properties" begin

    @testset "Monotonic A with Ca" begin
        c3 = PH.C3()
        A_vals = [PH.photosynthesis(c3, Ca = ca).A for ca in 50.0:50.0:2000.0]
        @test all(diff(A_vals) .>= 0)
    end

    @testset "Monotonic A with PAR" begin
        c3 = PH.C3()
        A_vals = [PH.photosynthesis(c3, PAR = par).A for par in 0.0:50.0:2000.0]
        @test all(diff(A_vals) .>= 0)
    end

    @testset "gs increases with RH" begin
        c3 = PH.C3()
        gs_vals = [PH.photosynthesis(c3, RH = rh).gs for rh in 0.10:0.05:0.95]
        @test all(diff(gs_vals) .>= 0)
    end

    @testset "Peaked temperature response (10-40 C)" begin
        c3 = PH.C3()
        A_vals = [PH.photosynthesis(c3, Tleaf = 273.15 + t).A for t in 10.0:1.0:40.0]
        peak_idx = argmax(A_vals)
        @test peak_idx > 1 && peak_idx < length(A_vals)
    end

    @testset "Dark respiration at PAR = 0" begin
        @test PH.photosynthesis(PH.C3(), PAR = 0.0).A < 0
    end

    @testset "CO2 compensation point" begin
        c3 = PH.C3()
        @test PH.photosynthesis(c3, Ca = 20.0).A < 0
        @test PH.photosynthesis(c3, Ca = 400.0).A > 0
    end

    @testset "gs >= gso" begin
        c3 = PH.C3()
        results = [PH.photosynthesis(c3, PAR = par, Ca = ca)
                   for par in 100.0:100.0:2000.0, ca in 200.0:100.0:2000.0]
        @test all(r.gs >= c3.gso for r in results)
    end

    @testset "O2 inhibition" begin
        c3 = PH.C3()
        A_vals = [PH.photosynthesis(c3, O2 = o2).A for o2 in 10e3:10e3:210e3]
        @test all(diff(A_vals) .<= 0)
    end

end
