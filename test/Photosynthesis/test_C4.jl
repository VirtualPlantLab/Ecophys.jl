
using Test
import Ecophys
import Unitful
import Unitful: K, J, K, mol, kPa, mbar, μbar, m, μmol, s, mmol
PH = Ecophys.Photosynthesis

let

    # C4 - Photosynthesis without units
    c4 = PH.C4()
    PAR_f = 1000 # μmol/m^2/s
    RH_f = 0.8
    Tleaf_f = 298.15 # K
    Ca_f = 400.0 # μmol/mol
    O2_f = 210e3 # μmol/mol
    gb_f = 1.0 # mol/m^2/s
    A_f, gs_f = PH.photosynthesis(c4,
        PAR = PAR_f,
        RH = RH_f,
        Tleaf = Tleaf_f,
        Ca = Ca_f,
        O2 = O2_f,
        gb = gb_f)
    @test abs(A_f - 30.5477) < 0.0001
    @test abs(gs_f - 0.0018) < 1e-4
    Ag_f, gs_f = PH.photosynthesis(c4,
        PAR = PAR_f,
        RH = RH_f,
        Tleaf = Tleaf_f,
        Ca = Ca_f,
        O2 = O2_f,
        gb = gb_f,
        net = false)
    @test Ag_f > A_f

    # C4Q - Photosynthesis without units
    c4_q = PH.C4Q()
    PAR_q = 1000μmol / m^2 / s
    RH_q = 0.8
    Tleaf_q = 298.15K
    Ca_q = 400.0μmol / mol
    O2_q = 210e3μmol / mol
    gb_q = 1.0mol / m^2 / s
    A_q, gs_q = PH.photosynthesis(c4_q,
        PAR = PAR_q,
        RH = RH_q,
        Tleaf = Tleaf_q,
        Ca = Ca_q,
        O2 = O2_q,
        gb = gb_q)
    Ag_q, gs_q = PH.photosynthesis(c4_q,
        PAR = PAR_q,
        RH = RH_q,
        Tleaf = Tleaf_q,
        Ca = Ca_q,
        O2 = O2_q,
        gb = gb_q,
        net = false)
    @test Ag_q > A_q

    # Check the same results with and without units
    @test A_q ≈ A_f * μmol / m^2 / s
    @test Ag_q ≈ Ag_f * μmol / m^2 / s
    @test gs_q ≈ gs_f * mol / m^2 / s

    # Temperature response curves
    Ag = [PH.photosynthesis(C4(), Tleaf = 273.15 + x).A for x in 0:1:40.0]
    @test all(Ag .> 0.0)

    # Light response curves
    Ag = [PH.photosynthesis(C4(), PAR = x, net = false).A for x in 50.0:50.0:2e3]
    @test all(Ag .> 0.0)

    # CO2 response curves
    Ag = [PH.photosynthesis(C4(), Ca = x, net = false).A for x in 50.0:50.0:2e3]
    @test all(Ag .> 0.0)
end
