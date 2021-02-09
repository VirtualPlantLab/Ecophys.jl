using Test
import Ecophys
import Unitful
import Unitful: K, J, K, mol, kPa, mbar, μbar, m, μmol, s, mmol
PH = Ecophys.Photosynthesis

# C3F - Photosynthesis without units
c3_f = PH.C3F()
PAR_f = 1000 # μmol/m^2/s
RH_f = 0.8
Tleaf_f = 298.15 # K
Ca_f = 400.0 # μmol/mol
O2_f = 210e3 # μmol/mol
gb_f = 1.0 # mol/m^2/s
A_f, gs_f  = PH.A_gs(c3_f, PAR_f, RH_f, Tleaf_f, Ca_f, O2_f, gb_f)
@test A_f ≈ 21.73243935392202
@test gs_f ≈ 0.0022225283032597217

# C3F - Photosynthesis with units
c3_q = PH.C3Q()
PAR_q = 1000μmol/m^2/s
RH_q = 0.8
Tleaf_q = 298.15K
Ca_q = 400.0μmol/mol
O2_q = 210.0mmol/mol
gb_q = 1.0mol/m^2/s
A_q, gs_q = PH.A_gs(c3_q, PAR_q, RH_q, Tleaf_q, Ca_q, O2_q, gb_q)

# Check the same results with and without units
@test A_q ≈ A_f*μmol/m^2/s
@test gs_q ≈ gs_f*mol/m^2/s

