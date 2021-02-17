using Test
import Ecophys
import Unitful
import Unitful: m, dimension, mol, μmol, s, W, K, mmol, Pa
PH = Ecophys.Photosynthesis
using Roots

#############################################################################################################
################################################# C3 leaves #################################################
#############################################################################################################

# Without units
opt = PH.SimpleOptical()
c3_f = PH.C3F()
gb_f = PH.simplegbF()
ws_f = 1.0
Tl_f = 300.0
Ta_f = 298.0
P_f  = 101325.0
PAR_f = 1000.0 # μmol/m^2/s -> 219 W/m^2
NIR_f = 219.0  # W/m^2
RH_f = 0.8
Ca_f = 400.0 # μmol/mol
O2_f = 210e3 # μmol/mol
res_f = PH.energybalance(Tl_f, gb_f, c3_f, opt, PAR_f, NIR_f, ws_f, RH_f, Ta_f, Ca_f, P_f, O2_f) # W/m^2
Tl_f = PH.solve_energy_balance(gb_f, c3_f, opt, PAR_f, NIR_f, ws_f, RH_f, Ta_f, Ca_f, P_f, O2_f, xatol = 0.1, order = Order2())

# With units
opt = PH.SimpleOptical()
c3_q = PH.C3Q()
gb_q = PH.simplegbQ()
ws_q = 1.0m/s
Tl_q = 300.0K
Ta_q = 298.0K
P_q = 101325.0Pa
PAR_q = 1000.0μmol/m^2/s
NIR_q = 219.0W/m^2
RH_q = 0.8
Ca_q = 400.0μmol/mol
O2_q = 210.0mmol/mol
res_q = PH.energybalance(Tl_q, gb_q, c3_q, opt, PAR_q, NIR_q, ws_q, RH_q, Ta_q, Ca_q, P_q, O2_q) # W/m^2
@test res_q == res_f*W/m^2
Tl_q = PH.solve_energy_balance(gb_q, c3_q, opt, PAR_q, NIR_q, ws_q, RH_q, Ta_q, Ca_q, P_q, O2_q, xatol = 0.01, order = Order2())
@test abs(Tl_q - Tl_f*K) < 1e-4K

# With effect of angles (horizontal)
gb_q = PH.gbAngleQ()
Tl_q = 300.0K
res_q = PH.energybalance(Tl_q, gb_q, c3_q, opt, PAR_q, NIR_q, ws_q, RH_q, Ta_q, Ca_q, P_q, O2_q) # W/m^2
@test res_q == res_f*W/m^2
Tl_q = PH.solve_energy_balance(gb_q, c3_q, opt, PAR_q, NIR_q, ws_q, RH_q, Ta_q, Ca_q, P_q, O2_q, xatol = 0.01, order = Order2())
@test abs(Tl_q - Tl_f*K) < 1e-4K

# With effect of angles (45°)
gb_q = PH.gbAngleQ(ang = 45.0)
Tl_q = 300.0K
res_q = PH.energybalance(Tl_q, gb_q, c3_q, opt, PAR_q, NIR_q, ws_q, RH_q, Ta_q, Ca_q, P_q, O2_q) # W/m^2
@test abs(res_q - 164.3657W/m^2) < 1e-4W/m^2
Tl_q = PH.solve_energy_balance(gb_q, c3_q, opt, PAR_q, NIR_q, ws_q, RH_q, Ta_q, Ca_q, P_q, O2_q, xatol = 0.01, order = Order2())
@test abs(Tl_q - 305.9689K) < 1e-4K

#############################################################################################################
################################################# C4 leaves #################################################
#############################################################################################################

# Without units
c4_f = PH.C4F()
Tl_f = 300.0
res_f = PH.energybalance(Tl_f, gb_f, c4_f, opt, PAR_f, NIR_f, ws_f, RH_f, Ta_f, Ca_f, P_f, O2_f) # W/m^2
Tl_f = PH.solve_energy_balance(gb_f, c4_f, opt, PAR_f, NIR_f, ws_f, RH_f, Ta_f, Ca_f, P_f, O2_f, xatol = 0.01, order = Order2())

# With units
c4_q = PH.C4Q()
Tl_q = 300.0K
gb_q = PH.simplegbQ()
res_q = PH.energybalance(Tl_q, gb_q, c4_q, opt, PAR_q, NIR_q,ws_q, RH_q, Ta_q, Ca_q, P_q, O2_q) # W/m^2
@test res_q == res_f*W/m^2
Tl_q = PH.solve_energy_balance(gb_q, c4_q, opt, PAR_q, NIR_q, ws_q, RH_q, Ta_q, Ca_q, P_q, O2_q, xatol = 0.01, order = Order2())
@test abs(Tl_q - Tl_f*K) < 1e-4K

# With effect of angles (horizontal)
gb_q = PH.gbAngleQ()
Tl_q = 300.0K
res_q = PH.energybalance(Tl_q, gb_q, c4_q, opt, PAR_q, NIR_q, ws_q, RH_q, Ta_q, Ca_q, P_q, O2_q) # W/m^2
@test res_q == res_f*W/m^2
Tl_q = PH.solve_energy_balance(gb_q, c4_q, opt, PAR_q, NIR_q, ws_q, RH_q, Ta_q, Ca_q, P_q, O2_q, xatol = 0.01, order = Order2())
@test abs(Tl_q - Tl_f*K) < 1e-4K

# With effect of angles (45°)
gb_q = PH.gbAngleQ(ang = 45.0)
Tl_q = 300.0K
res_q = PH.energybalance(Tl_q, gb_q, c4_q, opt, PAR_q, NIR_q, ws_q, RH_q, Ta_q, Ca_q, P_q, O2_q) # W/m^2
@test abs(res_q - 164.5309W/m^2) < 1e-4W/m^2
Tl_q = PH.solve_energy_balance(gb_q, c4_q, opt, PAR_q, NIR_q, ws_q, RH_q, Ta_q, Ca_q, P_q, O2_q, xatol = 0.01, order = Order2())
@test abs(Tl_q - 305.9749K) < 1e-4K