using Test
import Ecophys
import Unitful
import Unitful: m, dimension, mol, μmol, s, W, K, mmol, Pa
PH = Ecophys.Photosynthesis
using Roots

#############################################################################################################
################################################# C3 leaves #################################################
#############################################################################################################

#let 

# Without units
opt = PH.SimpleOptical()
c3 = PH.C3()
gb = PH.simplegb()
ws = 1.0
Tl = 300.0
Tair = 298.0
P = 101325.0
PAR = 1000.0 # μmol/m^2/s -> 219 W/m^2
NIR = 219.0  # W/m^2
RH = 0.8
Ca = 400.0 # μmol/mol
O2 = 210e3 # μmol/mol
res = PH.energybalance(Tl, gb, c3, opt, PAR, NIR, ws, RH, Tair, Ca, P, O2) # W/m^2
Tl = PH.solve_energy_balance(c3, gb = gb, opt = opt, PAR = PAR, NIR = NIR,
    ws = ws, RH = RH, Tair = Tair, Ca = Ca, P = P, O2 = O2,
    xatol = 0.01, order = Order2())

# With units
opt = PH.SimpleOptical()
c3_q = PH.C3Q()
gb_q = PH.simplegbQ()
ws_q = ws * m / s
Tl_q = 300.0K
Ta_q = Tair * K
P_q = P * Pa
PAR_q = PAR * μmol / m^2 / s
NIR_q = NIR * W / m^2
RH_q = RH
Ca_q = Ca * μmol / mol
O2_q = O2 * mmol / mol
res_q = PH.energybalance(Tl_q,
    gb_q,
    c3_q,
    opt,
    PAR_q,
    NIR_q,
    ws_q,
    RH_q,
    Ta_q,
    Ca_q,
    P_q,
    O2_q) # W/m^2
@test_broken res_q == res * W / m^2
Tl_q = PH.solve_energy_balance(c3_q, gb = gb_q, opt = opt, PAR = PAR_q,
    NIR = NIR_q, ws = ws_q, RH = RH_q, Tair = Ta_q,
    Ca = Ca_q, P = P_q, O2 = O2_q, xatol = 0.01,
    order = Order2())
@test_broken abs(Tl_q - Tl * K) < 1e-4K

# With effect of angles (horizontal)
gb_q = PH.gbAngleQ()
Tl_q = 300.0K
res_q = PH.energybalance(Tl_q,
    gb_q,
    c3_q,
    opt,
    PAR_q,
    NIR_q,
    ws_q,
    RH_q,
    Ta_q,
    Ca_q,
    P_q,
    O2_q) # W/m^2
@test_broken res_q == res * W / m^2
Tl_q = PH.solve_energy_balance(c3_q, gb = gb_q, opt = opt, PAR = PAR_q,
    NIR = NIR_q, ws = ws_q, RH = RH_q, Tair = Ta_q,
    Ca = Ca_q, P = P_q, O2 = O2_q, xatol = 0.01,
    order = Order2())
@test_broken abs(Tl_q - Tl * K) < 1e-4K

# With effect of angles (45°)
gb_q = PH.gbAngleQ(ang = 45.0)
Tl_q = 300.0K
res_q = PH.energybalance(Tl_q,
    gb_q,
    c3_q,
    opt,
    PAR_q,
    NIR_q,
    ws_q,
    RH_q,
    Ta_q,
    Ca_q,
    P_q,
    O2_q) # W/m^2
@test_broken abs(res_q - 175.3157W / m^2) < 1e-4W / m^2
Tl_q = PH.solve_energy_balance(c3_q, gb = gb_q, opt = opt, PAR = PAR_q,
    NIR = NIR_q, ws = ws_q, RH = RH_q, Tair = Ta_q,
    Ca = Ca_q, P = P_q, O2 = O2_q, xatol = 0.01,
    order = Order2())
@test_broken abs(Tl_q - 306.3603K) < 1e-4K

#############################################################################################################
################################################# C4 leaves #################################################
#############################################################################################################

# Without units
c4 = PH.C4()
Tl = 300.0
res = PH.energybalance(Tl, gb, c4, opt, PAR, NIR, ws, RH, Tair, Ca, P, O2) # W/m^2
Tl = PH.solve_energy_balance(c4, gb = gb, opt = opt, PAR = PAR, NIR = NIR,
    ws = ws, RH = RH, Tair = Tair, Ca = Ca, P = P, O2 = O2,
    xatol = 0.01, order = Order2())

# With units
c4_q = PH.C4Q()
Tl_q = 300.0K
gb_q = PH.simplegbQ()
res_q = PH.energybalance(Tl_q,
    gb_q,
    c4_q,
    opt,
    PAR_q,
    NIR_q,
    ws_q,
    RH_q,
    Ta_q,
    Ca_q,
    P_q,
    O2_q) # W/m^2
@test_broken res_q == res * W / m^2
Tl_q = PH.solve_energy_balance(c4_q, gb = gb_q, opt = opt, PAR = PAR_q,
    NIR = NIR_q, ws = ws_q, RH = RH_q, Tair = Ta_q,
    Ca = Ca_q, P = P_q, O2 = O2_q, xatol = 0.01,
    order = Order2())
@test_broken abs(Tl_q - Tl * K) < 1e-4K

# With effect of angles (horizontal)
gb_q = PH.gbAngleQ()
Tl_q = 300.0K
res_q = PH.energybalance(Tl_q,
    gb_q,
    c4_q,
    opt,
    PAR_q,
    NIR_q,
    ws_q,
    RH_q,
    Ta_q,
    Ca_q,
    P_q,
    O2_q) # W/m^2
@test_broken res_q == res * W / m^2
Tl_q = PH.solve_energy_balance(c4_q, gb = gb_q, opt = opt, PAR = PAR_q,
    NIR = NIR_q, ws = ws_q, RH = RH_q, Tair = Ta_q,
    Ca = Ca_q, P = P_q, O2 = O2_q, xatol = 0.01,
    order = Order2())
@test_broken abs(Tl_q - Tl * K) < 1e-4K

# With effect of angles (45°)
gb_q = PH.gbAngleQ(ang = 45.0)
Tl_q = 300.0K
res_q = PH.energybalance(Tl_q,
    gb_q,
    c4_q,
    opt,
    PAR_q,
    NIR_q,
    ws_q,
    RH_q,
    Ta_q,
    Ca_q,
    P_q,
    O2_q) # W/m^2
@test_broken abs(res_q - 175.4809W / m^2) < 1e-4W / m^2
Tl_q = PH.solve_energy_balance(c4_q, gb = gb_q, opt = opt, PAR = PAR_q,
    NIR = NIR_q, ws = ws_q, RH = RH_q, Tair = Ta_q,
    Ca = Ca_q, P = P_q, O2 = O2_q, xatol = 0.01,
    order = Order2())
@test_broken abs(Tl_q - 306.3663K) < 1e-4K

#end
