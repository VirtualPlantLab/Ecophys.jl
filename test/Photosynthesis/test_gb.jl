using Test
import Ecophys
import Unitful
import Unitful: m, dimension
PH = Ecophys.Photosynthesis

# Conditions for testing
ws_f = 1.0
Tl_f = 300.0
Ta_f = 298.0
P_f = 101325.0

ws_q = 1.0m/s
Tl_q = 300.0K
Ta_q = 298.0K
P_q = 101325.0Pa

# Classic boundary layer conductance - no units
gb_f = PH.simplegbF()
@test gb_f.d == 0.01
gbh_f, gbw_f, gbc_f = PH.gb(gb_f, ws_f, Tl_f, Ta_f, P_f)
@test gbh_f ≈ 0.44748020949010225
@test gbw_f ≈ 0.4811615155807551
@test gbc_f ≈ 0.6591912763456346

# Classic boundary layer conductance - units
gb_q = PH.simplegbQ()
@test gb_q.d == 0.01m
@test PH.enhancegb(gb_f) == (1.0, 1.0)
@test PH.enhancegb(gb_q) == (1.0, 1.0)
gbh_q, gbw_q, gbc_q = PH.gb(gb_q, ws_q, Tl_q, Ta_q, P_q)
@test gbh_q ≈ 0.44748020961388435mol/m^2/s
@test gbw_q ≈ 0.4811615155807551mol/m^2/s
@test gbc_q ≈ 0.6591912763456346mol/m^2/s

# Novel boundary layer conductance model - no units
gba_f = PH.gbAngleF()
@test gba_f.d == 0.01
gbh_f, gbw_f, gbc_f = PH.gb(gb_f, ws_f, Tl_f, Ta_f, P_f)
@test gbh_f ≈ 0.44748020949010225
@test gbw_f ≈ 0.4811615155807551
@test gbc_f ≈ 0.6591912763456346

gba_q = PH.gbAngleQ()
@test gba_q.d == 0.01m
@test PH.enhancegb(gba_f) == (1.0, 1.0)
@test PH.enhancegb(gba_q) == (1.0, 1.0)
gbh_q, gbw_q, gbc_q = PH.gb(gba_q, ws_q, Tl_q, Ta_q, P_q)
@test gbh_q ≈ gbh_f*mol/m^2/s
@test gbw_q ≈ gbw_f*mol/m^2/s
@test gbc_q ≈ gbc_f*mol/m^2/s

gba_q = PH.gbAngleQ(ang = 45.0)
@test gba_q.ang == 45.0
@test all(PH.enhancegb(gba_q) .≈ (1.2984999107526014, 1.0432376629376556))
gbh_q, gbw_q, gbc_q = PH.gb(gba_q, ws_q, Tl_q, Ta_q, P_q)
@test abs(gbh_q - 0.5222mol/m^2/s) < 1e-4mol/m^2/s
@test abs(gbw_q - 0.5615mol/m^2/s) < 1e-4mol/m^2/s
@test abs(gbc_q - 0.7692mol/m^2/s) < 1e-4mol/m^2/s
