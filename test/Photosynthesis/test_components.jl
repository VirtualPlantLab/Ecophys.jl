using Test
import Ecophys
import Unitful: K, J, K, mol, kPa, Quantity
import Unitful
PH = Ecophys.Photosynthesis

# Constants used in components
@test PH.GasConstant(Float64) == 8.314
@test PH.GasConstant(Float32) == 8.314f0
@test PH.GasConstant(Quantity) == 8.314J/K/mol
@test PH.Tref(Float64) == 298.15
@test PH.Tref(Float32) == 298.15f0
@test PH.Tref(Quantity) == 298.15K

# Arrhenius temperature response
@test PH.arrhenius(1.0, 65330.0, 298.15) == 1.0
@test PH.arrhenius(1.0, 65330.0J/mol, 298.15K) == 1.0
@test PH.arrhenius(1.0mol, 65330.0J/mol, 298.15K) == 1.0mol
@test PH.arrhenius(1.0, 65330.0, 298.15 + 10.0) ≈ 2.3520205168012107
@test PH.arrhenius(1.0, 65330.0, 298.15 - 10.0) ≈ 0.40066167345438386

# Peaked temperature response with optimum temperature
@test PH.peaked(1.0, 26900.0, 2e5, 298.15, 298.15) == 1.0
@test PH.peaked(1.0, 26900.0J/mol, 2e5J/mol, 298.15K, 298.15K) == 1.0
@test PH.peaked(1.0mol, 26900.0J/mol, 2e5J/mol, 298.15K, 298.15K) == 1.0mol
@test PH.peaked(1.0, 26900.0, 2e5, 298.15, 298.15 + 10.0) ≈ 0.524803439463303
@test PH.peaked(1.0, 26900.0, 2e5, 298.15, 298.15 - 10.0) ≈ 0.785398095155038

# Compute temperature optimum
@test PH.Topt(26900.0, 2e5, 650.0) ≈ 300.53561483043364
@test PH.Topt(26900.0J/mol, 2e5J/mol, 650.0J/mol/K) ≈ 300.53561483043364K

# Saturated vapour pressure
@test PH.es(298.15) ≈ 3.167688063635486
@test PH.es(298.15K) ≈ 3.167688063635486kPa
@test PH.es(35.86) == 0.0
@test PH.es(35.86K) == 0.0kPa

# VPD
@test PH.VPD(298.15, 1.0) == 0.0
@test PH.VPD(298.15, 0.0) == PH.es(298.15)
