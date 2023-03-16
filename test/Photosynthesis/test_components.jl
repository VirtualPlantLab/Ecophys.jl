
using Test
import Ecophys
import Unitful: K, J, K, mol, kPa, Quantity, μmol, m, s, W, Pa, N, kg
import Unitful
PH = Ecophys.Photosynthesis

let

# Constants used in different modules
@test PH.GasConstant(Float64) == 8.31446261815324
@test PH.GasConstant(Float32) == 8.31446261815324f0
@test PH.GasConstant(Quantity) == 8.31446261815324J/K/mol

@test PH.Tref(Float64) == 298.15
@test PH.Tref(Float32) == 298.15f0
@test PH.Tref(Quantity) == 298.15K

@test PH.StefanBoltzmann(Float64) == 5.670374419e-8
@test PH.StefanBoltzmann(Float32) == 5.670374419f-8
@test PH.StefanBoltzmann(Quantity) == 5.670374419e-8W/m^2/K^4

@test PH.Gravity(Float64) == 9.7803274
@test PH.Gravity(Float32) == 9.7803274f0
@test PH.Gravity(Quantity) == 9.7803274m/s^2

@test PH.PAREnergy(Float64) == 1/4.56
@test PH.PAREnergy(Float32) == 1/4.56f0
@test PH.PAREnergy(Quantity) == 1/4.56*J/μmol

# Arrhenius temperature response
@test PH.arrhenius(1.0, 65330.0, 298.15) == 1.0
@test PH.arrhenius(1.0, 65330.0J/mol, 298.15K) == 1.0
@test PH.arrhenius(1.0mol, 65330.0J/mol, 298.15K) == 1.0mol

# Peaked temperature response with optimum temperature
@test PH.peaked(1.0, 26900.0, 2e5, 298.15, 298.15) == 1.0
@test PH.peaked(1.0, 26900.0J/mol, 2e5J/mol, 298.15K, 298.15K) == 1.0
@test PH.peaked(1.0mol, 26900.0J/mol, 2e5J/mol, 298.15K, 298.15K) == 1.0mol

# Compute temperature optimum
@test (PH.Topt(26900.0, 2e5, 650.0) - 300.5) < 0.1
@test abs(PH.Topt(26900.0J/mol, 2e5J/mol, 650.0J/mol/K) - 300.5K) < 0.1K

# Saturated vapour pressure
@test abs(PH.es(298.15) - 3167.69) < 0.01
@test abs(PH.es(298.15K) - 3.16kPa) < 0.01kPa
@test PH.es(35.86) == 0.0
@test PH.es(35.86K) == 0.0kPa

# VPD
@test PH.VPD(298.15, 1.0) == 0.0
@test PH.VPD(298.15, 0.0) == PH.es(298.15)

# zeros with units
@test PH.zeroflux(Float64) == 0.0
@test PH.zeroflux(Float32) == 0f0
@test PH.zeroflux(Quantity) == 0.0μmol/m^2/s

# molar fractions with units
@test PH.zeroconc(Float64) == 0.0
@test PH.zeroconc(Float32) == 0f0
@test PH.zeroconc(Quantity) == 0.0mol/mol

# Thermal expansion
@test (PH.ThermalExpansion(298.15) - 0.0033) < 1e-4
@test abs(PH.ThermalExpansion(298.15K) - 0.0033/K) < 1e-4/K

# Kinematic viscosity
@test abs(PH.KinematicViscosity(298.15) - 1.5698e-5) < 1e-8
@test abs(PH.KinematicViscosity(298.15K) - 1.5698e-5m^2/s) < 1e-8m^2/s

# Thermal diffusivity
@test abs(PH.ThermalDiffusivity(298.15) - 1.942e-5) < 1e-8
@test abs(PH.ThermalDiffusivity(298.15K) - 1.942e-5m^2/s) < 1e-8m^2/s

# Molar volume 
@test abs(PH.MolarVolume(298.15, 101325) - 0.0244) < 1e-4
@test abs(PH.MolarVolume(298.15K, 101325Pa) - 0.0244m^3/mol) < 1e-4m^3/mol

# Air density
@test abs(PH.AirDensity(298.15)  - 1.1975) < 1e-4
@test abs(PH.AirDensity(298.15K) - 1.1975kg/m^3) < 1e-4kg/m^3

# Specific heat of the air
@test abs(PH.SpecificHeat(298.15) - 1006.7273) < 1e-4
@test abs(PH.SpecificHeat(298.15K) - 1006.7273J/kg/K) < 1e-4J/kg/K

# Molar heat of vaporization
@test PH.MolarHeatVapour(298.15) == 44076.55
@test PH.MolarHeatVapour(298.15K) == 44076.55J/mol

end