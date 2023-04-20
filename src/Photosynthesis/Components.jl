
# Constants are takemn from NIST Reference on Constants, Units and Uncertainty
# https://physics.nist.gov/cuu/Constants/index.html

# Air properties that vary time are fitted to The Engineering Toolbox
# https://www.engineeringtoolbox.com/

# This needs to be revised to include the effecty of humidity with the equations
# by Tsilingris (2007, Thermophysical and transport properties of humid air
# at temperature range between 0 °C and 100 °C)

###############################################################################
################################## Constants ##################################
###############################################################################

# Gas constant
GasConstant(::Type{T}) where T <: Real = T(8.31446261815324) # J/K/mol
GasConstant(::Type{T}) where T <: Quantity = 8.31446261815324J/K/mol

# Stefan-Boltzmann constant for thermal radiation emmission
StefanBoltzmann(::Type{T}) where T <: Real = T(5.670374419e-8) # W/m2/K4
StefanBoltzmann(::Type{T}) where T <: Quantity = 5.670374419e-8W/m^2/K^4

# Earth gravity
Gravity(::Type{T}) where T <: Real = T(9.7803274) # m/s^2
Gravity(::Type{T}) where T <: Quantity = 9.7803274m/s^2

# Energy of a umol of PAR
PAREnergy(::Type{T}) where T <: Real = T(1/4.56) # J/umol
PAREnergy(::Type{T}) where T <: Quantity = 1/4.56*J/μmol

# 25 Celsius in Kelvin
Tref(::Type{T}) where T <: Real = T(298.15)
Tref(::Type{T}) where T <: Quantity = 298.15K

# 0 Celsius in Kelvin
T0(::Type{T}) where T <: Real = T(273.15)
T0(::Type{T}) where T <: Quantity = 273.15K

# Temperature in Murray's model for es(T)
Tb(::Type{T}) where T <: Real = T(35.86)
Tb(::Type{T}) where T <: Quantity = 35.86K

# es(T = 35.86 K) in Murray's model
es0(::Type{T}) where T <: Real = T(610.78)
es0(::Type{T}) where T <: Quantity = 610.78Pa

# 0 umol/m2/s
zeroflux(::Type{T}) where T <: Real = T(0.0) # μmol/m^2/s
zeroflux(::Type{T}) where T <: Quantity = 0.0μmol/m^2/s

# 0 umol/mol
zeroconc(::Type{T}) where T <: Real = T(0.0) # μmol/mol
zeroconc(::Type{T}) where T <: Quantity = 0.0μmol/mol

###############################################################################
########################### Temperature responses #############################
###############################################################################

# Based on Medlyn et al. (2002, Plant, Cell & Environment)

# Exponential increase with temperature, no optimum, normalized at 25 °C (Equation 10.3)
function arrhenius(p25, Ha, Tk)
    mode = typeof(Ha)
    R = GasConstant(mode)
    T₂₅ = Tref(mode)
    p25*exp((Tk - T₂₅)*Ha/(T₂₅*R*Tk))
end

# Increase with temperature with an optimum, normalized at optimum (Equation 10.5)
function peaked(pₒₚₜ, Ha, Hd, Topt, Tk)
    mode = typeof(Ha)
    R = GasConstant(mode) 
    eterm1 = exp(Ha/R*(inv(Topt) - inv(Tk)))
    eterm2 = exp(Hd/R*(inv(Topt) - inv(Tk)))
    pₒₚₜ*Hd*eterm1/(Hd - Ha*(1.0 - eterm2))
end

# Calculate temperature optimum from entropy and other parameters
function Topt(Ha, Hd, S)
    mode = typeof(Ha)
    R = GasConstant(mode) 
    Hd/(S - R*log(Ha/(Hd - Ha)))
end

###############################################################################
################################ VPD responses ################################
###############################################################################

# Equation by Murray (1967, Internal Report at the RAND corporation)
# Compute saturate air vapour pressure as a function of absolute temperature
function es(T)
    mode = typeof(T)
    T₀ = T0(mode)
    b  = Tb(mode)
    es₀ = es0(mode)
    es₀*exp(17.2693882*(T - T₀)/(T - b))
end

# Compute VPD as a function of temperature and relative humidity
function VPD(T, RH)
    es(T)*(1.0 - RH)
end

###############################################################################
################################ Air properties ###############################
###############################################################################

# Air thermal expansion coefficient
function ThermalExpansion(Tavg::T) where T <: Real
    1e-3/(0.001Tavg + 0.0047)
end

function ThermalExpansion(Tavg::T) where T <: Quantity
    1e-3/(0.001Tavg + 0.0047K)
end

# Kinematic viscosity
function KinematicViscosity(Tavg::T) where T <: Real
    T((1e-5Tavg*Tavg + 0.0034Tavg - 0.3328)*1e-5) # m^2/s
end
function KinematicViscosity(Tavg::T) where T <: Quantity
    (1e-5/K^2*Tavg*Tavg + 0.0034/K*Tavg - 0.3328)*1e-5m^2/s # m^2/s
end

# Thermal diffusivity
function ThermalDiffusivity(Tavg::T) where T <: Real
    T((0.0001Tavg*Tavg + 0.0561Tavg - 6.1952)*1e-6) # m^2/s
end
function ThermalDiffusivity(Tavg::T) where T <: Quantity
    (0.0001Tavg*Tavg/K^2 + 0.0561*Tavg/K - 6.1952)*1e-6m^2/s # m^2/s
end

# Molar volume 
MolarVolume(Tavg, P) = GasConstant(typeof(Tavg))*Tavg/P

# Air density
AirDensity(Tavg::T) where T <: Real = 1/(0.0028Tavg + 0.0002) # kg/m^3
AirDensity(Tavg::T) where T <: Quantity = kg/m^3/(0.0028/K*Tavg + 0.0002)

# Specific heat of the air
SpecificHeat(Tavg::T) where T <: Real = (4e-7Tavg^2 - 0.0002Tavg + 1.0308)*1e3 # J/kg/K
SpecificHeat(Tavg::T) where T <: Quantity = (4e-7Tavg*Tavg/(K^2) - 0.0002/K*Tavg + 1.0308)*1e3J/kg/K

# Molar heat of vaporization
MolarHeatVapour(Tavg::T) where T <: Real = 56.897e3 - 0.043e3*Tavg # J/mol
MolarHeatVapour(Tavg::T) where T <: Quantity = (56.897e3 - 0.043e3/K*Tavg)*J/mol

