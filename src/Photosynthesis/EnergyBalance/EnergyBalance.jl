### This file contains public API ###
# SimpleOptical
# energybalance
# solve_energy_balance

import Roots: find_zero, Order2

abstract type Optical end

"""
    SimpleOptical(; αPAR = 0.85, αNIR = 0.20, ϵ = 0.95)

Simple optical properties of a leaf.

# Arguments

- `αPAR`: Absorption coefficient of PAR
- `αNIR`: Absorption coefficient of NIR
- `ϵ`: Emissivity for thermal radiation
"""
Base.@kwdef struct SimpleOptical{T <: Real} <: Optical
    αPAR::T = 0.85
    αNIR::T = 0.20
    ϵ::T = 0.95
end

"""
    energybalance(pgb, pAgs, pEb, PAR, NIR, ws, RH, Tair, Ca, P, O2)

Calculate the energy balance of a leaf.

# Arguments

- `pgb`: Boundary layer conductance model
- `pAgs`: Photosynthesis and stomatal conductance model
- `pEb`: Optical properties of the leaf
- `PAR`: Photosynthetically active radiation (umol/m2/s)
- `NIR`: Near-infrared radiation (W/m2)
- `ws`: Wind speed (m/s)
- `RH`: Relative humidity
- `Tair`: Air temperature (K)
- `Ca`: Atmospheric CO2 concentration (μmol/mol)
- `P`: Air pressure (kPa)
- `O2`: Atmospheric O2 concentration (μmol/mol)

# Details

Inputs maybe be either `Real` or `Quantity` types (i.e., with physical units). 
If `Quantity` types are used, the output will be a `Quantity` type.
"""
function energybalance(Tleaf, pgb, pAgs, pEb, PAR, NIR, ws, RH, Tair, Ca, P, O2)
    mode = typeof(PAR)
    Tavg = (Tleaf + Tair) / 2
    # Boundary layer conductance
    gbh, gbw, gbc = gb(pgb, ws, Tleaf, Tair, P)
    # Photosynthesis and stomnatal conductance
    A, gsc = photosynthesis(pAgs, PAR = PAR, RH = RH, Tleaf = Tleaf, Ca = Ca, O2 = O2, gb = gbc)
    # Solar radiation absorbed by the leaf
    pe = PAREnergy(mode) # Conversion of PAR from umol/m2/s to W/m2
    Q = PAR * pe * pEb.αPAR + NIR * pEb.αNIR 
    # Sensible heat
    Mv = MolarVolume(Tavg, P)
    ρ  = AirDensity(Tavg)
    Cp = SpecificHeat(Tavg)
    H = (Tleaf - Tair) * ρ * Cp * gbh * Mv
    # Latent heat
    vap_heat = MolarHeatVapour(Tavg)
    gsw = gsc * 1.56
    vpd = VPD(Tavg, RH)
    Tr = 1 / (1 / gsw + 1 / gbw) * vpd / P
    LE = vap_heat * Tr
    # Simple model of thermal radiation
    Rb = thermal(Tair, pEb) - thermal(Tleaf, pEb)
    # Residue of energy balance
    res = Q + Rb - H - LE    
    res
end

# Two-side emmission of thermal radiation
function thermal(Tleaf, p::Optical)
    mode = typeof(Tleaf)
    σ = StefanBoltzmann(mode)
    2 * p.ϵ * σ * Tleaf^4
end

"""
    solve_energy_balance(Ags::Union{C3Q, C4Q}; gb = simplegbQ(), 
                         opt = SimpleOptical(), PAR = 1000.0μmol/m^2/s, 
                         NIR = 250.0W/m^2, ws = 1.0m/s, RH = 0.75, 
                         Tair = 298.0K, Ca = 400.0μmol/mol, P = 101.0kPa, 
                         O2 = 210.0mmol/mol, order = Order2(), xatol = 0.01, 
                         maxfnevals = 100)
    solve_energy_balance(Ags::Union{C3, C4}; gb = simplegb(), 
                         opt = SimpleOptical(), PAR = 1000.0, NIR = 250.0, 
                         ws = 1.0, RH = 0.75, Tair = 298.0, Ca = 400.0, 
                         P = 101.0e3, O2 = 210.0e3, order = Order2(), xatol = 0.01, 
                         maxfnevals = 100)

Solve the leaf energy balance coupled to photosynthesis and transpiration.

# Arguments

- `Ags`: Photosynthesis and stomatal conductance model
- `gb`: Boundary layer conductance model
- `opt`: Optical properties of the leaf
- `PAR`: Photosynthetically active radiation (umol/m2/s)
- `NIR`: Near-infrared radiation (W/m2)
- `ws`: Wind speed (m/s)
- `RH`: Relative humidity
- `Tair`: Air temperature (K)
- `Ca`: Atmospheric CO2 concentration (μmol/mol)
- `P`: Air pressure (Pa)
- `O2`: Atmospheric O2 concentration (μmol/mol)
- `order`: Order of the root solving algorithm that finds leaf temperature
          (see Roots.jl package for more information).
- `xatol`: Absolute tolerance of the root solving algorithm (see Roots.jl package for more information),
- `maxfnevals`: Maximum number of function evaluations of the root solving algorithm (see Roots.jl package for more information).

# Details

Inputs maybe be either `Real` or `Quantity` types from Unitful.jl (i.e., with 
physical units). If `Quantity` types are used, the output will be a `Quantity` 
type.

# Returns

A named tuple with net CO2 assimilation (`A`, μmol/m^2/s), transpiration (`Tr`, 
mol/m^2/s) and leaf temperature (`Tleaf`, K).
"""
function solve_energy_balance(Ags::Union{C3Q, C4Q}; gb = simplegbQ(), 
                              opt = SimpleOptical(), PAR = 1000.0μmol/m^2/s, 
                              NIR = 250.0W/m^2, ws = 1.0m/s, RH = 0.75, 
                              Tair = 298.0K, Ca = 400.0μmol/mol, P = 101.0kPa, 
                              O2 = 210.0mmol/mol, order = Order2(), xatol = 0.01, 
                              maxfnevals = 100)
    solve_energy_balance(Ags, gb, opt, PAR, NIR, ws, RH, Tair, Ca, P, O2,
                         order, xatol, maxfnevals)
end
function solve_energy_balance(Ags::Union{C3, C4}; gb = simplegb(), 
                              opt = SimpleOptical(), PAR = 1000.0, NIR = 250.0, 
                              ws = 1.0, RH = 0.75, Tair = 298.0, Ca = 400.0, 
                              P = 101.0e3, O2 = 210.0e3, order = Order2(), xatol = 0.01, 
                              maxfnevals = 100)
    solve_energy_balance(Ags, gb, opt, PAR, NIR, ws, RH, Tair, Ca, P, O2,
                         order, xatol, maxfnevals)
end

function solve_energy_balance(pAgs, pgb, pEb, PAR, NIR, ws, RH, Tair, Ca, P, O2,
                              order = Order2(), xatol = 0.01, maxfnevals = 100)
    # Find the temperature
    Tleaf = find_zero(x -> energybalance(x, pgb, pAgs, pEb, PAR, NIR, ws, RH, Tair, Ca, P, O2), 
                                  (Tair - 10, Tair + 10), order, xatol = xatol, maxfnevals = maxfnevals)
    # Boundary layer conductances
    gbh, gbw, gbc = gb(pgb, ws, Tleaf, Tair, P)
    # A and gsc
    A, gsc = photosynthesis(pAgs, PAR = PAR, Tleaf = Tleaf, RH = RH, Ca = Ca, O2 = O2, gb = gbc)
    # Transpiration
    Tr  = transpiration(gsw = gsc*1.6, gbw = gbw, Tleaf = Tleaf, Tair = Tair, P = P, RH = RH)
    # Return the temperature and two fluxes
    return (Tleaf = Tleaf, A = A, Tr = Tr)
end
function solve_energy_balance(pAgs, pgb, pEb, PAR::Quantity, NIR::Quantity, 
                              ws::Quantity, RH, Tair::Quantity, Ca::Quantity, 
                              P::Quantity, O2::Quantity, order = Order2(), 
                              xatol = 0.01, maxfnevals = 100)
    Tau = Tair/1.0K
    Tleaf = find_zero(x -> energybalance(x*K, pgb, pAgs, pEb, PAR, NIR, ws, RH, Tair, Ca, P, O2), 
                                  (Tau - 10, Tau + 10), order, xatol = xatol, maxfnevals = maxfnevals)
    Tleaf *= 1K
    # Boundary layer conductances
    gbh, gbw, gbc = gb(pgb, ws, Tleaf, Tair, P)
    # A and gsc
    A, gsc = photosynthesis(pAgs, PAR = PAR, Tleaf = Tleaf, RH = RH, Ca = Ca, O2 = O2, gb = gbc)
    # Transpiration
    Tr  = transpiration(gsw = gsc*1.6, gbw = gbw, Tleaf = Tleaf, Tair = Tair, P = P, RH = RH)
    # Return the temperature and two fluxes
    return (Tleaf = Tleaf, A = A, Tr = Tr)
end


"""
    transpiration(;gsw = 0.1, gbw = 1.0, Tleaf = 300.0, Tair = 298.0, P = 101e3,
                   RH = 0.75)

Compute transpiration rate (mol/m^2/s) from conductance to water vapor and 
environmental variables.

# Arguments
- `gsw`: Stomatal conductance to water vapor (mol/m^2/s)
- `gbw`: Boundary layer conductance to water vapor (mol/m^2/s)
- `Tleaf`: Leaf temperature (K)
- `Tair`: Air temperature (K)
- `P`: Air pressure (Pa)
- `RH`: Relative humidity
"""
function transpiration(;gsw = 0.1, gbw = 1.0, Tleaf = 300.0, Tair = 298.0, P = 101e3,
                        RH = 0.75)
    Tavg = (Tleaf + Tair) / 2
    vpd = VPD(Tavg, RH)
    1/(1/gsw + 1/gbw)*vpd/P
end