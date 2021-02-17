

import Roots: find_zero, Order2

abstract type Optical end

struct SimpleOptical{T <: Real} <: Optical
    αPAR::T
    αNIR::T
    ϵ::T
end

function SimpleOptical(; αPAR=0.85, αNIR=0.15, ϵ=0.95)
    SimpleOptical(αPAR, αNIR, ϵ)
end

function energybalance(Tl, pgb, pAgs, pEb, PAR, NIR, ws, RH, Ta, Ca, P, O2)
    mode = typeof(PAR)
    Tavg = (Tl + Ta) / 2
    # Boundary layer conductance
    gbh, gbw, gbc = gb(pgb, ws, Tl, Ta, P)
    # Photosynthesis and stomnatal conductance
    A, gsc = A_gs(pAgs, PAR, RH, Tl, Ca, O2, gbc)
    # Solar radiation absorbed by the leaf
    pe = PAREnergy(mode) # Conversion of PAR from umol/m2/s to W/m2
    Q = PAR * pe * pEb.αPAR + NIR * pEb.αNIR 
    # Sensible heat
    Mv = MolarVolume(Tavg, P)
    ρ  = AirDensity(Tavg)
    Cp = SpecificHeat(Tavg)
    H = (Tl - Ta) * ρ * Cp * gbh * Mv
    # Latent heat
    vap_heat = MolarHeatVapour(Tavg)
    gsw = gsc * 1.56
    vpd = VPD(Tavg, RH)
    Tr = 1 / (1 / gsw + 1 / gbw) * vpd / P
    LE = vap_heat * Tr
    # Simple model of thermal radiation
    Rb = thermal(Ta, pEb) - thermal(Tl, pEb)
    # Residue of energy balance
    res = Q + Rb - H - LE    
    res
end

# Two-side emmission of thermal radiation
function thermal(Tl, p::Optical)
    mode = typeof(Tl)
    σ = StefanBoltzmann(mode)
    2 * p.ϵ * σ * Tl^4
end

# Solve the energy balance and return fluxes
function solve_energy_balance(pgb, pAgs, pEb, PAR, NIR, ws, RH, Ta::Real, Ca, P, O2;
                              order = Order2(), xatol = 0.01, maxfnevals = 100)
    Tl = find_zero(x -> energybalance(x, pgb, pAgs, pEb, PAR, NIR, ws, RH, Ta, Ca, P, O2), 
                                  (Ta - 10, Ta + 10), order, xatol = xatol, maxfnevals = maxfnevals)
end
function solve_energy_balance(pgb, pAgs, pEb, PAR, NIR, ws, RH, Ta::Quantity, Ca, P, O2;
                              order = Order2(), xatol = 0.01, maxfnevals = 100)
    Tau = Ta/1.0K
    Tl = find_zero(x -> energybalance(x*K, pgb, pAgs, pEb, PAR, NIR, ws, RH, Ta, Ca, P, O2), 
                                  (Tau - 10, Tau + 10), order, xatol = xatol, maxfnevals = maxfnevals)
    Tl*K
end