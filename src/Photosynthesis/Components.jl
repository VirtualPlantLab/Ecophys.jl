
###############################################################################
################################## Constants ##################################
###############################################################################

# Gas constant
GasConstant(::Type{T}) where T <: Real = T(8.314) # J/K/mol
GasConstant(::Type{T}) where T <: Unitful.Quantity = 8.314J/K/mol

# 25 Celsius in Kelvin
Tref(::Type{T}) where T <: Real = T(298.15)
Tref(::Type{T}) where T <: Unitful.Quantity = 298.15K
# 0 Celsius in Kelvin
T0(::Type{T}) where T <: Real = T(273.15)
T0(::Type{T}) where T <: Unitful.Quantity = 273.15K
# Temperature in Murray's model for es(T)
Tb(::Type{T}) where T <: Real = T(35.86)
Tb(::Type{T}) where T <: Unitful.Quantity = 35.86K
# es(T = 35.86 K) in Murray's model
es0(::Type{T}) where T <: Real = T(0.61078)
es0(::Type{T}) where T <: Unitful.Quantity = 0.61078kPa


# 0 umol/m2/s
zeroflux(::Type{T}) where T <: Real = T(0.0) # μmol/m^2/s
zeroflux(::Type{T}) where T <: Unitful.Quantity = 0.0μmol/m^2/s

# 0 umol/mol
zeroconc(::Type{T}) where T <: Real = T(0.0) # μmol/mol
zeroconc(::Type{T}) where T <: Unitful.Quantity = 0.0μmol/mol


###############################################################################
########################### Temperature responses #############################
###############################################################################

# Based on Medlyn et al. (2002, Plant, Cell & Environment)

# Exponential increase with temperature, no optimum, normalized at 25 °C (Equation 10.3)
function arrhenius(p₂₅, Hₐ, Tₖ)
    mode = typeof(Hₐ)
    R = GasConstant(mode)
    T₂₅ = Tref(mode)
    p₂₅*exp((Tₖ - T₂₅)*Hₐ/(T₂₅*R*Tₖ))
end

# Increase with temperature with an optimum, normalized at optimum (Equation 10.5)
function peaked(pₒₚₜ, Hₐ, Hd, Tₒₚₜ, Tₖ)
    mode = typeof(Hₐ)
    R = GasConstant(mode) 
    eterm1 = exp(Hₐ/R*(inv(Tₒₚₜ) - inv(Tₖ)))
    eterm2 = exp(Hd/R*(inv(Tₒₚₜ) - inv(Tₖ)))
    pₒₚₜ*Hd*eterm1/(Hd - Hₐ*(1.0 - eterm2))
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