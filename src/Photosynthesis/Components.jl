
###############################################################################
################################## Constants ##################################
###############################################################################

# Gas constant
GasConstant(::Type{T}) where T <: Real = T(8.314) # J/K/mol
GasConstant(::Type{T}) where T <: Unitful.Quantity = 8.314J/K/mol

# Base temperature for responses
Tref(::Type{T}) where T <: Real = T(298.15)
Tref(::Type{T}) where T <: Unitful.Quantity = 298.15K


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

# Increase with temperature with an optimum, normalized at 25 °C (Equation 10.5)
function peaked(pₒₚₜ, Hₐ, Hd, Tₒₚₜ, Tₖ)
    mode = typeof(Hₐ)
    R = GasConstant(mode) 
    eterm1 = exp(Hₐ/R*(inv(Tₒₚₜ) - inv(Tₖ)))
    eterm2 = exp(Hd/R*(inv(Tₒₚₜ) - inv(Tₖ)))
    pₒₚₜ*Hd*eterm1/(Hd - Hₐ*(1.0 - eterm2))
end
