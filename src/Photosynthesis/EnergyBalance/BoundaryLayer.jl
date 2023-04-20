### This file contains public API ###
# gb
# simplegb
# simplegbQ
# gbAngle
# gbAngleQ

import SpecialFunctions: beta

abstract type gbType end

"""
    gb(p::gbType, ws, Tleaf, Tair, P)

Compute boundary layer conductance for heat, water vapor and CO2.

# Arguments
- `p`: Model of boundary layer conductance
- `ws`: Wind speed (m/s)
- `Tleaf`: Leaf temperature (K)
- `Tair`: Air temperature (K)
- `P`: Air pressure (Pa)

# Returns
- `gbh`: Boundary layer conductance for heat (W/m²/K)
- `gbw`: Boundary layer conductance for water vapor (mol/m²/s)
- `gbc`: Boundary layer conductance for CO2 (mol/m²/s)
"""
function gb(p::gbType, ws, Tleaf, Tair, P)
    # Compute constants with the proper units
    mode = typeof(ws)
    Tavg = (Tleaf + Tair)/2
    a    = ThermalExpansion(Tavg)
    g    = Gravity(mode)
    ν    = KinematicViscosity(Tavg)
    k    = ThermalDiffusivity(Tavg)
    Mv   = MolarVolume(Tavg, P)
    R    = GasConstant(mode)

    # Nusselt number under free convection
    # DO NOT USE EXPONENTIATION BECAUSE IT MESSES UP INTEGRATION WITH UNITFUL.JL    
    Gr::Float64 = a*g*p.d*p.d*p.d*abs(Tleaf - Tair)/(ν*ν) # Grashof number
    Nufree::Float64 = 0.5*Gr^0.25 # Horizontal flat plate - laminar flow

    # Forced convection with effect of aspect ratio and inclination
    Re::Float64 = p.d*ws/ν
    Nuforced::Float64 = Nuforc(Re) # Flat plate - laminar & turbulent flow
    front::Float64, back::Float64 = enhancegb(p) # Empirical model based on data reviewed by Schuepp(1993)
    Nuforced = Nuforced/2*(front + back)

    # Mixed convection formula proposed by Schuepp (1993)
    Nu::Float64 = (Nufree^3.5 + Nuforced^3.5)^(1/3.5)

    # Boundary layer conductances
    gbh = Nu*k/p.d/Mv
    gbw = gbh/0.93
    gbc = gbw*1.37

    (gbh, gbw, gbc)
end

# Transition from laminar to turbulent flow
transition() = (log(999)/5e3, 2e4)

# Compute Nusselt number under forced convection for a horizontal flat plate
Nulam(Re)  = 0.6*Re^0.5
Nuturb(Re) = 0.032*Re^0.8
Nulog(Re, k, Reo)  =  1/(1 + exp(-k*(Re - Reo)))
function Nuforc(Re)
    k, Reo = transition()
    Nulam(Re)*Nulog(Re, k, Reo) + (1 - Nulog(Re, k, Reo))*Nuturb(Re)
end


##### Classic model of boundary layer conductance #####

abstract type simplegbType <: gbType end

"""
    simplegb(; d = 0.01)

Simple model of boundary layer conductance.

# Arguments
- `d`: Characteristic leaf length (m)
"""
Base.@kwdef struct simplegb{T <: Real} <: simplegbType
    d::T = 0.01 # Characteristic length (m)
end

"""
    simplegbQ(; d = 0.01m)

Simple model of boundary layer conductance using `Quantity` from Unitful.jl.

# Arguments
- `d`: Characteristic leaf length (m)
"""
Base.@kwdef struct simplegbQ{T <: Real} <: simplegbType
    d::Quantity{T, dimension(m)} = 0.01m # Characteristic length (m)
end

function enhancegb(p::simplegbType)
    return (1.0, 1.0)
end

##### Effect of leaf angle and aspect ratio on boundary layer conductance #####


abstract type gbAngleType <: gbType end
"""
    gbAngle(; d = 0.01, ang = 0.0, ar = 1.0, fangm = 1.381, fangk = 0.034, 
    α = 2.738, b0_0 = 0.455, d_b0 = 2.625, b0_n = 0.373, b0_KAR = 28.125,  db_0 = 0.085, 
    d_db = 0.437, db_n = 5.175, db_KAR = 0.884,  β_0 = 3.362, d_β = 17.664, β_n = 4.727, 
    β_KAR = 0.677)

Model of boundary layer conductance that accounts for inclination angle and leaf
aspect ratio (see documentation for details).

# Arguments
- `d`: Characteristic leaf length (m)
- `ang`: Leaf inclination angle (°)
- `ar`: Leaf aspect ratio (length/width)
- `fangm`: Maximum enhancement factor due to inclination angle
- `fangk`: Exponent in response to inclination angle
- `α`: Effect on back boundary layer conductance due to leaf inclination angle and aspect ratio
- `b0_0`: Parameter in the effect of aspect ratio (see documentation)
- `d_b0`: Parameter in the effect of aspect ratio (see documentation)
- `b0_n`: Parameter in the effect of aspect ratio (see documentation)
- `b0_KAR`: Parameter in the effect of aspect ratio (see documentation)
- `db_0`: Parameter in the effect of aspect ratio (see documentation)
- `d_db`: Parameter in the effect of aspect ratio (see documentation)
- `db_n`: Parameter in the effect of aspect ratio (see documentation)
- `db_KAR`: Parameter in the effect of aspect ratio (see documentation)
- `β_0`: Parameter in the effect of aspect ratio (see documentation)
- `d_β`: Parameter in the effect of aspect ratio (see documentation)
- `β_n`: Parameter in the effect of aspect ratio (see documentation)
- `β_KAR`: Parameter in the effect of aspect ratio (see documentation)
"""
Base.@kwdef struct gbAngle{T <: Real} <: gbAngleType
    d::T     = 0.01   # Characteristic length (m)
    ang::T   = 0.0  # Leaf inclination angle (rad)
    ar::T    = 1.0  # Leaf aspect ratio (length/width)
    # Enhancement of front boundary layer conductance due to leaf inclination angle
    fangm::T = 1.381 # Maximum enhancement factor
    fangk::T = 0.034 # Exponent in response
    # Effect on back boundary layer conductance due to leaf inclination angle and aspect ratio
    α::T     = 2.738
    # 1. Effect of aspect ratio on beta parameter of beta density function
    b0_0::T   = 0.455
    d_b0::T   = 2.625
    b0_n::T   = 0.373
    b0_KAR::T = 28.125
    # 2. Effect of aspect ratio on height factor of beta density function
    db_0::T   = 0.085
    d_db::T   = 0.437
    db_n::T   = 5.175
    db_KAR::T = 0.884
    # 3. Effect of aspect ratio on vertical displacement of beta density function
    β_0::T    = 3.362
    d_β::T    = 17.664
    β_n::T    = 4.727
    β_KAR::T  = 0.677
end

"""
    gbAngleQ(; d = 0.01m, ang = 0.0, ar = 1.0, fangm = 1.381, fangk = 0.034, 
    α = 2.738, b0_0 = 0.455, d_b0 = 2.625, b0_n = 0.373, b0_KAR = 28.125,  db_0 = 0.085, 
    d_db = 0.437, db_n = 5.175, db_KAR = 0.884,  β_0 = 3.362, d_β = 17.664, β_n = 4.727, 
    β_KAR = 0.677)

Model of boundary layer conductance that accounts for inclination angle and leaf
aspect ratio (see documentation for details) using `Quantity` for Unitful.jl.

# Arguments
- `d`: Characteristic leaf length (m)
- `ang`: Leaf inclination angle (°)
- `ar`: Leaf aspect ratio (length/width)
- `fangm`: Maximum enhancement factor due to inclination angle
- `fangk`: Exponent in response to inclination angle
- `α`: Effect on back boundary layer conductance due to leaf inclination angle and aspect ratio
- `b0_0`: Parameter in the effect of aspect ratio (see documentation)
- `d_b0`: Parameter in the effect of aspect ratio (see documentation)
- `b0_n`: Parameter in the effect of aspect ratio (see documentation)
- `b0_KAR`: Parameter in the effect of aspect ratio (see documentation)
- `db_0`: Parameter in the effect of aspect ratio (see documentation)
- `d_db`: Parameter in the effect of aspect ratio (see documentation)
- `db_n`: Parameter in the effect of aspect ratio (see documentation)
- `db_KAR`: Parameter in the effect of aspect ratio (see documentation)
- `β_0`: Parameter in the effect of aspect ratio (see documentation)
- `d_β`: Parameter in the effect of aspect ratio (see documentation)
- `β_n`: Parameter in the effect of aspect ratio (see documentation)
- `β_KAR`: Parameter in the effect of aspect ratio (see documentation)
"""
Base.@kwdef struct gbAngleQ{T <: Real, dT} <: gbAngleType
    d::dT = 0.01m # Characteristic length (m)
    ang::T = 0.0 # Leaf inclination angle (°)
    ar::T  = 1.0 # Leaf aspect ratio (length/width)
    # Enhancement of front boundary layer conductance due to leaf inclination angle
    fangm::T = 1.381 # Maximum enhancement factor
    fangk::T = 0.034 # Exponent in response
    # Effect on back boundary layer conductance due to leaf inclination angle and aspect ratio
    α::T = 2.738
    # 1. Effect of aspect ratio on beta parameter of beta density function
    b0_0::T   = 0.455
    d_b0::T   = 2.625
    b0_n::T   = 0.373
    b0_KAR::T = 28.125
    # 2. Effect of aspect ratio on height factor of beta density function
    db_0::T   = 0.085
    d_db::T   = 0.437
    db_n::T   = 5.175
    db_KAR::T = 0.884
    # 3. Effect of aspect ratio on vertical displacement of beta density function
    β_0::T    = 3.362
    d_β::T    = 17.664
    β_n::T    = 4.727
    β_KAR::T  = 0.677
end

# Compute the effects of leaf angle and aspect ratio on the front and back boundary layer conductance
function enhancegb(p::gbAngleType)
    # Enhancement on the front side
    front = (p.fangm - (p.fangm - 1.0)*exp(-p.ang*p.fangk))
    # Efect on the back side
    b0 = p.b0_0 + p.d_b0*p.ar^p.b0_n/(p.b0_KAR^p.b0_n + p.ar^p.b0_n)
    db = p.db_0 + p.d_db*p.ar^p.db_n/(p.db_KAR^p.db_n + p.ar^p.db_n)
    β = p.β_0 + p.d_β*p.ar^p.β_n/(p.β_KAR^p.β_n + p.ar^p.β_n)
    thres = b0 + db*((3.0/90.0)^(p.α - 1)*(1 - 3.0/90.0)^(β - 1))/beta(p.α, β)
    back = ifelse(p.ang > 3.0, b0 + db*((p.ang/90)^(p.α - 1)*(1 - p.ang/90)^(β - 1))/beta(p.α, β),
                               1 + (thres - 1)/3*p.ang)
    (front, back)
end