import SpecialFunctions: beta

abstract type gbtype end

function gb(p::gbtype, ws, Tl, Ta, P)
    # Compute constants with the proper units
    mode = typeof(ws)
    Tavg = (Tl + Ta)/2
    a    = ThermalExpansion(Tavg)
    g    = Gravity(mode)
    ν    = KinematicViscosity(Tavg)
    k    = ThermalDiffusivity(Tavg)
    Mv   = MolarVolume(Tavg, P)
    R    = GasConstant(mode)

    # Nusselt number under dree convection
    Gr = a*g*p.d^3*abs(Tl - Ta)/ν^2 # Grashof number
    Nufree = 0.5*Gr^0.25 # Horizontal flat plate - laminar flow

    # Forced convection with effect of aspect ratio and inclination
    Re = p.d*ws/ν
    Nuforced = Nuforc(Re) # Flat plate - laminar & turbulent flow
    front, back = enhancegb(p) # Empirical model based on data reviewer by Schuepp(1993)
    Nuforced = Nuforced/2*(front + back)

    # Mixed convection formula proposed by Schuepp (1993)
    Nu = (Nufree^3.5 + Nuforced^3.5)^(1/3.5)

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

abstract type simplegb <: gbtype end

struct simplegbF{T <: Real} <: simplegb
    d::T # Characteristic length (m)
end

struct simplegbQ{T <: Real} <: simplegb
    d::Quantity{T, dimension(m)} # Characteristic length (m)
end

function simplegbF(; d = 0.01)
    simplegbF(d)
end

function simplegbQ(; d = 0.01m)
    simplegbQ(d)
end

function enhancegb(p::simplegb)
    return (1.0, 1.0)
end

##### Effect of leaf angle and aspect ratio on boundary layer conductance #####


abstract type gbAngle <: gbtype end

struct gbAngleF{T <: Real} <: gbAngle
    d::T     # Characteristic length (m)
    ang::T  # Leaf inclination angle (rad)
    ar::T    # Leaf aspect ratio (length/width)
    # Enhancement of front boundary layer conductance due to leaf inclination angle
    fangm::T # Maximum enhancement factor
    fangk::T # Exponent in response
    # Effect on back boundary layer conductance due to leaf inclination angle and aspect ratio
    α::T
    # 1. Effect of aspect ratio on beta parameter of beta density function
    b0_0::T
    d_b0::T
    b0_n::T
    b0_KAR::T
    # 2. Effect of aspect ratio on height factor of beta density function
    db_0::T
    d_db::T
    db_n::T
    db_KAR::T
    # 3. Effect of aspect ratio on vertical displacement of beta density function
    β_0::T
    d_β::T
    β_n::T
    β_KAR::T 
end

struct gbAngleQ{T <: Real} <: gbAngle
    d::Quantity{T, dimension(m)}     # Characteristic length (m)
    ang::T  # Leaf inclination angle (°)
    ar::T    # Leaf aspect ratio (length/width)
    # Enhancement of front boundary layer conductance due to leaf inclination angle
    fangm::T # Maximum enhancement factor
    fangk::T # Exponent in response
    # Effect on back boundary layer conductance due to leaf inclination angle and aspect ratio
    α::T
    # 1. Effect of aspect ratio on beta parameter of beta density function
    b0_0::T
    d_b0::T
    b0_n::T
    b0_KAR::T
    # 2. Effect of aspect ratio on height factor of beta density function
    db_0::T
    d_db::T
    db_n::T
    db_KAR::T
    # 3. Effect of aspect ratio on vertical displacement of beta density function
    β_0::T
    d_β::T
    β_n::T
    β_KAR::T 
end

function gbAngleF(; d = 0.01, ang = 0.0, ar = 1.0, fangm = 1.381, fangk = 0.034, 
     α = 2.738, b0_0 = 0.455, d_b0 = 2.625, b0_n = 0.373, b0_KAR = 28.125,  db_0 = 0.085, 
     d_db = 0.437, db_n = 5.175, db_KAR = 0.884,  β_0 = 3.362, d_β = 17.664, β_n = 4.727, 
     β_KAR = 0.677)
  gbAngleF(d, ang , ar, fangm, fangk, α, b0_0, d_b0, b0_n, b0_KAR, 
    db_0, d_db, db_n, db_KAR, β_0, d_β, β_n, β_KAR)
end

function gbAngleQ(; d = 0.01m, ang = 0.0, ar = 1.0, fangm = 1.381, fangk = 0.034, 
    α = 2.738, b0_0 = 0.455, d_b0 = 2.625, b0_n = 0.373, b0_KAR = 28.125,  db_0 = 0.085, 
    d_db = 0.437, db_n = 5.175, db_KAR = 0.884,  β_0 = 3.362, d_β = 17.664, β_n = 4.727, 
    β_KAR = 0.677)
 gbAngleQ(d, ang , ar, fangm, fangk, α, b0_0, d_b0, b0_n, b0_KAR, 
   db_0, d_db, db_n, db_KAR, β_0, d_β, β_n, β_KAR)
end

# Compute the effects of leaf angle and aspect ratio on the front and back boundary layer conductance
function enhancegb(p::gbAngle)
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

