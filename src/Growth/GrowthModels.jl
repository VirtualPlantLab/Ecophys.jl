# In development. If you find any issue, please report it to the repository. Thank you for your collaboration!

### This file contains public API ###
# Organ, OrganQ
# compute_potential_GR

abstract type OrganType end

"""
    Organ(f = [0.53, 0.25, 0.05, 0.05, 0.06, 0.06],
          te = 100.0, tm = 50.0, tb = 10.0, wmax = 0.1)

Data structure to store all parameters for the growth model of a plant organ.

# Arguments
- `f`: fraction of dry matter allocated to different carbon pools of newly formed material, [carbohydrates, proteins, lipids, lignin, organic acids, minerals]
- `te`: time at which growth rate extinguishes, [DD]
- `tm`: time at which growth rate is maximum, [DD]
- `tb`: base temperature, [°C]
- `wmax`: maximum weight, [g] 
"""
Base.@kwdef mutable struct Organ <: OrganType
    f::Vector{Float64} = [0.53, 0.25, 0.05, 0.05, 0.06, 0.06] # Fractions of dry matter allocated to different carbon pools of newly formed material
    te::Float64 = 100.0 # Time at which growth rate extinguishes, [DD]
    tm::Float64 = 50.0 # Time at which growth rate is maximum, [DD]
    tb::Float64 = 10.0 # Moment at which growth starts, [DD]
    wmax::Float64 = 0.1 # Maximum organ weight, [g]
    tt::Float64 = 20.0 # Thermal time, [DD]
end

Base.@kwdef mutable struct OrganQ{T <: Real} <: OrganType
    f::Vector{T} = [0.53, 0.25, 0.05, 0.05, 0.06, 0.06] # Fractions of dry matter allocated to different carbon pools of newly formed material
    te::Quantity{T, dimension(°C)} = 100.0°C # Time at which growth rate extinguishes, [DD, °C]
    tm::Quantity{T, dimension(°C)} = 50.0°C # Time at which growth rate is maximum, [DD, °C]
    tb::Quantity{T, dimension(°C)} = 10.0°C # Moment at which growth starts, [DD, °C]
    wmax::Quantity{T, dimension(g)} = 0.1g # Maximum organ weight, [g]
    tt::Quantity{T, dimension(°C)} = 20.0°C # Thermal time, [DD, °C]
end

"""
    compute_potential_GR(par::OrganType)

Compute potential growth rate of the leaf according to thermal time (Yin et al 2003)
The potential growth rate is computed as a function of thermal time (tt, DD), and the parameters of the growth model: 
    - time at which growth rate extinguishes (te, DD),
    - time at which growth rate is maximum (tm, DD)
    - base thermal time at which growth starts (tb, DD)
    - maximum weight (wmax, g)

Citations:
Yin, Xinyou & Goudriaan, Jan & Lantinga, Egbert & Vos, Jan & Spiertz, Huub. (2003). A Flexible Sigmoid Function of Determinate Growth. Annals of botany. 91. 361-71. 10.1093/aob/mcg029. 
"""

function compute_potential_GR(t = 20.0, te = 100.0, tm = 50.0, tb = 10.0, wmax = 0.1)
    #Compute maximum growth rate
    cm = wmax * ((2 * te - tm) / (te * (te - tm))) * (tm / te) ^ (tm / (tm - te)) 
    #Compute growth rate
    exp = (tm - tb) / (te - tm)
    c = (t > te || t < tb) ? 0.0 : cm * ((te - t) / (te - tm)) * (((t - tb) / (tm - tb)) ^ exp)
end

function compute_potential_GR(par::Organ)
    #Extract parameters from organ
    t, te, tm, tb, wmax = par.tt, par.te, par.tm, par.tb, par.wmax
    return compute_potential_GR(t, te, tm, tb, wmax)
end

function compute_potential_GR(par::OrganQ)
    #Extract parameters from organ and unstrip from units
    t, te, tm, tb, wmax = ustrip(par.tt), ustrip(par.te), ustrip(par.tm), ustrip(par.tb), ustrip(par.wmax)
    return compute_potential_GR(t, te, tm, tb, wmax)*u"g/d"
end

#Example of compute_c using different data types as input:
compute_potential_GR()
compute_potential_GR(Organ())
compute_potential_GR(OrganQ())

