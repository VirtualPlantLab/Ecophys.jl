module Ecophys

export C3, C3Q, C4, C4Q, SimpleOptical, energybalance, solve_energy_balance, gb,
    simplegb, simplegbQ, gbAngle, gbAngleQ, photosynthesis, transpiration

# Photosynthesis models
include("Photosynthesis/Photosynthesis.jl")

# Photosynthesis
const photosynthesis = Photosynthesis.photosynthesis
const C3 = Photosynthesis.C3
const C3Q = Photosynthesis.C3Q
const C4 = Photosynthesis.C4
const C4Q = Photosynthesis.C4Q
const SimpleOptical = Photosynthesis.SimpleOptical
const simplegb = Photosynthesis.simplegb
const simplegbQ = Photosynthesis.simplegbQ
const gbAngle = Photosynthesis.gbAngle
const gbAngleQ = Photosynthesis.gbAngleQ
const gb = Photosynthesis.gb
const energybalance = Photosynthesis.energybalance
const solve_energy_balance = Photosynthesis.solve_energy_balance
const transpiration = Photosynthesis.transpiration
end
