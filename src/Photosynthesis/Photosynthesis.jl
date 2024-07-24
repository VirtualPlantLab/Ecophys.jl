module Photosynthesis

import Unitful
import Unitful: K, J, K, mol, kPa, Pa, m, μmol, mmol, s, W, kg, dimension, Quantity

export C3, C3Q, C4, C4Q, SimpleOptical, energybalance, solve_energy_balance, gb,
    simplegb, simplegbQ, gbAngle, gbAngleQ, photosynthesis, transpiration

# Most general abstract types
abstract type Ags end

# Functions and data structures shared by different photosynthesis models
include("Components.jl")
include("FvCB/C3.jl")
include("FvCB/C4.jl")
include("EnergyBalance/BoundaryLayer.jl")
include("EnergyBalance/EnergyBalance.jl")

end
