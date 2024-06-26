module Growth

#General abstract type for growth models
abstract type Organ end

#Functions and data structures shared by growth models
include("GrowthModels.jl")
include("GrowthComponents.jl")

end