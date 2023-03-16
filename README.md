# Ecophys

**This package is still in development**

This package contains modules describing different ecophysiological functions of 
plants, including processes such as photosynthesis, respiration, transpiration 
or phenology. They may be used as standalone or as a component of a plant growth
model.

## Installation

The package Ecophys is not in the Julia registry. To install it, you can use the 
following command:

```julia
Pkg.add(url = "https://github.com/AleMorales/Ecophys.jl.git")
```
## Photosynthesis

The module Photosynthesis contains functions to calculate leaf CO2 assimilation
and stomatal conductance for C3 and C4 species, based on the work by [Yin & Struik (2009, NJAS)](https://www.tandfonline.com/doi/full/10.1016/j.njas.2009.07.001).
To create a model, use the corresponding function (`C3()` or `C4()`) and pass the
parameters as keyword arguments (they all have default values that correspond to Tables 2 the original publication):

```julia
using Ecophys
c3 = C3(Vcmax25 = 140.0)
c4 = C4(Vcmax25 = 140.0)
```
To compute CO2 assimilation and stomatal conductance, use the `photosynthesis()` function, 
passing the photosynthesis model and the environmental conditions as inputs (with
defaults):

```julia
A_c3, gs_c3  = photosynthesis(c3, PAR = 100.0)
A_c4, gs_c4  = photosynthesis(c4, PAR = 100.0)
```

It is also possible to work with physical units using the Unitful.jl package. In 
such case, the functions `C3Q()` and `C4Q` should be used to create the model but
now the parameters are stored as `Quantity` objects:

```julia
using Unitful.DefaultSymbols # import symbols for units
c3Q = C3Q(Vcmax25 = 140.0μmol/m^2/s)
c4Q = C4Q(Vcmax25 = 140.0μmol/m^2/s)
```

And the environmental conditions should be passed as `Quantity` objects (defaults
are updated accordingly, see Unitful.jl documentation for details on how to
create `Quantity` objects):

```julia
A_c3, gs_c3  = photosynthesis(c3Q, PAR = 100.0μmol/m^2/s)
A_c4, gs_c4  = photosynthesis(c4Q, PAR = 100.0μmol/m^2/s)
```

## Leaf Energy Balance

Ecophys may also compute the leaf energy balance to couple photosynthesis, 
transpiration and leaf temperature. In addition to the models of photosynthesis
and stomatal conductance mentioned in the above, additional models of boundary
layer conductance and leaf optical properties are required. 

Currently, only a simple model of optical properties is avaiable that defines
the leaf absorptance in PAR and NIR and its emmisivity in the thermal domain.
This model is created using the `SimpleOptical()` function (defaults are provided):

```julia
using Ecophys
opt = SimpleOptical(αPAR = 0.80)
```

Two models to compute the boundary layer conductance are available. They differ
in the amount of information used regarding the geometry of the leaf. A simple
model only accounts for the leaf characteristic length and is the most common
approach (as before, a version that supports `Quantity` objects is also available):

```julia
gb = simplegb(d = 0.1)
gbQ = simplegbQ(d = 0.1m)
```

The second model is more complex as it takes into account the aspect ratio (length/width) of
the leaf as well as its inclination angle. It will also distinguish between the
boundary layer conductance of the front and back side of the leaf. This model 
relies on unpublished equations fitted to the data reviewed by 
[Schuepp (1993, New Phyto)](https://nph.onlinelibrary.wiley.com/doi/10.1111/j.1469-8137.1993.tb03898.x):

```julia
gbang = gbAngle(d = 0.1, ang = π/4, ar = 0.1)
gbangQ = gbAngleQ(d = 0.1m, ang = π/4, ar = 0.1)
```

The leaf energy balance is then computed using `solve_energy_balance()` which
will compute the leaf temperature that closes the energy balance as well as the
corresponding CO2 assimilationa and transpiration:

```julia
Tleaf, A, Tr = solve_energy_balance(c3; gb = gb, opt = opt, PAR = 100.0, ws = 5.0)
TleafangQ, AangQ, TrangQ = solve_energy_balance(c3Q; gb = gbangQ, opt = opt, PAR = 100.0μmol/m^2/s, ws = 5.0m/s)
```
