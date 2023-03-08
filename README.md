# BigLeaf

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bgctw.github.io/BigLeaf.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://bgctw.github.io/BigLeaf.jl/dev)
[![Build Status](https://github.com/bgctw/BigLeaf.jl/workflows/CI/badge.svg)](https://github.com/bgctw/BigLeaf.jl/actions)
[![Coverage](https://codecov.io/gh/bgctw/BigLeaf.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/bgctw/BigLeaf.jl)


**BigLeaf.jl** is a partial Julia port of Jürgen Knauer's 
[**bigleaf** R package](https://bitbucket.org/juergenknauer/BigLeaf) 
for the calculation of physical (e.g. aerodynamic conductance, surface temperature) 
and physiological (e.g. canopy conductance, water-use efficiency) ecosystem properties 
from eddy covariance data and accompanying meteorological measurements. 
All calculations are based on a 'big-leaf' representation of the vegetation and 
return representative bulk ecosystem/canopy variables.

Earlier, the package was named **Bigleaf.jl** (lowercase leaf). With the renaming
to **BigLeaf.jl**, a new UUID has been generated. From user's perspective, this
is a new package. From developers perspective its the same project with a breaking
change documented in the version number. The lowercase versions can 
still be installed via the General Julia registry and link to older releases
at this repository.

## Citation
Knauer J, El-Madany TS, Zaehle S, Migliavacca M (2018) BigLeaf—An R package for the calculation of physical and physiological ecosystem properties from eddy covariance data.
PLoS ONE 13(8): e0201114. [https://doi.org/10.1371/journal.pone.0201114](https://doi.org/10.1371/journal.pone.0201114)


## Installation and Loading

The `BigLeaf.jl` package can be installed with the usual command once:

```julia
using Pkg
Pkg.add("BigLeaf")
```

And then importet to the every Julia session by:
```julia
using BigLeaf
```

## Usage
See the [Documentation](https://bgctw.github.io/BigLeaf.jl/dev/) 
that includes a [walkthrough](https://bgctw.github.io/BigLeaf.jl/dev/walkthrough/).

[Please report bugs or issues here](https://github.com/bgctw/BigLeaf.jl/issues)

## Package content 
We are porting functionality of the [R package](https://bitbucket.org/juergenknauer/BigLeaf) as needed. Please
file an [issue](https://github.com/bgctw/BigLeaf.jl/issues) if you need a specific feature.

At the current state we ported
- Meteorological variables
- Boundary layer and Aerodynamic conductance
- Surface conductance
- Evapotranspiration
- Potential radiation
- Unit conversions
- Filtering data

