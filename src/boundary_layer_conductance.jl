"""
    compute_Gb!(df::AbstractDataFrame, approach; kwargs...)

Estimate boundary layer conductance.

# Arguments  
For `Val(:Thom_1972)`
- `df`       : DataFrame with required variables (depend on approach)
- `approach` : one of
  - `Val(:Thom_1972)`: see [`Gb_Thom`](@ref)
  - `Val(:Choudhury_1988)`: see [`Gb_Choudhury`](@ref)
  - `Val(:Su_2001)`: see [`Gb_Su`](@ref)
  - `Val(:constant_kB1)`: see [`Gb_constant_kB1`](@ref)
 
The different approaches required different variables present in `df` and
different keyword arguments.

# Value
updated DataFrame `df` with the following columns:
- `Gb_h`: Boundary layer conductance for heat transfer (m s-1)

To subsequently derived quantities see
- [`compute_Gb_quantities`](@ref) for Resistance, kB-1 constant, and CO2 conductance
- [`add_Gb!`](@ref) for conductances of other species given their Schmidt numbers.

# See also
[`Gb_Thom`](@ref), `Gb_Choudhury`](@ref), [`Gb_Su`](@ref), `Gb_Choudhury`](@ref), 
[`aerodynamic_conductance!`](@ref)

```jldoctest; output = false
using DataFrames
df = DataFrame(Tair=25,pressure=100,wind=[3,4,5],
  ustar=[0.5,0.6,0.65],H=[200,230,250]) 
leafwidth=0.01;LAI=5;zh=25;zr=40
compute_Gb!(df, Val(:Choudhury_1988); leafwidth, LAI, zh, zr)
≈(df.Gb_h[1], 0.379, rtol=1e-3)
# output
true
``` 
"""
function compute_Gb!(df::AbstractDataFrame, approach::Val{:Thom_1972}; kwargs...)
  ft(ustar) = Gb_Thom(ustar; kwargs...)
  transform!(df, :ustar => ByRow(ft) => :Gb_h)
end
function compute_Gb!(df::AbstractDataFrame, approach::Val{:constant_kB1}; kB_h, kwargs...)
  # do not use ByRow because kb_H can be a vector
  ft(ustar) = Gb_constant_kB1.(ustar, kB_h; kwargs...)
  transform!(df, :ustar => ft => :Gb_h)
end
function compute_Gb!(df::AbstractDataFrame, approach::Val{:Choudhury_1988}; kwargs...)
  Gb_Choudhury!(df; kwargs...)
end
function compute_Gb!(df::AbstractDataFrame, approach::Val{:Su_2001}; kwargs...)
  Gb_Su!(df; kwargs...)
end

"""
    compute_Gb_quantities(Gb_h)
    compute_Gb_quantities!(df:::AbstractDataFrame)

Based on boundary layer conductance for heat, compute derived quantities.

# Arguments
- `Gb_h` : Boundary layer conductance for heat transfer (m s-1)
- `df`   : DataFrame with above columns
- `constants=`[`bigleaf_constants`](@ref)`()`: entries `Sc_CO2` and `Pr`

# Value
NamedTuple with entries
- `Gb_h`: Boundary layer conductance for heat transfer (m s-1)
- `Rb_h`: Boundary layer resistance for heat transfer (s m-1)
- `kB_h`: kB-1 parameter for heat transfer
- `Gb_CO2`: Boundary layer conductance for CO2 (m s-1). 
"""
function compute_Gb_quantities(Gb_h, ustar; constants=bigleaf_constants())
  Rb_h = 1/Gb_h
  kB_h = Rb_h*constants[:k]*ustar
  Gb_CO2 = Gb_h / (constants[:Sc_CO2]/constants[:Pr])^0.67
  (;Rb_h, Gb_h, kB_h, Gb_CO2)
end
function compute_Gb_quantities!(df::AbstractDataFrame; constants=bigleaf_constants())
  ft(args...) = compute_Gb_quantities.(args...; constants)
  transform!(df, SA[:Gb_h, :ustar] => ft => AsTable)
end


"""
    add_Gb(Gb_h::Union{Missing,Number}, Sc::Vararg{Pair,N}; constants)
    add_Gb!(df::AbstractDataFrame, Sc::Vararg{Pair,N}; Gb_h = df.Gb_h, kwargs...) 

compute additional boundary layer conductance quantities for given Schmidt-numbers

# Arguments  
- `Gb_h`     : Boundary layer conductance for heat transfer (m s-1)
- `Sc`       : several `Pair{Symbol,Number}` Output name and Schmidt number of 
               additional conductances to be calculated
- `df`       : DataFrame to add output columns               
optional
- `constants=`[`bigleaf_constants`](@ref)`()`: Dictionary with entries 
  - `Pr` - Prandtl number 

# Details
Boundary layer conductance for other quantities x is calculated 
based on boundary layer for heat transfer as (Hicks et al. 1987):
 
  ``Gb_x = Gb_h / (Sc_x / Pr)^{0.67}``
   
where `Sc_x` is the Schmidt number of quantity x, and Pr is the Prandtl number (0.71).

# Value
a NameTuple or `df` with keys `Gb_x` where `x` are the keys in `Sc` and 
corresponding boundary layer conductances (m s-1).

# Examples
```jldoctest; output=false
using DataFrames
df = DataFrame(Gb_h=[0.02, missing, 0.055])
add_Gb!(df, :O2 => 0.84, :CH4 => 0.99)
propertynames(df)[2:3] == [:Gb_O2, :Gb_CH4]
# output
true
```
"""
function add_Gb!(df::AbstractDataFrame, Sc::Vararg{Pair,N}; Gb_h = df.Gb_h, kwargs...) where N
  N == 0 && return(df)
  Scn, Scv = get_names_and_values("Gb_", Sc...)
  ft() = add_Gb_.(Gb_h, Ref(Scn), Ref(Scv); kwargs...)
  transform!(df, [] => ft => AsTable)
end
function add_Gb(Gb_h::Union{Missing,Number}, Sc::Vararg{Pair,N}; kwargs...) where N
  Scn, Scv = get_names_and_values("Gb_", Sc...)
  add_Gb_(Gb_h, Scn, Scv; kwargs...)
end
function add_Gb_(Gb_h::Union{Missing,Number}, Scn::NTuple{N,Symbol}, Scv::NTuple{N}; 
  constants=bigleaf_constants()) where N 
  Gbxv = @. Gb_h / (Scv/constants[:Pr])^0.67
  Gbx = NamedTuple{Scn}(Gbxv)
end
function get_names_and_values(prefix::AbstractString, Sc::Vararg{Pair,N}) where N
  Scn = ntuple(i -> Symbol(prefix * string(Sc[i].first)), N)
  Scv = ntuple(i -> Sc[i].second, N)
  Scn, Scv
end



"""
    Gb_Thom(ustar; constants)
    compute_Gb!(df, Val{:Thom_1972})

Boundary Layer Conductance according to Thom 1972, an empirical formulation for the 
for heat transfer based on a simple ustar (friction velocity) dependency.

# Arguments  
- `ustar`     : Friction velocity (m s-1)
- `df`        : DataFrame with above variables
- `constants=`[`bigleaf_constants`](@ref)`()`
 
# Details
The empirical equation for Rb suggested by Thom 1972 is:
 
``Rb = 6.2 ustar^{-0.67}``
 
Gb (=1/Rb) for water vapor and heat are assumed to be equal in this package.

# Value
see [`compute_Gb!`](@ref)
 
# References
- Thom, A., 1972: Momentum, mass and heat exchange of vegetation.
  Quarterly Journal of the Royal Meteorological Society 98, 124-134.
- Hicks, B_B., Baldocchi, D_D., Meyers, T_P., Hosker, J_R., Matt, D_R., 1987:
  A preliminary multiple resistance routine for deriving dry deposition velocities
  from measured quantities. Water, Air, and Soil Pollution 36, 311-330.

```@example; output = false
using DataFrames
df = DataFrame(ustar = SA[0.1,missing,0.3])
compute_Gb!(df, Val(:Thom_1972))
propertynames(df) == [:ustar, :Rb_h, :Gb_h, :kB_h, :Gb_CO2]
``` 
"""
function Gb_Thom(ustar::Union{Missing,Number}; constants=bigleaf_constants())
  Rb_h = 6.2*ustar^-0.667
  Gb_h = 1/Rb_h
end

"""
    Gb_constant_kB1(ustar, kB_h; constants)
    compute_Gb!(df, Val{:constant_kB1})

Boundary Layer Conductance using constant kB-1 value for heat transfer.

# Arguments  
- `ustar`     : Friction velocity (m s-1)
- `df`        : DataFrame with above variables
- `kB_h`      : kB-1 value for heat transfer
- `constants=`[`bigleaf_constants`](@ref)`()`
 
# Details
Rb_h computed by ``kB_h/(k * ustar)``, where k is the von Karman constant.
"""
function Gb_constant_kB1(ustar, kB_h; constants=bigleaf_constants())
  Rb_h = kB_h/(constants[:k] * ustar)
  Gb_h = 1/Rb_h
  # Gb_CO2 = Gb_h / (constants[:Sc_CO2]/constants[:Pr])^0.67
  # (;Rb_h, Gb_h, kB_h, Gb_CO2)
end
  
"""
    Gb_Choudhury(ustar; leafwidth, LAI, wind_zh, constants)
      alpha   = 4.39 - 3.97*exp(-0.258*LAI)
    Gb_Choudhury!(df::AbstractDataFrame; leafwidth, LAI, wind_zh=nothing, zh, zr, d=0.7*zh,
      z0m = nothing, stab_formulation=Val(:Dyer_1970), constants)

Estimate the canopy boundary layer conductance 
for heat transfer according to Choudhury & Monteith 1988.

# Arguments  
- `Tair`             : Air temperature (degC)
- `pressure`         : Atmospheric pressure (kPa)
- `wind`             : Wind speed at sensor height (m s-1)
- `ustar`            : Friction velocity (m s-1)
- `H`                : Sensible heat flux (W m-2)
- `df`               : DataFrame or matrix containing the above variables
- `leafwidth`        : Leaf width (m)
- `LAI`              : One-sided leaf area index
- `wind_zh`          : Wind speed at canopy heihgt (m s-1), if not given then computed by
                       [`wind_profile`](@ref) using the arguments below
- `zh`               : Canopy height (m)
- `zr`               : Instrument (reference) height (m)
- `d`                : Zero-plane displacement height (-), can be calculated using 
                      `roughness_parameters`
- `z0m`              : Roughness length for momentum (m). If not provided, calculated 
                        from `roughness_parameters`  within `wind_profile` using zh and zr
- `stab_formulation` : Stability correction function used, 
                       see [`stability_correction`](@ref).
- `constants=`[`bigleaf_constants`](@ref)`()`
                        
# Value
see [`compute_Gb!`](@ref)
 
# Details
Boundary layer conductance according to Choudhury & Monteith 1988 is
given by:

``Gb_h = LAI((2a/\\alpha)*\\sqrt(u(zh)/w)*(1-exp(-\\alpha/2)))``

where ``\\alpha`` is modeled as an empirical relation to LAI (McNaughton & van den Hurk 1995):

``\\alpha = 4.39 - 3.97*exp(-0.258*LAI)``

``w`` is leafwidth and ``u(zh)`` is the wind speed at the canopy surface.

It can be approximated from measured wind speed at sensor height zr and a wind extinction 
coefficient ``\\alpha``: ``u(zh) = u(zr) / (exp(\\alpha(zr/zh -1)))``.
However, here (if not explicitly given) it is estimated by [`wind_profile`](@ref)

Gb (=1/Rb) for water vapor and heat are assumed to be equal in this package.
Gb for other quantities x is calculated as (Hicks et al. 1987):
``Gb_x = Gb / (Sc_x / Pr)^0.67``
where Sc_x is the Schmidt number of quantity x, and Pr is the Prandtl number (0.71).
         
# References
- Choudhury, B. J., Monteith J_L., 1988: A four-layer model for the heat
  budget of homogeneous land surfaces. Q. J. R. Meteorol. Soc. 114, 373-398.
- McNaughton, K. G., Van den Hurk, B_J_J_M., 1995: A 'Lagrangian' revision of
  the resistors in the two-layer model for calculating the energy budget of a
  plant canopy. Boundary-Layer Meteorology 74, 261-288.
- Hicks, B_B., Baldocchi, D_D., Meyers, T_P., Hosker, J_R., Matt, D_R., 1987:
  A preliminary multiple resistance routine for deriving dry deposition velocities
  from measured quantities. Water, Air, and Soil Pollution 36, 311-330.
"""
function Gb_Choudhury(; leafwidth, LAI, wind_zh, constants=bigleaf_constants())
  alpha   = 4.39 - 3.97*exp(-0.258*LAI)
  # (ismissing(wind_zh) || isnothing(wind_zh)) && return(
  #   (Rb_h=missing, Gb_h=missing, kB_h=missing, Gb_CO2=missing))
  wind_zh = max(0.01, wind_zh) ## avoid zero windspeed
  Gb_h = LAI*((0.02/alpha)*sqrt(wind_zh/leafwidth)*(1-exp(-alpha/2)))
  # Rb_h = 1/Gb_h
  # kB_h = Rb_h*constants[:k]*ustar
  # Gb_CO2 = Gb_h / (constants[:Sc_CO2]/constants[:Pr])^0.67
  # (;Rb_h, Gb_h, kB_h, Gb_CO2)
end
function Gb_Choudhury!(df::AbstractDataFrame; leafwidth, LAI, wind_zh, 
  #zh, zr, d=0.7*zh, z0m = nothing, stab_formulation=Val(:Dyer_1970),
  constants=bigleaf_constants()
  )
  # if isnothing(wind_zh)
  #   wind_zh = wind_profile(zh, df, d, z0m; zh, zr, stab_formulation, constants)
  # end
  # Broadcasting does not work over keyword arguments, need to pass as positional
  fwind(wind_zh, leafwidth, LAI; kwargs...) = Gb_Choudhury(;wind_zh, leafwidth, LAI, kwargs...)
  ft() = fwind.(wind_zh, leafwidth, LAI; constants)
  transform!(df, [] => ft => :Gb_h)
end

"""
    Gb_Su(Tair,pressure,ustar; zh, wind_zh, Dl, fc, N=2, Cd=0.2, hs=0.01, constants)
    Gb_Su!(df; wind_zh=nothing, Dl, fc=nothing, N=2, Cd=0.2, hs=0.01,
      z0m = nothing, zh, zr, d = 0.7*zh, LAI, stab_formulation=Val(:Dyer_1970), constants)

Estimate Boundary Layer Conductance to heat transfer using the physically based 
formulation according to Su et al. 2001. 

# Arguments
- `Tair`      : Air temperature (degC)
- `pressure`  : Atmospheric pressure (kPa)
- `ustar`     : Friction velocity (m s-1)
- `H`         : Sensible heat flux (W m-2)
- `wind`      : Wind speed at sensor height (m s-1)
- `df`        : DataFrame or matrix containing the above variables
- `Dl`        : Leaf characteristic dimension (m)
- `fc`        : Fractional vegetation cover [0-1] (if not provided, calculated from LAI)
- `LAI`       : One-sided leaf area index (-)
- `N`         : Number of leaf sides participating in heat exchange (defaults to 2)
- `Cd`        : Foliage drag coefficient (-)
- `hs`        : Roughness height of the soil (m)
- `LAI`       : One-sided leaf area index
- `wind_zh`   : Wind speed at canopy heihgt (m s-1) (if not given then computed by
                       [`wind_profile`](@ref) using the arguments below)
- `zh`        : Canopy height (m)
- `zr`        : Instrument (reference) height (m)
- `d`         : Zero-plane displacement height (-), can be calculated using 
         `roughness_parameters`
- `z0m`       : Roughness length for momentum (m). If not provided, calculated 
           from `roughness_parameters`  within `wind_profile` using zh and zr
- `stab_formulation` : Stability correction function used, 
                       see [`stability_correction`](@ref).
- `constants=`[`bigleaf_constants`](@ref)`()`

# Value
see [`compute_Gb!`](@ref)
    
# Details
The formulation is based on the kB-1 model developed by Massman 1999. 
Su et al. 2001 derived the following approximation:
 
``kB-1 = (k Cd fc^2) / (4Ct ustar/u(zh)) + kBs-1(1 - fc)^2``

If fc (fractional vegetation cover) is missing, it is estimated from LAI:
``fc = 1 - exp(-LAI/2)``

The wind speed at the top of the canopy is calculated using function
[`wind_profile`](@ref).

Ct is the heat transfer coefficient of the leaf (Massman 1999):

``Ct = Pr^-2/3 Reh^-1/2 N``

where Pr is the Prandtl number (set to 0.71), and Reh is the Reynolds number for leaves:

``Reh = Dl wind(zh) / v``
 
kBs-1, the kB-1 value for bare soil surface, is calculated according 
to Su et al. 2001:

``kBs^-1 = 2.46(Re)^0.25 - ln(7.4)``

Gb (=1/Rb) for water vapor and heat are assumed to be equal in this package.
Gb for other quantities x is calculated as (Hicks et al. 1987):
  ``Gb_x = Gb / (Sc_x / Pr)^0.67``
where Sc_x is the Schmidt number of quantity x, and Pr is the Prandtl number (0.71).

# References
- Su, Z., Schmugge, T., Kustas, W. & Massman, W., 2001: An evaluation of
  two models for estimation of the roughness height for heat transfer between
  the land surface and the atmosphere. Journal of Applied Meteorology 40, 1933-1951.
-  Massman, W., 1999: A model study of kB H- 1 for vegetated surfaces using
  localized near-field' Lagrangian theory. Journal of Hydrology 223, 27-43.
-  Hicks, B_B., Baldocchi, D_D., Meyers, T_P., Hosker, J_R., Matt, D_R., 1987:
   A preliminary multiple resistance routine for deriving dry deposition velocities
  from measured quantities. Water, Air, and Soil Pollution 36, 311-330.

```jldoctest; output = false
using DataFrames
df = DataFrame(Tair=25,pressure=100,wind=[3,4,5],ustar=[0.5,0.6,0.65],H=[200,230,250]) 
compute_Gb!(df,Val(:Su_2001), zh=25,zr=40,Dl=0.01,LAI=5)
# the same meteorological conditions, but larger leaves
compute_Gb!(df,Val(:Su_2001), zh=25,zr=40,Dl=0.1,LAI=5)
# same conditions, large leaves, and sparse canopy cover (LAI = 1.5)
compute_Gb!(df,Val(:Su_2001), zh=25,zr=40,Dl=0.1,LAI=1.5)
≈(df.Gb_h[1], 0.0638, rtol=1e-3)
# output
true
``` 
"""
function Gb_Su(Tair,pressure,ustar; wind_zh, Dl, fc, N=2, Cd=0.2, hs=0.01,
  constants=bigleaf_constants()
  )
  wind_zh = max(0.01,wind_zh) ## avoid zero windspeed
  v   = kinematic_viscosity(Tair,pressure;constants)
  Re  = Reynolds_Number(Tair,pressure,ustar,hs; constants)
  kBs = 2.46 * (Re)^0.25 - log(7.4)
  Reh = Dl * wind_zh / v
  Ct  = 1*constants[:Pr]^-0.6667*Reh^-0.5*N
  #
  kB_h = (constants[:k]*Cd)/(4*Ct*ustar/wind_zh)*fc^2 + kBs*(1 - fc)^2
  Rb_h = kB_h/(constants[:k]*ustar)
  Gb_h = 1/Rb_h
end
function Gb_Su!(df::AbstractDataFrame; wind_zh, Dl, fc=nothing, 
  N=2, Cd=0.2, hs=0.01, LAI,
#  zh, zr, d = 0.7*zh, stab_formulation=Val(:Dyer_1970), 
  constants=bigleaf_constants()
  )
  # if isnothing(wind_zh)
  #   wind_zh = wind_profile(df, zh, d, z0m; zh, zr, stab_formulation, constants)
  # end
  if isnothing(fc)
    isnothing(LAI) && error("one of 'fc' or 'LAI' must be provided")
    fc = length(LAI) == 1 ? (1-exp(-LAI/2)) : @. (1-exp(-LAI/2))
  end
  isnothing(Dl) && error(
    "need to provide keyword argument Dl with :Su_2001 method")
  inputcols = SA[:Tair,:pressure,:ustar]
  # Broadcasting does not work over keyword arguments, need to pass as positional
  fwind(wind_zh, Dl, fc, N, Cd, hs, args...; kwargs...) = Gb_Su(args...; wind_zh, Dl, fc, N, Cd, hs, kwargs...)
  ft(args...) = fwind.(wind_zh, Dl, fc, N, Cd, hs, args...; constants)
  transform!(df, inputcols => ft => :Gb_h)
end
