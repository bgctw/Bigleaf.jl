"""
    Monin_Obukhov_length(Tair, pressure, ustar, H; constants)
    Monin_Obukhov_length!(df;constants=BigleafConstants())
    Monin_Obukhov_length(df;constants=BigleafConstants())

calculates the Monin-Obukhov length.

# Arguments
- `Tair`      : Air temperature (degC)
- `pressure`  : Atmospheric pressure (kPa)
- `ustar`     : Friction velocity (m s-1)
- `H`         : Sensible heat flux (W m-2)
- `df`        : DataFrame containing the above variables
optional
- `constants=`[`BigleafConstants`](@ref)`()`: Dictionary with entries 
  - `Kelvin` - conversion degree Celsius to Kelvin 
  - `cp` - specific heat of air for constant pressure (J K-1 1) 
  - `k` - von Karman constant (-) 
  - `g` - gravitational acceleration (m s-2)

# Details
The Monin-Obukhov length (L) is given by:

``L = - (\\rho * cp * ustar^3 * Tair) / (k * g * H)``
 
where ``\\rho`` is air density (kg m-3).

# Note
Note that L gets very small for very low ustar values with implications
for subsequent functions using L as input. It is recommended to filter
data and exclude low ustar values (ustar < ~0.2) beforehand. 

# Value
Monin-Obukhov length L (m). The non-mutating DataFrame variant returns a vector,
the mutating variant add or modifies column `:MOL`.


# References
Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany. 

# See also
[`stability_parameter`](@ref)

```@example; output = false
Monin_Obukhov_length(
  Tair=25,pressure=100,
  ustar=seq(0.2,1,0.1),H=seq(40,200,20))
``` 
"""
function Monin_Obukhov_length(Tair::FT, pressure, ustar, H; constants=BigleafConstants()
  ) where FT
  rho  = air_density(Tair, pressure; constants)
  TairK = Tair + oftype(Tair, constants.Kelvin)
  #MOL  = (-rho*constants.cp*ustar^3*TairK) / (constants.k*constants.g*H)
  MOL  = (-rho*FT(constants.cp)*ustar^3*TairK) / (FT(constants.k)*FT(constants.g)*H)
end
function Monin_Obukhov_length!(df;constants=BigleafConstants())
  fr = (args...) -> Monin_Obukhov_length(args...; constants)
  transform!(df, SA[:Tair, :pressure, :ustar, :H] => ByRow(fr) => :MOL)
end
# function Monin_Obukhov_length(df::DFTable; kwargs...)
#     tmp = Monin_Obukhov_length.(df.Tair, df.pressure, df.ustar, df.H; kwargs...)
# end

"""
    stability_parameter(z,d,MOL)
    stability_parameter(z,d,Tair, pressure, ustar, H; constants)
    stability_parameter!(df::AbstractDataFrame; z,d, MOL=nothing, constants)

calculates stability parameter "zeta", a parameter characterizing stratification in the 
lower atmosphere.

# Arguments             
- `z`   : height (m)
- `d`   : Zero-plane displacement height (m)
- `MOL` : Monin-Obukhov-length L (m)
- `df`  : DataFrame containting the variables required by [`Monin_Obukhov_length`](@ref)
optional
- `constants=`[`BigleafConstants`](@ref)`()`

In the second variant and if `MOL=nothing` in the DataFrame variants, MOL is computed by 
[`Monin_Obukhov_length`](@ref).

# Details
The stability parameter ``\\zeta`` is given by:

``\\zeta = (z - d) / L``

where L is the Monin-Obukhov length (m), calculated by
[`Monin_Obukhov_length`](@ref). The displacement height can 
be estimated from the function [`roughness_parameters`](@ref).
         
# Value
``\\zeta``: stability parameter (-). The nonmutainting DataFrame variant returns a vector.
The mutating variant modifies or adds column [`:zeta`].

```jldoctest; output = false
using DataFrames
df = DataFrame(Tair=25.0, pressure=100.0, ustar=0.2:0.1:1.0, H=40:20.0:200)
z=40;d=15
zeta = stability_parameter.(z,d, df.Tair, df.pressure, df.ustar, df.H)
all(zeta .< 0)
# output
true
``` 
""" 
function stability_parameter(z,d,MOL)
  zeta = (z - d) / MOL
end
function stability_parameter(z,d,Tair, pressure, ustar::FT, H; 
  constants=BigleafConstants()) where FT
  MOL = Monin_Obukhov_length(Tair, pressure, ustar, H; constants);
  stability_parameter(FT(z),FT(d),MOL)
end
function stability_parameter!(df::AbstractDataFrame; z,d, MOL=nothing,
  constants=BigleafConstants())
  if isnothing(MOL)
    ft = (args...) -> stability_parameter.(z,d,args...; constants)
    transform!(df, SA[:Tair, :pressure, :ustar, :H] => ft => :zeta)
  else
    ft = () -> stability_parameter.(z,d,MOL)
    transform!(df, [] => ft => :zeta)
  end
end
# function stability_parameter(df::DFTable; z,d, MOL=nothing,
#   constants=BigleafConstants())
#   if isnothing(MOL) MOL = Monin_Obukhov_length(df; constants); end
#   stability_parameter.(z,d,MOL)
# end

abstract type StabilityCorrectionMethod end
struct Dyer1970 <: StabilityCorrectionMethod end
struct Businger1971 <: StabilityCorrectionMethod end
struct NoStabilityCorrection <: StabilityCorrectionMethod end

"""
    stability_correction(zeta; 
      stab_formulation=Dyer1970())
    stability_correction(z,d, Tair,pressure,ustar,H; constants,
      stab_formulation=Dyer1970())
    stability_correction!(df; zeta, z, d; 
      stab_formulation=Dyer1970(), constants =BigleafConstants())
    
Integrated Stability Correction Functions for Heat and Momentum

# Arguments
- `zeta`             : Stability parameter zeta (-)
- `Tair`,`pressure`,`ustar`,`H` : see [`Monin_Obukhov_length`](@ref)
- `z`,`d`            : see [`stability_parameter`](@ref)
- `df`  : DataFrame containting the variables required by [`Monin_Obukhov_length`](@ref)
- `stab_formulation` : Formulation for the stability function. Either 
            `Dyer1970()`, or `Businger1971()` or `NoStabilityCorrection()`

In the second and third form computes `zeta` by [`stability_parameter`](@ref) and
[`Monin_Obukhov_length`](@ref) and requires respective arguments.

# Details
These dimensionless values are needed to correct deviations
from the exponential wind profile under non-neutral conditions.
The functions give the integrated form of the universal functions. They
depend on the value of the stability parameter ``\\zeta``,
a function of heigh `z`,
which can be calculated from the function [`stability_parameter`](@ref).
The integration of the universal functions is:

``\\psi = -x * \\zeta`` 

for stable atmospheric conditions (``\\zeta`` >= 0), and

``\\psi = 2 * log( (1 + y(\\zeta)) / 2)``

for unstable atmospheric conditions (``\\zeta`` < 0).

The different formulations differ in their value of x and y.
  
# Value
a NamedTuple with the following columns:
- `psi_h`: the value of the stability function for heat and water vapor (-)
- `psi_m`: the value of the stability function for momentum (-)

# References
- Dyer, A_J., 1974: A review of flux-profile relationships. 
  Boundary-Layer Meteorology 7, 363-372.
- Dyer, A. J., Hicks, B_B., 1970: Flux-Gradient relationships in the
  constant flux layer. Quart. J. R. Meteorol. Soc. 96, 715-721.
- Businger, J_A., Wyngaard, J. C., Izumi, I., Bradley, E. F., 1971:
  Flux-Profile relationships in the atmospheric surface layer. 
  J. Atmospheric Sci. 28, 181-189.
- Paulson, C_A., 1970: The mathematical representation of wind speed
  and temperature profiles in the unstable atmospheric surface layer.
  Journal of Applied Meteorology 9, 857-861.
  Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany.

# Examples  
```jldoctest; output = false
using DataFrames
zeta = -2:0.5:0.5
df2 = DataFrame(stability_correction.(zeta; stab_formulation=Businger1971()))                         
propertynames(df2) == [:psi_h, :psi_m]
# output
true
``` 
"""            
function stability_correction(zeta::FT; stab_formulation=Dyer1970()) where FT
  # integration of universal functions (after Paulson1970 and Foken 2008)
  stab_formulation isa NoStabilityCorrection && return(
    (psi_h = zero(zeta), psi_m = zero(zeta)))
  ismissing(zeta) && return((psi_h = missing, psi_m = missing))
  is_stable = zeta >= 0 
  if is_stable
    x_h::FT, x_m::FT = get_stability_coefs_stable(stab_formulation)
    psi_h = x_h * zeta
    psi_m = x_m * zeta
  else
    y_h::FT, y_m::FT = get_stability_coefs_unstable(stab_formulation, zeta)
    psi_h = 2 * log( (1 + y_h ) / 2)
    psi_m = 2 * log( (1 + y_m ) / 2) +
                      log( ( 1 + y_m^2 ) / 2)
                      -2 * atan(y_m) + pi/2
  end
  (;psi_h, psi_m)
end 

get_stability_coefs_stable(::Businger1971) = (x_h = -7.8, x_m = -6)
function get_stability_coefs_unstable(::Businger1971, zeta)
  y_h = 0.95 * ( 1 - 11.6 * zeta)^0.5
  y_m = (1 - 19.3*zeta)^0.25
  (;y_h, y_m)
end

get_stability_coefs_stable(::Dyer1970) = (x_h = -5, x_m = -5)
function get_stability_coefs_unstable(::Dyer1970, zeta)
  y_h       = (1 - 16 * zeta)^0.5
  y_m       = (1 - 16 * zeta)^0.25
  (;y_h, y_m)
end

function stability_correction(z,d, Tair::Union{Missing,Number},pressure,ustar::FT,H; 
  stab_formulation=Dyer1970(), constants=BigleafConstants()) where FT
  stab_formulation isa NoStabilityCorrection && return(
    (psi_h = zero(FT), psi_m = zero(FT)))
  FT <: Missing && return(missing)
  MOL = Monin_Obukhov_length(Tair,pressure,ustar,H; constants)
  zeta  = stability_parameter(FT(z),FT(d),MOL)
  psis = stability_correction(zeta; stab_formulation)
  psis
end
# function stability_correction(Tair,pressure,ustar,H, z,d; kwargs...) 
#   Tables.columns(
#     stability_correction.(Tair,pressure,ustar,H, z,d; kwargs...))
# end
# function stability_correction(df::DFTable, z,d; kwargs...)   
#   tmp1 = stability_correction.(df.Tair, df.pressure, df.ustar, df.H, z, d; kwargs...)
#   tmp2 = Tables.columns(tmp1)
# end
function stability_correction!(df; zeta=nothing, z=nothing, d=nothing, 
  stab_formulation=Dyer1970(), constants =BigleafConstants())
  # cannot dispatch on keyword argument, hence need if-clause
  if stab_formulation isa NoStabilityCorrection
    zeroT = :ustar in propertynames(df) ? zero(eltype(skipmissing(df.ustar))) : zero(z)
    df[!,:psi_h] .= zeroT
    df[!,:psi_m] .= zeroT
    return(df)
  end
  if !isnothing(zeta)
    ft = () -> stability_correction.(zeta; stab_formulation)
    transform!(df, [] => ft => AsTable)
  else
    ft = (args...) -> stability_correction.(z, d, args...; stab_formulation, constants)
    transform!(df, SA[:Tair,:pressure,:ustar,:H] => ft => AsTable)
  end
end

function stability_correction(df::DFTable; z, d, 
  stab_formulation=Dyer1970(), constants =BigleafConstants()
  )
  # do not provide zeta, because can simply invoke stability_correction.(zeta)
  if stab_formulation isa NoStabilityCorrection
    zeroT = zero(eltype(skipmissing(df.ustar)))
    rows = map(Tables.rows(df)) do row
      (psi_h = zeroT, psi_m = zeroT)
    end
    return(Tables.columns(rows))
  end
  Tables.columns(stability_correction.(
    z, d, df.Tair, df.pressure, df.ustar, df.H; stab_formulation, constants))
end


