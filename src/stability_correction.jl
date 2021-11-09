"""
    Monin_Obukhov_length(Tair, pressure, ustar, H; constants)
    Monin_Obukhov_length!(df;constants=bigleaf_constants())
    Monin_Obukhov_length(df;constants=bigleaf_constants())

calculates the Monin-Obukhov length.

# Arguments
- `Tair`      : Air temperature (degC)
- `pressure`  : Atmospheric pressure (kPa)
- `ustar`     : Friction velocity (m s-1)
- `H`         : Sensible heat flux (W m-2)
- `df`        : DataFrame containing the above variables
optional
- `constants=`[`bigleaf_constants`](@ref)`()`: Dictionary with entries 
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
function Monin_Obukhov_length(Tair, pressure, ustar, H; constants=bigleaf_constants())
  rho  = air_density(Tair, pressure; constants)
  TairK = Tair + constants[:Kelvin]
  MOL  = (-rho*constants[:cp]*ustar^3*TairK) / (constants[:k]*constants[:g]*H)
end
function Monin_Obukhov_length!(df;constants=bigleaf_constants())
  ft(args...) = Monin_Obukhov_length(args...; constants)
  transform!(df, SA[:Tair, :pressure, :ustar, :H] => ByRow(ft) => :MOL)
end
function Monin_Obukhov_length(df::DFTable; kwargs...)
    tmp = Monin_Obukhov_length.(df.Tair, df.pressure, df.ustar, df.H; kwargs...)
    # eltype_agg = Union{eltype(df.Tair), eltype(df.pressure), eltype(df.ustar), eltype(df.H)}
    # convert(Vector{eltype_agg}, tmp)
end

"""
    stability_parameter(zr,d,MOL)
    stability_parameter!(df::AbstractDataFrame; zr,d, MOL=nothing, constants)
    stability_parameter(df::AbstractDataFrame; zr,d, MOL=nothing, constants)

calculates stability parameter "zeta", a parameter characterizing stratification in the 
lower atmosphere.

# Arguments             
- `zr`  : Instrument (reference) height (m)
- `d`   : Zero-plane displacement height (m)
- `MOL` : Monin-Obukhov-length L (m)
- `df`  : DataFrame containting the variables required by [`Monin_Obukhov_length`](@ref)
optional
- `constants=`[`bigleaf_constants`](@ref)`()`

If `MOL=nothing` in the DataFrame variants, it is computed by 
[`Monin_Obukhov_length`](@ref) using `df` and `constants`.

# Details
The stability parameter ``\\zeta`` is given by:

``\\zeta = (zr - d) / L``

where L is the Monin-Obukhov length (m), calculated from the ction
[`Monin_Obukhov_length`](@ref). The displacement height can 
be estimated from the function [`roughness_parameters`](@ref).
         
# Value
``\\zeta``: stability parameter (-). The nonmutainting DataFrame variant returns a vector.
The mutating variant modifies or adds column [`:zeta`].

```jldoctest; output = false
using DataFrames
df = DataFrame(Tair=25, pressure=100, ustar=0.2:0.1:1.0, H=40:20:200)
zeta = stability_parameter(df;zr=40,d=15)
all(zeta .< 0)
# output
true
``` 
""" 
function stability_parameter(zr,d,MOL)
  zeta = (zr - d) / MOL
end
function stability_parameter!(df::AbstractDataFrame; zr,d, MOL=nothing,
  constants=bigleaf_constants())
  df[!,:zeta] .= stability_parameter(df; zr,d,MOL,constants)
  df
end
function stability_parameter(df::DFTable; zr,d, MOL=nothing,
  constants=bigleaf_constants())
  if isnothing(MOL) MOL = Monin_Obukhov_length(df; constants); end
  stability_parameter.(zr,d,MOL)
end

"""
    stability_correction(zeta; 
      stab_formulation=Val(:Dyer_1970))
    stability_correction(Tair,pressure,ustar,H, z,d; constants,
      stab_formulation=Val(:Dyer_1970))
    stability_correction!(df, z, d; 
      stab_formulation=Val(:Dyer_1970), constants = bigleaf_constants())
    
Integrated Stability Correction Functions for Heat and Momentum

# Arguments
- `zeta`             : Stability parameter zeta (-)
- `stab_formulation` : Formulation for the stability function. Either 
            `Val(:Dyer_1970)`, or `Val(:Businger_1971)` or `Val(:no_stability_correction)`
- `Tair`,`pressure`,`ustar`,`H` : see [`Monin_Obukhov_length`](@ref)
- `z`,`d`            : see [`stability_parameter`](@ref)
- `df`  : DataFrame containting the variables required by [`Monin_Obukhov_length`](@ref)

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

``\\psi = -x * zeta`` 

for stable atmospheric conditions (``\\zeta`` >= 0), and

``\\psi = 2 * log( (1 + y(zeta)) / 2)``

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
df2 = DataFrame(stability_correction.(zeta; stab_formulation=Val(:Businger_1971)))                         
propertynames(df2) == [:psi_h, :psi_m]
# output
true
``` 
"""            
function stability_correction(zeta; stab_formulation=Val(:Dyer_1970))
  # integration of universal functions (after Paulson_1970 and Foken 2008)
  stab_formulation isa Val{:no_stability_correction} && return(
    (psi_h = zero(zeta), psi_m = zero(zeta)))
  ismissing(zeta) && return((psi_h = missing, psi_m = missing))
  is_stable = zeta >= 0 
  if is_stable
    x_h, x_m = get_stability_coefs_stable(stab_formulation)
    psi_h = x_h * zeta
    psi_m = x_m * zeta
  else
    y_h, y_m = get_stability_coefs_unstable(stab_formulation, zeta)
    psi_h = 2 * log( (1 + y_h ) / 2)
    psi_m = 2 * log( (1 + y_m ) / 2) +
                      log( ( 1 + y_m^2 ) / 2)
                      -2 * atan(y_m) + pi/2
  end
  (;psi_h, psi_m)
end 

get_stability_coefs_stable(::Val{:Businger_1971}) = (x_h = -7.8, x_m = -6)
function get_stability_coefs_unstable(::Val{:Businger_1971}, zeta)
  y_h = 0.95 * ( 1 - 11.6 * zeta)^0.5
  y_m = (1 - 19.3*zeta)^0.25
  (;y_h, y_m)
end

get_stability_coefs_stable(::Val{:Dyer_1970}) = (x_h = -5, x_m = -5)
function get_stability_coefs_unstable(::Val{:Dyer_1970}, zeta)
  y_h       = (1 - 16 * zeta)^0.5
  y_m       = (1 - 16 * zeta)^0.25
  (;y_h, y_m)
end

#TODO z or zr here?
function stability_correction(Tair::Union{Missing,Number},pressure,ustar,H, z,d; 
  stab_formulation=Val(:Dyer_1970), constants = bigleaf_constants())
  stab_formulation isa Val{:no_stability_correction} && return(
    (psi_h = zero(z), psi_m = zero(z)))
  MOL = Monin_Obukhov_length(Tair,pressure,ustar,H; constants)
  zeta  = stability_parameter(z,d,MOL)
  stability_correction(zeta; stab_formulation)
end
function stability_correction(Tair,pressure,ustar,H, z,d; kwargs...) 
  Tables.columns(
    stability_correction.(Tair,pressure,ustar,H, z,d; kwargs...))
end
function stability_correction(df::DFTable, z,d; kwargs...)   
  tmp1 = stability_correction.(df.Tair, df.pressure, df.ustar, df.H, z, d; kwargs...)
  tmp2 = Tables.columns(tmp1)
end

function stability_correction!(df; zeta=df.zeta, 
  stab_formulation=Val(:Dyer_1970), constants = bigleaf_constants())
  # cannot dispatch on keyword argument, hence need if-clause
  if stab_formulation isa Val{:no_stability_correction}
    df[!,:psi_h] .= 0.0
    df[!,:psi_m] .= 0.0
    return(df)
  end
  ft() = stability_correction.(zeta; stab_formulation)
  transform!(df, [] => ft => AsTable)
end

function stability_correction!(df, z, d; 
  stab_formulation=Val(:Dyer_1970), constants = bigleaf_constants())
  ft(args...) = stability_correction.(args..., z, d; stab_formulation, constants)
  transform!(df, SA[:Tair,:pressure,:ustar,:H] => ft => AsTable)
end
