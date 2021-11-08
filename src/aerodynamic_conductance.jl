"""
Aerodynamic Conductance

Bulk aerodynamic conductance, including options for the boundary layer conductance
formulation and stability correction functions.

# Arguments
- `data`              : DataFrame or matrix containing all required variables
- `Tair`              : Air temperature (deg C)
- `pressure`          : Atmospheric pressure (kPa)
- `wind`              : Wind speed (m s-1)
- `ustar`             : Friction velocity (m s-1)
- `H`                 : Sensible heat flux (W m-2)
- `zr`                : Instrument (reference) height (m)
- `zh`                : Canopy height (m)
- `d`                 : Zero-plane displacement height (m)
- `z0m`               : Roughness length for momentum (m), optional; if not provided, 
                      it is estimated from `roughness_parameters`
                      TODO : (method=Val(:wind_profile)). Only used if `wind_profile = true` and/or 
                      `Rb_model` = `Val(:Su_2001)` or
                      : `Val(:Choudhury_1988)`.
- `Dl`                : Characteristic leaf dimension (m) (if `Rb_model` = `Val(:Su_2001)`) 
                      :    or leaf width (if `Rb_model` = `Val(:Choudhury_1988)`); 
                      ignored otherwise.
- `N`                 : Number of leaf sides participating in heat exchange (1 or 2); 
                      only used if `Rb_model = Val(:Su_2001)`.
                      :    Defaults to 2.
- `fc`                : Fractional vegetation cover (-); only used if `Rb_model = Val(:Su_2001)`. 
                      See Details.
- `LAI`               : One-sided leaf area index (m2 m-2); only used if 
                      `Rb_model` = `Val(:Choudhury_1988)` or `Val(:Su_2001)`.
- `Cd`                : Foliage drag coefficient (-); only used if `Rb_model = Val(:Su_2001)`. 
- `hs`                : Roughness length of bare soil (m); only used if `Rb_model = Val(:Su_2001)`.
- `wind_profile`      : Should Ga for momentum be calculated based on the logarithmic wind 
                      profile equation?  Defaults to `false`.
- `stab_correction`   : Should stability correction be applied? Defaults to `true`. 
                      Ignored if `wind_profile = false`.                         
- `stab_formulation`  : Stability correction function. Either `Val(:Dyer_1970)` (default) or
                      :    `Val(:Businger_1971)`. Ignored if `wind_profile = false` or if 
                      `stab_correction = false`.
- `Rb_model`          : Boundary layer resistance formulation. 
                      One of `Val(:Thom_1972),Val(:Choudhury_1988),Val(:Su_2001),Val(:constant_kB1)`.
- `kB_h`              : kB-1 value for heat transfer; only used if `Rb_model = Val(:constant_kB1)`
- `Sc`                : Optional: Schmidt number of additional quantities to be calculated
- `Sc_name`           : Optional: Name of the additonal quantities, has to be of same length than 
                         `Sc_name`
- `constants=`[`bigleaf_constants`](@ref)`()`: Dictionary with entries 
  - `k` - von Karman constant 
  - `cp` - specific heat of air for constant pressure (J K-1 kg-1) 
  - `Kelvin` - conversion degree Celsius to Kelvin 
  - `g` - gravitational acceleration (m s-2) 
  - `pressure0` - reference atmospheric pressure at sea level (Pa) 
  - `Tair0` - reference air temperature (K) 
  - `Sc_CO2` - Schmidt number for CO2  
  - `Pr` - Prandtl number (if `Sc` is provided)

# Details

Aerodynamic conductance for heat (Ga_h) is calculated as:

``Ga_h = 1 / (Ra_m + Rb_h)``
 
where ``Ra_m`` is the aerodynamic resistance for momentum and ``Rb`` the (quasi-laminar) 
canopy boundary layer resistance ('excess resistance').
 
The aerodynamic resistance for momentum ``Ra_m`` is given by:

``Ra_m = u/ustar^2``

Note that this formulation accounts for changes in atmospheric stability, and does not 
require an additional stability correction function. 

An alternative method to calculate ``Ra_m`` is provided
(calculated if `wind_profile = true`):

``Ra_m = (ln((zr - d)/z0m) - psi_h) / (k ustar)``

If the roughness parameters z0m and d are unknown, they can be estimated using
[`roughness_parameters`](@ref). The argument `stab_formulation` determines the stability 
correction function used to account for the effect of atmospheric stability on Ra_m 
(Ra_m is lower for unstable and higher for stable stratification). Stratification is based 
on a stability parameter zeta (z-d/L), where z = reference height, d the zero-plane 
displacement height, and L the Monin-Obukhov length, calculated with 
[`Monin_Obukhov_length`](@ref)
The stability correction function is chosen by the argument `stab_formulation`. 
Options are `Val(:Dyer_1970)` and `Val(:Businger_1971)`.

The model used to determine the canopy boundary layer resistance for heat (``Rb_h``) is 
specified by the argument `Rb_model`. The following options are implemented:
`Val(:Thom_1972)` is an empirical formulation based on the 
friction velocity (ustar) (Thom 1972):

``Rb_h = 6.2ustar^-0.667``
  
The model by Choudhury & Monteith 1988 (`Rb_model = Val(:Choudhury_1988)`),
calculates ``Rb_h`` based on leaf width, LAI and ustar (Note that function argument `Dl`
represents leaf width (w) and not characteristic leaf dimension (Dl)
if `Rb_model` = `Val(:Choudhury_1988)`):
 
``Gb_h = LAI((0.02/\\alpha)*\\sqrt(u(zh)/w)*(1-exp(-\\alpha/2)))``
   
where ``\\alpha`` is a canopy attenuation coefficient modeled in dependence on LAI,
u(zh) is wind speed at canopy height (calculated from [`wind_profile`](@ref)),
and w is leaf width (m). See [`Gb_Choudhury`](@ref) for further details.
   
The option `Rb_model = Val(:Su_2001)` calculates ``Rb_h`` based on the physically-based 
Rb model by Su et al. 2001, a simplification of the model developed by Massman 1999:
 
``kB-1 = (k Cd fc^2) / (4Ct ustar/u(zh)) + kBs-1(1 - fc)^2``
    
where Cd is a foliage drag coefficient (defaults to 0.2), fc is fractional
vegetation cover, Bs-1 is the inverse Stanton number for bare soil surface,
and Ct is a heat transfer coefficient. See [`Gb_Su`](@ref) for 
details on the model.

The models calculate the parameter kB-1, which is related to ``Rb_h``:

``kB-1 = Rb_h * (k * ustar)``

Rb (and Gb) for water vapor and heat are assumed to be equal in this package.
Gb for other quantities x is calculated as (Hicks et al. 1987):

``Gb_x = Gb / (Sc_x / Pr)^0.67``

where Sc_x is the Schmidt number of quantity x, and Pr is the Prandtl number (0.71).
 
# Value
a dataframe with the following columns:
- `Ga_m`: Aerodynamic conductance for momentum transfer (m s-1)
- `Ra_m`: Aerodynamic resistance for momentum transfer (s m-1)
- `Ga_h`: Aerodynamic conductance for heat transfer (m s-1)
- `Ra_h`: Aerodynamic resistance for heat transfer (s m-1)
- `Gb_h`: Canopy boundary layer conductance for heat transfer (m s-1)
- `Rb_h`: Canopy boundary layer resistance for heat transfer (s m-1)
- `kB_h`: kB-1 parameter for heat transfer
- `zeta`: Stability parameter 'zeta' (missing if `wind_profile = false`)
- `psi_h`: Integrated stability correction function (missing if `wind_profile = false`)
- `Ra_CO2`: Aerodynamic resistance for CO2 transfer (s m-1)
- `Ga_CO2`: Aerodynamic conductance for CO2 transfer (m s-1)
- `Gb_CO2`: Canopy boundary layer conductance for CO2 transfer (m s-1)
- `Ga_Sc_name`: Aerodynamic conductance for `Sc_name` (m s-1). Only added if `Sc_name` and 
                  the respective `Sc` are provided.
- `Gb_Sc_name`: Boundary layer conductance for `Sc_name` (m s-1). Only added if `Sc_name` 
  and the respective `Sc` are provided.
       
# Note
The roughness length for water and heat (z0h) is not returned by this function, but 
it can be calculated from the following relationship (e.g. Verma 1989):

``kB-1 = ln(z0m/z0h)`` 

it follows:

``z0h = z0m / exp(kB-1)``

`kB-1` is an output of this function.

Input variables such as LAI, Dl, or zh can be either constants, or
vary with time (i.e. vectors of the same length as `data`).

Note that boundary layer conductance to water vapor transfer (Gb_w) is often 
assumed to equal Gb_h. This assumption is also made in this R package, for 
example in the function [`surface_conductance`](@ref).

If the roughness length for momentum (`z0m`) is not provided as input, it is estimated 
from the function `roughness_parameters` within `wind_profile` if `wind_profile = true` 
and/or `Rb_model` = `Val(:Su_2001)` or `Val(:Choudhury_1988)` The `roughness_parameters`
function estimates a single `z0m` value for the entire time period! If a varying `z0m` value 
(e.g. across seasons or years) is required, `z0m` should be provided as input argument.
     

TODO
For adding aerodynamic conductance for other species see [`add_Ga!`](@ref).
        
# References
- Verma, S., 1989: Aerodynamic resistances to transfers of heat, mass and momentum.
  In: Estimation of areal evapotranspiration, IAHS Pub, 177, 13-20.
- Verhoef, A., De Bruin, H., Van Den Hurk, B., 1997: Some practical notes on the parameter kB-1
  for sparse vegetation. Journal of Applied Meteorology, 36, 560-572.
- Hicks, B_B., Baldocchi, D_D., Meyers, T_P., Hosker, J_R., Matt, D_R., 1987:
  A preliminary multiple resistance routine for deriving dry deposition velocities
  from measured quantities. Water, Air, and Soil Pollution 36, 311-330.
- Monteith, J_L., Unsworth, M_H., 2008: Principles of environmental physics.
  Third Edition. Elsevier Academic Press, Burlington, USA. 

# See also
[`Gb_Thom`](@ref), [`Gb_Choudhury`](@ref), [`Gb_Su`](@ref) for calculations of Rb / Gb only

```@example; output = false
# df = DataFrame(Tair=25,pressure=100,wind=c(3,4,5),ustar=c(0.5,0.6,0.65),H=c(200,230,250))   
 
# # simple calculation of Ga  
# aerodynamic_conductance(df,Rb_model=Val(:Thom_1972)) 

# # calculation of Ga using a model derived from the logarithmic wind profile
# aerodynamic_conductance(df,Rb_model=Val(:Thom_1972),zr=40,zh=25,d=17.5,z0m=2,wind_profile=true) 

# # simple calculation of Ga, but a physically based canopy boundary layer model
# aerodynamic_conductance(df,Rb_model=Val(:Su_2001),zr=40,zh=25,d=17.5,Dl=0.05,N=2,fc=0.8)
``` 
"""
function aerodynamic_conductance!(df, Rb_model = Val(:Thom_1972);
  zr,zh=nothing,d = 0.7*zh ,z0m=nothing,Dl=nothing,N=2,fc=nothing,LAI=nothing,Cd=0.2,hs=0.01,
  leafwidth=nothing,
  wind_profile=false,
  stab_formulation=Val(:Dyer_1970),
  kB_h=nothing,constants=bigleaf_constants()
  )
  # add zeta, psi_h and compute z0m
  if !(stab_formulation isa Val{:no_stability_correction})
    stability_parameter!(df::AbstractDataFrame; zr,d, constants)
  else
    df[!,:zeta] .= missing
  end
  # adds columns psi_m and psi_h, (:no_stability_correction: just add 0s without using zeta)
  stability_correction!(df; stab_formulation, constants)
  if isnothing(z0m)
    z0m = roughness_parameters(Val(:wind_profile), df, zh, zr; psi_m = df.psi_m).z0m
  end
  # calculate canopy boundary layer conductance (Gb)
  Rb_model isa Val{:Thom_1972} && compute_Gb!(df, Rb_model; constants)
  Rb_model isa Val{:constant_kB1} && compute_Gb!(df, Rb_model; kB_h, constants)
  Rb_model isa Val{:Choudhury_1988} && compute_Gb!(df, Rb_model; 
    leafwidth, LAI, zh, zr, d, z0m, stab_formulation, constants)
  Rb_model isa Val{:Su_2001} && compute_Gb!(df, Rb_model; 
    Dl, fc, N, Cd, hs, z0m, zh, zr, d, LAI, stab_formulation, constants)
  # calculate aerodynamic conductance for momentum (Ga_m)
  # TODO factor out to own function
  df[!,:Ra_m] .= @. max((log((zr - d)/z0m) - df.psi_h),0) / (constants[:k]*ustar)
  df[!,:Ga_m] .= 1/df.Ra_m
  df[!,:Ra_h] .= df.Ra_m + df.Rb_h
  df[!,:Ga_h] .= 1/df.Ra_h
  df[!,:GA_CO2] .=  1/(df.Ra_m + 1/df.Gb_CO2)
  df[!,:Ra_CO2] = 1/df.Ga_CO2
  # TODO add Ga for other species:
  # Ga_x   <- 1/(Ra_m + 1/Gb_x)
  df
#   return(DataFrame(Ga_m,Ra_m,Ga_h,Ra_h,Gb_h,Rb_h,kB_h,zeta,psi_h,Ra_CO2,Gab_x))
end

"""
    add_Ga(Gb_h, Ga_m, Sc::Vararg{Pair,N}; constants)
    add_Ga!(df::AbstractDataFrame, Sc; Gb_h = df.Gb_h, Ga_m = df.Ga.m, kwargs...) 

compute additional aerodynamic conductance quantities for given Schmidt-numbers

# Arguments  
- `Gb_h`     : Boundary layer conductance for heat transfer (m s-1)
- `Ga_m`     : Aerodynamic conductance for momentum (m s-1)
- `Sc`       : several `Pair{Symbol,Number}` Output name and Schmidt number of 
               additional conductances to be calculated
- `df`       : DataFrame to add output columns               
optional
- `constants=`[`bigleaf_constants`](@ref)`()`: Dictionary with entries 
  - `Pr` - Prandtl number 

# Details
Aerodynamic conductance is calculated as

``Ga_x = 1/(1/Ga_m + 1/Gb_x)``

where Gb_x is the Boundary layer conductance for other quantities x is calculated 
based on boundary layer for heat transfer, Schmidt-Number, and  Prantl number,
as documented in [`add_Gb`](@ref).

# Value
a NameTuple or `df` with keys `Ga_x` where `x` are the keys in `Sc` and 
corresponding aerodynamic conductances (m s-1).

# Examples
```jldoctest; output=false
using DataFrames
df = DataFrame(Gb_h=[0.02, missing, 0.055], Ga_m = [0.03, 0.03, 0.03])
add_Ga!(df, :O2 => 0.84, :CH4 => 0.99)
propertynames(df)[3:4] == [:Ga_O2, :Ga_CH4]
# output
true
```
"""
function add_Ga(Gb_h::Union{Missing,Number}, Ga_m::Union{Missing,Number}, 
  Sc::Vararg{Pair,N}; kwargs...) where N
  Scn, Scv = get_names_and_values("Ga_", Sc...)
  add_Ga_(Gb_h, Ga_m, Scn, Scv; kwargs...)
end
function add_Ga_(Gb_h::Union{Missing,Number}, Ga_m::Union{Missing,Number}, 
  Scn::NTuple{N,Symbol}, Scv::NTuple{N};  constants=bigleaf_constants()) where N 
  Gbx = add_Gb_(Gb_h, Scn, Scv; constants)
  Gaxv = ntuple(i -> 1/(1/Ga_m + 1/Gbx[i]), length(Gbx))
  Gax = NamedTuple{Scn}(Gaxv)
end
function add_Ga!(df::AbstractDataFrame, Sc::Vararg{Pair,N}; Gb_h = df.Gb_h, Ga_m = df.Ga_m, kwargs...) where N
  N == 0 && return(df)
  Scn, Scv = get_names_and_values("Ga_", Sc...)
  ft() = add_Ga_.(Gb_h, Ga_m, Ref(Scn), Ref(Scv); kwargs...)
  transform!(df, [] => ft => AsTable)
end




