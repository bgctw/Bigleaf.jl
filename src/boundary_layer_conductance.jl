"""
Boundary Layer Conductance according to Thom 1972

An empirical formulation for the canopy boundary layer conductance 
for heat transfer based on a simple ustar (friction velocity) dependency.

# Arguments  
For `Val(:Thom_1972)`
- `ustar`     : Friction velocity (m s-1)
optional
- `constants=`[`bigleaf_constants`](@ref)`()`: Dictionary with entries 
 
# Details
The empirical equation for Rb suggested by Thom 1972 is:
 
``Rb = 6.2 ustar^{-0.67}``
 
Gb (=1/Rb) for water vapor and heat are assumed to be equal in this package.
Gb for other quantities x is calculated as (Hicks et al. 1987):
 
``Gb_x = Gb / (Sc_x / Pr)^{0.67}``
 
where `Sc_x` is the Schmidt number of quantity x, and Pr is the Prandtl number (0.71).
 
# Value
a NameTuple or DataFrame with the following columns:
- `Gb_h`: Boundary layer conductance for heat transfer (m s-1)
- `Rb_h`: Boundary layer resistance for heat transfer (s m-1)
- `kB_h`: kB-1 parameter for heat transfer
- `Gb_CO2`: Boundary layer conductance for CO2 (m s-1). 
- ...: Boundary layer conductance for keys provided with argument `Sc` (m s-1). 
 
# References
- Thom, A., 1972: Momentum, mass and heat exchange of vegetation.
  Quarterly Journal of the Royal Meteorological Society 98, 124-134.
- Hicks, B_B., Baldocchi, D_D., Meyers, T_P., Hosker, J_R., Matt, D_R., 1987:
  A preliminary multiple resistance routine for deriving dry deposition velocities
  from measured quantities. Water, Air, and Soil Pollution 36, 311-330.

# See also
[`Gb_Choudhury`](@ref), [`Gb_Su`](@ref), [`aerodynamic_conductance`](@ref)

```@example; output = false
# Gb_Thom(seq(0.1,1.4,0.1))

# ## calculate Gb for SO2 as well
# Gb_Thom(seq(0.1,1.4,0.1),Sc=1.25,Sc_name="SO2")
``` 
"""
function compute_Gb!(df::AbstractDataFrame, approach::Val{:Thom_1972}; kwargs...)
  compute_Gb_!(df, approach, :ustar; kwargs...) # inputcols
end
function compute_Gb!(df::AbstractDataFrame, approach::Val{:Choudhury_1988}; kwargs...)
  Gb_Choudhury!(df; kwargs...)
end
function compute_Gb!(df::AbstractDataFrame, approach::Val{:Su_2001}; kwargs...)
  Gb_Su!(df; kwargs...)
end

function compute_Gb_!(df::AbstractDataFrame, approach, inputcols; 
  Sc::AbstractVector=SA[], kwargs...
  )
  fGb(args...) = compute_Gb(approach, args...; kwargs...)
  transform!(df, :ustar => ByRow(fGb) => SA[:Rb_h, :Gb_h, :kB_h, :Gb_CO2])
end
compute_Gb(::Val{:Thom_1972}, args...; kwargs...) = Gb_Thom(args...; kwargs...)

"""
Boundary Layer Conductance according to Thom 1972

An empirical formulation for the canopy boundary layer conductance 
for heat transfer based on a simple ustar (friction velocity) dependency.

# Arguments  
- `ustar`     : Friction velocity (m s-1)
optional
- `Sc`        : `Pair{Symbol,Number}` Output name and Schmidt number of 
                additional quantities to be calculated
- `constants=`[`bigleaf_constants`](@ref)`()`: Dictionary with entries 
  - `k` - von-Karman constant 
  - `Sc_CO2` - Schmidt number for CO2  
  - `Pr` - Prandtl number (if `Sc` is provided)
 
# Details
The empirical equation for Rb suggested by Thom 1972 is:
 
``Rb = 6.2 ustar^{-0.67}``
 
Gb (=1/Rb) for water vapor and heat are assumed to be equal in this package.
Gb for other quantities x is calculated as (Hicks et al. 1987):
 
``Gb_x = Gb / (Sc_x / Pr)^{0.67}``
 
where `Sc_x` is the Schmidt number of quantity x, and Pr is the Prandtl number (0.71).
 
# Value
a NameTuple or DataFrame with the following columns:
- `Gb_h`: Boundary layer conductance for heat transfer (m s-1)
- `Rb_h`: Boundary layer resistance for heat transfer (s m-1)
- `kB_h`: kB-1 parameter for heat transfer
- `Gb_CO2`: Boundary layer conductance for CO2 (m s-1). 
- ...: Boundary layer conductance for keys provided with argument `Sc` (m s-1). 
 
# References
- Thom, A., 1972: Momentum, mass and heat exchange of vegetation.
  Quarterly Journal of the Royal Meteorological Society 98, 124-134.
- Hicks, B_B., Baldocchi, D_D., Meyers, T_P., Hosker, J_R., Matt, D_R., 1987:
  A preliminary multiple resistance routine for deriving dry deposition velocities
  from measured quantities. Water, Air, and Soil Pollution 36, 311-330.

# See also
[`Gb_Choudhury`](@ref), [`Gb_Su`](@ref), [`aerodynamic_conductance`](@ref)

```@example; output = false
# Gb_Thom(seq(0.1,1.4,0.1))

# ## calculate Gb for SO2 as well
# Gb_Thom(seq(0.1,1.4,0.1),Sc=1.25,Sc_name="SO2")
``` 
"""
function Gb_Thom(ustar::Union{Missing,Number}; constants=bigleaf_constants())
  Rb_h = 6.2*ustar^-0.667
  Gb_h = 1/Rb_h
  kB_h = Rb_h*constants[:k]*ustar
  Gb_CO2 = Gb_h / (constants[:Sc_CO2]/constants[:Pr])^0.67
  (;Rb_h, Gb_h, kB_h, Gb_CO2)
end
# function Gb_Thom(ustar::Missing; constants=bigleaf_constants()) 
#   (Rb_h = missing, Gb_h = missing, kB_h = missing, Gb_CO2 = missing)
# end





"""
    add_Gb(Gb_h::Union{Missing,Number}, Sc::Vararg{Pair,N}; constants=bigleaf_constants())
    add_Gb(df::AbstractDataFrame, Sc::Vararg{Pair,N}; Gb_h = df.Gb_h, kwargs...) 

compute additional boundary layer conductance quantities for given Schmidt-numbers

# Arguments  
- `Gb_h`     : Boundary layer conductance for heat transfer (m s-1)
- `Sc`       : several `Pair{Symbol,Number}` Output name and Schmidt number of 
                additional conductances to be calculated
- `df`       : DataFrame to add output columns               
optional
- `constants=`[`bigleaf_constants`](@ref)`()`: Dictionary with entries 
  - `Pr` - Prandtl number 

# Value
a NameTuple or `df` with keys as in `Sc` and corresponding conductances (m s-1)
"""
function add_Gb(Gb_h::Union{Missing,Number}, Sc::Vararg{Pair,N}; 
  constants=bigleaf_constants()) where N
   Scn = ntuple(i -> Sc[i].first, N)
   Scv = ntuple(i -> Sc[i].second, N)
   Gbxv = @. Gb_h / (Scv/constants[:Pr])^0.67
   Gbx = NamedTuple{Scn}(Gbxv)
end
function add_Gb!(df::AbstractDataFrame, Sc::Vararg{Pair,N}; Gb_h = df.Gb_h, kwargs...) where N
  N == 0 && return(df)
  ft(Gb_h) = add_Gb(Gb_h, Sc...; kwargs...)
  transform!(df, :Gb_h => ByRow(ft) => AsTable)
end



  
"""
Boundary Layer Conductance according to Choudhury & Monteith 1988

A formulation for the canopy boundary layer conductance 
for heat transfer according to Choudhury & Monteith 1988.

# Arguments  
- `df`             : DataFrame or matrix containing all required variables
- `Tair`             : Air temperature (degC)
- `pressure`         : Atmospheric pressure (kPa)
- `wind`             : Wind speed at sensor height (m s-1)
- `ustar`            : Friction velocity (m s-1)
- `H`                : Sensible heat flux (W m-2)
- `leafwidth`        : Leaf width (m)
- `LAI`              : One-sided leaf area index
- `zh`               : Canopy height (m)
- `zr`               : Instrument (reference) height (m)
- `d`                : Zero-plane displacement height (-), can be calculated using `roughness_parameters`
- `z0m`              : Roughness length for momentum (m). If not provided, calculated from `roughness_parameters` 
                         within `wind_profile`
- `stab_formulation` : Stability correction function used (If `stab_correction = true`).
                         Either `Val(:Dyer_1970)` or `Val(:Businger_1971)`.
- `Sc`               : Optional: Schmidt number of additional quantities to be calculated
- `Sc_name`          : Optional: Name of the additonal quantities, has to be of same length than 
                        `Sc_name`
- `constants=`[`bigleaf_constants`](@ref)`()`: Dictionary with entries 
  _ k - von-Karman constant 
  _ Sc_CO2 - Schmidt number for CO2  
  _ Pr - Prandtl number (if `Sc` is provided)
                        
# Value
 A data frame with the following columns:
 - Gb_h: Boundary layer conductance for heat transfer (m s-1)
 - Rb_h: Boundary layer resistance for heat transfer (s m-1)
 - kB_h: kB-1 parameter for heat transfer
 
# Details
Boundary layer conductance according to Choudhury & Monteith 1988 is
given by:

``Gb_h = LAI((2a/\\alpha)*\\sqrt(u(zh)/w)*(1-exp(-\\alpha/2)))``

where ``\\alpha`` is modeled as an empirical relation to LAI (McNaughton & van den Hurk 1995):

``\\alpha = 4.39 - 3.97*exp(-0.258*LAI)``

``w`` is leafwidth and ``u(zh)`` is the wind speed at the canopy surface.

It can be approximated from measured wind speed at sensor height zr and a wind extinction 
coefficient ``\\alpha``: ``u(zh) = u(zr) / (exp(\\alpha(zr/zh -1)))``.
However, here it is computed from [`wind_profile`](@ref)

Gb (=1/Rb) for water vapor and heat are assumed to be equal in this package.
Gb for other quantities x is calculated as (Hicks et al. 1987):
``Gb_x = Gb / (Sc_x / Pr)^0.67``
where Sc_x is the Schmidt number of quantity x, and Pr is the Prandtl number (0.71).
         
# Note
If the roughness length for momentum (`z0m`) is not provided as input, it is estimated 
in the DataFrame variant from the function [`roughness_parameters`](@ref) within 
[`wind_profile`(@ref)]. 
This function
estimates a single `z0m` value for the entire time period! If a varying `z0m` value 
(e.g. across seasons or years) is required, `z0m` should be provided as input argument.
         
# References
- Choudhury, B. J., Monteith J_L., 1988: A four-layer model for the heat
  budget of homogeneous land surfaces. Q. J. R. Meteorol. Soc. 114, 373-398.
- McNaughton, K. G., Van den Hurk, B_J_J_M., 1995: A 'Lagrangian' revision of
  the resistors in the two-layer model for calculating the energy budget of a
  plant canopy. Boundary-Layer Meteorology 74, 261-288.
- Hicks, B_B., Baldocchi, D_D., Meyers, T_P., Hosker, J_R., Matt, D_R., 1987:
  A preliminary multiple resistance routine for deriving dry deposition velocities
  from measured quantities. Water, Air, and Soil Pollution 36, 311-330.
            
# See also
[`Gb_Thom`](@ref), [`Gb_Su`](@ref), [`aerodynamic_conductance`](@ref)
   
```@example; output = false
# ## bulk canopy boundary layer resistance for a closed canopy (LAI=5) 
# ## with large leaves (leafwdith=0.1)            
# df = DataFrame(Tair=25,pressure=100,wind=c(3,4,5),ustar=c(0.5,0.6,0.65),H=c(200,230,250))    
# Gb_Choudhury(data=df,leafwidth=0.1,LAI=5,zh=25,d=17.5,zr=40)

# ## same conditions, but smaller leaves (leafwidth=0.01)
# Gb_Choudhury(data=df,leafwidth=0.01,LAI=5,zh=25,d=17.5,zr=40) 
``` 
"""
function Gb_Choudhury(ustar; leafwidth, LAI, wind_zh, constants=bigleaf_constants())
  alpha   = 4.39 - 3.97*exp(-0.258*LAI)
  wind_zh = max(0.01, wind_zh) ## avoid zero windspeed
  Gb_h = LAI*((0.02/alpha)*sqrt(wind_zh/leafwidth)*(1-exp(-alpha/2)))
  # TODO facture out derving 3 others form Gb_h?
  Rb_h = 1/Gb_h
  kB_h = Rb_h*constants[:k]*ustar
  Gb_CO2 = Gb_h / (constants[:Sc_CO2]/constants[:Pr])^0.67
  (;Rb_h, Gb_h, kB_h, Gb_CO2)
end

function Gb_Choudhury!(df::AbstractDataFrame; leafwidth, LAI, zh, zr, d=0.7*zh,
  z0m = nothing,
  stab_formulation=Val(:Dyer_1970),
  constants=bigleaf_constants()
  )
  wind_zh = wind_profile(df, zh, d, z0m; zh, zr, stab_formulation, constants)
  # Broadcasting does not work over keyword arguments, need to pass as positional
  fwind(ustar, wind_zh; kwargs...) = Gb_Choudhury(ustar;wind_zh, kwargs...)
  ft(ustar) = fwind.(ustar, wind_zh; leafwidth, LAI, constants)
  transform!(df, :ustar => ft => AsTable)
end

"""
Boundary Layer Conductance according to Su et al. 2001

A physically based formulation for the canopy boundary layer conductance
to heat transfer according to Su et al. 2001. 

# Arguments
- data      DataFrame or matrix containing all required variables
- Tair      Air temperature (degC)
- pressure  Atmospheric pressure (kPa)
- ustar     Friction velocity (m s-1)
- wind      Wind speed (m s-1)
- H         Sensible heat flux (W m-2)
- zh        Canopy height (m)
- zr        Reference height (m)
- d         Zero-plane displacement height (-), can be calculated using `roughness_parameters`
- z0m       Roughness length for momentum (m). If not provided, calculated from `roughness_parameters` 
                 within `wind_profile`
- Dl        Leaf characteristic dimension (m)
- fc        Fractional vegetation cover [0-1] (if not provided, calculated from LAI)
- LAI       One-sided leaf area index (-)
- N         Number of leaf sides participating in heat exchange (defaults to 2)
- Cd        Foliage drag coefficient (-)
- hs        Roughness height of the soil (m)
- stab_formulation Stability correction function used (If `stab_correction = true`).
                        Either `Val(:Dyer_1970)` or `Val(:Businger_1971)`.
- Sc        Optional: Schmidt number of additional quantities to be calculated
- Sc_name   Optional: Name of the additional quantities, has to be of same length than 
                 `Sc_name`
- constants Kelvin - conversion degree Celsius to Kelvin 
                 pressure0 - reference atmospheric pressure at sea level (Pa) 
                 Tair0 - reference air temperature (K) 
                 Sc_CO2 - Schmidt number for CO2  
                 Pr - Prandtl number (if `Sc` is provided)

# Value
 A DataFrame with the following columns:
 - Gb_h: Boundary layer conductance for heat transfer (m s-1)
 - Rb_h: Boundary layer resistance for heat transfer (s m-1)
 - kB_h: kB-1 parameter for heat transfer
    
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

# Note
If the roughness length for momentum (`z0m`) is not provided as input, it is estimated 
from the function `roughness_parameters` within `wind_profile`. This function
estimates a single `z0m` value for the entire time period! If a varying `z0m` value 
(e.g. across seasons or years) is required, `z0m` should be provided as input argument.


# References
- Su, Z., Schmugge, T., Kustas, W. & Massman, W., 2001: An evaluation of
  two models for estimation of the roughness height for heat transfer between
  the land surface and the atmosphere. Journal of Applied Meteorology 40, 1933-1951.
-  Massman, W., 1999: A model study of kB H- 1 for vegetated surfaces using
  localized near-field' Lagrangian theory. Journal of Hydrology 223, 27-43.
-  Hicks, B_B., Baldocchi, D_D., Meyers, T_P., Hosker, J_R., Matt, D_R., 1987:
   A preliminary multiple resistance routine for deriving dry deposition velocities
  from measured quantities. Water, Air, and Soil Pollution 36, 311-330.

# See also
[`Gb_Thom`](@ref), [`Gb_Choudhury`](@ref), [`aerodynamic_conductance`](@ref)

```@example; output = false
# # Canopy boundary layer resistance (and kB-1 parameter) for a set of meteorological conditions,
# # a leaf characteristic dimension of 1cm, and an LAI of 5
# df = DataFrame(Tair=25,pressure=100,wind=c(3,4,5),ustar=c(0.5,0.6,0.65),H=c(200,230,250)) 
# Gb_Su(data=df,zh=25,zr=40,d=17.5,Dl=0.01,LAI=5)

# # the same meteorological conditions, but larger leaves
# Gb_Su(data=df,zh=25,zr=40,d=17.5,Dl=0.1,LAI=5)

# # same conditions, large leaves, and sparse canopy cover (LAI = 1.5)
# Gb_Su(data=df,zh=25,zr=40,d=17.5,Dl=0.1,LAI=1.5)
``` 
"""
function Gb_Su(Tair,pressure,ustar; zh, wind_zh, Dl, fc, N=2, Cd=0.2, hs=0.01,
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
  Gb_CO2 = Gb_h / (constants[:Sc_CO2]/constants[:Pr])^0.67
  (;Rb_h, Gb_h, kB_h, Gb_CO2)
end

function Gb_Su!(df::AbstractDataFrame; Dl, fc=nothing, N=2, Cd=0.2, hs=0.01, 
  z0m = nothing, zh, zr, d = 0.7*zh, LAI,
  stab_formulation=Val(:Dyer_1970), constants=bigleaf_constants()
  )
  wind_zh = wind_profile(df, zh, d, z0m; zh, zr, stab_formulation, constants)
  if isnothing(fc)
    isnothing(LAI) && error("one of 'fc' or 'LAI' must be provided")
    fc = (1-exp(-LAI/2)) 
  end
  inputcols = SA[:Tair,:pressure,:ustar]
  # Broadcasting does not work over keyword arguments, need to pass as positional
  fwind(wind_zh, args...; kwargs...) = Gb_Su(args...; wind_zh, kwargs...)
  ft(args...) = fwind.(wind_zh, args...; zh, Dl, fc, N=2, Cd=0.2, hs=0.01, constants)
  transform!(df, inputcols => ft => AsTable)
end
