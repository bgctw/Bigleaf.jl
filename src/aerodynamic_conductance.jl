"""
Aerodynamic Conductance

Bulk aerodynamic conductance, including options for the boundary layer conductance
formulation and stability correction functions.

# Arguments
- `data`              : Data_frame or matrix containing all required variables
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
- `zeta`: Stability parameter 'zeta' (NA if `wind_profile = false`)
- `psi_h`: Integrated stability correction function (NA if `wind_profile = false`)
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
function aerodynamic_conductance(df;
  zr,zh,d,z0m=nothing,Dl,N=2,fc=nothing,LAI,Cd=0.2,hs=0.01,wind_profile=false,
  stab_correction=true,stab_formulation=Val(:Dyer_1970),
  Rb_model=c(Val(:Thom_1972),Val(:Choudhury_1988),Val(:Su_2001),Val(:constant_kB1)),
  kB_h=nothing,Sc=nothing,Sc_name=nothing,constants=bigleaf_constants()
  )

  ## calculate canopy boundary layer conductance (Gb)
  if Rb_model in SA[Val(:Thom_1972),Val(:Choudhury_1988),Val(:Su_2001)]
    if Rb_model == Val(:Thom_1972)
      Gb_mod = Gb_Thom(ustar=ustar,Sc=Sc,Sc_name=Sc_name,constants=constants)
    elseif Rb_model == Val(:Choudhury_1988)
      Gb_mod = Gb_Choudhury(df,Tair=Tair,pressure=pressure,wind=wind,ustar=ustar,
                             H=H,leafwidth=Dl,LAI=LAI,zh=zh,zr=zr,d=d,z0m=z0m,
                             stab_formulation=stab_formulation,Sc=Sc,Sc_name=Sc_name,
                             constants=constants)
    elseif Rb_model == Val(:Su_2001)
      Gb_mod = Gb_Su(data=df,Tair=Tair,pressure=pressure,ustar=ustar,wind=wind,
                      H=H,zh=zh,zr=zr,d=d,z0m=z0m,Dl=Dl,N=N,fc=fc,LAI=LAI,Cd=Cd,hs=hs,
                      stab_formulation=stab_formulation,Sc=Sc,Sc_name=Sc_name,
                      constants=constants)  
    end
    kB_h = Gb_mod.kB_h
    Rb_h = Gb_mod.Rb_h
    Gb_h = Gb_mod.Gb_h
    # TODO
    # Gb_x = DataFrame(Gb_mod[,grep(names(Gb_mod),pattern="Gb_")[-1]])
    # colnames(Gb_x) = grep(colnames(Gb_mod),pattern="Gb_",value=true)[-1]
  elseif Rb_model == Val(:constant_kB1)
    isnothing(kB_h) && error(
      "value of kB-1 has to be specified if Rb_model is set to 'constant_kB-1'!")
    Rb_h = kB_h/(constants[:k] * ustar)
    Gb_h = 1/Rb_h
    if (!isnothing(Sc) || !isnothing(Sc_name))
      length(Sc) != length(Sc_name) && error(
        "arguments 'Sc' and 'Sc_name' must have the same length")
      !is_numeric(Sc) && error("argument 'Sc' must be numeric")
      Sc   = SA[constants[:Sc_CO2], Sc]
      Gb_x = DataFrame(lapply(Sc,function(x) Gb_h / (x/constants[:Pr])^0.67))
      colnames(Gb_x) = paste0("Gb_",c("CO2",Sc_name))
    end
  end 
  
#   ## calculate aerodynamic conductance for momentum (Ga_m)
#   if (wind_profile)
#     if (isnothing(z0m) & Rb_model in c(Val(:constant_kB1),Val(:Thom_1972)))
#       stop("z0m must be provided if wind_profile=true!")
# elseif (isnothing(z0m) & Rb_model in c(Val(:Choudhury_1988),Val(:Su_2001)))
#       # z0m estimated as in Choudhury_1988 or Su_2001
#       z0m = roughness_parameters(method=Val(:wind_profile),zh=zh,zr=zr,d=d,data=df,
#                                   Tair=Tair,pressure=pressure,wind=wind,ustar=ustar,H=H,
#                                   stab_roughness=true,stab_formulation=stab_formulation,
#                                   constants=constants)[,"z0m"]
# end
    
#     if (stab_correction)
      
#       zeta  =  stability_parameter(data=df,Tair=Tair,pressure=pressure,ustar=ustar,
#                                     H=H,zr=zr,d=d,constants=constants)
      
#       if (stab_formulation in c(Val(:Dyer_1970),Val(:Businger_1971)))
        
#         psi_h = stability_correction(zeta,formulation=stab_formulation)[,"psi_h"]
        
# else 
#         stop("'stab_formulation' has to be one of 'Dyer_1970' or 'Businger_1971'.
#              Choose 'stab_correction = false' if no stability correction should be applied.")
# end
      
#       Ra_m  = pmax((log((zr - d)/z0m) - psi_h),0) / (constants[:k]*ustar)
      
# else 
        
#         Ra_m  = pmax((log((zr - d)/z0m)),0) / (constants[:k]*ustar)
#         zeta = psi_h = rep(NA_integer_,length=length(Ra_m))
        
# end
    
# else 
    
#     if ((!missing(zr) | !missing(d) | !missing(z0m)) & Rb_model in c(Val(:constant_kB1),Val(:Thom_1972)))
#       @warn"Provided roughness length parameters (zr,d,z0m) are not used if 'wind_profile = false' (the default). Ra_m is calculated as wind / ustar^2")
# end
    
#     Ra_m = wind / ustar^2
#     zeta = psi_h = rep(NA_integer_,length=length(Ra_m))
    
# end
  
#   Ga_m   = 1/Ra_m
#   Ra_h   = Ra_m + Rb_h
#   Ga_h   = 1/Ra_h
#   Ga_x   = 1/(Ra_m + 1/Gb_x)
#   Ra_CO2 = 1/Ga_x[,1]
#   colnames(Ga_x) = paste0("Ga_",c("CO2",Sc_name))
  
#   Gab_x = cbind(Ga_x,Gb_x)
#   Gab_x = Gab_x[rep(c(1,ncol(Gab_x)-(ncol(Gab_x)/2-1)),ncol(Gab_x)/2) + sort(rep(0:(ncol(Gab_x)/2-1),2))] # reorder columns

#   return(DataFrame(Ga_m,Ra_m,Ga_h,Ra_h,Gb_h,Rb_h,kB_h,zeta,psi_h,Ra_CO2,Gab_x))
end