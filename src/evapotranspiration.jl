"""
    potential_ET(Tair, pressure, Rn, ::Val{:PriestleyTaylor}; ...)
    potential_ET(Tair, pressure, Rn, G, S, ::Val{:PriestleyTaylor}; ...)

    potential_ET(df, approach; ...)
    potential_ET!(df, approach; ...)

Potential evapotranspiration according to Priestley & Taylor 1972 or
the Penman-Monteith equation with a prescribed surface conductance.

# Arguments             
- Tair:      Air temperature (degC)
- pressure:  Atmospheric pressure (kPa)
- Rn:        Net radiation (W m-2)
- VPD:       Vapor pressure deficit (kPa)
- Ga:        Aerodynamic conductance to heat/water vapor (m s-1)
- approach:  Approach used. 
  Either `Val(:PriestleyTaylor)` (default), or `Val(:PenmanMonteith)`.
- df:      Data_frame or matrix containing all required variables; optional
optional:
- G:         Ground heat flux (W m-2). Defaults to zero.
- S:         Sum of all storage fluxes (W m-2) . Defaults to zero.
- `Esat_formula`: formula used in [`Esat_from_Tair`](@ref)
- `constants=`[`bigleaf_constants`](@ref)`()`: Dictionary with entries 
  - cp - specific heat of air for constant pressure (J K-1 kg-1) 
  - eps - ratio of the molecular weight of water vapor to dry air 
  - Pa2kPa - conversion pascal (Pa) to kilopascal (kPa) 
  - for :PenmanMonteith:
    - Rd - gas constant of dry air (J kg-1 K-1)
    - Rgas - universal gas constant (J mol-1 K-1)
    - Kelvin - conversion degree Celsius to Kelvin
additional optional arguments if G and S are not specified by postition
and in the non-mutating DataFrame variant.
- G and S: as in the positional forms, but checked for missing or containing missings
- missing_G_as_NA:  if `true`, missing optional G are treated as `NA`s, 
  otherwise set to 0. 
- missing_S_as_NA:  if `true`, missing optional S are treated as `NA`s, 
  otherwise set to 0. 
additional optional for PriestleyTaylor:
- alpha:     Priestley-Taylor coefficient
additional optional for PenmanMonteith:
- Gs_pot:    Potential/maximum surface conductance (mol m-2 s-1); 
  defaults to 0.6 mol m-2 s-1;

# Details
Potential evapotranspiration is calculated according to Priestley & Taylor, 1972
if `approach = Val(:PriestleyTaylor)` (the defau
``LE_{pot,PT} = (\\alpha * \\Delta * (Rn - G - S)) / (\\Delta + \\gamma)``

``\\alpha`` is the Priestley-Taylor coefficient, ``\\Delta`` is the slope 
of the saturation vapor pressure curve (kPa K-1), and ``\\gamma`` is the 
psychrometric constant (kPa K-1).
if `approach = Val(:PenmanMonteith)`, potential evapotranspiration is calculated according
to the Penman-Monteith equat

``LE_{pot,PM} = (\\Delta * (Rn - G - S) + \\rho * cp * VPD * Ga) / 
(\\Delta + \\gamma * (1 + Ga/Gs_{pot})``

where ``\\Delta`` is the slope of the saturation vapor pressure curve (kPa K-1),
``\\rho`` is the air density (kg m-3), 
and ``\\gamma`` is the psychrometric constant (kPa K-1).
The value of `Gs_pot` is typically a maximum value of Gs observed at the site, e.g. the 90th
percentile of Gs within the growing season.

If keyword arguments `G` or `S` in the are set to `missing`, ground heat flux and storage flux
are assumed zero respectively. Note, that these keyword argument are not available in the
mutating form, where the caller has to add tree these columns before 
(see [fill_GS_missings](@ref)).
         
Both methods are provided with several forms:
- all required inputs as positional arguments
- provideing G and S as positional arguments with handing of missing values
- providing a DataFrame with columns corresponding to required inputs


# Value
NamedTuple or DataFrame with the following entries:
- `ET_pot`: Potential evapotranspiration (kg m-2 s-1)
- `LE_pot`: Potential latent heat flux (W m-2)
For the mutating form, the original df with columns `ET_pot`, `LE_pot`, `G`, and `S` 
updated or added.


# References 
- Priestley, C_H_B., Taylor, R_J., 1972: On the assessment of surface heat flux
and evaporation using large-scale parameters. Monthly Weather Review 100, 81-92.  
- Allen, R_G., Pereira L_S., Raes D., Smith M., 1998: Crop evapotranspiration -
Guidelines for computing crop water requirements - FAO Irrigation and drainage
paper 56.
- Novick, K_A., et al. 2016: The increasing importance of atmospheric demand
for ecosystem water and carbon fluxes. Nature Climate Change 6, 1023 - 1027.
            
# See also 
[`surface_conductance`](@ref)

# Examples

```jldoctest; output=false
# Calculate potential ET of a surface that receives a net radiation of 500 Wm-2
# using Priestley-Taylor:
Tair,pressure,Rn = 30.0,100.0,500.0
ET_pot, LE_pot = potential_ET(Tair,pressure,Rn, Val(:PriestleyTaylor))    
≈(ET_pot, 0.0002035969; rtol = 1e-5)
# output
true
``` 
```@jldoctest; output=false
# Calculate potential ET for a surface with known Gs (0.5 mol m-2 s-1) and Ga (0.1 m s-1)
# using Penman-Monteith:
Tair,pressure,Rn = 30.0,100.0,500.0
VPD,Ga = 2.0, 0.1
ET_pot,LE_pot = potential_ET(Tair,pressure,Rn,VPD, Ga, Val(:PenmanMonteith); 
  Gs_pot=0.5, infoGS=false)    
# now cross-check with the inverted equation
#Ga2 = surface_conductance(Tair=20,pressure=100,VPD=2,Ga=0.1,Rn=400,LE=LE_pot_PM)
#Ga2 ≈ GA
true
# output
true
``` 
"""
function potential_ET(Tair, pressure, Rn, approach::Val{:PriestleyTaylor};
  G=zero(Tair),S=zero(Tair), kwargs...)
  #
  potential_ET(Tair, pressure, Rn, G, S, approach; kwargs...)
end,
function potential_ET(Tair, pressure, Rn, G, S, ::Val{:PriestleyTaylor};
  alpha=1.26,
  Esat_formula=Val(:Sonntag_1990),
  constants=bigleaf_constants())
  #
  gamma  = psychrometric_constant(Tair,pressure;constants)
  Delta  = Esat_from_Tair_deriv(Tair; formula = Esat_formula,constants)
  LE_pot = (alpha * Delta * (Rn - G - S)) / (Delta + gamma)
  ET_pot = LE_to_ET(LE_pot,Tair)
  (ET_pot = ET_pot, LE_pot = LE_pot)
end,
function potential_ET(Tair, pressure, Rn, VPD, Ga, approach::Val{:PenmanMonteith};
  G=zero(Tair),S=zero(Tair), kwargs...)
  #
  potential_ET(Tair, pressure, Rn, VPD, Ga, G, S, approach; kwargs...)
end,
function potential_ET(Tair, pressure, Rn, VPD, Ga, G, S, ::Val{:PenmanMonteith};
  Gs_pot=0.6,
  Esat_formula=Val(:Sonntag_1990),
  constants=bigleaf_constants())
  #
  gamma  = psychrometric_constant(Tair,pressure;constants)
  Delta  = Esat_from_Tair_deriv(Tair; formula = Esat_formula,constants)
  Gs_pot = mol_to_ms(Gs_pot,Tair,pressure;constants)
  rho    = air_density(Tair,pressure;constants)
  LE_pot = (Delta * (Rn - G - S) + rho * constants[:cp] * VPD * Ga) / 
    (Delta + gamma * (1 + Ga / Gs_pot))
  ET_pot = LE_to_ET(LE_pot,Tair)
  (ET_pot = ET_pot, LE_pot = LE_pot)
end,
function potential_ET(df, approach::Val{:PriestleyTaylor}; 
  G=missing,S=missing, infoGS=true, kwargs...) 
  #
  dfGS = get_df_GS(df, G,S; infoGS) 
  f(args...) = potential_ET(args..., approach; kwargs...)
  select(hcat(select(df,:Tair, :pressure, :Rn), dfGS; copycols = false),
    All() => ByRow(f) => AsTable
  )
end,
function potential_ET(df, approach::Val{:PenmanMonteith}; 
  G=missing,S=missing, missing_G_as_NA=false, missing_S_as_NA=false, infoGS=true,
  kwargs...)
  #
  dfGS = get_df_GS(df, G,S; infoGS) 
  f(args...) = potential_ET(args..., approach; kwargs...)
  select(hcat(select(df,:Tair, :pressure, :Rn, :VPD, :Ga), dfGS; copycols = false),
    All() => ByRow(f) => AsTable
  )
end,
function potential_ET!(df, approach::Val{:PriestleyTaylor}; 
  G=missing,S=missing, infoGS=true, kwargs...) 
  dfGS = get_df_GS(df, G,S; infoGS) 
  # temporarily add G and S to the DataFrame to mutate
  df._tmp_G .= dfGS.G
  df._tmp_S .= dfGS.S
  f(args...) = potential_ET(args..., approach; kwargs...)
  transform!(df,
    [:Tair, :pressure, :Rn, :_tmp_G, :_tmp_S] => ByRow(f) => AsTable
   )
   select!(df, Not([:_tmp_G, :_tmp_S]))
end,
function potential_ET!(df, approach::Val{:PenmanMonteith}; 
  G=missing,S=missing, infoGS=true, kwargs...) 
  dfGS = get_df_GS(df, G,S; infoGS) 
  df._tmp_G .= dfGS.G
  df._tmp_S .= dfGS.S
  f(args...) = potential_ET(args..., approach; kwargs...)
  transform!(df, 
    [:Tair, :pressure, :Rn, :VPD, :Ga, :_tmp_G, :_tmp_S] => ByRow(f) => AsTable
  )
  select!(df, Not([:_tmp_G, :_tmp_S]))
end

function get_df_GS(df, G, S; infoGS=true)
  nout = nrow(df)
  G_ = if ismissing(G)
    infoGS && @info("Ground heat flux G is not provided and set to 0.")
    Zeros(nout)
  else 
    G
  end
  S_ = if ismissing(S)
    infoGS && @info("Storage heat flux S is not provided and set to 0.")
    Zeros(nout)
  else
    S
  end
  DataFrame(G  = G_, S = S_)
end

# function fill_GS_missings(df::DataFrame, G,S, missing_G_as_NA, missing_S_as_NA; infoGS=true)
#   nout = nrow(df)
#   G_ = ifelse(
#     ismissing(G),
#     (infoGS && @info("Ground heat flux G is not provided and set to 0."); Zeros(nout)),
#     @. ifelse(missing_G_as_NA, G, coalesce(G, 0.0))
#     )
#   S_ = ifelse(
#   ismissing(S),
#   (infoGS && @info("Storage heat flux S is not provided and set to 0."); Zeros(nout)),
#   @. ifelse(missing_S_as_NA, G, coalesce(S, 0.0))
#   )
#   DataFrame(G  = G_, S = S_)
# end

function fill_vec(G; is_replace_missing = true, fillvalue = zero(G))
  @. ifelse(is_replace_missing, coalesce(G, fillvalue), G)
end



"""
    TODO

Evapotranspiration (ET) split up into imposed ET and equilibrium ET.

- Tair:      Air temperature (deg C)
- pressure:  Atmospheric pressure (kPa)
- VPD:       Air vapor pressure deficit (kPa)
- Gs:        surface conductance to water vapor (m s-1)
- Rn:        Net radiation (W m-2)
- G:         Ground heat flux (W m-2); optional
- S:         Sum of all storage fluxes (W m-2); optional
- missing_G_as_NA:  if `TRUE`, missing G are treated as `NA`s, otherwise set 0. 
- missing_S_as_NA:  if `TRUE`, missing S are treated as `NA`s, otherwise set 0.
- Esat_formula:  Optional: formula to be used for the calculation of esat and  slope of esat.
                     One of `"Sonntag_1990"` (Default), `"Alduchov_1996"`, `"Allen_1998"`.
                     See [`Esat_slope`](@ref). 
- `constants=`[`bigleaf_constants`](@ref)`()`: Dictionary with entries 

- constants 
  - cp - specific heat of air for constant pressure (J K-1 kg-1) 
  - eps - ratio of the molecular weight of water vapor to dry  (-) 
  - Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
                 
# Details
Total evapotranspiration can be written in the form (Jarvis & McNaughton 6):

``ET = \\Omega ET_{eq} + (1 - \\Omega)ET_{imp}``

where ``\\Omega`` is the decoupling coefficient as calculated from
[`decoupling`](@ref). `ET_eq` is the equilibrium evapotranspiration e,
the ET rate that would occur under uncoupled conditions, where tbudget
is dominated by radiation (when Ga -> 0):

  ``ET_eq = (\\Delta * (Rn - G - S) * \\lambda) / (\\Del\\gamma)         

where ``\\Delta`` is the slope of the saturation vapor pressur(kPa K-1),
``\\lambda`` is the latent heat of vaporization (J kg-1), and \\gamma``
is the psychrometric constant (kPa K-1).
`ET_imp` is the imposed evapotranspiration rate, the ET rate
that would occur under fully coupled conditions (when Ga -> inf):

  ``ET_imp = (\\rho * cp * VPD * Gs * \\lambda) / \\gamma``

where ``\\rho`` is the air density (kg m-3).

# Note
Surface conductance (Gs) can be calculated with [`surface_conductance`](@ref)      
Aerodynamic conductance (Ga) can be calculated using erodynamic_conductance`](@ref).
      
# Value
A DataFrame with the following columns:
- ET_eq: Equilibrium ET (kg m-2 s-1)
- ET_imp: Imposed ET (kg m-2 s-1)
- LE_eq: Equilibrium LE (W m-2)
- LE_imp: Imposed LE (W m-2)      

#References
- Jarvis, P_G., McNaughton, K_G., 1986: Stomatal control of transpiration:
scaling up from leaf to region. Advances in Ecological Rese1-49.
- Monteith, J_L., Unsworth, M_H., 2008: Principles of ironmPhysics.
3rd edition. Academic Press, London. 
            
#See also
[`decoupling`](@ref)            
            
```@example; output = false
df = DataFrame(Tair=20,pressure=100,VPD=seq(0.5,4,0.5),
                 Gs_ms=seq(0.01,0.002,length_out=8),Rn=seq(50,400,50)          
equilibrium_imposed_ET(df)            
``` 
"""
function equilibrium_imposed_ET(Tair,pressure,VPD,Gs, Rn,
  G=NULL,S=NULL,missing_G_as_NA=false,missing_S_as_NA=false,
  Esat_formula=Val(:Sonntag_1990),
  constants=bigleaf_constants())
  # 
  #check_input(data,list(Tair,pressure,VPD,Rn,Gs,G,S))
  G,S = fill_GS_missings(G,S,missing_G_as_NA, missing_S_as_NA)
  rho    = air_density(Tair,pressure,constants)
  gamma  = psychrometric_constant(Tair,pressure,constants)
  Delta  = Esat_from_Tair_deriv(Tair,Esat_formula,constants)
  LE_eq  = (Delta * (Rn - G - S)) / (gamma + Delta)
  LE_imp = (rho * constants[:cp] * Gs * VPD) / gamma
  #
  ET_imp = LE_to_ET(LE_imp,Tair)
  ET_eq  = LE_to_ET(LE_eq,Tair)
  (ET_pot = ET_pot, LE_pot = LE_pot, LE_imp = LE_imp)
end


"""
    fill_GS_missings!(df::DataFrame, G=df.G, S=df.S, 
      missing_G_as_NA=false, missing_S_as_NA=false; infoGS=true)

Fill missing values in ground heat flux and storage flux.

# Arguments
- df: DataFrame where columns G and S should be updated
- G: Ground heat flux
- S: Storage heat flux
- missing_G_as_NA: set to true to not fill NAs in G to propagate to computations
- missing_S_as_NA: set to true to not fill NAs in S to propagate to computations
- infoGS: set to false to avoid log info message 

```@example
true
# using DataFrames
# df = DataFrame(Tair = 20.0:1.0:30.0,pressure = 100.0, Rn = 500.0)
# fill_GS_missings!(df, missing, missing)
# all( df.S .== df.G .== 0.0)
```
"""
function fill_GS_missings!(df::DataFrame, G=df.G, S=df.S, 
  missing_G_as_NA=false, missing_S_as_NA=false; infoGS=true
  )
  nout = nrow(df)
  df.G .= ifelse(
      ismissing(G),
      (infoGS && @info("Ground heat flux G is not provided and set to 0."); Zeros(nout)),
      @. ifelse(missing_G_as_NA, G, coalesce(G, 0.0))
      )
  df.S .= ifelse(
      ismissing(S),
      (infoGS && @info("Storage heat flux S is not provided and set to 0."); Zeros(nout)),
      @. ifelse(missing_S_as_NA, G, coalesce(S, 0.0))
    )
  df
end

"""
    TODO; implement surface_conductance.

This stub is there to satisfy links im Help-pages.
"""
function surface_conductance()
  error("not yet implemented.")
end



