###########################
#### Evapotranspiration ###
###########################

#' Potential Evapotranspiration
#' 
#' Potential evapotranspiration according to Priestley & Taylor 1972 or
#'              the Penman-Monteith equation with a prescribed surface conductance.
#' 
#' - data      Data_frame or matrix containing all required variables; optional
#' - Tair      Air temperature (degC)
#' - pressure  Atmospheric pressure (kPa)
#' - Rn        Net radiation (W m-2)
#' - G         Ground heat flux (W m-2); optional
#' - S         Sum of all storage fluxes (W m-2); optional
#' - VPD       Vapor pressure deficit (kPa); only used if `approach = Val(:PenmanMonteith)`.
#' - Ga        Aerodynamic conductance to heat/water vapor (m s-1); only used if `approach = Val(:PenmanMonteith)`.
#' - approach  Approach used. Either `Val(:PriestleyTaylor)` (default), or `Val(:PenmanMonteith)`.
#' - alpha     Priestley-Taylor coefficient; only used if `approach = Val(:PriestleyTaylor)`.
#' - Gs_pot    Potential/maximum surface conductance (mol m-2 s-1); defaults to 0.6 mol m-2 s-1;
#'                  only used if `approach = Val(:PenmanMonteith)`.
#' - missing_G_as_NA  if `true`, missing G are treated as `NA`s, otherwise set to 0. 
#' - missing_S_as_NA  if `true`, missing S are treated as `NA`s, otherwise set to 0. 
#' - Esat_formula  Optional: formula to be used for the calculation of esat and the slope of esat.
#'                      One of `"Sonntag_1990"` (Default), `"Alduchov_1996"`, or `"Allen_1998"`.
#'                      See [`Esat_slope`](@ref). 
#' - constants cp - specific heat of air for constant pressure (J K-1 kg-1) 
#'                  eps - ratio of the molecular weight of water vapor to dry air 
#'                  Pa2kPa - conversion pascal (Pa) to kilopascal (kPa) 
#'                  Rd - gas constant of dry air (J kg-1 K-1) (only used if `approach = Val(:PenmanMonteith)`) 
#'                  Rgas - universal gas constant (J mol-1 K-1) (only used if `approach = Val(:PenmanMonteith)`) 
#'                  Kelvin - conversion degree Celsius to Kelvin (only used if `approach = Val(:PenmanMonteith)`) 
#' 
#' # Details
#'  Potential evapotranspiration is calculated according to Priestley & Taylor, 1972
#'          if `approach = Val(:PriestleyTaylor)` (the default):
#' 
#'            ``LE_pot,PT = (\\alpha * \\Delta * (Rn - G - S)) / (\\Delta + \\gamma)``
#'
#'          ``\\alpha`` is the Priestley-Taylor coefficient, ``\\Delta`` is the slope 
#'          of the saturation vapor pressure curve (kPa K-1), and ``\\gamma`` is the 
#'          psychrometric constant (kPa K-1).
#'          if `approach = Val(:PenmanMonteith)`, potential evapotranspiration is calculated according
#'          to the Penman-Monteith equation:
#' 
#'          ``LE_pot,PM = (\\Delta * (Rn - G - S) + \\rho * cp * VPD * Ga) / (\\Delta + \\gamma * (1 + Ga/Gs_pot)``
#'          
#'          where ``\\Delta`` is the slope of the saturation vapor pressure curve (kPa K-1),
#'          ``\\rho`` is the air density (kg m-3), and ``\\gamma`` is the psychrometric constant (kPa K-1).
#'          The value of `Gs_pot` is typically a maximum value of Gs observed at the site, e.g. the 90th
#'          percentile of Gs within the growing season.
#'          
#' # Value
#'  a DataFrame with the following columns:
#'         - ET_pot: Potential evapotranspiration (kg m-2 s-1)
#'         - LE_pot: Potential latent heat flux (W m-2)
#'         
#' # Note
#' If the first argument `data` is provided (either a matrix or a DataFrame),
#'       the following variables can be provided as character (in which case they are interpreted as
#'       the column name of `data`) or as numeric vectors, in which case they are taken
#'       directly for the calculations. If `data` is not provided, all input variables have to be
#'       numeric vectors.        
#'   
#' #References
#' Priestley, C_H_B., Taylor, R_J., 1972: On the assessment of surface heat flux
#'             and evaporation using large-scale parameters. Monthly Weather Review 100, 81-92.  
#'             
#'             Allen, R_G., Pereira L_S., Raes D., Smith M., 1998: Crop evapotranspiration -
#'             Guidelines for computing crop water requirements - FAO Irrigation and drainage
#'             paper 56.
#'              
#'             Novick, K_A., et al. 2016: The increasing importance of atmospheric demand
#'             for ecosystem water and carbon fluxes. Nature Climate Change 6, 1023 - 1027.
#'             
#' #See also
#' [`surface_conductance`](@ref)
#'                                
#' ```@example; output = false
#' ``` 
#' # Calculate potential ET of a surface that receives a net radiation of 500 Wm-2
#' # using Priestley-Taylor:
#' potential_ET(Tair=30,pressure=100,Rn=500,alpha=1.26,approach=Val(:PriestleyTaylor))    
#' 
#' # Calculate potential ET for a surface with known Gs (0.5 mol m-2 s-1) and Ga (0.1 m s-1)
#' # using Penman-Monteith:
#' LE_pot_PM = potential_ET(Gs_pot=0.5,Tair=20,pressure=100,VPD=2,Ga=0.1,Rn=400,
#'                           approach=Val(:PenmanMonteith))[,"LE_pot"]
#' LE_pot_PM
#' 
#' # now cross-check with the inverted equation
#' surface_conductance(Tair=20,pressure=100,VPD=2,Ga=0.1,Rn=400,LE=LE_pot_PM)
"""
"""
function potential_ET(data,Tair="Tair",pressure="pressure",Rn="Rn",G=nothing,S=nothing,
                         VPD="VPD",Ga="Ga_h",approach=c(Val(:PriestleyTaylor),Val(:PenmanMonteith)),
                         alpha=1.26,Gs_pot=0.6,missing_G_as_NA=false,missing_S_as_NA=false,
                         Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                         constants=bigleaf_constants())
  
  approach = match_arg(approach)
  
  check_input(data,list(Tair,pressure,Rn,G,S))
  
  if(!isnothing(G))
    if (!missing_G_as_NA)
      G[ismissing(G)] = 0
end
else 
    cat("Ground heat flux G is not provided and set to 0.",fill=true)
    G = 0
end
  
  if(!isnothing(S))
    if(!missing_S_as_NA)
      S[ismissing(S)] = 0 
end
else 
    cat("Energy storage fluxes S are not provided and set to 0.",fill=true)
    S = 0
end
  
  gamma  = psychrometric_constant(Tair,pressure,constants)
  Delta  = Esat_slope(Tair,Esat_formula,constants)[,"Delta"]
  
  
  if (approach == Val(:PriestleyTaylor))
    
    LE_pot = (alpha * Delta * (Rn - G - S)) / (Delta + gamma)
    ET_pot = LE_to_ET(LE_pot,Tair)
    
elseif (approach == Val(:PenmanMonteith))
    
    check_input(data,list(Gs_pot,VPD,Ga))
    
    Gs_pot = mol_to_ms(Gs_pot,Tair=Tair,pressure=pressure,constants=constants)
    rho    = air_density(Tair,pressure,constants)
    
    LE_pot = (Delta * (Rn - G - S) + rho * constants[:cp] * VPD * Ga) / 
      (Delta + gamma * (1 + Ga / Gs_pot))
    ET_pot = LE_to_ET(LE_pot,Tair)
end
  
  return(DataFrame(ET_pot,LE_pot))
  
end




#' Reference Evapotranspiration
#' 
#' Reference evapotranspiration calculated from the Penman-Monteith
#'              equation with a prescribed surface conductance.
#'              This function is deprecated. Use potential_ET(...,approach=Val(:PenmanMonteith)) instead.
#' 
#' - data      Data_frame or matrix containing all required variables; optional
#' - Gs_ref    Reference surface conductance (m s-1); defaults to 0.0143 m s-1.
#' - Tair      Air temperature (degC)
#' - pressure  Atmospheric pressure (kPa)
#' - VPD       Vapor pressure deficit (kPa)
#' - Ga        Aerodynamic conductance to heat/water vapor (m s-1)
#' - Rn        Net radiation (W m-2)
#' - G         Ground heat flux (W m-2); optional
#' - S         Sum of all storage fluxes (W m-2); optional
#' - missing_G_as_NA  if `true`, missing G are treated as `NA`s, otherwise set to 0. 
#' - missing_S_as_NA  if `true`, missing S are treated as `NA`s, otherwise set to 0. 
#' - Esat_formula  Optional: formula to be used for the calculation of esat and the slope of esat.
#'                      One of `"Sonntag_1990"` (Default), `"Alduchov_1996"`, or `"Allen_1998"`.
#'                      See [`Esat_slope`](@ref). 
#' - constants cp - specific heat of air for constant pressure (J K-1 kg-1) 
#'                  eps - ratio of the molecular weight of water vapor to dry air 
#'                  Rd - gas constant of dry air (J kg-1 K-1) (only if `approach = Val(:PenmanMonteith)`) 
#'                  Rgas - universal gas constant (J mol-1 K-1) (only if `approach = Val(:PenmanMonteith)`) 
#'                  Kelvin - conversion degree Celsius to Kelvin (only if `approach = Val(:PenmanMonteith)`) 
#' 
#' @export                            
function reference_ET(data,Gs_ref=0.0143,Tair="Tair",pressure="pressure",VPD="VPD",Rn="Rn",Ga="Ga_h",
                         G=nothing,S=nothing,missing_G_as_NA=false,missing_S_as_NA=false,
                         Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                         constants=bigleaf_constants())
  
  stop("this function is deprecated (since bigleaf version 0.6.0). For the calculation of potential ET from the Penman-Monteith equation (as formerly calculated with this function), use function potential_ET() with the argument approach='Penman-Monteith'. Note that the default value for argument 'Gs_pot' is expressed now in mol m-2 s-1 for simplicity (0.6 mol m-2 s-1).")
  
end





#' Equilibrium and Imposed Evapotranspiration
#' 
#' Evapotranspiration (ET) split up into imposed ET and equilibrium ET.
#' 
#' - data      Data_frame or matrix containing all required input variables
#' - Tair      Air temperature (deg C)
#' - pressure  Atmospheric pressure (kPa)
#' - VPD       Air vapor pressure deficit (kPa)
#' - Gs        surface conductance to water vapor (m s-1)
#' - Rn        Net radiation (W m-2)
#' - G         Ground heat flux (W m-2); optional
#' - S         Sum of all storage fluxes (W m-2); optional
#' - missing_G_as_NA  if `true`, missing G are treated as `NA`s, otherwise set to 0. 
#' - missing_S_as_NA  if `true`, missing S are treated as `NA`s, otherwise set to 0.
#' - Esat_formula  Optional: formula to be used for the calculation of esat and the slope of esat.
#'                      One of `"Sonntag_1990"` (Default), `"Alduchov_1996"`, or `"Allen_1998"`.
#'                      See [`Esat_slope`](@ref). 
#' - constants cp - specific heat of air for constant pressure (J K-1 kg-1) 
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) 
#'                  Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#'                  
#' # Details
 Total evapotranspiration can be written in the form (Jarvis & McNaughton 1986):
#' 
#'            ``ET = \\Omega ET_eq + (1 - \\Omega)ET_imp``
#'          
#'          where ``\\Omega`` is the decoupling coefficient as calculated from
#'          [`decoupling`](@ref). `ET_eq` is the equilibrium evapotranspiration rate,
#'          the ET rate that would occur under uncoupled conditions, where the heat budget
#'          is dominated by radiation (when Ga -> 0):
#'          
#'            ``ET_eq = (\\Delta * (Rn - G - S) * \\lambda) / (\\Delta + \\gamma)``
#'          
#'          where ``\\Delta`` is the slope of the saturation vapor pressure curve (kPa K-1),
#'          ``\\lambda`` is the latent heat of vaporization (J kg-1), and ``\\gamma``
#'          is the psychrometric constant (kPa K-1).
#'          `ET_imp` is the imposed evapotranspiration rate, the ET rate
#'          that would occur under fully coupled conditions (when Ga -> inf):
#'          
#'            ``ET_imp = (\\rho * cp * VPD * Gs * \\lambda) / \\gamma``
#'          
#'          where ``\\rho`` is the air density (kg m-3).
#' 
#' # Note
#' Surface conductance (Gs) can be calculated with [`surface_conductance`](@ref).
#'       Aerodynamic conductance (Ga) can be calculated using [`aerodynamic_conductance`](@ref).
#'       
#' # Value
 A DataFrame with the following columns:
#'         - ET_eq: Equilibrium ET (kg m-2 s-1)
#'         - ET_imp: Imposed ET (kg m-2 s-1)
#'         - LE_eq: Equilibrium LE (W m-2)
#'         - LE_imp: Imposed LE (W m-2)      
#' 
#' #References
#' Jarvis, P_G., McNaughton, K_G., 1986: Stomatal control of transpiration:
#'             scaling up from leaf to region. Advances in Ecological Research 15, 1-49.
#'             
#'             Monteith, J_L., Unsworth, M_H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London. 
#'             
#' #See also
#' [`decoupling`](@ref)            
#'             
#' ```@example; output = false
#' ``` 
#' df = DataFrame(Tair=20,pressure=100,VPD=seq(0.5,4,0.5),
#'                  Gs_ms=seq(0.01,0.002,length_out=8),Rn=seq(50,400,50))            
#' equilibrium_imposed_ET(df)            
#'             
"""
"""
function equilibrium_imposed_ET(data,Tair="Tair",pressure="pressure",VPD="VPD",Gs="Gs_ms",
                                   Rn="Rn",G=nothing,S=nothing,missing_G_as_NA=false,missing_S_as_NA=false,
                                   Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                                   constants=bigleaf_constants())
  
  check_input(data,list(Tair,pressure,VPD,Rn,Gs,G,S))
  
  if(!isnothing(G))
    if (!missing_G_as_NA){G[ismissing(G)] = 0}
else 
    cat("Ground heat flux G is not provided and set to 0.",fill=true)
    G = 0
end
  
  if(!isnothing(S))
    if(!missing_S_as_NA){S[ismissing(S)] = 0 }
else 
    cat("Energy storage fluxes S are not provided and set to 0.",fill=true)
    S = 0
end
  
  rho    = air_density(Tair,pressure,constants)
  gamma  = psychrometric_constant(Tair,pressure,constants)
  Delta  = Esat_slope(Tair,Esat_formula,constants)[,"Delta"]
  
  LE_eq  = (Delta * (Rn - G - S)) / (gamma + Delta)
  LE_imp = (rho * constants[:cp] * Gs * VPD) / gamma
  
  ET_imp = LE_to_ET(LE_imp,Tair)
  ET_eq  = LE_to_ET(LE_eq,Tair)
  
  return(DataFrame(ET_eq,ET_imp,LE_eq,LE_imp))
end
