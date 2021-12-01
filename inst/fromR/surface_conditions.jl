#############################
#### Surface conditions  ####
#############################

#' Big-Leaf Surface Conditions
#' 
#' Calculates meteorological conditions at the big-leaf surface
#'              by inverting bulk transfer equations for water, energy, and carbon
#'              fluxes.
#' 
#' - data             DataFrame or matrix containing all required input variables
#' - Tair             Air temperature (deg C)
#' - pressure         Atmospheric pressure (kPa)
#' - H                Sensible heat flux (W m-2)
#' - LE               Latent heat flux (W m-2)
#' - VPD              Vapor pressure deficit (kPa)
#' - Ga               Aerodynamic conductance for heat/water vapor (m s-1)
#' - calc_surface_CO2 Calculate surface CO2 concentration? Defaults to `false`.
#' - Ca               Atmospheric CO2 concentration (mol mol-1). Required if `calc_surface_CO2 = true`.
#' - NEE              Net ecosystem exchange (umol m-2 s-1). Required if `calc_surface_CO2 = true`.
#' - Ga_CO2           Aerodynamic conductance for CO2 (m s-1). Required if `calc_surface_CO2 = true`.
#' - Esat_formula     Optional: Esat_formula to be used for the calculation of esat and the slope of esat.
#'                         One of `"Sonntag_1990"` (Default), `"Alduchov_1996"`, or `"Allen_1998"`.
#'                         See [`Esat_slope`](@ref).           
#' - constants        cp - specific heat of air for constant pressure (J K-1 kg-1)  
#'                         eps - ratio of the molecular weight of water vapor to dry air (-) 
#'                         Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' # Details
 Canopy surface temperature and humidity are calculated by inverting bulk transfer equations of
#'          sensible and latent heat, respectively. 'Canopy surface' in this case refers to 
#'          the surface of the big-leaf (i.e. at height d + z0h; the apparent sink of sensible heat and water vapor).
#'          Aerodynamic canopy surface temperature is given by:
#'          
#'            ``Tsurf = Tair + H / (\\rho * cp * Ga)``
#'          
#'          where ``\\rho`` is air density (kg m-3). 
#'          Vapor pressure at the canopy surface is:
#'          
#'            ``esurf = e + (LE * \\gamma)/(Ga * \\rho * cp)``
#'          
#'          where ``\\gamma`` is the psychrometric constant (kPa K-1).
#'          Vapor pressure deficit (VPD) at the canopy surface is calculated as:
#'          
#'            ``VPD_surf = Esat_surf - esurf``
#'            
#'          CO2 concentration at the canopy surface is given by:
#'          
#'            ``Ca_surf = Ca + NEE / Ga_CO2``
#'          
#'          Note that Ga is assumed to be equal for water vapor and sensible heat.
#'          Ga is further assumed to be the inverse of the sum of the turbulent part
#'          and the canopy boundary layer conductance (1/Ga = 1/Ga_m + 1/Gb; 
#'          see [`aerodynamic_conductance!`](@ref)). Ga_CO2, the aerodynamic conductance
#'          for CO2 is also calculated by [`aerodynamic_conductance!`](@ref).
#'          If Ga is replaced by Ga_m (i.e. only the turbulent conductance part), 
#'          the results of the functions represent conditions outside the canopy
#'          boundary layer, i.e. in the canopy airspace.
#' 
#' # Note
#' The following sign convention for NEE is employed (relevant if 
#'       `calc_surface_CO2 = true`): 
#'       negative values of NEE denote net CO2 uptake by the ecosystem.
#' 
#' # Value
 a DataFrame with the following columns:
#'         - Tsurf: Surface temperature (deg C) 
#'         - esat_surf: Saturation vapor pressure at the surface (kPa) 
#'         - esurf: vapor pressure at the surface (kPa) 
#'         - VPD_surf: vapor pressure deficit at the surface (kPa) 
#'         - qsurf: specific humidity at the surface (kg kg-1) 
#'         - rH_surf: relative humidity at the surface (-) 
#'         - Ca_surf: CO2 concentration at the surface (umol mol-1)             
#'         
#' ```@example; output = false
#' ```
#' # calculate surface temperature, water vapor, VPD etc. at the surface
#' # for a given temperature and turbulent fluxes, and under different 
#' # aerodynamic conductance.
#' surface_conditions(Tair=25,pressure=100,LE=100,H=200,VPD=1.2,Ga=c(0.02,0.05,0.1)) 
#'          
#' # now calculate also surface CO2 concentration
#' surface_conditions(Tair=25,pressure=100,LE=100,H=200,VPD=1.2,Ga=c(0.02,0.05,0.1),
#'                    Ca=400,Ga_CO2=c(0.02,0.05,0.1),NEE=-20,calc_surface_CO2=true)
#'                    
#' #References
#' Knauer, J. et al., 2018: Towards physiologically meaningful water-use efficiency estimates
#'             from eddy covariance data. Global Change Biology 24, 694-710.
#'             
#'             Blanken, P_D. & Black, T_A., 2004: The canopy conductance of a boreal aspen forest,
#'             Prince Albert National Park, Canada. Hydrological Processes 18, 1561-1578.
#'             
#'             Shuttleworth, W. J., Wallace, J_S., 1985: Evaporation from sparse crops-
#'             an energy combination theory. Quart. J. R. Met. Soc. 111, 839-855.
#'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
#' @export 
function surface_conditions(data,Tair="Tair",pressure="pressure",LE="LE",H="H",
                               VPD="VPD",Ga="Ga_h",calc_surface_CO2=false,Ca="Ca",Ga_CO2="Ga_CO2",
                               NEE="NEE",Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                               constants=BigleafConstants())
  
  check_input(data,list(Tair,pressure,LE,H,VPD,Ga))
  
  rho   = air_density(Tair,pressure,constants)
  gamma = psychrometric_constant(Tair,pressure,constants)
  
  # 1) Temperature
  Tsurf = Tair + H / (rho * constants.cp * Ga)
  
  # 2) Humidity
  esat      = Esat_slope(Tair,Esat_formula,constants)[,"Esat"]
  e         = esat - VPD
  esat_surf = Esat_slope(Tsurf,Esat_formula,constants)[,"Esat"]
  esurf     = e + (LE * gamma)/(Ga * rho * constants.cp)
  VPD_surf  = pmax(esat_surf - esurf,0)
  qsurf     = VPD_to_q(VPD_surf,Tsurf,pressure,Esat_formula,constants)
  rH_surf   = VPD_to_rH(VPD_surf,Tsurf,Esat_formula)
  
  # 3) CO2 concentration
  if (calc_surface_CO2)
    check_input(data,Ca,NEE,Ga_CO2)
    Ca_surf = surface_CO2(Ca,NEE,Ga_CO2,Tair,pressure)
else 
    Ca_surf = rep(NA_integer_,length(Tair))
end
  
  return(DataFrame(Tsurf,esat_surf,esurf,VPD_surf,qsurf,rH_surf,Ca_surf))
end



#' CO2 Concentration at the Canopy Surface
#'
#' the CO2 concentration at the canopy surface derived from net ecosystem
#'              CO2 exchange and measured atmospheric CO2 concentration.
#'              
#' - Ca       Atmospheric CO2 concentration (umol mol-1)
#' - NEE      Net ecosystem exchange (umol CO2 m-2 s-1)
#' - Ga_CO2   Aerodynamic conductance for CO2 (m s-1)
#' - Tair     Air temperature (degC)
#' - pressure Atmospheric pressure (kPa)
#' 
#' # Details
 CO2 concentration at the canopy surface is calculated as:
#' 
#'          ``Ca_surf = Ca + NEE / Ga_CO2``
#'        
#'        Note that this equation can be used for any gas measured (with NEE
#'        replaced by the net exchange of the respective gas and Ga_CO2 by the Ga of 
#'        that gas).
#' 
#' # Note
#' the following sign convention is employed: negative values of NEE denote
#'       net CO2 uptake by the ecosystem.
#' 
#' # Value
 - Ca_surf -: CO2 concentration at the canopy surface (umol mol-1)
#' 
#' ```@example; output = false
#' ``` 
#' surface_CO2(Ca=400,NEE=-30,Ga_CO2=0.05,Tair=25,pressure=100)
#' 
"""
"""
function surface_CO2(Ca,NEE,Ga_CO2,Tair,pressure)
  
  Ga_CO2 = ms_to_mol(Ga_CO2,Tair,pressure)
  
  Ca_surf = Ca + NEE/Ga_CO2
  
  return(Ca_surf)
  
end



#' Radiometric Surface Temperature
#' 
#' Radiometric surface temperature from longwave radiation
#'              measurements.
#'              
#' - data        DataFrame or matrix containing all required input variables            
#' - LW_up       Longwave upward radiation (W m-2)
#' - LW_down     Longwave downward radiation (W m-2)
#' - emissivity  Emissivity of the surface (-)
#' - constants   sigma - Stefan-Boltzmann constant (W m-2 K-4) 
#'                    Kelvin - conversion degree Celsius to Kelvin 
#' 
#' # Details
 Radiometric surface temperature (Trad) is calculated as:
#' 
#'            ``Trad = ((LW_up - (1 - \epsilon)*LW_down) / (\sigma \epsilon))^(1/4)``   
#' 
#' # Value
 a DataFrame with the following columns:
#'         - Trad_K: Radiometric surface temperature (K) 
#'         - Trad_degC: Radiometric surface temperature (degC) 
#' 
#' ```@example; output = false
#' ``` 
#' # determine radiometric surface temperature for the site DE-Tha in June 2014 
#' # assuming an emissivity of 0.98.
#' # (Note that variable 'LW_down' was only included for the DE-Tha example dataset
#' # and not for the others due restrictions on file size) 
#' Trad = radiometric_surface_temp(DE_Tha_Jun_2014,emissivity=0.98)
#' summary(Trad)
#' 
#' #References
#' Wang, W., Liang, S., Meyers, T. 2008: Validating MODIS land surface
#'             temperature products using long-term nighttime ground measurements.
#'             Remote Sensing of Environment 112, 623-635.
#' 
"""
"""
function radiometric_surface_temp(data,LW_up="LW_up",LW_down="LW_down",
                                     emissivity,constants=BigleafConstants())
  
  check_input(data,list(LW_up,LW_down))
  
  Trad_K    = ((LW_up - (1 - emissivity)*LW_down) / (constants.sigma * emissivity))^(1/4)
  Trad_degC = Trad_K - constants.Kelvin
  
  return(DataFrame(Trad_K,Trad_degC))
end


### can be added at some point, but the emissivity values seem to be fairly low
function # emissivity_from_albedo(albedo)
#   
#   emissivity = -0.16*albedo + 0.99
#   
#   return(emissivity)
# }

