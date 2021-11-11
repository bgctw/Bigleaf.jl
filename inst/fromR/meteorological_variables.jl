#################################
### Meteorological variables ####
#################################

#' Air Density
#' 
#' Air density of moist air from air temperature and pressure.
#' 
#' - Tair      Air temperature (deg C)
#' - pressure  Atmospheric pressure (kPa)
#' - constants Kelvin - conversion degC to Kelvin 
#'                  Rd - gas constant of dry air (J kg-1 K-1) 
#'                  kPa2Pa - conversion kilopascal (kPa) to pascal (Pa)
#' 
#' # Details
 Air density (``\\rho``) is calculated as:
#' 
#'   ``\\rho = pressure / (Rd * Tair)``
#' 
#' # Value
 - ``\\rho``: air density (kg m-3)
#' 
#' ```@example; output = false
#' ``` 
#' # air density at 25degC and standard pressure (101.325kPa)
#' air_density(25,101.325)
#' 
#' #References
#' Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany. 
#' 
"""
"""
function air_density(Tair,pressure,constants=bigleaf_constants())
  
  Tair     = Tair + constants[:Kelvin]
  pressure = pressure * constants[:kPa2Pa]
  
  rho = pressure / (constants[:Rd] * Tair) 
  
  return(rho)
end



#' Atmospheric Pressure from Hypsometric Equation
#' 
#' An estimate of mean pressure at a given elevation as predicted by the
#'              hypsometric equation.
#'              
#' - elev      Elevation a_s_l. (m)
#' - Tair      Air temperature (deg C)
#' - VPD       Vapor pressure deficit (kPa); optional
#' - constants Kelvin- conversion degC to Kelvin 
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) 
#'                  Rd - gas constant of dry air (J kg-1 K-1) 
#'                  g -  gravitational acceleration (m s-2) 
#'                  Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' # Details
 Atmospheric pressure is approximated by the hypsometric equation:
#' 
#'          ``pressure = pressure_0 / (exp(g * elevation / (Rd Temp)))``
#'       
#' # Note
#' The hypsometric equation gives an estimate of the standard pressure
#'       at a given altitude. 
#'       If VPD is provided, humidity correction is applied and the
#'       virtual temperature instead of air temperature is used. VPD is 
#'       internally converted to specific humidity.
#'
#' # Value
 - pressure -: Atmospheric pressure (kPa)
#'                            
#' #References
#' Stull B., 1988: An Introduction to Boundary Layer Meteorology.
#'             Kluwer Academic Publishers, Dordrecht, Netherlands.
#' 
#' ```@example; output = false
#' ```
#' # mean pressure at 500m altitude at 25 deg C and VPD of 1 kPa
#' pressure_from_elevation(500,Tair=25,VPD=1)
#' 
#' @export                           
function pressure_from_elevation(elev,Tair,VPD=nothing,constants=bigleaf_constants())
  
  Tair     = Tair + constants[:Kelvin]
  
  if(isnothing(VPD))
    
    pressure = constants[:pressure0] / exp(constants[:g] * elev / (constants[:Rd]*Tair))
    
else 
    
    pressure1   = constants[:pressure0] / exp(constants[:g] * elev / (constants[:Rd]*Tair))
    Tv          = virtual_temp(Tair - constants[:Kelvin],pressure1 * constants[:Pa2kPa],
                                VPD,Esat_formula="Sonntag_1990",constants) + constants[:Kelvin]
    
    pressure    = constants[:pressure0] / exp(constants[:g] * elev / (constants[:Rd]*Tv))
end
  
  pressure = pressure * constants[:Pa2kPa]
  return(pressure)
end 




#' Saturation Vapor Pressure (Esat) and Slope of the Esat Curve
#' 
#' Calculates saturation vapor pressure (Esat) over water and the
#'              corresponding slope of the saturation vapor pressure curve.
#' 
#' - Tair      Air temperature (deg C)
#' - Esat_formula   Esat_formula to be used. Either `"Sonntag_1990"` (Default), 
#'                  `"Alduchov_1996"`, or `"Allen_1998"`.
#' - constants Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' # Details
 Esat (kPa) is calculated using the Magnus equation:
#' 
#'    ``Esat = a * exp((b * Tair) / (c + Tair)) / 1000``
#'  
#'  where the coefficients a, b, c take different values depending on the formula used.
#'  The default values are from Sonntag 1990 (a=611.2, b=17.62, c=243.12). This version
#'  of the Magnus equation is recommended by the WMO (WMO 2008; p1.4-29). Alternatively,
#'  parameter values determined by Alduchov & Eskridge 1996 or Allen et al. 1998 can be 
#'  used (see references).
#'  The slope of the Esat curve (``\\Delta``) is calculated as the first derivative of the function:
#'  
#'    ``\\Delta = dEsat / dTair``
#'  
#'  which is solved using `\link[stats]{D`}.
#' 
#' # Value
 A dataframe with the following columns: 
#'  - Esat: Saturation vapor pressure (kPa)
#'  - Delta: Slope of the saturation vapor pressure curve (kPa K-1)
#'    
#' #References
#' Sonntag D. 1990: Important new values of the physical constants of 1986, vapor 
#'             pressure formulations based on the ITS-90 and psychrometric formulae. 
#'             Zeitschrift fuer Meteorologie 70, 340-344.
#'             
#'             World Meteorological Organization 2008: Guide to Meteorological Instruments
#'             and Methods of Observation (WMO-No.8). World Meteorological Organization,
#'             Geneva. 7th Edition.
#'             
#'             Alduchov, O. A. & Eskridge, R. E., 1996: Improved Magnus form approximation of 
#'             saturation vapor pressure. Journal of Applied Meteorology, 35, 601-609
#'             
#'             Allen, R_G., Pereira, L_S., Raes, D., Smith, M., 1998: Crop evapotranspiration -
#'             Guidelines for computing crop water requirements - FAO irrigation and drainage
#'             paper 56, FAO, Rome.
#' 
#' ```@example; output = false
#' ``` 
#' Esat_slope(seq(0,45,5))[,"Esat"]  # Esat in kPa
#' Esat_slope(seq(0,45,5))[,"Delta"] # the corresponding slope of the Esat curve (Delta) in kPa K-1
#'        
#' @importFrom stats D                  
"""
"""
function Esat_slope(Tair,Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                       constants=bigleaf_constants())
  
  Esat_formula = match_arg(Esat_formula)
  
  if (Esat_formula == "Sonntag_1990")
    a = 611.2
    b = 17.62
    c = 243.12
elseif (Esat_formula == "Alduchov_1996")
    a = 610.94
    b = 17.625
    c = 243.04
elseif (Esat_formula == "Allen_1998")
    a = 610.8
    b = 17.27
    c = 237.3
end
  
  # saturation vapor pressure
  Esat = a * exp((b * Tair) / (c + Tair))
  Esat = Esat * constants[:Pa2kPa]
  
  # slope of the saturation vapor pressure curve
  Delta = eval(D(expression(a * exp((b * Tair) / (c + Tair))),name="Tair"))
  Delta = Delta * constants[:Pa2kPa]
  
  return(DataFrame(Esat,Delta))
end



#' Psychrometric Constant
#' 
#' Calculates the psychrometric 'constant'.
#' 
#' - Tair      Air temperature (deg C)
#' - pressure  Atmospheric pressure (kPa)
#' - constants cp - specific heat of air for constant pressure (J K-1 kg-1) 
#'                  eps - ratio of the molecular weight of water vapor to dry air (-)
#'                  
#' # Details
 The psychrometric constant (``\\gamma``) is given as:
#' 
#'    ``\\gamma = cp * pressure / (eps * \\lambda)``
#'  
#'  where ``\\lambda`` is the latent heat of vaporization (J kg-1), 
#'  as calculated from [`latent_heat_vaporization`](@ref).
#'  
#' # Value
 - ``\\gamma`` -: the psychrometric constant (kPa K-1)
#'  
#' #References
#' Monteith J_L., Unsworth M_H., 2008: Principles of Environmental Physics.
#'             3rd Edition. Academic Press, London. 
#' 
#' ```@example; output = false
#' ``` 
#' psychrometric_constant(seq(5,45,5),100)
#' 
"""
"""
function psychrometric_constant(Tair,pressure,constants=bigleaf_constants())
  
  lambda = latent_heat_vaporization(Tair)
  gamma  = (constants[:cp] * pressure) / (constants[:eps] * lambda)
  
  return(gamma)
end



#' Latent Heat of Vaporization
#' 
#' Latent heat of vaporization as a function of air temperature.
#' 
#' - Tair  Air temperature (deg C)
#' 
#' # Details
 The following Esat_formula is used:
#' 
#'   ``\\lambda = (2.501 - 0.00237*Tair)10^6``
#' 
#' # Value
 - ``\\lambda`` -: Latent heat of vaporization (J kg-1) 
#' 
#' #References
#' Stull, B., 1988: An Introduction to Boundary Layer Meteorology (p.641)
#'             Kluwer Academic Publishers, Dordrecht, Netherlands
#'             
#'             Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany. 
#' 
#' ```@example; output = false
#' ``` 
#' latent_heat_vaporization(seq(5,45,5))             
#'             
"""
"""
function latent_heat_vaporization(Tair) 
  
  k1 = 2.501
  k2 = 0.00237
  lambda = ( k1 - k2 * Tair ) * 1e+06
  
  return(lambda)
end






#' Solver Function for Wet-Bulb Temperature
#' 
#' Solver function used in wetbulb_temp()
#' 
#' - ea           Air vapor pressure (kPa)
#' - Tair         Air temperature (degC)
#' - gamma        Psychrometric constant (kPa K-1)
#' - accuracy     Accuracy of the result (degC)
#' - Esat_formula Optional: Esat_formula to be used for the calculation of esat and the slope of esat.
#'                     One of `"Sonntag_1990"` (Default), `"Alduchov_1996"`, or `"Allen_1998"`.
#'                     See [`Esat_slope`](@ref). 
#' - constants    Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' # Note
#' Arguments `accuracy` and `Esat_formula` are passed to this function by wetbulb_temp().
#' 
#' @importFrom stats optimize 
#' 
#' @keywords internal
function wetbulb_solver(ea,Tair,gamma,accuracy,Esat_formula,constants=bigleaf_constants())
  wetbulb_optim = optimize(function(Tw){abs(ea - c((Esat_slope(Tw,Esat_formula,constants)[,"Esat"] - 0.93*gamma*(Tair - Tw))))},
                            interval=c(-100,100),tol=accuracy)
  return(wetbulb_optim)
end



#' Wet-Bulb Temperature
#' 
#' calculates the wet bulb temperature, i.e. the temperature
#'              that the air would have if it was saturated.
#' 
#' - Tair      Air temperature (deg C)
#' - pressure  Atmospheric pressure (kPa)
#' - VPD       Vapor pressure deficit (kPa)
#' - accuracy  Accuracy of the result (deg C)
#' - Esat_formula  Optional: Esat_formula to be used for the calculation of esat and the slope of esat.
#'                      One of `"Sonntag_1990"` (Default), `"Alduchov_1996"`, or `"Allen_1998"`.
#'                      See [`Esat_slope`](@ref). 
#' - constants cp - specific heat of air for constant pressure (J K-1 kg-1) 
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) 
#'                  Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' # Details
 Wet-bulb temperature (Tw) is calculated from the following expression:
#'          
#'            ``e = Esat(Tw) - gamma* (Tair - Tw)``
#'          
#'          The equation is solved for Tw using `\link[stats]{optimize`}.
#'          Actual vapor pressure e (kPa) is calculated from VPD using the function [`VPD_to_e`](@ref).
#'          The psychrometric constant gamma (kPa K-1) is calculated from [`psychrometric_constant`](@ref).
#'          
#' # Value
 - Tw -: wet-bulb temperature (degC)      
#'              
#' #References
#' Monteith J_L., Unsworth M_H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London.
#'             
#' ```@example; output = false
#' ``` 
#' wetbulb_temp(Tair=c(20,25),pressure=100,VPD=c(1,1.6))             
#'        
#' @importFrom stats optimize                  
"""
"""
function wetbulb_temp(Tair,pressure,VPD,accuracy=1e-03,Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                         constants=bigleaf_constants())
  
  if (!is_numeric(accuracy))
    stop("'accuracy' must be numeric!")
end
  
  if (accuracy > 1)
    print("'accuracy' is set to 1 degC")
    accuracy = 1
end
  
  # determine number of digits to print
  ndigits = as_numeric(strsplit(format(accuracy,scientific = true),"-")[[1]][2])
  ndigits = ifelse(ismissing(ndigits),0,ndigits)
  
  
  gamma  = psychrometric_constant(Tair,pressure)
  ea     = VPD_to_e(VPD,Tair,Esat_formula)
  
  Tw = sapply(seq_along(ea),function(i) round(wetbulb_solver(ea[i],Tair[i],gamma[i],
                                                              accuracy=accuracy,Esat_formula,constants)$minimum,ndigits))
  
  return(Tw)
  
end











#' Solver function for dew point temperature
#' 
#' Solver function used in dew_point()
#' 
#' - ea           Air vapor pressure (kPa)
#' - accuracy     Accuracy of the result (degC)
#' - Esat_formula Optional: Esat_formula to be used for the calculation of esat and the slope of esat.
#'                     One of `"Sonntag_1990"` (Default), `"Alduchov_1996"`, or `"Allen_1998"`.
#'                     See [`Esat_slope`](@ref). 
#' - constants    Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' # Note
#' Arguments `accuracy` and `Esat_formula` are passed to this function by dew_point().
#' 
#' @importFrom stats optimize 
#' 
#' @keywords internal
function dew_point_solver(ea,accuracy,Esat_formula,constants=bigleaf_constants())
  
  Td_optim = optimize(function(Td){abs(ea - Esat_slope(Td,Esat_formula,constants)[,"Esat"])},
                       interval=c(-100,100),tol=accuracy)
  return(Td_optim)
end



#' Dew Point
#' 
#' calculates the dew point, the temperature to which air must be 
#'              cooled to become saturated (i.e. e = Esat(Td))
#'
#' - Tair     Air temperature (degC)
#' - VPD      Vapor pressure deficit (kPa)
#' - accuracy Accuracy of the result (deg C)
#' - Esat_formula  Optional: Esat_formula to be used for the calculation of esat and the slope of esat.
#'                      One of `"Sonntag_1990"` (Default), `"Alduchov_1996"`, or `"Allen_1998"`.
#'                      See [`Esat_slope`](@ref). 
#' - constants Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' # Details
 Dew point temperature (Td) is defined by:
#' 
#'           ``e = Esat(Td)``
#'    
#'          where e is vapor pressure of the air and Esat is the vapor pressure deficit.
#'          This equation is solved for Td using `\link[stats]{optimize`}.
#'          
#' # Value
 - Td -: dew point temperature (degC)
#' 
#' #References
#' Monteith J_L., Unsworth M_H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London.
#'             
#' ```@example; output = false
#' ```
#' dew_point(c(25,30),1.5)                
#' 
#' @importFrom stats optimize 
#' @export              
function dew_point(Tair,VPD,accuracy=1e-03,Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                      constants=bigleaf_constants())
  
  if (!is_numeric(accuracy))
    stop("'accuracy' must be numeric!")
end
  
  if (accuracy > 1)
    print("'accuracy' is set to 1 degC")
    accuracy = 1
end
  
  # determine number of digits to print
  ndigits = as_numeric(strsplit(format(accuracy,scientific = true),"-")[[1]][2])
  ndigits = ifelse(ismissing(ndigits),0,ndigits)
  
  ea = VPD_to_e(VPD,Tair,Esat_formula)
  Td = sapply(seq_along(ea),function(i) round(dew_point_solver(ea[i],accuracy=accuracy,
                                                                Esat_formula,constants)$minimum,ndigits))
  
  return(Td)
end







#' Virtual Temperature
#' 
#' Virtual temperature, defined as the temperature at which dry air would have the same
#'              density as moist air at its actual temperature.
#' 
#' - Tair      Air temperature (deg C)
#' - pressure  Atmospheric pressure (kPa)
#' - VPD       Vapor pressure deficit (kPa)
#' - Esat_formula  Optional: Esat_formula to be used for the calculation of esat and the slope of esat. 
#'                      One of `"Sonntag_1990"` (Default), `"Alduchov_1996"`, or `"Allen_1998"`.
#'                      See [`Esat_slope`](@ref). 
#' - constants Kelvin - conversion degree Celsius to Kelvin 
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) 
#' 
#' # Details
 the virtual temperature is given by:
#'  
#'    ``Tv = Tair / (1 - (1 - eps) e/pressure)``
#' 
#'  where Tair is in Kelvin (converted internally). Likewise, VPD is converted 
#'  to actual vapor pressure (e in kPa) with [`VPD_to_e`](@ref) internally.
#' 
#' # Value
 - Tv -: virtual temperature (deg C)
#' 
#' #References
#' Monteith J_L., Unsworth M_H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London.
#'  
#' ```@example; output = false
#' ``` 
#' virtual_temp(25,100,1.5)                        
#'               
"""
"""
function virtual_temp(Tair,pressure,VPD,Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                         constants=bigleaf_constants())
  
  e    = VPD_to_e(VPD,Tair,Esat_formula)
  Tair = Tair + constants[:Kelvin]
  
  Tv = Tair / (1 - (1 - constants[:eps]) * e/pressure) 
  Tv = Tv - constants[:Kelvin]
  
  return(Tv)
end



#' Kinematic Viscosity of Air
#' 
#' calculates the kinematic viscosity of air.
#' 
#' - Tair      Air temperature (deg C)
#' - pressure  Atmospheric pressure (kPa)
#' - constants Kelvin - conversion degree Celsius to Kelvin 
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) 
#'                  Tair0 - reference air temperature (K) 
#'                  kPa2Pa - conversion kilopascal (kPa) to pascal (Pa)
#' 
#' # Details
 where v is the kinematic viscosity of the air (m2 s-1), 
#'          given by (Massman 1999b):
#'          
#'            ``v = 1.327 * 10^-5(pressure0/pressure)(Tair/Tair0)^1.81``
#'          
#' # Value
 - v -: kinematic viscosity of air (m2 s-1)
#' 
#' #References
#' Massman, W_J., 1999b: Molecular diffusivities of Hg vapor in air, 
#'             O2 and N2 near STP and the kinematic viscosity and thermal diffusivity
#'             of air near STP. Atmospheric Environment 33, 453-457.      
#'             
#' ```@example; output = false
#' ``` 
#' kinematic_viscosity(25,100)    
#' 
#' @export         
function kinematic_viscosity(Tair,pressure,constants=bigleaf_constants())
  
  Tair     = Tair + constants[:Kelvin]
  pressure = pressure * constants[:kPa2Pa]
  
  v  = 1.327e-05*(constants[:pressure0]/pressure)*(Tair/constants[:Tair0])^1.81
  return(v)
end
