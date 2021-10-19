#' Air Density
#' 
#' Air density of moist air from air temperature and pressure_
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param constants Kelvin - conversion degC to Kelvin \cr
#'                  Rd - gas constant of dry air (J kg-1 K-1) \cr
#'                  kPa2Pa - conversion kilopascal (kPa) to pascal (Pa)
#' 
#' @details Air density (\eqn{\rho}) is calculated as:
#' 
#'   \deqn{\rho = pressure / (Rd * Tair)}
#' 
#' @return \item{\eqn{\rho}}{air density (kg m-3)}
#' 
#' @examples 
#' # air density at 25degC and standard pressure (101.325kPa)
#' air_density(25,101.325)
#' 
#' @references Foken, T, 2008: Micrometeorology_ Springer, Berlin, Germany_ 
#' 
#' @export
function air_density(Tair,pressure; constants=bigleaf_constants())
  Tair_K     = Tair + constants[:Kelvin]
  pressure_Pa = pressure * constants[:kPa2Pa]
  rho = pressure_Pa / (constants[:Rd] * Tair_K) 
end



#' Atmospheric Pressure from Hypsometric Equation
#' 
#' An estimate of mean pressure at a given elevation as predicted by the
#'              hypsometric equation_
#'              
#' @param elev      Elevation a_s_l_ (m)
#' @param Tair      Air temperature (deg C)
#' @param VPD       Vapor pressure deficit (kPa); optional
#' @param constants Kelvin- conversion degC to Kelvin \cr
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) \cr
#'                  Rd - gas constant of dry air (J kg-1 K-1) \cr
#'                  g -  gravitational acceleration (m s-2) \cr
#'                  Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' @details Atmospheric pressure is approximated by the hypsometric equation:
#' 
#'          \deqn{pressure = pressure_0 / (exp(g * elevation / (Rd Temp)))}
#'       
#' @note The hypsometric equation gives an estimate of the standard pressure
#'       at a given altitude_ 
#'       If VPD is provided, humidity correction is applied and the
#'       virtual temperature instead of air temperature is used_ VPD is 
#'       internally converted to specific humidity_
#'
#' @return \item{pressure -}{Atmospheric pressure (kPa)}
#'                            
#' @references Stull B_, 1988: An Introduction to Boundary Layer Meteorology_
#'             Kluwer Academic Publishers, Dordrecht, Netherlands_
#' 
#' @examples
#' # mean pressure at 500m altitude at 25 deg C and VPD of 1 kPa
#' pressure_from_elevation(500,Tair=25,VPD=1)
#' 
#' @export                           
function pressure_from_elevation(elev,Tair,VPD=missing; constants=bigleaf_constants())
  Tair_K     = Tair + constants[:Kelvin]
  if is_missing(VPD)
    pressure = constants[:pressure0] / exp(constants[:g] * elev / (constants[:Rd]*Tair_K))
  else 
    pressure1   = constants[:pressure0] / exp(constants[:g] * elev / (constants[:Rd]*Tair_K))
    Tv          = virtual_temp(Tair_K - constants[:Kelvin],pressure1 * constants[:Pa2kPa],
                                VPD,Esat_formula="Sonntag_1990",constants) + constants[:Kelvin]
    pressure    = constants[:pressure0] / exp(constants[:g] * elev / (constants[:Rd]*Tv))
  end
  pressure = pressure * constants[:Pa2kPa]
end 

#' Psychrometric Constant
#' 
#' Calculates the psychrometric 'constant'_
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-)
#'                  
#' @details The psychrometric constant (\eqn{\gamma}) is given as:
#' 
#'    \deqn{\gamma = cp * pressure / (eps * \lambda)}
#'  
#'  where \eqn{\lambda} is the latent heat of vaporization (J kg-1), 
#'  as calculated from \code{\link{latent_heat_vaporization}}_
#'  
#' @return \item{\eqn{\gamma} -}{the psychrometric constant (kPa K-1)}
#'  
#' @references Monteith J_L_, Unsworth M_H_, 2008: Principles of Environmental Physics_
#'             3rd Edition_ Academic Press, London_ 
#' 
#' @examples 
#' psychrometric_constant(seq(5,45,5),100)
#' 
#' @export
function psychrometric_constant(Tair,pressure; constants=bigleaf_constants())
  lambda = latent_heat_vaporization(Tair)
  gamma  = (constants[:cp] * pressure) / (constants[:eps] * lambda)
end

#' Latent Heat of Vaporization
#' 
#' Latent heat of vaporization as a function of air temperature_
#' 
#' @param Tair  Air temperature (deg C)
#' 
#' @details The following formula is used:
#' 
#'   \deqn{\lambda = (2.501 - 0.00237*Tair)10^6}
#' 
#' @return \item{\eqn{\lambda} -}{Latent heat of vaporization (J kg-1)} 
#' 
#' @references Stull, B_, 1988: An Introduction to Boundary Layer Meteorology (p_641)
#'             Kluwer Academic Publishers, Dordrecht, Netherlands
#'             
#'             Foken, T, 2008: Micrometeorology_ Springer, Berlin, Germany_ 
#' 
#' @examples 
#' latent_heat_vaporization(seq(5,45,5))             
#'             
#' @export
function latent_heat_vaporization(Tair) 
  k1 = 2.501
  k2 = 0.00237
  lambda = ( k1 - k2 * Tair ) * 1e+06
end

#' Solver Function for Wet-Bulb Temperature
#' 
#' Solver function used in wetbulb_temp()
#' 
#' @param ea           Air vapor pressure (kPa)
#' @param Tair         Air temperature (degC)
#' @param gamma        Psychrometric constant (kPa K-1)
#' @param accuracy     Accuracy of the result (degC)
#' @param Esat_formula Optional: formula to be used for the calculation of esat and the slope of esat_
#'                     One of \code{"Sonntag_1990"} (Default), \code{"Alduchov_1996"}, or \code{"Allen_1998"}_
#'                     See \code{\link{Esat_slope}}_ 
#' @param constants    Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' @note Arguments \code{accuracy} and \code{Esat_formula} are passed to this function by wetbulb_temp()_
#' 
#' @importFrom stats optimize 
#' 
#' @keywords internal
# function wetbulb_solver(ea,Tair,gamma,accuracy,Esat_formula; constants=bigleaf_constants())
#   wetbulb_optim = optimize(function(Tw)
#   abs(ea - c((Esat_slope(Tw,Esat_formula,constants)[,"Esat"] - 0.93*gamma*(Tair - Tw))))end,
#                             interval=c(-100,100),tol=accuracy)
#   return(wetbulb_optim)
# end



#' Wet-Bulb Temperature
#' 
#' calculates the wet bulb temperature, i_e_ the temperature
#'              that the air would have if it was saturated_
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param accuracy  Accuracy of the result (deg C)
#' @param Esat_formula  Optional: formula to be used for the calculation of esat and the slope of esat_
#'                      One of \code{"Sonntag_1990"} (Default), \code{"Alduchov_1996"}, or \code{"Allen_1998"}_
#'                      See \code{\link{Esat_slope}}_ 
#' @param constants cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) \cr
#'                  Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' @details Wet-bulb temperature (Tw) is calculated from the following expression:
#'          
#'            \deqn{e = Esat(Tw) - gamma* (Tair - Tw)}
#'          
#'          The equation is solved for Tw using \code{\link[stats]{optimize}}_
#'          Actual vapor pressure e (kPa) is calculated from VPD using the function \code{\link{VPD_to_e}}_
#'          The psychrometric constant gamma (kPa K-1) is calculated from \code{\link{psychrometric_constant}}_
#'          
#' @return \item{Tw -}{wet-bulb temperature (degC)}      
#'              
#' @references Monteith J_L_, Unsworth M_H_, 2008: Principles of Environmental Physics_
#'             3rd edition_ Academic Press, London_
#'             
#' @examples 
#' wetbulb_temp(Tair=c(20,25),pressure=100,VPD=c(1,1.6))             
#'        
#' @importFrom stats optimize                  
#' @export
# function wetbulb_temp(Tair,pressure,VPD;accuracy=1e-03,Esat_formula="Sonntag_1990",
#                          constants=bigleaf_constants()){
#   if (!is_numeric(accuracy)){
#     stop("'accuracy' must be numeric!")
#   end
#   if (accuracy > 1){
#     print("'accuracy' is set to 1 degC")
#     accuracy = 1
#   end
#   # determine number of digits to print
#   ndigits = as_numeric(strsplit(format(accuracy,scientific = TRUE),"-")[[1]][2])
#   ndigits = ifelse(is_na(ndigits),0,ndigits)
#   gamma  = psychrometric_constant(Tair,pressure)
#   ea     = VPD_to_e(VPD,Tair,Esat_formula)
#   Tw = sapply(seq_along(ea),function(i) round(wetbulb_solver(ea[i],Tair[i],gamma[i],
#                                                               accuracy=accuracy,Esat_formula,constants)$minimum,ndigits))
#   return(Tw)
# end


#' Solver function for dew point temperature
#' 
#' Solver function used in dew_point()
#' 
#' @param ea           Air vapor pressure (kPa)
#' @param accuracy     Accuracy of the result (degC)
#' @param Esat_formula Optional: formula to be used for the calculation of esat and the slope of esat_
#'                     One of \code{"Sonntag_1990"} (Default), \code{"Alduchov_1996"}, or \code{"Allen_1998"}_
#'                     See \code{\link{Esat_slope}}_ 
#' @param constants    Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' @note Arguments \code{accuracy} and \code{Esat_formula} are passed to this function by dew_point()_
#' 
#' @importFrom stats optimize 
#' 
#' @keywords internal
# function dew_point_solver(ea,accuracy,Esat_formula; constants=bigleaf_constants()){
  
#   Td_optim = optimize(function(Td){abs(ea - Esat_slope(Td,Esat_formula,constants)[,"Esat"])end,
#                        interval=c(-100,100),tol=accuracy)
#   return(Td_optim)
# end



#' Dew Point
#' 
#' calculates the dew point, the temperature to which air must be 
#'              cooled to become saturated (i_e_ e = Esat(Td))
#'
#' @param Tair     Air temperature (degC)
#' @param VPD      Vapor pressure deficit (kPa)
#' @param accuracy Accuracy of the result (deg C)
#' @param Esat_formula  Optional: formula to be used for the calculation of esat and the slope of esat_
#'                      One of \code{"Sonntag_1990"} (Default), \code{"Alduchov_1996"}, or \code{"Allen_1998"}_
#'                      See \code{\link{Esat_slope}}_ 
#' @param constants Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' @details Dew point temperature (Td) is defined by:
#' 
#'           \deqn{e = Esat(Td)}
#'    
#'          where e is vapor pressure of the air and Esat is the vapor pressure deficit_
#'          This equation is solved for Td using \code{\link[stats]{optimize}}_
#'          
#' @return \item{Td -}{dew point temperature (degC)}
#' 
#' @references Monteith J_L_, Unsworth M_H_, 2008: Principles of Environmental Physics_
#'             3rd edition_ Academic Press, London_
#'             
#' @examples
#' dew_point(c(25,30),1.5)                
#' 
#' @importFrom stats optimize 
#' @export              
# function dew_point(Tair,VPD,accuracy=1e-03,Esat_formula="Sonntag_1990",
#                       constants=bigleaf_constants()){
  
#   if (!is_numeric(accuracy)){
#     stop("'accuracy' must be numeric!")
#   }
  
#   if (accuracy > 1){
#     print("'accuracy' is set to 1 degC")
#     accuracy = 1
#   }
  
#   # determine number of digits to print
#   ndigits = as_numeric(strsplit(format(accuracy,scientific = TRUE),"-")[[1]][2])
#   ndigits = ifelse(is_na(ndigits),0,ndigits)
  
#   ea = VPD_to_e(VPD,Tair,Esat_formula)
#   Td = sapply(seq_along(ea),function(i) round(dew_point_solver(ea[i],accuracy=accuracy,
#                                                                 Esat_formula,constants)$minimum,ndigits))
  
#   return(Td)
# }







#' Virtual Temperature
#' 
#' Virtual temperature, defined as the temperature at which dry air would have the same
#'              density as moist air at its actual temperature_
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param VPD       Vapor pressure deficit (kPa)
#' @param Esat_formula  Optional: formula to be used for the calculation of esat and the slope of esat_ 
#'                      One of \code{"Sonntag_1990"} (Default), \code{"Alduchov_1996"}, or \code{"Allen_1998"}_
#'                      See \code{\link{Esat_slope}}_ 
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  eps - ratio of the molecular weight of water vapor to dry air (-) 
#' 
#' @details the virtual temperature is given by:
#'  
#'    \deqn{Tv = Tair / (1 - (1 - eps) e/pressure)}
#' 
#'  where Tair is in Kelvin (converted internally)_ Likewise, VPD is converted 
#'  to actual vapor pressure (e in kPa) with \code{\link{VPD_to_e}} internally_
#' 
#' @return \item{Tv -}{virtual temperature (deg C)}
#' 
#' @references Monteith J_L_, Unsworth M_H_, 2008: Principles of Environmental Physics_
#'             3rd edition_ Academic Press, London_
#'  
#' @examples 
#' virtual_temp(25,100,1.5)                        
#'               
#' @export
function virtual_temp(Tair,pressure,VPD;Esat_formula="Sonntag_1990",
                         constants=bigleaf_constants())
  e    = VPD_to_e(VPD,Tair,Esat_formula)
  Tair = Tair + constants[:Kelvin]
  Tv = Tair / (1 - (1 - constants[:eps]) * e/pressure) 
  Tv = Tv - constants[:Kelvin]
  return(Tv)
end



#' Kinematic Viscosity of Air
#' 
#' calculates the kinematic viscosity of air_
#' 
#' @param Tair      Air temperature (deg C)
#' @param pressure  Atmospheric pressure (kPa)
#' @param constants Kelvin - conversion degree Celsius to Kelvin \cr
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) \cr
#'                  Tair0 - reference air temperature (K) \cr
#'                  kPa2Pa - conversion kilopascal (kPa) to pascal (Pa)
#' 
#' @details where v is the kinematic viscosity of the air (m2 s-1), 
#'          given by (Massman 1999b):
#'          
#'            \deqn{v = 1.327 * 10^-5(pressure0/pressure)(Tair/Tair0)^1.81}
#'          
#' @return \item{v -}{kinematic viscosity of air (m2 s-1)}
#' 
#' @references Massman, W_J_, 1999b: Molecular diffusivities of Hg vapor in air, 
#'             O2 and N2 near STP and the kinematic viscosity and thermal diffusivity
#'             of air near STP_ Atmospheric Environment 33, 453-457_      
#'             
#' @examples 
#' kinematic_viscosity(25,100)    
#' 
#' @export         
function kinematic_viscosity(Tair,pressure; constants=bigleaf_constants())
  
  Tair     = Tair + constants[:Kelvin]
  pressure = pressure * constants[:kPa2Pa]
  
  v  = 1.327e-05*(constants[:pressure0]/pressure)*(Tair/constants[:Tair0])^1.81
  return(v)
end
