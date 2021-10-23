#######################################
#### Canopy-Atmosphere Decoupling  ####
#######################################

#' Canopy-Atmosphere Decoupling Coefficient
#' 
#' The canopy-atmosphere decoupling coefficient 'Omega'. 
#' 
#' - data        Data_frame or matrix containing all required input variables
#' - Tair        Air temperature (deg C)
#' - pressure    Atmospheric pressure (kPa)
#' - Ga          Aerodynamic conductance to heat/water vapor (m s-1)
#' - Gs          Surface conductance (m s-1)
#' - approach    Approach used to calculate omega. Either `"Jarvis&McNaughton_1986"` (default)
#'                    or `"Martin_1989"`.
#' - LAI         Leaf area index (m2 m-2), only used if `approach = "Martin_1989"`.
#' - Esat_formula  Optional: formula to be used for the calculation of esat and the slope of esat.
#'                      One of `"Sonntag_1990"` (Default), `"Alduchov_1996"`, or `"Allen_1998"`.
#'                      See [`Esat_slope`](@ref). 
#' - constants   Kelvin - conversion degree Celsius to Kelvin 
#'                    cp - specific heat of air for constant pressure (J K-1 kg-1) 
#'                    eps - ratio of the molecular weight of water vapor to dry air (-) 
#'                    sigma - Stefan-Boltzmann constant (W m-2 K-4) 
#'                    Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' # Details
 The decoupling coefficient Omega ranges from 0 to 1 and quantifies the
#'          linkage of the conditions (foremost humidity and temperature) at the canopy surface
#'          to the ambient air. Values close to 0 indicate well coupled conditions
#'          characterized by high physiological (i.e. stomatal) control on transpiration
#'          and similar conditions at the canopy surface compared to the atmosphere above
#'          the canopy. Values close to 1 indicate the opposite, i.e. decoupled conditions and 
#'          a low stomatal control on transpiration (Jarvis & McNaughton 1986). 
#'          The `"Jarvis&McNaughton_1986"` approach (default option) is the original
#'          formulation for the decoupling coefficient, given by (for an amphistomatous 
#'          canopy):
#'          
#'          ``\\Omega = \frac{\epsilon + 1``{\epsilon + 1 + \frac{Ga}{Gc}}}{%
#'          \\Omega = (\epsilon + 1) / ( \epsilon + 1 + Ga/Gc)}
#'          
#'          where ``\epsilon = \frac{s``{\\gamma}}{\epsilon = s/\\gamma} is a dimensionless coefficient
#'          with s being the slope of the saturation vapor pressure curve (Pa K-1), and ``\\gamma`` the 
#'          psychrometric constant (Pa K-1).
#'          
#'          The approach `"Martin_1989"` by Martin 1989 additionally takes radiative coupling
#'          into account:
#'          
#'          ``\\Omega = \frac{\epsilon + 1 + \frac{Gr``{Ga}}{\epsilon + (1 + \frac{Ga}{Gs}) (1 + \frac{Gr}{Ga})}}{%
#'          \\Omega = (\epsilon + 1 + Gr/Ga) / (\epsilon + (1 + Ga/Gs) (1 + Gr/Ga))}
#' 
#' # Value
 - ``\\Omega`` -: the decoupling coefficient Omega (-)
#' 
#' #References
#' Jarvis P_G., McNaughton K_G., 1986: Stomatal control of transpiration:
#'             scaling up from leaf to region. Advances in Ecological Research 15, 1-49. 
#'             
#'             Martin P., 1989: The significance of radiative coupling between
#'             vegetation and the atmosphere. Agricultural and Forest Meteorology 49, 45-53.
#' 
#' #See also
#' [`aerodynamic_conductance`](@ref), [`surface_conductance`](@ref),
#'          [`equilibrium_imposed_ET`](@ref)
#' 
#' ```@example; output = false
#' ``` 
#' # Omega calculated following Jarvis & McNaughton 1986
#' set_seed(3)
#' df = DataFrame(Tair=rnorm(20,25,1),pressure=100,Ga_h=rnorm(20,0.06,0.01),
#'                  Gs_ms=rnorm(20,0.005,0.001))
#' decoupling(df,approach="Jarvis&McNaughton_1986")
#' 
#' # Omega calculated following Martin 1989 (requires LAI)
#' decoupling(df,approach="Martin_1989",LAI=4)
#' 
"""
"""
function decoupling(data,Tair="Tair",pressure="pressure",Ga="Ga_h",Gs="Gs_ms",
                       approach=c("Jarvis&McNaughton_1986","Martin_1989"),
                       LAI,Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                       constants=bigleaf_constants())
  
  approach    = match_arg(approach)
  
  check_input(data,list(Tair,pressure,Ga,Gs))
  
  Delta   = Esat_slope(Tair,Esat_formula,constants)[,"Delta"]
  gamma   = psychrometric_constant(Tair,pressure,constants)
  epsilon = Delta/gamma
  
  if (approach == "Jarvis&McNaughton_1986")
    
    Omega = (epsilon + 1) / (epsilon + 1 + Ga/Gs)
    
else if (approach == "Martin_1989") 
    
    if (is_null(LAI))
      
      stop("LAI is not provided!")
      
else 
      
      Gr    = longwave_conductance(Tair,LAI,constants)
      Omega = (epsilon + 1 + Gr/Ga) / (epsilon + 1 + Ga/Gs + Gr/Gs + Gr/Ga)
      
end
end
  
  return(Omega)
  
end



#' Longwave Radiative Transfer Conductance of the Canopy
#' 
#' - Tair      Air temperature (deg C)
#' - LAI       Leaf area index (m2 m-2)
#' - constants Kelvin - conversion degree Celsius to Kelvin 
#'                  sigma - Stefan-Boltzmann constant (W m-2 K-4) 
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1)
#'                  
#' # Details
 the following formula is used (Martin, 1989):
#' 
#'          ``Gr = 4 \sigma Tair^3 LAI / cp``                             
#'                                       
#' # Value
 - Gr -: longwave radiative transfer conductance of the canopy (m s-1)
#'                  
#' #References
#' Martin P., 1989: The significance of radiative coupling between
#'             vegetation and the atmosphere. Agricultural and Forest Meteorology 49, 45-53.
#'          
#' ```@example; output = false
#' ``` 
#' longwave_conductance(25,seq(1,8,1))            
#'                  
#' @export             
function longwave_conductance(Tair,LAI,constants=bigleaf_constants())
  
  Tair = Tair + constants[:Kelvin]
  
  Gr = 4 * constants[:sigma] * Tair^3 * LAI / constants[:cp]
  
  return(Gr)
end
