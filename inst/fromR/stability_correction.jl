#####################################################
### Stability parameters and stability correction ### ---------------------------------------------------
#####################################################

#' Monin-Obukhov Length
#' 
#' calculates the Monin-Obukhov length.
#' 
#' - data      DataFrame or matrix containing all required variables
#' - Tair      Air temperature (deg C)
#' - pressure  Atmospheric pressure (kPa)
#' - ustar     Friction velocity (m s-1)
#' - H         Sensible heat flux (W m-2)
#' - constants Kelvin - conversion degree Celsius to Kelvin 
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1) 
#'                  k - von Karman constant (-) 
#'                  g - gravitational acceleration (m s-2)
#' 
#' # Details
 The Monin-Obukhov length (L) is given by:
#' 
#'              ``L = - (\\rho * cp * ustar^3 * Tair) / (k * g * H)``
#'              
#'              where ``rho`` is air density (kg m-3).
#' 
#' # Value
 - L -: Monin-Obukhov length (m)
#' 
#' # Note
#' Note that L gets very small for very low ustar values with implications
#'       for subsequent functions using L as input. It is recommended to filter
#'       data and exclude low ustar values (ustar < ~0.2) beforehand. 
#' 
#' #References
#' Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany. 
#' 
#' #See also
#' [`stability_parameter`](@ref)
#' 
#' ```@example; output = false
#' ``` 
#' Monin_Obukhov_length(Tair=25,pressure=100,ustar=seq(0.2,1,0.1),H=seq(40,200,20))
#' 
"""
"""
function Monin_Obukhov_length(data,Tair="Tair",pressure="pressure",ustar="ustar",
                                 H="H",constants=BigleafConstants())
  
  check_input(data,list(Tair,pressure,ustar,H))
  
  rho  = air_density(Tair,pressure,constants)
  Tair = Tair + constants.Kelvin
  MOL  = (-rho*constants.cp]*ustar^3*Tair) / (constants[:k]*constants[:g*H)
  
  return(MOL)
end



#' Stability Parameter "zeta"
#' 
#' calculates "zeta", a parameter characterizing stratification in 
#'              the lower atmosphere.
#' 
#' - data      DataFrame or matrix containing all required variables
#' - Tair      Air temperature (degC)
#' - pressure  Atmospheric pressure (kPa)
#' - ustar     Friction velocity (m s-1)
#' - H         Sensible heat flux (W m-2)
#' - zr        Instrument (reference) height (m)
#' - d         Zero-plane displacement height (m)
#' - constants Kelvin - conversion degree Celsius to Kelvin 
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1) 
#'                  k - von Karman constant (-) 
#'                  g - gravitational acceleration (m s-2)
#' 
#' # Details
 The stability parameter ``\\zeta`` is given by:
#' 
#'            ``\\zeta = (zr - d) / L``
#'          
#'          where L is the Monin-Obukhov length (m), calculated from the function
#'          [`Monin_Obukhov_length`](@ref). The displacement height d can 
#'          be estimated from the function [`roughness_parameters`](@ref).
#'          
#' # Value
 - ``\\zeta`` - : stability parameter (-)
#' 
#' ```@example; output = false
#' ``` 
#' df = DataFrame(Tair=25,pressure=100,ustar=seq(0.2,1,0.1),H=seq(40,200,20))
#' stability_parameter(df,zr=40,d=15)
#' 
#' @export           
function stability_parameter(data,Tair="Tair",pressure="pressure",ustar="ustar",
                                H="H",zr,d,constants=BigleafConstants())
  
  check_input(data,list(Tair,pressure,ustar,H))
  
  MOL  = Monin_Obukhov_length(data,Tair,pressure,ustar,H,constants)
  zeta = (zr - d) / MOL
  
  return(zeta)
  
end



#' Integrated Stability Correction Functions for Heat and Momentum
#' 
#' dimensionless stability functions needed to correct deviations
#'              from the exponential wind profile under non-neutral conditions.
#'              
#' - zeta         Stability parameter zeta (-)
#' - formulation  Formulation for the stability function. Either `Val(:Dyer_1970)`, 
#'                     or `Val(:Businger_1971)`
#'
#' # Details
 The functions give the integrated form of the universal functions. They
#'          depend on the value of the stability parameter ``\\zeta``,
#'          which can be calculated from the function [`stability_parameter`](@ref).
#'          The integration of the universal functions is:
#'          
#'            ``\\psi = -x * zeta`` 
#'          
#'          for stable atmospheric conditions (``\\zeta`` >= 0), and
#'          
#'            ``\\psi = 2 * log( (1 + y) / 2) ``
#'          
#'          for unstable atmospheric conditions (``\\zeta`` < 0).
#'          
#'          The different formulations differ in their value of x and y.
#'   
#' # Value
 a DataFrame with the following columns:
#'          - psi_h: the value of the stability function for heat and water vapor (-)
#'          - psi_m: the value of the stability function for momentum (-)
#' 
#' #References
#' Dyer, A_J., 1974: A review of flux-profile relationships. 
#'             Boundary-Layer Meteorology 7, 363-372.
#'             
#'             Dyer, A. J., Hicks, B_B., 1970: Flux-Gradient relationships in the
#'             constant flux layer. Quart. J. R. Meteorol. Soc. 96, 715-721.
#'             
#'             Businger, J_A., Wyngaard, J. C., Izumi, I., Bradley, E. F., 1971:
#'             Flux-Profile relationships in the atmospheric surface layer. 
#'             J. Atmospheric Sci. 28, 181-189.
#'             
#'             Paulson, C_A., 1970: The mathematical representation of wind speed
#'             and temperature profiles in the unstable atmospheric surface layer.
#'             Journal of Applied Meteorology 9, 857-861.
#' 
#'             Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany.
#'
#' ```@example; output = false
#' ``` 
#' zeta = seq(-2,0.5,0.05)
#' stability_correction(zeta)
#' stability_correction(zeta,formulation=Val(:Businger_1971))                          
#'             
#' @export   
function stability_correction(zeta,formulation=c(Val(:Dyer_1970),Val(:Businger_1971)))
  
  formulation  = match_arg(formulation)
  
  check_input(nothing,zeta)
  
  psi_h = psi_m = numeric()
  
  # universal functions
  if (formulation == Val(:Businger_1971))
    x_h = -7.8
    x_m = -6
    y_h = 0.95 * ( 1 - 11.6 * zeta)^0.5
    y_m = (1 - 19.3*zeta)^0.25
elseif (formulation == Val(:Dyer_1970))
    x_h = x_m = -5
    y_h       = (1 - 16 * zeta)^0.5
    y_m       = (1 - 16 * zeta)^0.25
end
  
  # integration of universal functions (after Paulson_1970 and Foken 2008)
  # stable
  stable = zeta >= 0 | ismissing(zeta)
  psi_h[stable] = x_h * zeta[stable]
  psi_m[stable] = x_m * zeta[stable]
  # unstable
  unstable = zeta < 0 | ismissing(zeta)
  psi_h[unstable] = 2 * log( (1 + y_h[unstable] ) / 2)
  psi_m[unstable] = 2 * log( (1 + y_m[unstable] ) / 2) +
                     log( ( 1 + y_m[unstable]^2 ) / 2)
                     -2 * atan(y_m[unstable]) + pi/2
  
  return(DataFrame(psi_h,psi_m))
  
end 
