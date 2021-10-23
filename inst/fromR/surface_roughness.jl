###########################
#### Surface roughness ####
###########################

#' Roughness Reynolds Number
#' 
#' calculates the Roughness Reynolds Number.
#' 
#' - Tair      Air temperature (deg C)
#' - pressure  Atmospheric pressure (kPa)
#' - ustar     Friction velocity (m s-1)
#' - z0m       Roughness length (m)
#' - constants Kelvin - conversion degree Celsius to Kelvin 
#'                  pressure0 - reference atmospheric pressure at sea level (Pa) 
#'                  Tair0 - reference air temperature (K)
#'                  
#' # Details
 The Roughness Reynolds Number is calculated as in Massman 1999a:
#'          
#'            ``Re = z0m * ustar / v``
#'          
#'          where `v` is the kinematic viscosity (m2 s-1).
#'          
#' # Value
 - Re -: Roughness Reynolds Number (-)
#' 
#' #References
#' Massman, W_J., 1999a: A model study of kB H- 1 for vegetated surfaces using
#'            'localized near-field' Lagrangian theory. Journal of Hydrology 223, 27-43.
#' 
#' ```@example; output = false
#' ``` 
#' Reynolds_Number(25,100,0.5,z0m=0.5)                             
#' 
"""
"""
function Reynolds_Number(Tair,pressure,ustar,z0m,constants=bigleaf_constants())
  
  v  = kinematic_viscosity(Tair,pressure,constants)
  Re = z0m*ustar/v
  
  return(Re)
end



#' Roughness Parameters
#' 
#' A simple approximation of the two roughness parameters displacement height (d)
#'              and roughness length for momentum (z0m).
#'              
#' - method    Method to use, one of `"canopy_height","canopy_height&LAI","wind_profile"` 
#'                  NOTE: if `method = "canopy_height"`, only the following three arguments
#'                  are used. If `method = "canopy_height&LAI"`, only `zh, LAI, cd`, 
#'                  and `hs` are required.     
#' - zh        Vegetation height (m)          
#' - frac_d    Fraction of displacement height on canopy height (-)
#' - frac_z0m  Fraction of roughness length on canopy height (-)
#' - LAI       Leaf area index (-) 
#' - zr        Instrument (reference) height (m)
#' - cd        Mean drag coefficient for individual leaves. Defaults to 0.2. 
#'                  Only needed if `method = "canopy_height&LAI"`.
#' - hs        roughness length of the soil surface (m). Only needed if `method = "canopy_height&LAI"`
#'                  The following arguments are only needed if `method = "wind_profile"`!
#' - data      Data_frame or matrix containing all required variables
#' - Tair      Air temperature (deg C)
#' - pressure  Atmospheric pressure (kPa)
#' - wind      Wind speed at height zr (m s-1)
#' - ustar     Friction velocity (m s-1)
#' - H         Sensible heat flux (W m-2)
#' - d         Zero-plane displacement height (m); optional
#' - z0m       Roughness length for momentum (m); optional
#' - stab_roughness   Should stability correction be considered? Default is `TRUE`.
#' - stab_formulation Stability correction function used (If `stab_correction = TRUE`).
#'                         Either `"Dyer_1970"` or `"Businger_1971"`.
#' - constants k - von-Karman constant (-) 
#'                  Kelvin - conversion degree Celsius to Kelvin 
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1) 
#'                  g - gravitational acceleration (m s-2) 
#'                  se_median - conversion standard error (SE) of the mean to SE of the median
#'                  
#' 
#' # Details
 The two main roughness parameters, the displacement height (d)
#'          and the roughness length for momentum (z0m) can be estimated from simple
#'          empirical relationships with canopy height (zh). If `method = "canopy_height"`,
#'          the following formulas are used:  
#'          
#'            ``d = frac_d * zh``
#'          
#'            ``z0m = frac_z0m * zh``
#'          
#'          where frac_d defaults to 0.7 and frac_z0m to 0.1.
#'          
#'          Alternatively, d and z0m can be estimated from both canopy height and LAI
#'          (If `method = "canopy_height&LAI"`).
#'          Based on data from Shaw & Pereira 1982, Choudhury & Monteith 1988 proposed 
#'          the following semi-empirical relations:
#'          
#'            ``X = cd * LAI``
#'          
#'            ``d = 1.1 * zh * ln(1 + X^(1/4))`` 
#'          
#'            ``z0m = hs + 0.3 * zh * X^(1/2)   for 0 <= X <= 0.2``
#'          
#'            ``z0m = hs * zh * (1 - d/zh)   for 0.2 < X`` 
#'          
#'          If `method = "wind_profile"`, z0m is estimated by solving
#'          the wind speed profile for z0m:
#'          
#'            ``z0m = median((zr - d) * exp(-k*wind / ustar - psi_m)``
#'                  
#'          By default, d in this equation is fixed to 0.7*zh, but can be set to any
#'          other value. psi_m is 0 if `stab_roughness = false`.       
#' 
#' # Value
 a DataFrame with the following columns:
#'         - d: Zero-plane displacement height (m)
#'         - z0m: Roughness length for momentum (m)
#'         - z0m_se: Only if `method = wind_profile`: Standard Error of the median for z0m (m)
#'
#' #References
#' Choudhury, B. J., Monteith J_L., 1988: A four-layer model for the heat
#'             budget of homogeneous land surfaces. Q. J. R. Meteorol. Soc. 114, 373-398.
#'             
#'             Shaw, R. H., Pereira, A., 1982: Aerodynamic roughness of a plant canopy: 
#'             a numerical experiment. Agricultural Meteorology, 26, 51-65.
#'    
#' #See also
#' [`wind_profile`](@ref)
#'     
#' ```@example; output = false
#' ``` 
#' # estimate d and z0m from canopy height for a dense (LAI=5) and open (LAI=2) canopy
#' roughness_parameters(method="canopy_height&LAI",zh=25,LAI=5)
#' roughness_parameters(method="canopy_height&LAI",zh=25,LAI=2)   
#'    
#' # fix d to 0.7*zh and estimate z0m from the wind profile
#' df = DataFrame(Tair=c(25,25,25),pressure=100,wind=c(3,4,5),ustar=c(0.5,0.6,0.65),H=200)
#' roughness_parameters(method="wind_profile",zh=25,zr=40,frac_d=0.7,data=df)
#' 
#' # assume d = 0.8*zh
#' roughness_parameters(method="wind_profile",zh=25,zr=40,frac_d=0.8,data=df) 
#' 
#' @importFrom stats median sd complete_cases 
#' @export                                  
function roughness_parameters(method=c("canopy_height","canopy_height&LAI","wind_profile"),zh,
                                 frac_d=0.7,frac_z0m=0.1,LAI,zr,cd=0.2,hs=0.01,data,Tair="Tair",pressure="pressure",
                                 wind="wind",ustar="ustar",H="H",d=NULL,z0m=NULL,
                                 stab_roughness=TRUE,stab_formulation=c("Dyer_1970","Businger_1971"),
                                 constants=bigleaf_constants())
  
  method           = match_arg(method)
  stab_formulation = match_arg(stab_formulation)
  
  if (method == "canopy_height")
    
    d      = frac_d*zh
    z0m    = frac_z0m*zh
    z0m_se = NA
    
else if (method == "canopy_height&LAI")
    
    X = cd * LAI
    d = 1.1 * zh * log(1 + X^(1/4))
    
    if (X >= 0 & X <= 0.2)
      z0m = hs + 0.3 * X^(1/2)
else 
      z0m = 0.3 * zh * (1 - d/zh)
end
    z0m_se = NA
    
else if (method == "wind_profile")
    
    check_input(data,Tair,pressure,wind,ustar,H)
    
    if (is_null(d))
      
      d = frac_d * zh
      
end
    
    if (stab_roughness)
      
      zeta  = stability_parameter(data=data,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                                   zr=zr,d=d,constants=constants)
      psi_m = stability_correction(zeta,formulation=stab_formulation)[,"psi_m"]
      z0m_all = (zr - d) * exp(-constants[:k]*wind / ustar - psi_m)
      
else 
      
      z0m_all = (zr - d) * exp(-constants[:k]*wind / ustar)
      
end
    
    z0m_all[z0m_all > zh] = NA
    
    z0m    = median(z0m_all,na_rm=TRUE)
    z0m_se = constants[:se_median] * (sd(z0m_all,na_rm=TRUE) / sqrt(length(z0m_all[complete_cases(z0m_all)])))
    
end
  
  return(DataFrame(d,z0m,z0m_se))
end



#' Wind Speed at Given Heights in the Surface Layer
#' 
#' Wind speed at a given height above the canopy estimated from single-level
#'              measurements of wind speed.
#'          
#' - data      Data_frame or matrix containing all required variables
#' - z         Heights above ground for which wind speed is calculated.
#'                  Needs to be of same length as `data` or of length 1
#' - Tair      Air temperature (deg C)
#' - pressure  Atmospheric pressure (kPa)                                                                                  
#' - ustar     Friction velocity (m s-1)
#' - H         Sensible heat flux (W m-2)
#' - wind      Wind speed at height zr (m s-1); only used if `stab_correction = TRUE`
#' - zr        Instrument (reference) height (m)
#' - zh        Canopy height (m)
#' - d         Zero-plane displacement height (-)
#' - frac_d    Fraction of displacement height on canopy height (-);
#'                  only used if `d` is not available
#' - z0m       Roughness length (m), optional; only used if `stab_correction = false` (default=0.1) 
#' - frac_z0m  Fraction of roughness length on canopy height (-), optional; only used if `z0m` is not provided.
#'                  Default is 0.1.
#' - estimate_z0m Should `z0m` be estimated from the logarithmic wind profile? If `TRUE` (the default),
#'                     arguments `z0m` and `frac_z0m` are ignored.
#'                     See [`roughness_parameters`](@ref) for details. 
#' - stab_correction Should stability correction be applied? Defaults to `TRUE`
#' - stab_formulation Stability correction function used (If `stab_correction = TRUE`).
#'                         Either `"Dyer_1970"` or `"Businger_1971"`.
#' - constants k - von-Karman constant (-) 
#'                  Kelvin - conversion degree Celsius to Kelvin 
#'                  cp - specific heat of air for constant pressure (J K-1 kg-1) 
#'                  g - gravitational acceleration (m s-2)
#'
#' # Details
 The underlying assumption is the existence of a logarithmic wind profile
#'          above the height d + z0m (the height at which wind speed mathematically reaches zero
#'          according to the Monin-Obhukov similarity theory).
#'          In this case, the wind speed at a given height z is given by:
#'          
#'            ``u(z) = (ustar/k) * (ln((z - d) / z0m) - \psi{m``}
#'          
#'          The roughness parameters zero-plane displacement height (d) and roughness length (z0m)
#'          can be approximated from [`roughness_parameters`](@ref). ``\psi{m``} is omitted
#'          if `stab_correction = false` (not recommended). If `estimate_z0m = TRUE`,
#'          z0m is first estimated from the wind profile equation and then used in the equation
#'          above for the calculation of `u(z)` (see e.g. Newman & Klein 2014).        
#'                                                             
#' #Note
#' Note that this equation is only valid for z >= d + z0m, and it is not 
#'       meaningful to calculate values closely above d + z0m. All values in `heights`
#'       smaller than d + z0m will return 0.                                 
#'                                  
#' # Value
 A vector of wind speed at heights `z`. 
#'         
#' #References
#' Monteith, J_L., Unsworth, M_H., 2008: Principles of Environmental Physics.
#'             3rd edition. Academic Press, London. 
#'             
#'             Newman, J_F., Klein, P_M., 2014: The impacts of atmospheric stability on
#'             the accuracy of wind speed extrapolation methods. Resources 3, 81-105.
#'         
#' #See also
#' [`roughness_parameters`](@ref)
#' 
#' ```@example; output = false
#' ``` 
#' heights = seq(18,40,2)  # heights above ground for which to calculate wind speed
#' df = DataFrame(Tair=25,pressure=100,wind=c(3,4,5),ustar=c(0.5,0.6,0.65),H=c(200,230,250)) 
#' ws = DataFrame(matrix(NA,ncol=length(heights),nrow=nrow(df)))
#' colnames(ws) = paste0(heights,"m")
#' for (i in seq_along(heights))
#'   ws[,i] = wind_profile(df,z=heights[i],zr=40,zh=25,d=16)
#' }
#' 
#' @export                                                                                                                          
function wind_profile(data,z,Tair="Tair",pressure="pressure",ustar="ustar",H="H",wind="wind",
                         zr,zh,d=NULL,frac_d=0.7,z0m=NULL,frac_z0m=NULL,estimate_z0m=TRUE,
                         stab_correction=TRUE,stab_formulation=c("Dyer_1970","Businger_1971"),
                         constants=bigleaf_constants())
  
  stab_formulation = match_arg(stab_formulation)
  
  check_input(data,ustar)
  
  ## determine roughness parameters
  if (is_null(d))
    if (is_null(frac_d))
      stop("Either 'd' or 'frac_d' must be specified")
end
    d = frac_d * zh
end
  
  if (is_null(z0m) & !estimate_z0m)
    if (is_null(frac_z0m))
      stop("Either 'z0m' or 'frac_z0m' must be specified if 'estimate_z0m' = false")
end
    z0m = frac_z0m * zh
end
  
  
  if (estimate_z0m)
    
    if (!is_null(z0m) | !is_null(frac_z0m))
      cat("Note that arguments 'z0m' and 'frac_z0m' are ignored if 'estimate_z0m' = TRUE. z0m is
           calculated from the logarithmic wind_profile equation.",fill=TRUE)
end
    
    check_input(data,Tair,pressure,wind,ustar,H)
    
    z0m = roughness_parameters(method="wind_profile",zh=zh,zr=zr,d=d,data=data,
                                Tair=Tair,pressure=pressure,wind=wind,ustar=ustar,H=H,
                                stab_roughness=TRUE,stab_formulation=stab_formulation,
                                constants=constants)[,"z0m"]
end
  
  if ( any(z < (d + z0m) & !is_na(d + z0m)) )
    warning("function is only valid for heights above d + z0m! Wind speed for heights below d + z0m will return 0!") 
end
  
  ## calculate wind speeds at given heights z
  if (stab_correction)
      
    zeta  = stability_parameter(data=data,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                                 zr=z,d=d,constants=constants)
    psi_m = stability_correction(zeta,formulation=stab_formulation)[,"psi_m"]
    wind_heights = pmax(0,(ustar / constants[:k]) * (log(pmax(0,(z - d)) / z0m) - psi_m))
      
else 
      
    wind_heights = pmax(0,(ustar / constants[:k]) * (log(pmax(0,(z - d)) / z0m)))
      
end
  
  return(wind_heights)
end

