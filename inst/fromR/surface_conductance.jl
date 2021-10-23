##############################
#### Surface conductance  ####
##############################

#' Surface Conductance to Water Vapor
#' 
#' Calculates surface conductance to water vapor from the inverted Penman-Monteith
#'              equation (by default) or from a simple flux-gradient approach.
#' 
#' - data      Data_frame or matrix containing all required input variables
#' - Tair      Air temperature (deg C)
#' - pressure  Atmospheric pressure (kPa)
#' - Rn        Net radiation (W m-2)
#' - G         Ground heat flux (W m-2); optional
#' - S         Sum of all storage fluxes (W m-2); optional
#' - LE        Latent heat flux (W m-2)
#' - VPD       Vapor pressure deficit (kPa)
#' - Ga        Aerodynamic conductance to heat/water vapor (m s-1)
#' - missing_G_as_NA  if `TRUE`, missing G are treated as `NA`s, otherwise they are set to 0.
#'                         Only used if `formulation = "Penman-Monteith"`.
#' - missing_S_as_NA  if `TRUE`, missing S are treated as `NA`s, otherwise they are set to 0. 
#'                          Only used if `formulation = "Penman-Monteith"`.
#' - formulation Formulation used. Either `"Penman-Monteith"` (the default) 
#'                    using the inverted Penman-Monteith equation, or `"Flux-Gradient"`,
#'                    for a simple flux-gradient approach requiring ET, pressure, and VPD only.
#' - Esat_formula  Optional: formula to be used for the calculation of esat and the slope of esat.
#'                      One of `"Sonntag_1990"` (Default), `"Alduchov_1996"`, or `"Allen_1998"`. 
#'                      Only used if `formulation = "Penman-Monteith"`. See `\link{Esat_slope`}.
#' - constants   cp - specific heat of air for constant pressure (J K-1 kg-1) \cr
#'                    eps - ratio of the molecular weight of water vapor to dry air (-) \cr
#'                    Rd - gas constant of dry air (J kg-1 K-1) \cr
#'                    Rgas - universal gas constant (J mol-1 K-1) \cr
#'                    Kelvin - conversion degree Celsius to Kelvin \cr
#'                    Mw - molar mass of water vapor (kg mol-1) \cr
#'                    Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' 
#' # Details
 If `formulation = "Penman-Monteith"` (the default), surface conductance (Gs) in m s-1 
#'          is calculated from the inverted Penman-Monteith equation:
#' 
#'     \deqn{Gs = ( LE * Ga * \gamma ) / ( \Delta * A + \rho * cp * Ga * VPD - LE * ( \Delta + \gamma ) )}
#'  
#'  Where \eqn{\gamma} is the psychrometric constant (kPa K-1), \eqn{\Delta} is the slope of the 
#'  saturation vapor pressure curve (kPa K-1), and \eqn{\rho} is air density (kg m-3).
#'  Available energy (A) is defined as A = Rn - G - S. If G and/or S are not provided, A = Rn.
#'  
#'  By default, any missing data in G and S are set to 0. If `missing_S_as_NA = TRUE`
#'  or `missing_S_as_NA = TRUE`, Gs will give `NA` for these timesteps.
#'  
#'  If `formulation="Flux-Gradient"`, Gs (in mol m-2 s-1) is calculated from VPD and ET only:
#'  
#'     \deqn{Gs = ET/pressure * VPD}
#'  
#'  where ET is in mol m-2 s-1. Note that this formulation assumes fully coupled conditions (i.e. Ga = inf).
#'  This formulation is equivalent to the inverted form of Eq.6 in McNaughton & Black 1973:
#'  
#'     \deqn{Gs = LE * \gamma / (\rho * cp * VPD)}
#'  
#'  which gives Gs in m s-1. Note that Gs > Gc (canopy conductance) under conditions 
#'  when a significant fraction of ET comes from interception or soil evaporation. 
#'  
#'  If `pressure` is not available, it can be approximated by elevation using the 
#'  function `\link{pressure_from_elevation`}
#'  
#' # Value
 a dataframe with the following columns: 
#'  \item{Gs_ms}{Surface conductance in m s-1}
#'  \item{Gs_mol}{Surface conductance in mol m-2 s-1}
#' 
#' ```@example; output = false
#' ``` 
#' ## filter data to ensure that Gs is a meaningful proxy to canopy conductance (Gc)
#' DE_Tha_Jun_2014_2 = filter_data(DE_Tha_Jun_2014,quality_control=FALSE,
#'                                  vars_qc=c("Tair","precip","VPD","H","LE"),
#'                                  filter_growseas=FALSE,filter_precip=TRUE,
#'                                  filter_vars=c("Tair","PPFD","ustar","LE"),
#'                                  filter_vals_min=c(5,200,0.2,0),
#'                                  filter_vals_max=c(NA,NA,NA,NA),NA_as_invalid=TRUE,
#'                                  quality_ext="_qc",good_quality=c(0,1),
#'                                  missing_qc_as_bad=TRUE,GPP="GPP",doy="doy",
#'                                  year="year",tGPP=0.5,ws=15,min_int=5,precip="precip",
#'                                  tprecip=0.1,precip_hours=24,records_per_hour=2)
#' 
#' # calculate Gs based on a simple gradient approach
#' Gs_gradient = surface_conductance(DE_Tha_Jun_2014_2,Tair="Tair",pressure="pressure",
#'                                    VPD="VPD",formulation="Flux-Gradient")
#' summary(Gs_gradient)
#' 
#' # calculate Gs from the the inverted PM equation (now Rn, and Ga are needed),
#' # using a simple estimate of Ga based on Thom 1972
#' Ga = aerodynamic_conductance(DE_Tha_Jun_2014_2,Rb_model="Thom_1972")[,"Ga_h"]
#' 
#' # if G and/or S are available, don't forget to indicate (they are ignored by default).
#' # Note that Ga is not added to the DataFrame 'DE_Tha_Jun_2014'
#' Gs_PM = surface_conductance(DE_Tha_Jun_2014_2,Tair="Tair",pressure="pressure",
#'                              Rn="Rn",G="G",S=NULL,VPD="VPD",Ga=Ga,
#'                              formulation="Penman-Monteith")
#' summary(Gs_PM)
#' 
#'                               
#' # now add Ga to the DataFrame 'DE_Tha_Jun_2014' and repeat
#' DE_Tha_Jun_2014_2$Ga = Ga
#' Gs_PM2 = surface_conductance(DE_Tha_Jun_2014_2,Tair="Tair",pressure="pressure",
#'                               Rn="Rn",G="G",S=NULL,VPD="VPD",Ga="Ga",
#'                               formulation="Penman-Monteith")
#' # note the difference to the previous version (Ga="Ga")
#' summary(Gs_PM2)
#' 
#' @references Monteith, J., 1965: Evaporation and environment. In Fogg, G. E. (Ed.),
#'             The state and movement of water in living organisms (pp.205-234).
#'             19th Symp. Soc. Exp. Biol., Cambridge University Press, Cambridge
#' 
#'             McNaughton, K_G., Black, T_A., 1973: A study of evapotranspiration
#'             from a Douglas Fir forest using the energy balance approach. Water
#'             Resources Research 9, 1579-1590.
#'             
"""
"""
function surface_conductance(data,Tair="Tair",pressure="pressure",Rn="Rn",G=NULL,S=NULL,
                                VPD="VPD",LE="LE",Ga="Ga_h",missing_G_as_NA=FALSE,missing_S_as_NA=FALSE,
                                formulation=c("Penman-Monteith","Flux-Gradient"),
                                Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                                constants=bigleaf_constants())
  
  formulation = match_arg(formulation)
  
  if (formulation == "Flux-Gradient")
  
    check_input(data,list(Tair,pressure,VPD,LE))
    
    Gs_mol = (LE_to_ET(LE,Tair)/constants[:Mw]) * pressure / VPD
    Gs_ms  = mol_to_ms(Gs_mol,Tair,pressure)
    
else if (formulation == "Penman-Monteith")
    
    check_input(data,list(Tair,pressure,VPD,LE,Rn,Ga,G,S))
    
    if(!is_null(G))
      if (!missing_G_as_NA){G[is_na(G)] = 0}
else 
      cat("Ground heat flux G is not provided and set to 0.",fill=TRUE)
      G = 0
end
    
    if(!is_null(S))
      if(!missing_S_as_NA){S[is_na(S)] = 0 }
else 
      cat("Energy storage fluxes S are not provided and set to 0.",fill=TRUE)
      S = 0
end
    
    Delta = Esat_slope(Tair,Esat_formula,constants)[,"Delta"]
    gamma = psychrometric_constant(Tair,pressure,constants)
    rho   = air_density(Tair,pressure,constants)
    
    Gs_ms  = ( LE * Ga * gamma ) / ( Delta * (Rn-G-S) + rho * constants[:cp] * Ga * VPD - LE * ( Delta + gamma ) )
    Gs_mol = ms_to_mol(Gs_ms,Tair,pressure)
    
end
  
  return(DataFrame(Gs_ms,Gs_mol))

end