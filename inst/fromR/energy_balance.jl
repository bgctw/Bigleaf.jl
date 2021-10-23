########################
#### Energy balance ####
########################

#' Biochemical Energy
#' 
#' Radiant energy absorbed in photosynthesis or heat release by respiration calculated
#'              from net ecosystem exchange of CO2 (NEE).  
#' 
#' - NEE     Net ecosystem exchange (umol CO2 m-2 s-1)
#' - alpha   Energy taken up/released by photosynthesis/respiration per mol CO2 fixed/respired (J umol-1)
#' 
#' # Details
 The following sign convention is employed: NEE is negative when carbon is taken up by the ecosystem.
#'          Positive values of the resulting biochemical energy mean that energy (heat) is taken up by the ecosystem, 
#'          negative ones that heat is released.
#'          The value of alpha is taken from Nobel 1974 (see Meyers & Hollinger 2004), but other values
#'          have been used (e.g. Blanken et al., 1997)
#' 
#' # Value
 - Sp -: biochemical energy (W m-2)
#' 
#' #References
#' Meyers, T_P., Hollinger, S_E. 2004: An assessment of storage terms in the surface energy
#'             balance of maize and soybean. Agricultural and Forest Meteorology 125, 105-115.
#'             
#'             Nobel, P_S., 1974: Introduction to Biophysical Plant Physiology.
#'             Freeman, New York.
#'             
#'             Blanken, P_D. et al., 1997: Energy balance and canopy conductance of a boreal aspen
#'             forest: Partitioning overstory and understory components. 
#'             Journal of Geophysical Research 102, 28915-28927. 
#'             
#' ```@example; output = false
#' ``` 
#' # Calculate biochemical energy taken up by the ecosystem with 
#' # a measured NEE of -30umol CO2 m-2 s-1             
#' biochemical_energy(NEE=-30)            
#'            
#' @export 
function biochemical_energy(NEE,alpha=0.422)
  Sp = alpha*-NEE
  return(Sp)
end




#' Energy-Use Efficiency (EUE)
#' 
#' Fraction of net radiation fixed by primary productivity.
#' 
#' - GPP     Gross primary productivity exchange (umol CO2 m-2 s-1)
#' - alpha   Energy taken up/released by photosynthesis/respiration (J umol-1)
#' - Rn      Net radiation (W m-2)
#' 
#' # Details
 Energy use efficiency is calculated as:
#' 
#'            ``EUE = sum(GPP)/sum(Rn)``
#'          
#'          where the sums are calculated for complete cases of GPP and Rn over
#'          the entire time period.
#' 
#' # Value
 - EUE -: Energy use efficiency (-)
#' 
#' #See also
#' [`light_use_efficiency`](@ref)
#' 
#' ```@example; output = false
#' ```
#' energy_use_efficiency(GPP=20,Rn=500)
#' 
#' @export 
function energy_use_efficiency(GPP,alpha=0.422,Rn)
  
  Sp = biochemical_energy(-GPP,alpha)
  
  comp  = complete_cases(Sp,Rn) 
  
  Sp_sum = sum(Sp[comp],na_rm=T)
  Rn_sum = sum(Rn[comp],na_rm=T)
  
  EUE = Sp_sum/Rn_sum
  
  return(c("EUE"=EUE))
end





#' Energy Balance Closure
#' 
#' Calculates the degree of the energy balance non-closure for the entire time span
#'              based on the ratio of two sums (energy balance ratio), and ordinary least squares (OLS).
#' 
#' - data  Data_frame or matrix containing all required variables.
#' - Rn    Net radiation (W m-2)
#' - G     Ground heat flux (W m-2); optional
#' - S     Sum of all storage fluxes (W m-2); optional
#' - LE    Latent heat flux (W m-2)
#' - H     Sensible heat flux (W m-2)
#' - instantaneous    should the energy balance be calculated at the time step 
#'                         of the observations (`TRUE`), or over the entire time period
#'                         provided as input (`false`)
#' - missing_G_as_NA  if `TRUE`, missing G are treated as `NA`s ,otherwise set to 0. 
#' - missing_S_as_NA  if `TRUE`, missing S are treated as `NA`s, otherwise set to 0.
#' 
#' 
#' # Details
 The energy balance ratio (EBR) is calculated as:
#'          
#'            ``EBR = sum(LE + H)/sum(Rn - G - S)``
#'          
#'          the sum is taken for all time steps with complete observations (i.e. where
#'          all energy balance terms are available).
#' 
#' # Value
 a named vector containing:
#'         - n: number of complete (all energy balance terms available) observations
#'         - intercept: intercept of the OLS regression
#'         - slope: slope of the OLS regression
#'         - r_squared: r^2 of the OLS regression
#'         - EBR: energy balance ratio
#'         
#'         if `instantaneous = TRUE`, only `EBR` is returned.
#' 
#' #References
#' Wilson K., et al. 2002: Energy balance closure at FLUXNET sites.
#'             Agricultural and Forest Meteorology 113, 223-243.
#'
#' ```@example; output = false
#' ``` 
#' ## characterize energy balance closure for DE-Tha in June 2014
#' energy_closure(DE_Tha_Jun_2014,instantaneous=false)
#' 
#' ## look at half-hourly closure 
#' EBR_inst = energy_closure(DE_Tha_Jun_2014,instantaneous=TRUE)
#' summary(EBR_inst)
#' 
#' @importFrom stats complete_cases lm
"""
"""
function energy_closure(data,Rn="Rn",G=NULL,S=NULL,LE="LE",H="H",instantaneous=false,
                           missing_G_as_NA=false,missing_S_as_NA=false)
  
  check_input(data,list(Rn,LE,H,G,S))
  
  if(!is_null(G))
    if (!missing_G_as_NA){G[is_na(G)] = 0}
else 
    cat("Ground heat flux G is not provided and set to 0.",fill=TRUE)
    G = rep(0,nrow(data))
end
  
  if(!is_null(S))
    if(!missing_S_as_NA){S[is_na(S)] = 0 }
else 
    cat("Energy storage fluxes S are not provided and set to 0.",fill=TRUE)
    S = rep(0,nrow(data))
end
  
  if (!instantaneous)
    comp = complete_cases(Rn,G,S,LE,H)
    n    = sum(comp)
    
    EBR = sum(LE[comp] + H[comp]) / sum(Rn[comp] - G[comp] - S[comp])
    
    emod = lm(c(LE + H) ~ c(Rn - G - S))
    intercept = summary(emod)$coef[1,1]
    slope     = summary(emod)$coef[2,1]
    r_squared = summary(emod)$r_squared
    
    return(c("n"=n,"intercept"=round(intercept,3),"slope"=round(slope,3),"r^2"=round(r_squared,3),"EBR"=round(EBR,3)))
    
else 
    
    EBR = (LE + H) /(Rn - G - S)
    
    return(EBR)
end
end   



#' Isothermal Net Radiation
#' 
#' Calculates the isothermal net radiation, i.e. the net radiation 
#'              that the surface would receive if it had the same temperature than
#'              the air.
#'              
#' - data       Data_frame or matrix containing all required variables
#' - Rn         Net radiation (W m-2)
#' - Tair       Air temperature (degC)
#' - Tsurf      Surface temperature (degC)
#' - emissivity Emissivity of the surface (-)
#' - constants  sigma - Stefan-Boltzmann constant (W m-2 K-4) 
#'                   Kelvin - conversion degree Celsius to Kelvin 
#'
#' # Details
 The isothermal net radiation (Rni) is given by:
#'          
#'            ``Rni = Rn + \epsilon * \sigma * (Tsurf^4 - Tair^4)``
#'          
#'          where ``\epsilon`` is the emissivity of the surface. Tsurf and Tair
#'          are in Kelvin.
#'          
#' # Value
 - Rni -: isothermal net radiation (W m-2)
#' 
#' #References
#' Jones, H. 2014: Plants and Microclimate. 3rd edition, Cambridge
#'             University Press.
#' 
#' ```@example; output = false
#' ``` 
#' # calculate isothermal net radiation of a surface that is 2degC warmer than the air.
#' isothermal_Rn(Rn=400,Tair=25,Tsurf=27,emissivity=0.98) 
#' 
"""
"""
function isothermal_Rn(data,Rn="Rn",Tair="Tair",Tsurf="Tsurf",emissivity,
                          constants=bigleaf_constants())
  
  check_input(data,list(Rn,Tair,Tsurf))
  
  Tair  = Tair + constants[:Kelvin]
  Tsurf = Tsurf + constants[:Kelvin]
  
  Rni = Rn + emissivity * constants[:sigma] * (Tsurf^4 - Tair^4)
  
  return(Rni)
  
end