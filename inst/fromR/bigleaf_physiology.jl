############################
#### Big-Leaf Physiology ###-------------------------------------------------------------------------------
############################

#' Bulk Intercellular CO2 Concentration
#' 
#' Bulk canopy intercellular CO2 concentration (Ci) calculated based on Fick's law
#'              given surface conductance (Gs), gross primary productivity (GPP) and 
#'              atmospheric CO2 concentration (Ca).
#'                            
#' - data             Data_Frame or matrix with all required columns                            
#' - Ca               Atmospheric or surface CO2 concentration (umol mol-1)              
#' - GPP              Gross primary productivity (umol CO2 m-2 s-1)
#' - Gs               Surface conductance to water vapor (mol m-2 s-1)
#' - Rleaf            Ecosystem respiration stemming from leaves (umol CO2 m-2 s-1); defaults to 0          
#' - missing_Rleaf_as_NA if Rleaf is provided, should missing values be treated as `NA` (`TRUE`)
#'                            or set to 0 (`false`, the default)?
#' - constants        DwDc - Ratio of the molecular diffusivities for water vapor and CO2 (-)
#' 
#' # Details
 Bulk intercellular CO2 concentration (Ci) is given by:
#' 
#'            ``Ci = Ca - (GPP - Rleaf)/(Gs/1.6)``
#'          
#'          where Gs/1.6 (mol m-2 s-1) represents the surface conductance to CO2.
#'          Note that Gs is required in mol m-2 s-1 for water vapor. Gs is converted to
#'          its value for CO2 internally.
#'          Ca can either be atmospheric CO2 concentration (as measured), or surface
#'          CO2 concentration as calculated from [`surface_CO2`](@ref).
#'          
#' #Note
#' The equation is based on Fick's law of diffusion and is equivalent to the
#'       often used equation at leaf level (ci = ca - An/gs).
#'       Note that GPP and Gs have a different interpretation than An and gs.
#'       Gs comprises non-physiological contributions (i.e. physical evaporation)
#'       and is confounded by physical factors (e.g. energy balance non-closure).
#'       GPP does not account for dark respiration and is further subject to uncertainties
#'       in the NEE partitioning algorithm used. Leaf respiration can be provided,
#'       but it is usually not known at ecosystem level (as a consequence, Ci is likely to be 
#'       slightly underestimated)
#'       This function should be used with care and the resulting Ci might not be
#'       readily comparable to its leaf-level analogue and/or physiological meaningful.          
#' 
#' # Value
 - Ci -: Bulk canopy intercellular CO2 concentration (umol mol-1)
#' 
#' #References
#' Kosugi Y. et al., 2013: Determination of the gas exchange phenology in an
#'             evergreen coniferous forest from 7 years of eddy covariance flux data using
#'             an extended big-leaf analysis. Ecol Res 28, 373-385.
#'             
#'             Keenan T., Sabate S., Gracia C., 2010: The importance of mesophyll conductance in
#'             regulating forest ecosystem productivity during drought periods. Global Change Biology
#'             16, 1019-1034.
#'             
#' ```@example; output = false
#' ``` 
#' # calculate bulk canopy Ci of a productive ecosystem
#' intercellular_CO2(Ca=400,GPP=40,Gs=0.7)
#'  
#' # note the sign convention for NEE
#' 
"""
"""
function intercellular_CO2(data,Ca="Ca",GPP="GPP",Gs="Gs_mol",Rleaf=NULL,
                              missing_Rleaf_as_NA=false,constants=bigleaf_constants())
  
  check_input(data,list(Ca,GPP,Gs))
  
  if(!is_null(Rleaf))
    if(!missing_Rleaf_as_NA){Rleaf[is_na(Rleaf)] = 0 }
else 
    cat("Respiration from the leaves is ignored and set to 0.",fill=TRUE)
    Rleaf = 0
end
  
  Ci = Ca - (GPP - Rleaf)/(Gs/constants[:DwDc])
  
  return(Ci)
  
end



#' Bulk Canopy Photosynthetic Capacity (Vcmax and Jmax)
#' 
#' Bulk canopy maximum carboxylation rate (Vcmax25), and maximum electron
#'              transport rate (Jmax25) at 25 degrees Celsius from bulk intercellular 
#'              CO2 concentration using the Farquhar et al. 1980 model for C3 photosynthesis.
#'           
#' - data      Data_Frame or matrix with all required columns   
#' - C3        C3 vegetation (`TRUE`, the default) or C4 vegetation (`false`)?              
#' - Temp      Surface (or air) temperature (degC) 
#' - GPP       Gross primary productivity (umol m-2 s-1)
#' - Ci        Bulk canopy intercellular CO2 concentration (umol mol-1)
#' - PPFD      Photosynthetic photon flux density (umol m-2 s-1) 
#' - PPFD_j    PPFD threshold, below which the canopy is considered to 
#'                  be RuBP regeneration limited. Defaults to 500 umol m-2 s-1.
#' - PPFD_c    PPFD threshold, above which the canopy is considered to 
#'                  be Rubisco limited. Defaults to 1000 umol m-2 s-1.
#' - Rleaf     Ecosystem respiration stemming from leaves (umol CO2 m-2 s-1); defaults to 0 
#' - Oi        Intercellular O2 concentration (mol mol-1)
#' - Kc25      Michaelis-Menten constant for CO2 at 25 degC (umol mol-1)
#' - Ko25      Michaelis-Menten constant for O2 at 25 degC (mmol mol-1)
#' - Gam25     Photorespiratory CO2 compensation point ('Gamma star') 
#'                  at 25 degC (umol mol-1)
#' - Kc_Ha     Activation energy for Kc (kJ mol-1)
#' - Ko_Ha     Activation energy for Ko (kJ mol-1)
#' - Gam_Ha    Activation energy for Gam (kJ mol-1)
#' - Vcmax_Ha  Activation energy for Vcmax (kJ mol-1)
#' - Vcmax_Hd  Deactivation energy for Vcmax (kJ mol-1)
#' - Vcmax_dS  Entropy term for Vcmax (kJ mol-1 K-1)
#' - Jmax_Ha   Activation energy for Jmax (kJ mol-1)
#' - Jmax_Hd   Deactivation energy for Jmax (kJ mol-1)
#' - Jmax_dS   Entropy term for Jmax (kJ mol-1 K-1)
#' - Theta     Curvature term in the light response function of J (-)
#' - alpha_canopy Canopy absorptance (-)
#' - missing_Rleaf_as_NA if Rleaf is provided, should missing values be treated as `NA` (`TRUE`)
#'                            or set to 0 (`false`, the default)?
#' - Ci_C4        intercellular CO2 concentration below which photosynthesis
#'                     is considered to be CO2-limited (umol mol-1), ignored
#'                     if `C3 = TRUE`. 
#' - constants    Kelvin - conversion degree Celsius to Kelvin 
#'                     Rgas - universal gas constant (J mol-1 K-1) 
#'                     kJ2J - conversion kilojoule (kJ) to joule (J) 
#'                     J2kJ - conversion joule (J) to kilojoule (kJ) 
#'                     se_median - conversion standard error (SE) of the mean to SE of the median
#'                  
#' # Details
 The maximum carboxylation rate at 25degC (Vcmax25) and the maximum electron
#'          transport rate at 25degC (Jmax25), which characterize photosynthetic capacity,
#'          are calculated as at leaf level. 
#'          The required variables Gs and Ci can be calculated from 
#'          [`surface_conductance`](@ref) and [`intercellular_CO2`](@ref), respectively.
#'          
#'          Gas exchange parameters are taken from Bernacchi et al. 2001 (apparent values, which
#'          assume an infinite mesophyll conductance). Negative and very low Ci values 
#'          (the threshold is set to Ci < 80umol mol-1 at the moment) are filtered out.
#'          
#'          Vcmax is calculated from the photosynthesis model by Farquhar et al. 1980.
#'          If net photosynthesis is Rubisco-limited (RuBP-saturated carboxylation
#'          rate, i.e. light has to be (near-)saturating):
#'         
#'            ``Vcmax = (GPP * (Ci + Kc*(1.0 + Oi/Ko))) / (Ci - Gam)``
#'          
#'          where Kc and Ko are the Michaelis-Menten constants for CO2 and O2 (mmol mol-1),
#'          respectively, Oi is the O2 concentration, and Gam is the photorespiratory CO2
#'          compensation point (umol mol-1).
#'          Under low-light conditions, the electron transport rate J is calculated from
#'          the RuBP regeneration-limited photosynthesis rate:
#'          
#'            ``J = (GPP * (4.0 * Ci + 8.0 * Gam) / (Ci - Gam)``
#'          
#'          In this function, bulk canopy photosynthesis is assumed to be Rubisco/RuBP-regeneration
#'          limited, if incoming PPFD is above/below a specified threshold or range. These ranges
#'          are determined by the parameters `PPFD_j` and `PPFD_c`. If, for example,
#'          `PPFD_j = c(100,400)`, all conditions with a PPFD between 100 and 400 are assumed
#'          to be in the RuBP-regeneration (i.e. light-limited) photosynthesis domain. The electron
#'          transport rate J is then only calculated for periods that meet this criterion.
#'          
#'          Jmax is calculated from J and absorbed irradiance:
#'          
#'            \deqn{J = (APPFD_PSII + Jmax - sqrt((APPFD_PSII + Jmax)^2 - 
#'                     4.0 * Theta * APPFD_PSII * Jmax)) / (2.0 * Theta)
#'                 }
#'               
#'          where APPFD_PSII is the absorbed PPFD by photosystem II (PS II), 
#'          and Theta is a curvature parameter. APPFD_PSII is calculated as
#'          
#'            ``PPFD * alpha_canopy * 0.85 * beta``
#'          
#'          where alpha_canopy is canopy-scale absorptance, 0.85 is a correction factor,
#'          and beta is the fraction of photons absorbed by PS II (assumed 0.5).
#'          alpha_canopy accounts for non-absorbing components of the ecosystem such as
#'          stems or soil, and is very likely ecosystem-specific. This parameter is relatively
#'          sensitive for the determination of Jmax25 at some sites.
#'          
#'          Vcmax and Jmax at canopy level are assumed to follow the same temperature response
#'          as at leaf level. Hence, the respective parameter k at 25degC (k25) is calculated as 
#'          (see e.g. Kattge & Knorr 2007):
#'          
#'            \deqn{k25 = k / 
#'                        ( exp(Ha * (Temp - Tref) / (Tref * Rgas * Temp)) *
#'                        (1 + exp((Tref * dS - Hd) / (Tref * Rgas))) /
#'                        (1 + exp((Temp * dS - Hd) / (Temp * Rgas)))
#'                        )
#'                  }
#'          
#'          where Ha is the activation energy (kJ mol-1), Hd is the deactivation energy (kJ mol-1),
#'          and dS is the entropy term (kJ mol-1 K-1) of the respective parameter. Tref is set
#'          to 298.15 K.
#'          
#'          For C4 photosynthesis, the simplified model by von Caemmerer 2000 is used.
#'          For light-saturated photosynthesis, Vcmax is given by:
#'          
#'            ``Vcmax = GPP``
#'          
#'          Note that in addition to the range `PPFD_c`, the range `Ci_C4`
#'          discards all periods with low Ci, in which photosynthesis is likely to
#'          be CO2-limited (see von Caemmerer 2000 for details).
#'          
#'          In the light-limited case, J is calculated as:
#'          
#'            ``J = 3 * GPPj / (1 - 0.5) ``
#'          
#'          The calculation of Jmax25 and Vcmax25 is identical to C3 photosynthesis
#'          as described above.
#'          
#' #Note
#'   The critical assumption is that bulk canopy photosynthesis is limited by
#'         one of the two limitation states. Incoming PPFD is assumed to determine
#'         the limitation states. Note however that the ranges (`PPFD_j` and `PPFD_c`)
#'         are likely ecosystem-specific. E_g. dense canopies presumably require higher
#'         `PPFD_c` thresholds than open canopies. A threshold of 500 umol m-2 s-1 PPFD
#'         for Rubisco-limited photosynthesis was assumed a reasonable working assumption (see Kosugi et al. 2013).
#'         Here, `PPFD_c` defaults to 1000 umol m-2 s-1. Note that even under very high/low irradiances,
#'         not all photosynthetically active plant material of an ecosystem will be in the same
#'         limitation state. Note that parameters describing bulk canopy photosynthetic capacity are not directly 
#'         comparable to their leaf-level counterparts, as the former integrate over the entire canopy
#'         depth (i.e. are given per ground area, and not per leaf area).
#'         In general, the function should be used with care!
#'          
#' # Value
 a DataFrame with the following columns:
#'         - Vcmax25: maximum bulk canopy carboxylation rate at 25degC (umol m-2 (ground) s-1)
#'         - Jmax25: maximum bulk canopy electron transport rate at 25degC (umol m-2 (ground) s-1)
#'        
#' #References
#' Lloyd J. et al., 1995: A simple calibrated model of Amazon rainforest productivity
#'             based on leaf biochemical properties. Plant, Cell and Environment 18, 1129-1145.
#' 
#'             Rayment M_B., Loustau D., Jarvis P_G., 2002: Photosynthesis and respiration
#'             of black spruce at three organizational scales: shoot, branch and canopy.
#'             Tree Physiology 22, 219-229.
#' 
#'             Kosugi Y. et al., 2013: Determination of the gas exchange phenology in an
#'             evergreen coniferous forest from 7 years of eddy covariance flux data using
#'             an extended big-leaf analysis. Ecol Res 28, 373-385. 
#'             
#'             Ueyama M. et al, 2016: Optimization of a biochemical model with eddy covariance
#'             measurements in black spruce forests of Alaska for estimating CO2 fertilization
#'             effects. Agricultural and Forest Meteorology 222, 98-111.
#'             
#'             Bernacchi C_J., Singsaas E_L., Pimentel C., Portis JR A_R., Long S_P., 2001:
#'             Improved temperature response functions for models of Rubisco-limited
#'             photosynthesis. Plant, Cell and Environment 24, 253-259. 
#'             
#'             Bernacchi C_J., Pimentel C., Long S_P., 2003: In vivo temperature response
#'             functions of parameters required to model RuBP-limited photosynthesis.
#'             Plant, Cell and Environment 26, 1419-1430.
#'             
#'             von Caemmerer, 2000: Biochemical models of leaf photosynthesis. Techniques
#'             in plant sciences No. 2. CSIRO Publishing, Collingwood VIC, Australia.
#'
#' #See also
#' [`intercellular_CO2`](@ref), [`Arrhenius_temp_response`](@ref)
#'
#' ```@example; output = false
#' ``` 
#' DE_Tha_Jun_2014_2 = filter_data(DE_Tha_Jun_2014,quality_control=false,
#'                                  vars_qc=c("Tair","precip","VPD","H","LE"),
#'                                  filter_growseas=false,filter_precip=TRUE,
#'                                  filter_vars=c("Tair","PPFD","ustar","LE"),
#'                                  filter_vals_min=c(5,200,0.2,0),
#'                                  filter_vals_max=c(NA,NA,NA,NA),NA_as_invalid=TRUE,
#'                                  quality_ext="_qc",good_quality=c(0,1),
#'                                  missing_qc_as_bad=TRUE,GPP="GPP",doy="doy",
#'                                  year="year",tGPP=0.5,ws=15,min_int=5,precip="precip",
#'                                  tprecip=0.1,precip_hours=24,records_per_hour=2)
#' 
#' # calculate Ga
#' Ga = aerodynamic_conductance(DE_Tha_Jun_2014_2,Rb_model="Thom_1972")[,"Ga_h"]
#' 
#' # calculate Gs from the the inverted PM equation
#' Gs_PM = surface_conductance(DE_Tha_Jun_2014_2,Tair="Tair",pressure="pressure",
#'                              Rn="Rn",G="G",S=NULL,VPD="VPD",Ga=Ga,
#'                              formulation=Val(:Penman-Monteith))[,"Gs_mol"]
#' 
#' # calculate Ci 
#' Ci = intercellular_CO2(DE_Tha_Jun_2014_2,Ca="Ca",GPP="GPP",Gs=Gs_PM) 
#' 
#' # calculate Vcmax25 and Jmax25
#' photosynthetic_capacity(DE_Tha_Jun_2014_2,Temp="Tair",Ci=Ci,PPFD_j=c(200,500),PPFD_c=1000)
#' 
#' @importFrom stats optimize
#'                                                                        
#' @export                  
function photosynthetic_capacity(data,C3=TRUE,Temp,GPP="GPP",Ci,PPFD="PPFD",PPFD_j=c(200,500),PPFD_c=1000,
                                    Rleaf=NULL,Oi=0.21,Kc25=404.9,Ko25=278.4,Gam25=42.75,
                                    Kc_Ha=79.43,Ko_Ha=36.38,Gam_Ha=37.83,Vcmax_Ha=65.33,Vcmax_Hd=200,
                                    Vcmax_dS=0.635,Jmax_Ha=43.9,Jmax_Hd=200,Jmax_dS=0.640,
                                    Theta=0.7,alpha_canopy=0.8,missing_Rleaf_as_NA=false,Ci_C4=100,
                                    constants=bigleaf_constants())
  
  check_input(data,list(Temp,GPP,Ci,PPFD))
  
  Temp = Temp + constants[:Kelvin]
  Tref = 25.0 + constants[:Kelvin]
  
  if (C3){  # C3 vegetation
    Kc_Ha    = Kc_Ha * constants[:kJ2J]
    Ko_Ha    = Ko_Ha * constants[:kJ2J]
    Gam_Ha   = Gam_Ha * constants[:kJ2J]
    
    # Temperature dependencies of photosynthetic parameters 
    Kc  = Kc25 * exp(Kc_Ha * (Temp - Tref) / (Tref*constants[:Rgas]*Temp))
    Ko  = Ko25 * exp(Ko_Ha * (Temp - Tref) / (Tref*constants[:Rgas]*Temp))
    Gam = Gam25 * exp(Gam_Ha * (Temp - Tref) / (Tref*constants[:Rgas]*Temp))
    Ko  = Ko * constants[:J2kJ]
    
    # basic filtering on Ci 
    Ci[Ci < 80 | is_na(Ci)] = NA
    
    # Presumed limitation states 
    GPPc = GPPj = GPP
    GPPj[PPFD < PPFD_j[1] | PPFD > PPFD_j[2] | is_na(PPFD)] = NA
    GPPc[PPFD < PPFD_c | is_na(PPFD)] = NA
    
    if(!is_null(Rleaf))
      if(!missing_Rleaf_as_NA){Rleaf[is_na(Rleaf)] = 0 }
else 
      cat("Respiration from the leaves is ignored and set to 0.",fill=TRUE)
      Rleaf = 0
end
    
    # calculate Vcmax and J (electron transport rate)
    Vcmax = (GPPc-Rleaf) * (Ci + Kc*(1.0 + Oi/Ko)) / (Ci - Gam)
    J     = (GPPj-Rleaf) * (4.0 * Ci + 8.0 * Gam) / (Ci - Gam)
    
    
else {  # C4 vegetation
    
    # Presumed limitation states (C4) 
    GPPc = GPPj = GPP
    GPPj[PPFD < PPFD_j[1] | PPFD > PPFD_j[2] | is_na(PPFD) | Ci < 0] = NA
    GPPc[PPFD < PPFD_c | Ci < Ci_C4 | is_na(PPFD)] = NA
    
    Vcmax = GPPc
    J     = 3 * GPPj / (1 - 0.5)
    
end  
  
  
  # calculate Jmax from J
  APPFD_PSII = PPFD * alpha_canopy * 0.85 * 0.5
  
  calcJmax  = which(complete_cases(J,APPFD_PSII))
  if (length(calcJmax) > 0)
    Jmax = sapply(calcJmax, function(i) tryCatch(optimize(function(Jmax){abs(J[i] - c((APPFD_PSII[i] + Jmax - 
                                                                                          sqrt((APPFD_PSII[i] + Jmax)^2 - 
                                                                                                 4.0 * Theta * APPFD_PSII[i] * Jmax)) /
                                                                                         (2.0 * Theta)))},
                                                           interval=c(0,1000),tol=1e-02)$minimum,
                                                  error=function(err){NA}
    )
    )
else 
    warning("Not enough observations to calculate Jmax!")
    Jmax = NA
end
  
  # calculate Vcmax25 and Jmax25
  Vcmax25 = Arrhenius_temp_response(Vcmax,Temp-constants[:Kelvin],Ha=Vcmax_Ha,
                                     Hd=Vcmax_Hd,dS=Vcmax_dS,constants=constants)
  
  Jmax25 = Arrhenius_temp_response(Jmax,Temp[calcJmax]-constants[:Kelvin],Ha=Jmax_Ha,
                                    Hd=Jmax_Hd,dS=Jmax_dS,constants=constants)
  
  
  # calculate medians and standard errors of the median
  Vcmax25_Median = median(Vcmax25,na_rm=TRUE)
  Vcmax25_SE     = constants[:se_median] * sd(Vcmax25,na_rm=TRUE)/sqrt((sum(!is_na(Vcmax25))))
  Jmax25_Median  = median(Jmax25,na_rm=TRUE)
  Jmax25_SE      = constants[:se_median] * sd(Jmax25,na_rm=TRUE)/sqrt((sum(!is_na(Jmax25))))
  
  return(c("Vcmax25"=round(Vcmax25_Median,2),"Vcmax25_SE"=round(Vcmax25_SE,2),
           "Jmax25"=round(Jmax25_Median,2),"Jmax25_SE"=round(Jmax25_SE,2)))
  
end




#' (Modified) Arrhenius Temperature Response Function
#' 
#' (Modified) Arrhenius function describing
#'              the temperature response of biochemical parameters.
#'              
#' - param Parameter measured at measurement temperature (umol m-2 s-1)            
#' - Temp  Measurement temperature (degC)
#' - Ha    Activation energy for param (kJ mol-1)
#' - Hd    Deactivation energy for param (kJ mol-1)
#' - dS    Entropy term for param (kJ mol-1 K-1)
#' - constants Kelvin - conversion degree Celsius to Kelvin 
#'                  Rgas - universal gas constant (J mol-1 K-1) 
#'                  kJ2J - conversion kilojoule (kJ) to joule (J)
#'                  
#' # Details
 The function returns the biochemical rate at a reference
#'          temperature of 25degC given a predefined temperature response function.
#'          This temperature response is given by a modified form of the Arrhenius
#'          function:
#' 
#'             \deqn{param25 = param / 
#'                         ( exp(Ha * (Temp - Tref) / (Tref*Rgas*Temp)) *
#'                         (1 + exp((Tref*dS - Hd) / (Tref * Rgas))) /
#'                         (1 + exp((Temp*dS - Hd) / (Temp * Rgas)))
#'                         )
#'                  }
#'                  
#'          where param is the value/rate of the parameter at measurement temperature,
#'          Temp is temperature in K, Tref is reference temperature (298.15K), and Rgas
#'          is the universal gas constant (8.314 J K-1 mol-1). Ha is the activation
#'          energy (kJ mol-1), Hd is the deactivation energy (kJ mol-1), and dS the
#'          entropy term (kJ mol-1 K-1) of the respective parameter.
#'          
#'          If either Hd or dS or both are not provided, the equation above reduces
#'          to the first term (i.e. the common Arrhenius equation without the deactivation
#'          term.)         
#'                                  
#' # Value
 param25 - value of the input parameter at the reference temperature of 25degC (umol m-2 s-1)
#'               
#' #References
#' Johnson F_H., Eyring H., Williams R_W. 1942: 
#'             The nature of enzyme inhibitions in bacterial luminescence: sulfanilamide,
#'             urethane, temperature and pressure. Journal of cellular and comparative
#'             physiology 20, 247-268.
#' 
#'             Kattge J., Knorr W., 2007: Temperature acclimation in a biochemical
#'             model of photosynthesis: a reanalysis of data from 36 species.
#'             Plant, Cell and Environment 30, 1176-1190.
#'             
#' @export             
function Arrhenius_temp_response(param,Temp,Ha,Hd,dS,constants=bigleaf_constants())
  
  Temp = Temp + constants[:Kelvin]
  Tref = 25.0 + constants[:Kelvin]
  
  Ha = ifelse(missing(Ha),NA,Ha*constants[:kJ2J])
  Hd = ifelse(missing(Hd),NA,Hd*constants[:kJ2J])
  dS = ifelse(missing(dS),NA,dS*constants[:kJ2J])
  
  if (is_na(Ha))
    
    stop("Activation energy (Ha) has to be provided!")
    
end
  
  if (is_na(Hd) & is_na(dS))
    
    param25 = param / exp(Ha * (Temp - Tref) / (Tref*constants[:Rgas]*Temp))
  
else if (!is_na(Hd) & !is_na(dS))
    
    param25 = param /
      ( exp(Ha * (Temp - Tref) / (Tref*constants[:Rgas]*Temp)) *
        (1 + exp((Tref*dS - Hd) / (Tref * constants[:Rgas]))) /
        (1 + exp((Temp*dS - Hd) / (Temp * constants[:Rgas])))
      )
    
else if ((!is_na(Hd) & is_na(dS)) | (is_na(Hd) & !is_na(dS)) )

    warning("Both Hd and dS have to be provided for a temperature response
             that considers a temperature optimum and a deactivation term!
             Continue considering activation energy (Ha) only...")
    
    param25 = param / exp(Ha * (Temp - Tref) / (Tref*constants[:Rgas]*Temp))
    
end

  return(param25)
  
end



#' Stomatal Slope Parameter "g1"
#' 
#' Estimation of the intrinsic WUE metric "g1" (stomatal slope) 
#'              from nonlinear regression.
#' 
#' - data       Data_frame or matrix containing all required columns
#' - Tair       Air (or surface) temperature (deg C)
#' - pressure   Atmospheric pressure (kPa)
#' - GPP        Gross primary productivity (umol CO2 m-2 s-1)
#' - Gs         Surface conductance to water vapor (mol m-2 s-1)
#' - VPD        Vapor pressure deficit (kPa)
#' - Ca         Atmospheric CO2 concentration (air or surface) (umol mol-1)
#' - Rleaf      Ecosystem respiration stemming from leaves (umol CO2 m-2 s-1); defaults to 0 
#' - model      Stomatal model used. One of `"USO","Ball&Berry","Leuning"`.
#' - robust_nls Use robust nonlinear regression (`\link[robustbase]{nlrob`})? Default is `false`.
#' - nmin       Minimum number of data required to perform the fit; defaults to 40.
#' - fitg0      Should g0 and g1 be fitted simultaneously? 
#' - g0         Minimum stomatal conductance (mol m-2 s-1); ignored if `fitg0 = TRUE`.
#' - fitD0      Should D0 be fitted along with g1 (and g0 if `fitg0 = TRUE`)?; only used if `model = "Leuning"`.
#' - D0         Stomatal sensitivity parameter to VPD; only used if `model = "Leuning"` and `fitD0 = false`.
#' - Gamma      Canopy CO2 compensation point (umol mol-1); only used if `model = "Leuning"`. 
#'                   Can be a constant or a variable. Defaults to 50 umol mol-1.
#' - constants  Kelvin - conversion degree Celsius to Kelvin 
#'                   Rgas - universal gas constant (J mol-1 K-1) 
#'                   DwDc - Ratio of the molecular diffusivities for water vapor and CO2
#' - missing_Rleaf_as_NA if Rleaf is provided, should missing values be treated as `NA` (`TRUE`)
#'                            or set to 0 (`false`, the default)?
#' - ...        Additional arguments to `\link[stats]{nls`} or `\link[robustbase]{nlrob`} if `robust_nls = TRUE`.
#' 
#' # Details
 All stomatal models were developed at leaf-level, but its parameters 
#'          can also be estimated at ecosystem level (but be aware of caveats).
#'          
#'          The unified stomatal optimization (USO) model is given by (Medlyn et al. 2011):
#'      
#'             ``gs = g0 + 1.6*(1.0 + g1/sqrt(VPD)) * An/ca``
#'          
#'          The semi-empirical model by Ball et al. 1987 is defined as:
#'          
#'             ``gs = g0 + g1* ((An * rH) / ca)``
#'          
#'          Leuning 1995 suggested a revised version of the Ball&Berry model:
#'          
#'             ``gs = g0 + g1*An / ((ca - \\gamma) * (1 + VPD/D0))``
#'          
#'          where ``\\gamma`` is by default assumed to be constant, but likely varies with temperature and among
#'          plant species. 
#'          The equations above are valid at leaf-level. At ecosystem level, An is replaced by GPP (or GPP - Rleaf,
#'          where Rleaf is leaf respiration), and gs (stomatal conductance) by Gs (surface conductance). 
#'          The parameters in the models are estimated using nonlinear regression (`\link[stats]{nls`}) if
#'          `robust_nls = false` and weighted nonlinear regression if `robust_nls = TRUE`.
#'          The weights are calculated from `\link[robustbase]{nlrob`}, and `\link[stats]{nls`}
#'          is used for the actual fitting.
#'          Alternatively to measured VPD and Ca (i.e. conditions at instrument height), conditions at 
#'          the big-leaf surface can be provided. Those can be calculated using [`surface_conditions`](@ref).
#'          
#' 
#' # Value
 A `nls` model object, containing information on the fitted parameters, their uncertainty range,
#'         model fit, etc.
#' 
#' #References
#' Medlyn B_E., et al., 2011: Reconciling the optimal and empirical approaches to
#'             modelling stomatal conductance. Global Change Biology 17, 2134-2144.
#'             
#'             Ball T_J., Woodrow I_E., Berry J_A. 1987: A model predicting stomatal conductance
#'             and its contribution to the control of photosynthesis under different environmental conditions.
#'             In: Progress in Photosynthesis Research, edited by J_Biggins, pp. 221-224, Martinus Nijhoff Publishers,
#'             Dordrecht, Netherlands.
#'             
#'             Leuning R., 1995: A critical appraisal of a combined stomatal-photosynthesis
#'             model for C3 plants. Plant, Cell and Environment 18, 339-355.
#'             
#'             Knauer, J. et al., 2018: Towards physiologically meaningful water-use efficiency estimates
#'             from eddy covariance data. Global Change Biology 24, 694-710.
#' 
#' #See also
#' [`surface_conductance`](@ref)
#' 
#' ```@example; output = false
#' ``` 
#' ## filter data to ensure that Gs is a meaningful proxy to canopy conductance (Gc)
#' DE_Tha_Jun_2014_2 = filter_data(DE_Tha_Jun_2014,quality_control=false,
#'                                  vars_qc=c("Tair","precip","VPD","H","LE"),
#'                                  filter_growseas=false,filter_precip=TRUE,
#'                                  filter_vars=c("Tair","PPFD","ustar","LE"),
#'                                  filter_vals_min=c(5,200,0.2,0),
#'                                  filter_vals_max=c(NA,NA,NA,NA),NA_as_invalid=TRUE,
#'                                  quality_ext="_qc",good_quality=c(0,1),
#'                                  missing_qc_as_bad=TRUE,GPP="GPP",doy="doy",
#'                                  year="year",tGPP=0.5,ws=15,min_int=5,precip="precip",
#'                                  tprecip=0.1,precip_hours=24,records_per_hour=2)
#' 
#' # calculate Gs from the the inverted PM equation
#' Ga = aerodynamic_conductance(DE_Tha_Jun_2014_2,Rb_model="Thom_1972")[,"Ga_h"]
#' 
#' # if G and/or S are available, don't forget to indicate (they are ignored by default).
#' Gs_PM = surface_conductance(DE_Tha_Jun_2014_2,Tair="Tair",pressure="pressure",
#'                              Rn="Rn",G="G",S=NULL,VPD="VPD",Ga=Ga,
#'                              formulation=Val(:Penman-Monteith))[,"Gs_mol"]
#'                              
#' ### Estimate the stomatal slope parameter g1 using the USO model
#' mod_USO = stomatal_slope(DE_Tha_Jun_2014_2,model="USO",GPP="GPP",Gs=Gs_PM,
#'                           robust_nls=false,nmin=40,fitg0=false)
#'                           
#' ### Use robust regression to minimize influence of outliers in Gs                           
#' mod_USO = stomatal_slope(DE_Tha_Jun_2014_2,model="USO",GPP="GPP",Gs=Gs_PM,
#'                           robust_nls=TRUE,nmin=40,fitg0=false)
#' 
#' ### Estimate the same parameter from the Ball&Berry model and prescribe g0
#' mod_BB = stomatal_slope(DE_Tha_Jun_2014_2,model="Ball&Berry",GPP="GPP",
#'                          robust_nls=false,Gs=Gs_PM,g0=0.01,nmin=40,fitg0=false)
#' 
#' ## same for the Leuning model, but this time estimate both g1 and g0 (but fix D0)
#' mod_Leu = stomatal_slope(DE_Tha_Jun_2014_2,model="Leuning",GPP="GPP",Gs=Gs_PM,
#'                           robust_nls=false,nmin=40,fitg0=false,D0=1.5,fitD0=false)
#' 
#' @importFrom stats nls na_exclude
#' @importFrom robustbase nlrob
#' 
#' @export 
function stomatal_slope(data,Tair="Tair",pressure="pressure",GPP="GPP",Gs="Gs_mol",
                           VPD="VPD",Ca="Ca",Rleaf=NULL,model=c("USO","Ball&Berry","Leuning"),
                           robust_nls=false,nmin=40,fitg0=false,g0=0,fitD0=false,
                           D0=1.5,Gamma=50,missing_Rleaf_as_NA=false,
                           constants=bigleaf_constants(),...)
  
  model = match_arg(model)
  
  check_input(data,list(Tair,pressure,GPP,Gs,VPD,Ca))

  df   = DataFrame(Tair,pressure,GPP,Gs,VPD,Ca)
  DwDc = constants[:DwDc]  # ...to work within nls()
  
  
  if (model == "Leuning")
    check_input(data,Gamma)
    df$Gamma = Gamma
end
  
  
  
  if(!is_null(Rleaf))
    if(!missing_Rleaf_as_NA){Rleaf[is_na(Rleaf)] = 0 }
else 
    cat("Respiration from the leaves is ignored and set to 0.",fill=TRUE)
    Rleaf = 0
end
  
  GPP = (GPP - Rleaf)
  
  
  if (model == "Leuning")
    nr_data = sum(!is_na(GPP) & !is_na(Gs) & !is_na(VPD) & !is_na(Ca) & !is_na(Gamma))
else 
    nr_data = sum(!is_na(GPP) & !is_na(Gs) & !is_na(VPD) & !is_na(Ca))
end
  
  
  if (nr_data < nmin)
    stop("number of data is less than 'nmin'. g1 is not fitted to the data.")
else 
    
    if (model == "USO")
      
      if (fitg0)
        if (robust_nls)
          df$DwDc = rep(DwDc,nrow(df))
          mod_weights = nlrob(Gs ~ g0 + DwDc*(1.0 + g1/sqrt(VPD))*GPP/Ca,data=df,start=list(g0=0,g1=3),
                               na_action=na_exclude,...)$w
          mod = nls(Gs ~ g0 + DwDc*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g0=0,g1=3),weights=mod_weights,...)
else 
          mod = nls(Gs ~ g0 + DwDc*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g0=0,g1=3),...)
end
else 
        if (robust_nls)
          df$g0   = rep(g0,nrow(df))    # g0 as constant does not work in the nlrob function...
          df$DwDc = rep(DwDc,nrow(df))  # same with constants[:DwDc]
          mod_weights = nlrob(Gs ~ g0 + DwDc*(1.0 + g1/sqrt(VPD))*GPP/Ca,data=df,start=list(g1=3),
                               na_action=na_exclude,...)$w
          mod = nls(Gs ~ g0 + DwDc*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g1=3),weights=mod_weights,...)
else 
          mod = nls(Gs ~ g0 + DwDc*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g1=3),...)
end
end
      
else if (model == "Leuning")
      
      if (fitg0)
        if (fitD0)
          if (robust_nls)
            mod_weights = nlrob(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),data=df,
                                 start=list(g0=0,g1=9,D0=1.5),na_action=na_exclude,...)$w
            mod = nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g0=0,g1=9,D0=1.5),
                       weights=mod_weights,...)
else 
            mod = nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g0=0,g1=9,D0=1.5),...)
end
else 
          if (robust_nls)
            df$D0  = rep(D0,nrow(df))
            mod_weights = nlrob(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),data=df,
                                 start=list(g0=0,g1=9),na_action=na_exclude,...)$w
            mod = nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g0=0,g1=9),
                       weights=mod_weights,...)
else 
            mod = nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g0=0,g1=9),...)
end
end
else 
        if (fitD0)
          if (robust_nls)
            df$g0    = rep(g0,nrow(df))
            mod_weights = nlrob(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),data=df,
                                 start=list(g1=9,D0=1.5),na_action=na_exclude,...)$w
            mod = nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g1=9,D0=1.5),
                       weights=mod_weights,...)
else 
            mod = nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g1=9,D0=1.5),...)
end
else 
          if (robust_nls)
            df$g0  = rep(g0,nrow(df))
            df$D0  = rep(D0,nrow(df))
            mod_weights = nlrob(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),data=df,
                                 start=list(g1=9),na_action=na_exclude,...)$w
            mod = nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g1=9),
                       weights=mod_weights,...)
else 
            mod = nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g1=9),...)
end
end
end
      
else if (model == "Ball&Berry")
      
      rH = VPD_to_rH(VPD,Tair)
      df$rH = rH
      
      if (fitg0)
        if (robust_nls)
          mod_weights = nlrob(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g0=0,g1=9),data=df,
                               na_action=na_exclude,...)$w
          mod = nls(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g0=0,g1=9),weights=mod_weights,...)
else 
          mod = nls(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g0=0,g1=9),...)
end
else 
        if (robust_nls)
          df$g0   = rep(g0,nrow(df))
          mod_weights = nlrob(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g1=9),data=df,
                               na_action=na_exclude,...)$w
          mod = nls(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g1=9),weights=mod_weights,...)
else 
          mod = nls(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g1=9),...)
end
end
      
end
    
end
  
  return(mod)
  
end






#' Ecosystem Light Response
#' 
#' calculates GPP_ref at a reference (usually saturating) PPFD and 
#'              ecosystem quantum yield (alpha) using a rectangular light response curve.
#' 
#' - data      Data_frame or matrix containing all required columns
#' - NEE       Net ecosystem exchange (umol CO2 m-2 s-1)
#' - Reco      Ecosystem respiration (umol CO2 m-2 s-1)
#' - PPFD      Photosynthetic photon flux density (umol m-2 s-1)
#' - PPFD_ref  Reference PPFD (umol m-2 s-1) for which GPP_ref is estimated.
#'                  Default is 2000 umol m-2 s-1.
#' - ...       Additional arguments to `\link[stats]{nls`}
#' 
#' # Details
 A rectangular light response curve is fitted to NEE data. The curve
#'          takes the form as described in Falge et al. 2001:
#'          
#'             \deqn{-NEE = \\alpha PPFD / (1 - (PPFD / PPFD_ref) + \\alpha 
#'                          PPFD / GPP_ref) - Reco}
#'                       
#'          where ``\\alpha`` is the ecosystem quantum yield (umol CO2 m-2 s-1) (umol quanta m-2 s-1)-1, 
#'          and GPP_ref is the GPP at the reference PPFD (usually at saturating light). ``\\alpha`` 
#'          represents the slope of the light response curve, and is a measure for the light use
#'          efficiency of the canopy. 
#'          
#'          The advantage of this equation over the standard rectangular light response
#'          curve is that GPP_ref at PPFD_ref is more readily interpretable
#'          as it constitutes a value observed in the ecosystem, in contrast to 
#'          GPP_ref (mostly named 'beta') in the standard model that occurs at infinite light.
#'          `PPFD_ref` defaults to 2000 umol m-2 s-1, but other values can be used. For 
#'          further details refer to Falge et al. 2001.
#' 
#' #Note
#'   Note the sign convention. Negative NEE indicates that carbon is taken up
#'         by the ecosystem. Reco has to be 0 or larger.
#' 
#' # Value
 A `nls` model object containing estimates (+/- SE) for alpha and GPP_ref.
#' 
#' #References
#' Falge E., et al. 2001: Gap filling strategies for defensible annual
#'             sums of net ecosystem exchange. Agricultural and Forest Meteorology 107,
#'             43-69.
#'             
#'             Gilmanov T_G., et al. 2003: Gross primary production and light response
#'             parameters of four Southern Plains ecosystems estimated using long-term
#'             CO2-flux tower measurements. Global Biogeochemical Cycles 17, 1071.
#'             
#'             Reichstein M., Stoy P_C., Desai A_R., Lasslop G., Richardson A. 2012: 
#'             Partitioning of net fluxes. In: Eddy Covariance. A practical guide to
#'             measurement and data analysis. Aubinet M., Vesala T., Papale D. (Eds.).
#'             Springer.
#' 
#' @importFrom stats nls
#' 
"""
"""
function light_response(data,NEE="NEE",Reco="Reco",PPFD="PPFD",PPFD_ref=2000,...)

  check_input(data,list(NEE,Reco,PPFD))
  
  mod = nls(-NEE ~ alpha * PPFD / (1 - (PPFD / PPFD_ref) + (alpha * PPFD / GPP_ref)) - Reco,
             start=list(alpha=0.05,GPP_ref=30),...)
  
  return(mod)
end  




  

#' Light-Use Efficiency (LUE)
#' 
#' Amount of carbon fixed (GPP) per incoming light.
#' 
#' - GPP     Gross ecosystem productivity (umol CO2 m-2 s-1)
#' - PPFD    Photosynthetic photon flux density (umol quanta m-2 s-1)
#' 
#' # Details
 Light use efficiency is calculated as
#'          
#'             ``LUE = sum(GPP)/sum(PPFD)``
#'          
#'          where both GPP and PPFD are in umol m-2 s-1. A more meaningful 
#'          (as directly comparable across ecosystems) approach is to take 
#'          absorbed PPFD rather than incoming PPFD as used here.
#' 
#' # Value
 - LUE -: Light use efficiency (-)
#' 
#' #See also
#' [`energy_use_efficiency`](@ref)
#' 
#' ```@example; output = false
#' ```
#' light_use_efficiency(GPP=20,PPFD=1500)
#' 
#' @export 
function light_use_efficiency(GPP,PPFD)
  
  comp = complete_cases(GPP,PPFD)
  
  LUE = sum(GPP[comp],na_rm=TRUE)/sum(PPFD[comp],na_rm=TRUE)
  
  return(c("LUE"=LUE))
end
  



#' Stomatal Sensitivity to VPD
#' 
#' Sensitivity of surface conductance to vapor pressure deficit.
#' 
#' - data  Data_frame or matrix containing all required columns
#' - Gs    Surface conductance to water vapor (mol m-2 s-1)
#' - VPD   Vapor pressure deficit (kPa)
#' - ...   Additional arguments to `\link[stats]{nls`}
#' 
#' # Details
 The function fits the following equation (Oren et al. 1999):
#' 
#'             ``Gs = -m ln(VPD) + b``
#'
#'          where b is the reference surface conductance (Gs) at VPD=1kPa (in mol m-2 s-1),
#'          and m is the sensitivity parameter of Gs to VPD (in mol m-2 s-1 log(kPa)-1).
#'          The two parameters b and m are fitted using `\link[stats]{nls`}.
#'          VPD can be the one directly measured at instrument height, or the
#'          one at the surface, as returned by [`surface_conditions`](@ref).
#'          
#' # Value
 A `nls` model object containing (amongst others) estimates for the mean
#'         and standard errors of the parameters m and b.
#' 
#' #References
#' Oren R., et al. 1999: Survey and synthesis of intra- and interspecific
#'             variation in stomatal sensitivity to vapour pressure deficit. Plant,
#'             Cell & Environment 22, 1515-1526. 
#'             
#'             Novick K_A., et al. 2016: The increasing importance of atmospheric demand
#'             for ecosystem water and carbon fluxes. Nature Climate Change 6, 1023 - 1027.
#'          
#' #See also
#' [`surface_conductance`](@ref)          
#' 
#' ```@example; output = false
#' ```
#' ## calculate Ga, Gs, and the stomatal sensitivity to VPD for the site FR-Pue in
#' ## May 2012. Data are filtered for daytime, sufficiently high ustar, etc.
#' FR_Pue_May_2012_2 = filter_data(FR_Pue_May_2012,quality_control=TRUE,
#'                                  vars_qc=c("Tair","precip","H","LE"),
#'                                  filter_growseas=false,filter_precip=TRUE,
#'                                  filter_vars=c("Tair","PPFD","ustar","VPD"),
#'                                  filter_vals_min=c(5,200,0.2,0.3),
#'                                  filter_vals_max=c(NA,NA,NA,NA),
#'                                  NA_as_invalid=TRUE,quality_ext="_qc",
#'                                  good_quality=c(0,1),missing_qc_as_bad=TRUE,
#'                                  precip="precip",tprecip=0.1,precip_hours=24,
#'                                  records_per_hour=2)
#' Ga = aerodynamic_conductance(FR_Pue_May_2012_2)
#' Gs = surface_conductance(FR_Pue_May_2012_2,Ga=Ga[,"Ga_h"])
#' stomatal_sensitivity(FR_Pue_May_2012_2,Gs=Gs[,"Gs_mol"],VPD="VPD")
#' 
#' @importFrom stats nls
#' 
"""
"""
function stomatal_sensitivity(data,Gs="Gs_mol",VPD="VPD",...)
  
  check_input(data,list(Gs,VPD))
  
  mod = nls(Gs ~ -m * log(VPD) + b,start=list(m=0.05,b=0.2),...)
  
  return(mod)
end



  
  
  