###############################
### WUE metrics calculation ###
###############################

#' Water-Use Efficiency Metrics
#' 
#' Calculation of various water use efficiency (WUE) metrics.
#' 
#' - data      Data_frame or matrix containing all required variables
#' - GPP       Gross primary productivity (umol CO2 m-2 s-1)
#' - NEE       Net ecosystem exchange (umol CO2 m-2 s-1)
#' - LE        Latent heat flux (W m-2)
#' - VPD       Vapor pressure deficit (kPa)
#' - Tair      Air temperature (deg C)
#' - constants Cmol - molar mass of carbon (kg mol-1) 
#'                  umol2mol - conversion micromole (umol) to mole (mol) 
#'                  kg2g - conversion kilogram (kg) to gram (g)
#'
#' # Details
 the following metrics are calculated:
#' 
#'          Water-use efficiency (WUE):
#'          
#'            ``WUE = GPP / ET``
#'          
#'          Water-use efficiency based on NEE (WUE_NEE):
#'          
#'            ``WUE_NEE = NEE / ET``
#'          
#'          Inherent water-use efficiency (IWUE; Beer et al. 2009):
#'          
#'            ``IWUE = (GPP * VPD) / ET``
#'          
#'          Underlying water-use efficiency (uWUE; Zhou et al. 2014):
#'          
#'            ``uWUE= (GPP * sqrt(VPD)) / ET``
#'          
#'          All metrics are calculated based on the median of all values. E_g.
#'          WUE = median(GPP/ET,na_rm=TRUE)
#' 
#' # Value
 a named vector with the following elements:
#'         - WUE: Water-use efficiency (gC (kg H20)-1)
#'         - WUE_NEE: Water-use efficiency based on NEE (gC (kg H20)-1)
#'         - IWUE: Inherent water-use efficiency (gC kPa (kg H20)-1)
#'         - uWUE: Underlying water-use efficiency (gC kPa^0.5 (kg H20)-1)
#' 
#' #Note
#' Units for VPD can also be hPa. Units change accordingly.
#'       WUE_NEE is calculated based on the absolute value of NEE (the sign convention does not matter here).
#' 
#' #References
#' Beer, C., et al., 2009: Temporal and among-site variability of inherent
#'             water use efficiency at the ecosystem level. Global Biogeochemical Cycles 23, GB2018.
#'             
#'             Zhou, S., et al., 2014: The effect of vapor pressure deficit on water
#'             use efficiency at the sub-daily time scale. Geophysical Research Letters 41.
#'     
#' #See also
#' [`stomatal_slope`](@ref) for a measure of intrinsic WUE          
#'                             
#' ```@example; output = false
#' ``` 
#' ## filter data for dry periods and daytime at DE-Tha in June 2014
#' DE_Tha_Jun_2014_2 = filter_data(DE_Tha_Jun_2014,quality_control=false,
#'                                  vars_qc=c("Tair","precip","VPD","H","LE"),
#'                                  filter_growseas=false,filter_precip=TRUE,
#'                                  filter_vars=c("Tair","PPFD","ustar"),
#'                                  filter_vals_min=c(5,200,0.2),
#'                                  filter_vals_max=c(NA,NA,NA),NA_as_invalid=TRUE,
#'                                  quality_ext="_qc",good_quality=c(0,1),
#'                                  missing_qc_as_bad=TRUE,GPP="GPP",doy="doy",
#'                                  year="year",tGPP=0.5,ws=15,min_int=5,precip="precip",
#'                                  tprecip=0.1,precip_hours=24,records_per_hour=2)
#' 
#' ## calculate WUE metrics in the filtered periods
#' WUE_metrics(DE_Tha_Jun_2014_2)
#'                         
#' @importFrom stats median                                     
"""
"""
function WUE_metrics(data,GPP="GPP",NEE="NEE",LE="LE",VPD="VPD",Tair="Tair",
                        constants=bigleaf_constants())
  
  check_input(data,list(GPP,NEE,LE,VPD,Tair))
  
  ET  = LE_to_ET(LE,Tair)                 # kg H2O m-2 s-1
  GPP = (GPP * constants[:umol2mol] * constants[:Cmol]) * constants[:kg2g]  # gC m-2 s-1
  NEE = (NEE * constants[:umol2mol] * constants[:Cmol]) * constants[:kg2g]  # gC m-2 s-1
  
  WUE     = median(GPP/ET,na_rm=TRUE)
  WUE_NEE = median(abs(NEE)/ET,na_rm=TRUE)
  IWUE    = median((GPP*VPD)/ET,na_rm=TRUE)
  uWUE    = median((GPP*sqrt(VPD))/ET,na_rm=TRUE)
  
  return(c(WUE=WUE,WUE_NEE=WUE_NEE,IWUE=IWUE,uWUE=uWUE))
end
