#' Eddy Covariance Data of AT-Neu (Neustift)
#' 
#' Halfhourly eddy covariance Data of the site AT-Neu,
#'              a mountain meadow in Austria.
#'              (\url{http://sites_fluxdata_org/AT-Neu}). 
#'              Data are from July 2010.
#'              
#' @format A data frame with 1488 observations and 31 columns:
#'  \describe
#'    - year: year of measurement
#'    - month: month of measurement
#'    - doy: day of year
#'    - hour: hour (0 - 23.5)
#'    - Tair: Air temperature (degC) [TA_F]
#'    - Tair_qc: Quality control of `Tair` [TA_F_QC]
#'    - PPFD: Photosynthetic photon flux density (umol m-2 s-1) [PPFD_IN]
#'    - PPFD_qc: Quality control of `PPFD` [PPFD_IN_QC]
#'    - VPD: Vapor pressure deficit (kPa) [VPD_F]
#'    - VPD_qc: Quality control of `VPD` [VPD_F_QC]
#'    - pressure: Atmospheric pressure (kPa) [PA_F]
#'    - precip: precipitation (mm) [P_F]
#'    - precip_qc: Quality control of `precip` [P_F_QC]
#'    - ustar: friction velocity (m s-1) [USTAR]
#'    - wind: horizontal wind velocity (m s-1) [WS_F]
#'    - wind_qc: Quality control of `wind` [WS_F_QC]
#'    - Ca: CO2 concentration (ppm) [CO2_F_MDS]
#'    - Ca_qc: Quality control of `Ca` [CO2_F_MDS_QC]
#'    - LW_up: upward longwave radiation (W m-2) [LW_OUT]
#'    - Rn: Net radiation (W m-2) [NETRAD]
#'    - LE: Latent heat flux (W m-2) [LE_F_MDS]
#'    - LE_qc: Quality control of `LE` [LE_F_MDS_QC]
#'    - H: Sensible heat flux (W m-2) [H_F_MDS]
#'    - H_qc: Quality control of `H` [H_F_MDS_QC]
#'    - G: Ground heat flux (W m-2) [G_F_MDS]
#'    - G_qc: Quality control of `G` [G_F_MDS_QC]
#'    - NEE: Net ecosystem exchange (umol m-2 s-1) [NEE_VUT_USTAR50]
#'    - NEE_qc: Quality control of `NEE` [NEE_VUT_USTAR50_QC]
#'    - GPP: Gross primary productivity from nighttime partitioning (umol m-2 s-1) [GPP_NT_VUT_USTAR50]
#'    - GPP_qc: Quality control of `GPP` [NEE_VUT_USTAR50_QC]
#'    - Reco: Ecosystem respiration from nighttime partitioning (umol m-2 s-1) [RECO_NT_VUT_USTAR50]
#'  }
#'  
#' # Note
#' The original variable names as provided by the FLUXNET2015 dataset are 
#'       given in squared brackets. Note that variable units have been converted
#'       in some cases (e.g. VPD from hPa to kPa).      
#'  
#' @source original data were downloaded from
#'         \url{https://fluxnet_fluxdata_org/} (accessed 09 November 2016)  
"AT_Neu_Jul_2010"



#' Eddy Covariance Data of DE-Tha (Tharandt)
#' 
#' Halfhourly eddy covariance Data of the site DE-Tha,
#'              a spruce forest in Eastern Germany 
#'              (\url{http://sites_fluxdata_org/DE-Tha}). 
#'              Data are from June 2014.
#'              
#' @format A data frame with 1440 observations and 32 columns:
#'  \describe
#'    - year: year of measurement
#'    - month: month of measurement
#'    - doy: day of year
#'    - hour: hour (0 - 23.5)
#'    - Tair: Air temperature (degC) [TA_F]
#'    - Tair_qc: Quality control of `Tair` [TA_F_QC]
#'    - PPFD: Photosynthetic photon flux density (umol m-2 s-1) [PPFD_IN]
#'    - PPFD_qc: Quality control of `PPFD` [PPFD_IN_QC]
#'    - VPD: Vapor pressure deficit (kPa) [VPD_F]
#'    - VPD_qc: Quality control of `VPD` [VPD_F_QC]
#'    - pressure: Atmospheric pressure (kPa) [PA_F]
#'    - precip: precipitation (mm) [P_F]
#'    - precip_qc: Quality control of `precip` [P_F_QC]
#'    - ustar: friction velocity (m s-1) [USTAR]
#'    - wind: horizontal wind velocity (m s-1) [WS_F]
#'    - wind_qc: Quality control of `wind` [WS_F_QC]
#'    - Ca: CO2 concentration (ppm) [CO2_F_MDS]
#'    - Ca_qc: Quality control of `Ca` [CO2_F_MDS_QC]
#'    - LW_up: upward longwave radiation (W m-2) [LW_OUT]
#'    - LW_down: downward longwave radiation (W m-2) [LW_IN_F]
#'    - Rn: Net radiation (W m-2) [NETRAD]
#'    - LE: Latent heat flux (W m-2) [LE_F_MDS]
#'    - LE_qc: Quality control of `LE` [LE_F_MDS_QC]
#'    - H: Sensible heat flux (W m-2) [H_F_MDS]
#'    - H_qc: Quality control of `H` [H_F_MDS_QC]
#'    - G: Ground heat flux (W m-2) [G_F_MDS]
#'    - G_qc: Quality control of `G` [G_F_MDS_QC]
#'    - NEE: Net ecosystem exchange (umol m-2 s-1) [NEE_VUT_USTAR50]
#'    - NEE_qc: Quality control of `NEE` [NEE_VUT_USTAR50_QC]
#'    - GPP: Gross primary productivity from nighttime partitioning (umol m-2 s-1) [GPP_NT_VUT_USTAR50]
#'    - GPP_qc: Quality control of `GPP` [NEE_VUT_USTAR50_QC]
#'    - Reco: Ecosystem respiration from nighttime partitioning (umol m-2 s-1) [RECO_NT_VUT_USTAR50]
#'  }
#'  
#' # Note
#' The original variable names as provided by the FLUXNET2015 dataset are 
#'       given in squared brackets. Note that variable units have been converted
#'       in some cases (e.g. VPD from hPa to kPa).    
#'  
#' @source original data were downloaded from
#'         \url{https://fluxnet_fluxdata_org/} (accessed 09 November 2016) 
"DE_Tha_Jun_2014"



#' Eddy Covariance Data of FR-Pue (Puechabon)
#' 
#' Halfhourly eddy covariance Data of the site FR-Pue,
#'              a Mediterranean evergreen oak forest in Southern France
#'              (\url{http://sites_fluxdata_org/FR-Pue}).
#'              Data are from May 2012.
#'              
#' @format A data frame with 1488 observations and 29 columns:
#'  \describe
#'    - year: year of measurement
#'    - month: month of measurement
#'    - doy: day of year
#'    - hour: hour (0 - 23.5)
#'    - Tair: Air temperature (degC) [TA_F]
#'    - Tair_qc: Quality control of `Tair` [TA_F_QC]
#'    - PPFD: Photosynthetic photon flux density (umol m-2 s-1) [PPFD_IN]
#'    - PPFD_qc: Quality control of `PPFD` [PPFD_IN_QC]
#'    - VPD: Vapor pressure deficit (kPa) [VPD_F]
#'    - VPD_qc: Quality control of `VPD` [VPD_F_QC]
#'    - pressure: Atmospheric pressure (kPa) [PA_F]
#'    - precip: precipitation (mm) [P_F]
#'    - precip_qc: Quality control of `precip` [P_F_QC]
#'    - ustar: friction velocity (m s-1) [USTAR]
#'    - wind: horizontal wind velocity (m s-1) [WS_F]
#'    - wind_qc: Quality control of `wind` [WS_F_QC]
#'    - Ca: CO2 concentration (ppm) [CO2_F_MDS]
#'    - Ca_qc: Quality control of `Ca` [CO2_F_MDS_QC]
#'    - LW_up: upward longwave radiation (W m-2) [LW_OUT]
#'    - Rn: Net radiation (W m-2) [NETRAD]
#'    - LE: Latent heat flux (W m-2) [LE_F_MDS]
#'    - LE_qc: Quality control of `LE` [LE_F_MDS_QC]
#'    - H: Sensible heat flux (W m-2) [H_F_MDS]
#'    - H_qc: Quality control of `H` [H_F_MDS_QC]
#'    - NEE: Net ecosystem exchange (umol m-2 s-1) [NEE_VUT_USTAR50]
#'    - NEE_qc: Quality control of `NEE` [NEE_VUT_USTAR50_QC]
#'    - GPP: Gross primary productivity from nighttime partitioning (umol m-2 s-1) [GPP_NT_VUT_USTAR50]
#'    - GPP_qc: Quality control of `GPP` [NEE_VUT_USTAR50_QC]
#'    - Reco: Ecosystem respiration from nighttime partitioning (umol m-2 s-1) [RECO_NT_VUT_USTAR50]
#'  }
#'  
#' # Note
#' The original variable names as provided by the FLUXNET2015 dataset are 
#'       given in squared brackets. Note that variable units have been converted
#'       in some cases (e.g. VPD from hPa to kPa).       
#'  
#' @source original data were downloaded from
#'         \url{https://fluxnet_fluxdata_org/} (accessed 09 November 2016) 
"FR_Pue_May_2012"