#' Eddy Covariance Data of AT-Neu (Neustift)
#' 
#' Halfhourly eddy covariance Data of the site AT-Neu,
#'              a mountain meadow in Austria.
#'              (\url{http://sites_fluxdata_org/AT-Neu}). 
#'              Data are from July 2010.
#'              
#' @format A data frame with 1488 observations and 31 columns:
#'  \describe
#'    \item{year}{year of measurement}
#'    \item{month}{month of measurement}
#'    \item{doy}{day of year}
#'    \item{hour}{hour (0 - 23.5)}
#'    \item{Tair}{Air temperature (degC) [TA_F]}
#'    \item{Tair_qc}{Quality control of \code{Tair} [TA_F_QC]}
#'    \item{PPFD}{Photosynthetic photon flux density (umol m-2 s-1) [PPFD_IN]}
#'    \item{PPFD_qc}{Quality control of \code{PPFD} [PPFD_IN_QC]}
#'    \item{VPD}{Vapor pressure deficit (kPa) [VPD_F]}
#'    \item{VPD_qc}{Quality control of \code{VPD} [VPD_F_QC]}
#'    \item{pressure}{Atmospheric pressure (kPa) [PA_F]}
#'    \item{precip}{precipitation (mm) [P_F]}
#'    \item{precip_qc}{Quality control of \code{precip} [P_F_QC]}
#'    \item{ustar}{friction velocity (m s-1) [USTAR]}
#'    \item{wind}{horizontal wind velocity (m s-1) [WS_F]}
#'    \item{wind_qc}{Quality control of \code{wind} [WS_F_QC]}
#'    \item{Ca}{CO2 concentration (ppm) [CO2_F_MDS]}
#'    \item{Ca_qc}{Quality control of \code{Ca} [CO2_F_MDS_QC]}
#'    \item{LW_up}{upward longwave radiation (W m-2) [LW_OUT]}
#'    \item{Rn}{Net radiation (W m-2) [NETRAD]}
#'    \item{LE}{Latent heat flux (W m-2) [LE_F_MDS]}
#'    \item{LE_qc}{Quality control of \code{LE} [LE_F_MDS_QC]}
#'    \item{H}{Sensible heat flux (W m-2) [H_F_MDS]}
#'    \item{H_qc}{Quality control of \code{H} [H_F_MDS_QC]}
#'    \item{G}{Ground heat flux (W m-2) [G_F_MDS]}
#'    \item{G_qc}{Quality control of \code{G} [G_F_MDS_QC]}
#'    \item{NEE}{Net ecosystem exchange (umol m-2 s-1) [NEE_VUT_USTAR50]}
#'    \item{NEE_qc}{Quality control of \code{NEE} [NEE_VUT_USTAR50_QC]}
#'    \item{GPP}{Gross primary productivity from nighttime partitioning (umol m-2 s-1) [GPP_NT_VUT_USTAR50]}
#'    \item{GPP_qc}{Quality control of \code{GPP} [NEE_VUT_USTAR50_QC]}
#'    \item{Reco}{Ecosystem respiration from nighttime partitioning (umol m-2 s-1) [RECO_NT_VUT_USTAR50]}
#'  }
#'  
#' @note The original variable names as provided by the FLUXNET2015 dataset are 
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
#'    \item{year}{year of measurement}
#'    \item{month}{month of measurement}
#'    \item{doy}{day of year}
#'    \item{hour}{hour (0 - 23.5)}
#'    \item{Tair}{Air temperature (degC) [TA_F]}
#'    \item{Tair_qc}{Quality control of \code{Tair} [TA_F_QC]}
#'    \item{PPFD}{Photosynthetic photon flux density (umol m-2 s-1) [PPFD_IN]}
#'    \item{PPFD_qc}{Quality control of \code{PPFD} [PPFD_IN_QC]}
#'    \item{VPD}{Vapor pressure deficit (kPa) [VPD_F]}
#'    \item{VPD_qc}{Quality control of \code{VPD} [VPD_F_QC]}
#'    \item{pressure}{Atmospheric pressure (kPa) [PA_F]}
#'    \item{precip}{precipitation (mm) [P_F]}
#'    \item{precip_qc}{Quality control of \code{precip} [P_F_QC]}
#'    \item{ustar}{friction velocity (m s-1) [USTAR]}
#'    \item{wind}{horizontal wind velocity (m s-1) [WS_F]}
#'    \item{wind_qc}{Quality control of \code{wind} [WS_F_QC]}
#'    \item{Ca}{CO2 concentration (ppm) [CO2_F_MDS]}
#'    \item{Ca_qc}{Quality control of \code{Ca} [CO2_F_MDS_QC]}
#'    \item{LW_up}{upward longwave radiation (W m-2) [LW_OUT]}
#'    \item{LW_down}{downward longwave radiation (W m-2) [LW_IN_F]}
#'    \item{Rn}{Net radiation (W m-2) [NETRAD]}
#'    \item{LE}{Latent heat flux (W m-2) [LE_F_MDS]}
#'    \item{LE_qc}{Quality control of \code{LE} [LE_F_MDS_QC]}
#'    \item{H}{Sensible heat flux (W m-2) [H_F_MDS]}
#'    \item{H_qc}{Quality control of \code{H} [H_F_MDS_QC]}
#'    \item{G}{Ground heat flux (W m-2) [G_F_MDS]}
#'    \item{G_qc}{Quality control of \code{G} [G_F_MDS_QC]}
#'    \item{NEE}{Net ecosystem exchange (umol m-2 s-1) [NEE_VUT_USTAR50]}
#'    \item{NEE_qc}{Quality control of \code{NEE} [NEE_VUT_USTAR50_QC]}
#'    \item{GPP}{Gross primary productivity from nighttime partitioning (umol m-2 s-1) [GPP_NT_VUT_USTAR50]}
#'    \item{GPP_qc}{Quality control of \code{GPP} [NEE_VUT_USTAR50_QC]}
#'    \item{Reco}{Ecosystem respiration from nighttime partitioning (umol m-2 s-1) [RECO_NT_VUT_USTAR50]}
#'  }
#'  
#' @note The original variable names as provided by the FLUXNET2015 dataset are 
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
#'    \item{year}{year of measurement}
#'    \item{month}{month of measurement}
#'    \item{doy}{day of year}
#'    \item{hour}{hour (0 - 23.5)}
#'    \item{Tair}{Air temperature (degC) [TA_F]}
#'    \item{Tair_qc}{Quality control of \code{Tair} [TA_F_QC]}
#'    \item{PPFD}{Photosynthetic photon flux density (umol m-2 s-1) [PPFD_IN]}
#'    \item{PPFD_qc}{Quality control of \code{PPFD} [PPFD_IN_QC]}
#'    \item{VPD}{Vapor pressure deficit (kPa) [VPD_F]}
#'    \item{VPD_qc}{Quality control of \code{VPD} [VPD_F_QC]}
#'    \item{pressure}{Atmospheric pressure (kPa) [PA_F]}
#'    \item{precip}{precipitation (mm) [P_F]}
#'    \item{precip_qc}{Quality control of \code{precip} [P_F_QC]}
#'    \item{ustar}{friction velocity (m s-1) [USTAR]}
#'    \item{wind}{horizontal wind velocity (m s-1) [WS_F]}
#'    \item{wind_qc}{Quality control of \code{wind} [WS_F_QC]}
#'    \item{Ca}{CO2 concentration (ppm) [CO2_F_MDS]}
#'    \item{Ca_qc}{Quality control of \code{Ca} [CO2_F_MDS_QC]}
#'    \item{LW_up}{upward longwave radiation (W m-2) [LW_OUT]}
#'    \item{Rn}{Net radiation (W m-2) [NETRAD]}
#'    \item{LE}{Latent heat flux (W m-2) [LE_F_MDS]}
#'    \item{LE_qc}{Quality control of \code{LE} [LE_F_MDS_QC]}
#'    \item{H}{Sensible heat flux (W m-2) [H_F_MDS]}
#'    \item{H_qc}{Quality control of \code{H} [H_F_MDS_QC]}
#'    \item{NEE}{Net ecosystem exchange (umol m-2 s-1) [NEE_VUT_USTAR50]}
#'    \item{NEE_qc}{Quality control of \code{NEE} [NEE_VUT_USTAR50_QC]}
#'    \item{GPP}{Gross primary productivity from nighttime partitioning (umol m-2 s-1) [GPP_NT_VUT_USTAR50]}
#'    \item{GPP_qc}{Quality control of \code{GPP} [NEE_VUT_USTAR50_QC]}
#'    \item{Reco}{Ecosystem respiration from nighttime partitioning (umol m-2 s-1) [RECO_NT_VUT_USTAR50]}
#'  }
#'  
#' @note The original variable names as provided by the FLUXNET2015 dataset are 
#'       given in squared brackets. Note that variable units have been converted
#'       in some cases (e.g. VPD from hPa to kPa).       
#'  
#' @source original data were downloaded from
#'         \url{https://fluxnet_fluxdata_org/} (accessed 09 November 2016) 
"FR_Pue_May_2012"