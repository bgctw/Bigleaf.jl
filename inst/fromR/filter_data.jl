########################
### Filter functions ###----------------------------------------------------------------------------
########################

#' Basic Eddy Covariance Data Filtering
#'
#' Filters time series of EC data for high-quality values and specified
#'              meteorological conditions.
#' 
#' - data            Data_frame or matrix containing all required input variables in 
#'                        half-hourly or hourly resolution. Including year, month, day information
#' - quality_control Should quality control be applied? Defaults to `TRUE`.
#' - filter_growseas Should data be filtered for growing season? Defaults to `FALSE`.
#' - filter_precip   Should precipitation filtering be applied? Defaults to `FALSE`.
#' - filter_vars     Additional variables to be filtered. Vector of type character.
#' - filter_vals_min Minimum values of the variables to be filtered. Numeric vector of 
#'                        the same length than `filter_vars`. Set to `NA` to be ignored.
#' - filter_vals_max Maximum values of the variables to be filtered. Numeric vector of 
#'                        the same length than `filter_vars`. Set to `NA` to be ignored.
#' - NA_as_invalid   If `TRUE` (the default) missing data are filtered out (applies to all variables).
#' - vars_qc         Character vector indicating the variables for which quality filter should 
#'                        be applied. Ignored if `quality_control = FALSE`.
#' - quality_ext     The extension to the variables' names that marks them as 
#'                        quality control variables. Ignored if `quality_control = FALSE`.                       
#' - good_quality    Which values indicate good quality (i.e. not to be filtered) 
#'                        in the quality control (qc) variables? Ignored if `quality_control = FALSE`.
#' - missing_qc_as_bad If quality control variable is `NA`, should the corresponding data point be
#'                          treated as bad quality? Defaults to `TRUE`. Ignored if `quality_control = FALSE`.                        
#' - precip          Precipitation (mm time-1)
#' - GPP             Gross primary productivity (umol m-2 s-1); Ignored if `filter_growseas = FALSE`.
#' - doy             Day of year; Ignored if `filter_growseas = FALSE`.
#' - year            Year; Ignored if `filter_growseas = FALSE`.
#' - tGPP            GPP threshold (fraction of 95th percentile of the GPP time series).
#'                        Must be between 0 and 1. Ignored if `filter_growseas` is `FALSE`.
#' - ws              Window size used for GPP time series smoothing. 
#'                        Ignored if `filter_growseas = FALSE`.
#' - min_int         Minimum time interval in days for a given state of growing season.
#'                        Ignored if `filter_growseas = FALSE`.
#' - tprecip         Precipitation threshold used to identify a precipitation event (mm). 
#'                        Ignored if `filter_precip = FALSE`.
#' - precip_hours    Number of hours removed following a precipitation event (h).
#'                        Ignored if `filter_precip = FALSE`.
#' - records_per_hour Number of observations per hour. I_e. 2 for half-hourly data.
#' - filtered_data_to_NA Logical. If `TRUE` (the default), all variables in the input
#'                              DataFrame/matrix are set to `NA` for the time step where ANY of the
#'                              `filter_vars` were beyond their acceptable range (as
#'                              determined by `filter_vals_min` and `filter_vals_max`).
#'                              If `FALSE`, values are not filtered, and an additional column 'valid'
#'                              is added to the DataFrame/matrix, indicating if any value of a row
#'                              did (1) or did not fulfill the filter criteria (0).
#' - constants frac2percent - conversion between fraction and percent
#' 
#' # Details
 This routine consists of two parts:
#' 
#'          1) Quality control: All variables included in `vars_qc` are filtered for 
#'             good quality data. For these variables a corresponding quality variable with 
#'             the same name as the variable plus the extension as specified in `quality_ext`
#'             must be provided. For time steps where the value of the quality indicator is not included
#'             in the argument `good_quality`, i.e. the quality is not considered as 'good', 
#'             its value is set to `NA`.
#'             
#'          2) Meteorological filtering. Under certain conditions (e.g. low ustar), the assumptions
#'             of the EC method are not fulfilled. Further, some data analysis require certain meteorological
#'             conditions, such as periods without rainfall, or active vegetation (growing season, daytime).
#'             The filter applied in this second step serves to exclude time periods that do not fulfill the criteria
#'             specified in the arguments. More specifically, time periods where one of the variables is higher
#'             or lower than the specified thresholds (`filter_vals_min` and `filter_vals_max`)
#'             are set to `NA` for all variables. If a threshold is set to `NA`, it will be ignored.
#'          
#' # Value
 If `filtered_data_to_NA = TRUE` (default), the input DataFrame/matrix with 
#'         observations which did not fulfill the filter criteria set to `NA`. 
#'         If `filtered_data_to_NA = FALSE`, the input DataFrame/matrix with an additional 
#'         column "valid", which indicates whether all the data of a time step fulfill the 
#'         filtering criteria (1) or not (0).
#'         
#' @note The thresholds set with `filter_vals_min` and `filter_vals_max` filter all data
#'       that are smaller than ("<"), or greater than (">") the specified thresholds. That means
#'       if a variable has exactly the same value as the threshold, it will not be filtered. Likewise,
#'       `tprecip` filters all data that are greater than `tprecip`. 
#' 
#'       Variables considered of bad quality (as specified by the corresponding quality control variables)      
#'       will be set to `NA` by this routine. Data that do not fulfill the filtering criteria are set to
#'       `NA` if `filtered_data_to_NA = TRUE`. Note that with this option *all* variables of the same
#'       time step are set to `NA`. Alternatively, if `filtered_data_to_NA = FALSE` data are not set to `NA`,
#'       and a new column "valid" is added to the DataFrame/matrix, indicating if any value of a row
#'       did (1) or did not fulfill the filter criteria (0).
#'       
#' 
#' ```@example; output = false
#' ``` 
#' # Example of data filtering; data are for a month within the growing season,
#' # hence growing season is not filtered.
#' # If filtered_data_to_NA=TRUE, all values of a row are set to NA if one filter
#' # variable is beyond its bounds. 
#' DE_Tha_Jun_2014_2 = filter_data(DE_Tha_Jun_2014,quality_control=FALSE,
#'                                  vars_qc=c("Tair","precip","H","LE"),
#'                                  filter_growseas=FALSE,filter_precip=TRUE,
#'                                  filter_vars=c("Tair","PPFD","ustar"),
#'                                  filter_vals_min=c(5,200,0.2),
#'                                  filter_vals_max=c(NA,NA,NA),NA_as_invalid=TRUE,
#'                                  quality_ext="_qc",good_quality=c(0,1),
#'                                  missing_qc_as_bad=TRUE,GPP="GPP",doy="doy",
#'                                  year="year",tGPP=0.5,ws=15,min_int=5,precip="precip",
#'                                  tprecip=0.1,precip_hours=24,records_per_hour=2,
#'                                  filtered_data_to_NA=TRUE)
#'
#'  ## same, but with filtered_data_to_NA=FALSE
#'  DE_Tha_Jun_2014_3 = filter_data(DE_Tha_Jun_2014,quality_control=FALSE,
#'                                  vars_qc=c("Tair","precip","H","LE"),
#'                                  filter_growseas=FALSE,filter_precip=TRUE,
#'                                  filter_vars=c("Tair","PPFD","ustar"),
#'                                  filter_vals_min=c(5,200,0.2),
#'                                  filter_vals_max=c(NA,NA,NA),NA_as_invalid=TRUE,
#'                                  quality_ext="_qc",good_quality=c(0,1),
#'                                  missing_qc_as_bad=TRUE,GPP="GPP",doy="doy",
#'                                  year="year",tGPP=0.5,ws=15,min_int=5,precip="precip",
#'                                  tprecip=0.1,precip_hours=24,records_per_hour=2,
#'                                  filtered_data_to_NA=FALSE)
#'                                  
#'  # note the additional column 'valid' in DE_Tha_Jun_2014_3.
#'  # To remove time steps marked as filtered out (i.e. 0 values in column 'valid'):
#'  DE_Tha_Jun_2014_3[DE_Tha_Jun_2014_3["valid"] == 0,] = NA
#'   
#'   
#' @importFrom stats aggregate
#' @export                     
function filter_data(data,quality_control=TRUE,filter_growseas=FALSE,
                        filter_precip=FALSE,filter_vars=NULL,
                        filter_vals_min,filter_vals_max,NA_as_invalid=TRUE,
                        vars_qc=NULL,quality_ext="_qc",good_quality=c(0,1),
                        missing_qc_as_bad=TRUE,GPP="GPP",doy="doy",
                        year="year",tGPP=0.5,ws=15,min_int=5,precip="precip",
                        tprecip=0.01,precip_hours=24,records_per_hour=2,
                        filtered_data_to_NA=TRUE,constants=bigleaf_constants())
  
  
  ### I) Quality control filter
  if (quality_control)
    
    if (is_null(vars_qc))
      stop("quality_control (qc) is TRUE, but no qc variables are provided!")
end
    
    if (any(!vars_qc %in% colnames(data)))
      
      missing_vars = vars_qc[which(!vars_qc %in% colnames(data))]
      stop(paste("Variable ",missing_vars," is included in 'vars_qc', but does not exist in the input data!"))
      
end
    
    vars_qc_qc = paste0(vars_qc,quality_ext)
    if (any(!vars_qc_qc %in% colnames(data)))
      
      missing_vars_qc = vars_qc_qc[which(!vars_qc_qc %in% colnames(data))]
      missing_vars2   = substr(missing_vars_qc,1,nchar(missing_vars_qc) - nchar(quality_ext))
      stop(paste("Quality control for variable ",missing_vars2,"(",missing_vars_qc,") does not exist in the input data!")) 
end
    
    ## data quality
    cat("Quality control:",fill=TRUE)
    for (var in vars_qc)
      var_qc = paste0(var,quality_ext)
      check_input(data,var)
      check_input(data,var_qc)

      if (missing_qc_as_bad)
        data[get(paste0(var,quality_ext)) > max(good_quality) | is_na(get(paste0(var,quality_ext))),var] = NA   # exclude bad quality data or those where qc flag is not available
        qc_invalid      = sum(get(paste0(var,quality_ext)) > max(good_quality) | is_na(get(paste0(var,quality_ext)))) # count & report
else { # same, but consider missing quality flag variables as good
        data[get(paste0(var,quality_ext)) > max(good_quality) & !is_na(get(paste0(var,quality_ext))),var] = NA
        qc_invalid      = sum(get(paste0(var,quality_ext)) > max(good_quality) & !is_na(get(paste0(var,quality_ext))))
end
      
      qc_invalid_perc = round((qc_invalid/nrow(data))*constants$frac2percent,2)
      
      cat(var,": ",qc_invalid," data points (",qc_invalid_perc,"%) set to NA",fill=TRUE,sep="")
end
end
  
  
  #### II) Data filter
  valid = rep(1L,nrow(data))
  
  # 1) GPP
  growseas_invalid = numeric()
  if(filter_growseas)
    check_input(data,doy,year,GPP)
    date             = strptime(paste0(year,"-",doy),format="%Y-%j")
    GPP_daily        = aggregate(GPP,by=list(strftime(date)),mean,na_rm=TRUE)[,2]
    growing_season   = filter_growing_season(GPP_daily,tGPP=tGPP,ws=ws,min_int=min_int)
    growseas_invalid = which(sapply(growing_season,rep,48) == 0)
end

  # 2) precipitation
  precip_invalid = numeric()
  if (filter_precip)
    check_input(data,precip)
    if (NA_as_invalid)
      precip_events = which(precip > tprecip | is_na(precip))
else 
      precip_events = which(precip > tprecip)
end
    precip_invalid = unique(as_numeric(unlist(sapply(precip_events, function(x) x:(min(x+precip_hours*records_per_hour,nrow(data),na_rm=TRUE))))))
end

  # 3) all other filter variables (as defined in filter_vars)
  invalids = list(growseas_invalid,precip_invalid)
  
  if (!is_null(filter_vars))
    for (var in filter_vars)
      v  = which(filter_vars == var)
      vf = v + 2
      check_input(data,var)
      if (NA_as_invalid)
        invalids[[vf]] = which(get(var) < filter_vals_min[v] | get(var) > filter_vals_max[v] | is_na(get(var)))
else 
        invalids[[vf]] = which(get(var) < filter_vals_min[v] | get(var) > filter_vals_max[v] & !is_na(get(var)))
end
end
end
  
  # 4) calculate number and percentage of filtered values
  invalids_perc = sapply(invalids, function(x) round((length(x)/nrow(data))*constants$frac2percent,2))
  
  additional_invalids = sapply(2:length(invalids), function(x) 
    length(setdiff(invalids[[x]],unique(unlist(invalids[1:(x-1)])))))
  
  additional_invalids_perc = round(additional_invalids/nrow(data)*constants$frac2percent,2)
  
  
  # 5) write to output
  if (filter_growseas | filter_precip | length(filter_vars) > 0)
    
    var_names = c("growing season","precipitation",filter_vars)
    
    if (quality_control)
      cat("-------------------------------------------------------------------",fill=TRUE)
end
      
    cat("Data filtering:",fill=TRUE)
  
    cat(length(growseas_invalid)," data points (",invalids_perc[1],"%) excluded by growing season filter",fill=TRUE,sep="")
  
    invisible(sapply(c(1:(length(invalids)-1)), function(x) cat(additional_invalids[x]," additional data points (",
                                                                additional_invalids_perc[x],"%) excluded by ",var_names[x+1],
                                                                " filter (",length(unlist(invalids[x+1]))," data points = ",
                                                                invalids_perc[x+1]," % in total)",fill=TRUE,sep="")))
  
  
    invalid        = unique(unlist(invalids))
    valid[invalid] = 0
  
    excl_perc = round((length(invalid)/nrow(data))*constants$frac2percent,2)
  
    cat(length(invalid)," data points (",excl_perc,"%) excluded in total",fill=TRUE,sep="")
    cat(nrow(data) - length(invalid)," valid data points (",constants$frac2percent-excl_perc,"%) remaining.",fill=TRUE,sep="")
  
  
    # 6) return input data frame with filtered time steps set to NA or an additional 'valid' column
    if (filtered_data_to_NA)
      data_filtered = data
      data_filtered[valid < 1,] = NA
else 
      data_filtered = DataFrame(data,valid)
end
  
else 
    
    data_filtered = data
    
end
  
  return(data_filtered)
end





#' GPP-based Growing Season Filter
#' 
#' Filters annual time series for growing season based on smoothed daily GPP data.
#' 
#' - GPPd    daily GPP (any unit) 
#' - tGPP    GPP threshold (fraction of 95th percentile of the GPP time series).
#'                Takes values between 0 and 1. 
#' - ws      window size used for GPP time series smoothing
#' - min_int minimum time interval in days for a given state of growing season
#' 
#' # Details
 The basic idea behind the growing season filter is that vegetation is 
#'          considered to be active when its carbon uptake (GPP) is above a specified 
#'          threshold, which is defined relative to the peak GPP (95th percentile) 
#'          observed in the year. 
#'          The GPP-threshold is calculated as:
#'          
#'          \deqn{GPP_threshold = quantile(GPPd,0.95)*tGPP}
#'          
#'          GPPd time series are smoothed with a moving average to avoid fluctuations 
#'          in the delineation of the growing season. The window size defaults to 15 
#'          days, but depending on the ecosystem, other values can be appropriate. 
#'          
#'          The argument `min_int` serves to avoid short fluctuations in the 
#'          status growing season vs. no growing season by defining a minimum length
#'          of the status. If a time interval shorter than `min_int` is labeled
#'          as growing season or non-growing season, it is changed to the status of 
#'          the neighboring values.
#'          
#' # Value
 a vector of type integer of the same length as the input GPPd in which 0 indicate
#'         no growing season (dormant season) and 1 indicate growing season.
#'                 
#' @importFrom stats quantile filter                                 
#' @export  
function filter_growing_season(GPPd,tGPP,ws=15,min_int=5)
  
  if(sum(is_na(GPPd)) < 0.5*length(GPPd))
    
    growseas      = rep(1,length(GPPd))
    GPP_threshold = quantile(GPPd,probs=0.95,na_rm=TRUE)*tGPP
    
    ## smooth GPP
    GPPd_smoothed = filter(GPPd,method="convolution",filter=rep(1/ws,ws))
    
    ## set values at the beginning and end of the time series to the mean of the original values
    wsd = floor(ws/2)
    GPPd_smoothed[1:wsd] = mean(GPPd[1:(2*wsd)],na_rm=TRUE)
    GPPd_smoothed[(length(GPPd)-(wsd-1)):length(GPPd)] = mean(GPPd[(length(GPPd)-(2*wsd-1)):length(GPPd)],na_rm=TRUE)
    
    # check for occurence of missing values and set them to mean of the values surrounding them
    missing = which(is_na(GPPd_smoothed))
    if (length(missing) > 0)
      if (length(missing) > 10){warning("Attention, there is a gap in 'GPPd' of length n = ",length(missing))}
      replace_val = mean(GPPd_smoothed[max(1,missing[1] - 4):min((missing[length(missing)] + 4),length(GPPd_smoothed))],na_rm=TRUE)
      GPPd_smoothed[missing] = replace_val
end
    
    # filter daily GPP
    growseas[GPPd_smoothed < GPP_threshold] = 0
    
    ## change short intervals to the surrounding values to avoid 'wrong' fluctuations
    intervals = rle(growseas)
    short_int = which(intervals$lengths <= min_int)
    
    if (length(short_int) > 0)
      start = numeric()
      end   = numeric()
      
      for (i in 1:length(short_int))
        
        start[i] = sum(intervals$lengths[1:short_int[i]-1]) + 1
        end[i]   = start[i]+intervals$lengths[short_int[i]] - 1
        
        val = unique(growseas[start[i]:end[i]])
        
        if (val == 0 & growseas[start[i]-1] == 1)
          growseas[start[i]:end[i]] = 1   
else if (val == 1 & growseas[start[i]-1] == 0)
          growseas[start[i]:end[i]] = 0
end
end
end
    
    growseas = as_integer(growseas)
    
else 
    
    warning("number of available GPPd data is less than half the total number of days per year. Filter is not applied!")
    growseas = as_integer(rep(1,length(GPPd)))
    
end
  
  return(growseas)
end
