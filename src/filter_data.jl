"""
GPP-based Growing Season Filter

Filters annual time series for growing season based on smoothed daily GPP data.

# Arguments
- GPPd:    daily GPP (any unit) 
- tGPP:    GPP threshold (fraction of 95th percentile of the GPP time series).
               Takes values between 0 and 1. 
- ws:      window size used for GPP time series smoothing
- min_int: minimum time interval in days for a given state of growing season

# Details
The basic idea behind the growing season filter is that vegetation is 
considered to be active when its carbon uptake (GPP) is above a specified 
threshold, which is defined relative to the peak GPP (95th percentile) 
observed in the year. 
The GPP-threshold is calculated as:

``GPP_threshold = quantile(GPPd,0.95)*tGPP``

GPPd time series are smoothed with a moving average to avoid fluctuations 
in the delineation of the growing season. The window size defaults to 15 
days, but depending on the ecosystem, other values can be appropriate. 

The argument `min_int` serves to avoid short fluctuations in the 
status growing season vs. no growing season by defining a minimum length
of the status. If a time interval shorter than `min_int` is labeled
as growing season or non-growing season, it is changed to the status of 
the neighboring values.
         
# Value
a vector of type integer of the same length as the input GPPd in which 0 indicate
        no growing season (dormant season) and 1 indicate growing season.
"""
function filter_growing_season(GPPd,tGPP,ws=15,min_int=5)
  nday = length(GPPd)
  if sum(ismissing(GPPd)) >= 0.5*nday 
    @warn "number of available GPPd data is less than half the total number of days " *
    "per year. Filter is not applied!"
    return(Trues(nday))
  end
  growseas = repeat([true], nday)
  GPP_threshold = quantile(skipmissing(GPPd), 0.95)*tGPP
  
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


