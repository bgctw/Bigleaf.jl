"""
    filter_growing_season(GPPd, tGPP; ws=15, min_int=5, warngap=true)

Filters annual time series for growing season based on smoothed daily GPP data.

# Arguments
- GPPd:    daily GPP (any unit) 
- tGPP:    GPP threshold (fraction of 95th percentile of the GPP time series).
               Takes values between 0 and 1. 
optional               
- ws:      window size used for GPP time series smoothing
- min_int: minimum time interval in days for a given state of growing season
- warngap: set to false to suppress warning on too few non-missing data

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
the neighboring values, i.e its opposite.
         
# Value
a `BitVector` of the same length as the input GPPd in which `false` indicate
no growing season (dormant season) and `true` indicate growing season.
"""
function filter_growing_season(GPPd,tGPP;ws=15,min_int=5,warngap=true)
  nday = length(GPPd)
  if sum(.!ismissing.(GPPd)) < 0.5*nday 
    error("Expected number of available GPPd data " * 
    "to be at least half the total number of days ($nday). " *
    "But was only ($(sum(.!ismissing.(GPPd)))).")
  end
  GPP_threshold = quantile(skipmissing(GPPd), 0.95)*tGPP
  # smooth GPP
  #fromR: GPPd_smoothed = filter(GPPd,method="convolution",filter=rep(1/ws,ws))
  GPPd_smoothed = moving_average(GPPd, ws)
  # set values at the beginning and end of the time series to the mean of the original values
  wsd = floor(Int, ws/2)
  GPPd_smoothed[1:wsd] .= mean(skipmissing(GPPd[1:(2*wsd)]))
  GPPd_smoothed[(length(GPPd)-(wsd-1)):length(GPPd)] .= mean(skipmissing(GPPd[(length(GPPd)-(2*wsd-1)):length(GPPd)]))
  # check for occurence of missing values and set them to mean of the values surrounding them
  imissing = findall(ismissing, GPPd_smoothed)
  if length(imissing) > 0
    warngap && length(imissing) > 10 && @warn "Attention, there is a gap in 'GPPd' of length n = $(length(imissing))"
    #TODO check and correct for several gaps, see Impute package
    replace_val = mean(skipmissing(GPPd_smoothed[max(1,imissing[1] - 4):min((imissing[length(imissing)] + 4),length(GPPd_smoothed))]))
    GPPd_smoothed = coalesce.(GPPd_smoothed, replace_val)
  end
  # filter daily GPP
  growseas = GPPd_smoothed .>= GPP_threshold
  # change short intervals to the surrounding values to avoid 'wrong' fluctuations
  # switch shortest interval successively and recompute interval lengths
  intervals = rle(growseas)[2]
  imin = argmin(intervals)
  while intervals[imin] < min_int
      end_int = cumsum(intervals)
      start_int = vcat(0, end_int[1:end-1]) .+ 1
      growseas[start_int[imin]:end_int[imin]] .= !growseas[start_int[imin]]
      intervals = rle(growseas)[2]
      imin = argmin(intervals)
  end
  return(growseas)
end

