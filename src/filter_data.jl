"""
    setinvalid_qualityflag!(df; 
      vars=["LE","H","NEE","Tair","VPD","wind"],
      qc_suffix="_qc",
      good_quality_threshold = 1.0,
      missing_qc_as_bad = true,
      setvalmissing = true, 
    )

Set records with quality flags indicating problems to false in :valid column. 

# Arguments
- df: DataFrame with column :GPP
optional
- `vars=["LE","H","NEE","Tair","VPD","wind"]`: columns to theck for quality
- `qc_suffix="_qc"`: naming of the corresponding quality-flag column
- `good_quality_threshold = 1.0`: threshold in quality flag up to which
  data is considered good quality
- `missing_qc_as_bad = true`: set to false to not mark records with missing 
  quality flag as invalid
- `setvalmissing = true`: set to false to prevent replacing values in value column 
  corresponding to problematic quality flag to missing.

# Value
df with modified :valid and value columns. 

# Example
```jldoctest; output = false
using DataFrames
df = DataFrame(
  NEE = 1:3, NEE_qc = [1,1,2],
  GPP = 10:10:30, GPP_qc = [1,missing,1])
setinvalid_qualityflag!(df; vars = ["NEE", "GPP"])
df.valid == [true, false, false]
ismissing(df.GPP[2]) && ismissing(df.NEE[3])
# output
true
```
"""
function setinvalid_qualityflag!(df; setvalmissing = true, kwargs...)
  if setvalmissing
    set_badquality_missing!(df; kwargs...)
  end
  setinvalid_qualityflag!(df, Val(false); kwargs...)
end

function set_badquality_missing!(df;
  vars=["LE","H","NEE","Tair","VPD","wind"],
  qc_suffix="_qc",
  good_quality_threshold = 1.0,
  missing_qc_as_bad = true
  )
  fqc(x,xqc) = if missing_qc_as_bad
     ifelse(!ismissing(xqc) && xqc <= good_quality_threshold, x, missing)
  else
    ifelse(ismissing(xqc) || xqc <= good_quality_threshold, x, missing)
  end
  tmp = map(vars) do var
      qcvar = var * qc_suffix
      [var, qcvar] => ByRow(fqc) => var
  end
  transform!(df, tmp...)
end

function setinvalid_qualityflag!(df, setvalsmissing::Val{false};
  vars=["LE","H","NEE","Tair","VPD","wind"],
  qc_suffix="_qc",
  good_quality_threshold = 1.0,
  missing_qc_as_bad=true
  )
  vars_qc = vars .* qc_suffix
  fvalid_var = missing_qc_as_bad ? valid_nonmissingvar_ : valid_missingvar_
  function fsel(valid, x...) 
     nvar = length(x) ÷ 2
     vs = x[1:nvar]
     vs_qc = x[(nvar+1):end]
     #valid = repeat(@MVector([true]), length(vs[1]))
     for i = 1:nvar
         # works only from 1.7 @. valid = valid && fvalid_var(vs[i], vs_qc[i], good_quality_threshold)
         valid .= valid .&& fvalid_var.(vs[i], vs_qc[i], good_quality_threshold)
     end
     valid
  end
  if !hasproperty(df, :valid) df[!,:valid] .= true; end
  select!(df, :, Cols(vcat("valid", vars, vars_qc)) => fsel => :valid)     
end

function valid_nonmissingvar_(x, xqc, good_quality_threshold) 
  !ismissing(x) && !ismissing(xqc) && xqc <= good_quality_threshold
end
function valid_missingvar_(x, xqc, good_quality_threshold) 
 !ismissing(x) && (ismissing(xqc) || xqc <= good_quality_threshold)
end

"""
    setinvalid_range!(df, var_ranges...; setvalmissing = true, ...)

Set records with values outside specified ranges to false in :valid column. 

If their is no limit towards small or
large values, supply `-Inf` or `Inf` as the minimum or maximum respectively.
If there were false values in the :value column before, they are kept false.
In addition, set values outside ranges to missing.

# Arguments
- `df`: DataFrame with column :GPP
- `var_ranges`: Pair `Varname_symbol => (min,max)`: closed valid interval for 
  respective column 
optional
- setvalmissing: set to false to prevent replacing values in value column outside ranges to missing.

# Value
df with modified :valid and value columns. 
```jldoctest; output = false
using DataFrames
df = DataFrame(NEE = [missing, 2,3], GPP = 10:10:30)
setinvalid_range!(df, :NEE => (-2.0,4.0), :GPP => (8.0,28.0))
df.valid == [false, true, false]
ismissing(df.GPP[3])
# output
true
```
"""
function setinvalid_range!(df, var_ranges::Vararg{Pair,N}; setvalmissing = true, kwargs...) where N
  if setvalmissing
    set_badrange_missing!(df, var_ranges...; kwargs...)
  end
  setinvalid_range!(df, Val(false), var_ranges...; kwargs...)
end

function set_badrange_missing!(df, var_ranges::Vararg{Pair,N}) where N
  tmp = map(var_ranges) do p
      var, (min, max) = p
      var => (x -> @.(ifelse(!ismissing(x) && min <= x <= max, x, missing))) => var
  end
  transform!(df, tmp...)
end

function setinvalid_range!(df, setvalsmissing::Val{false}, var_ranges::Vararg{Pair,N}) where N
  function fval(valid, x...)
    for (p,xi) in zip(var_ranges, x)
      var, (min, max) = p
      @. valid = valid && !ismissing(xi) && (min <= xi <= max)
    end
    valid
  end
  vars = map(x -> x.first, var_ranges)
  if !hasproperty(df, :valid) df[!,:valid] .= true; end
  select!(df, :, SA[:valid, vars...] => fval => :valid) 
end

"""
    setinvalid_nongrowingseason!(df, tGPP; kwargs...)

Set non-growseason to false in :valid column.

# Arguments
- df: DataFrame with column :GPP
- tGPP: scalar threshold of daily GPP (see [`get_growingseason`](@ref))
optional:
- `update_GPPd`: set to true additionally update `:GPPd_smoothed` column to
  results from [`get_growingseason`](@ref)
- and others passed to [`get_growingseason`](@ref)

# Value
df with modified columns :valid and if  `:GPPd_smoothed`, 
where all non-growing season records are set to false.
"""
function setinvalid_nongrowingseason!(df, tGPP; update_GPPd_smoothed = false, kwargs...)
  if !hasproperty(df, :valid) df[!,:valid] .= true; end
  # non-copying dataframe where we can add grouping column __day
  dft = transform(df, :time => ByRow(Date) => :__day, copycols = false)
  function mean_nonmissing(x)
      mx = mean(skipmissing(x))
      # if there is no non-missing record in a group, mean returns nan
      ifelse(isnan(mx), missing, mx)
  end
  dfday = combine(groupby(dft, :__day), :GPP => mean_nonmissing => :GPPd, nrow)
  transform!(dfday, 
    :GPPd => (x -> get_growingseason(x, tGPP; kwargs...)) => SA[:valid, :GPPd_smoothed])
  # rep(dfd.valid, each = df.nrow)
  df[!,:valid] .= df.valid .&& vcat(fill.(dfday.valid, dfday.nrow)...)
  if update_GPPd_smoothed
      df[!,:GPPd_smoothed] .= vcat(fill.(dfday.GPPd_smoothed, dfday.nrow)...)
  end
  df
end

"""
    get_growingseason(GPPd, tGPP; ws=15, min_int=5, warngap=true)

Filters annual time series for growing season based on smoothed daily GPP data.

# Arguments
- `GPPd`:    daily GPP (any unit) 
- `tGPP`:    GPP threshold (fraction of 95th percentile of the GPP time series).
               Takes values between 0 and 1. 
optional               
- `ws`:      window size used for GPP time series smoothing
- `min_int`: minimum time interval in days for a given state of growing season
- `warngap`: set to false to suppress warning on too few non-missing data

# Details
The basic idea behind the growing season filter is that vegetation is 
considered to be active when its carbon uptake (GPP) is above a specified 
threshold, which is defined relative to the peak GPP (95th percentile) 
observed in the year. 
The GPP-threshold is calculated as:

``GPP_{threshold} = quantile(GPPd,0.95)*tGPP``

GPPd time series are smoothed with a moving average to avoid fluctuations 
in the delineation of the growing season. The window size defaults to 15 
days, but depending on the ecosystem, other values can be appropriate. 

The argument `min_int` serves to avoid short fluctuations in the 
status growing season vs. no growing season by defining a minimum length
of the status. If a time interval shorter than `min_int` is labeled
as growing season or non-growing season, it is changed to the status of 
the neighboring values, i.e its opposite.
         
# Value
A NamedTuple with entries
- `is_growingseason`: a `BitVector` of the same length as the input GPPd in which `false` 
  indicate no growing season (dormant season) and `true` indicate growing season.
- `GPPd_smoothed`: smoothed GPPd
"""
function get_growingseason(GPPd,tGPP;ws=15,min_int=5,warngap=true)
  nday = length(GPPd)
  if sum(.!ismissing.(GPPd)) < 0.5*nday 
    error("Expected number of available GPPd data " * 
    "to be at least half the total number of days ($nday). " *
    "But was only ($(sum(.!ismissing.(GPPd)))).")
  end
  GPP_threshold = quantile(skipmissing(GPPd), 0.95)*tGPP
  #@show GPP_threshold
  # smooth GPP
  #fromR: GPPd_smoothed = filter(GPPd,method="convolution",filter=rep(1/ws,ws))
  GPPd_smoothed = moving_average(GPPd, ws)
  # set values at the beginning and end of the time series to the mean of original values
  wsd = floor(Int, ws/2)
  GPPd_smoothed[1:wsd] .= mean(skipmissing(GPPd[1:(2*wsd)]))
  GPPd_smoothed[(length(GPPd)-(wsd-1)):length(GPPd)] .= 
    mean(skipmissing(GPPd[(length(GPPd)-(2*wsd-1)):length(GPPd)]))
  # check for occurence of missing values and set to mean of the values surrounding them
  imissing = findall(ismissing, GPPd_smoothed)
  if length(imissing) > 0
    warngap && length(imissing) > 10 && 
      @warn "Attention, there is a gap in 'GPPd' of length n = $(length(imissing))"
    #TODO check and correct for several gaps, see Impute package
    replace_val = mean(skipmissing(GPPd_smoothed[max(1,imissing[1] - 4):min(
      (imissing[length(imissing)] + 4),length(GPPd_smoothed))]))
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
  (is_growingseason = growseas, GPPd_smoothed = GPPd_smoothed)
end



