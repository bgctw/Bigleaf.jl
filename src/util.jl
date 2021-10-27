# """
#     toDataFrame(x::AbstractVector{<:SLArray}) 

# Convert an Vector{SLVector} to DataFrame.

# Such objects are returned when broadcasting on a function
# that retuns an SLVector.
# """
# function toDataFrame(x::AbstractVector{<:SLArray}) 
#     names = collect(keys(first(x)))
#     df = DataFrame(Tables.table(VectorOfArray(x)'))
#     rename!(df, names)
# end

# function vec_tuple_to_DataFrame(data; names = string.(keys(first(data))))
#     DataFrame(collect(map(idx -> getindex.(data, idx), eachindex(collect(first(data))))), collect(names))
# end

"""
    frac_hour(float::AbstractFloat)
    frac_hour(period::Type{<:Period}, float::AbstractFloat)

Create a period in given type (defaults to `Nanosecond`) from
fractional hours.

```jldoctest; output = false
using Dates
frac_hour(1+1/60) == Hour(1) + Minute(1)
# output
true
```
"""
function frac_hour(period::Type{<:Period}, float::AbstractFloat)
    #adapted from https://stackoverflow.com/a/51448499
    full_hour, Δ = divrem(float, 1)
    partial = period(round(Dates.value(period(Hour(1))) * Δ))
    Hour(full_hour) + partial
end
frac_hour(float::AbstractFloat) = frac_hour(Nanosecond, float)
  

"""
    moving_average(vs,n; nmin = n ÷ 2) 

Compute the moving average over a vector allowing for missings.
        
# Arguments    
- vs: numeric vector to average over
- n: window size: number of items to average
- nmin: minimum number of non-missing records

Values outside the edges are assumed missing.    

If the number of non-missing records within a window is smaller than nmin
then the averaged value is assumed missing. This avoids average at edges of 
periods with many missings to be very sensitive to the edge values.
"""
function moving_average(vs,n; nmin = n ÷ 2) 
    kernel = -trunc(Int,(n-1)/2):ceil(Int,(n-1)/2)
    #[mean(skipmissing(@view vs[i:(i+n-1)])) for i in 1:(length(vs)-(n-1))]    
    vsp = PaddedView(missing, vs, (-n:length(vs)+n,))
    fagg = function(x)
        sum(.!ismissing.(x)) < nmin && return(missing)
        mean(skipmissing(x))
    end
    #ans = [fagg(@view vsp[i.+kernel]) for i in 100:101]    
    ans = [fagg(@view vsp[i.+kernel]) for i in 1:length(vs)]    
end


  
