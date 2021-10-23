export toDataFrame

"""
    toDataFrame(x::AbstractVector{<:SLArray}) 

Convert an Vector{SLVector} to DataFrame.

Such objects are returned when broadcasting on a function
that retuns an SLVector.
"""
function toDataFrame(x::AbstractVector{<:SLArray}) 
    names = collect(keys(first(x)))
    df = DataFrame(Tables.table(VectorOfArray(x)'))
    rename!(df, names)
end

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
  
  
