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

