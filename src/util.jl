export toDataFrame

function toDataFrame(x::AbstractVector{<:SLArray}) 
    names = collect(keys(first(x)))
    df = DataFrame(Tables.table(VectorOfArray(x)'))
    rename!(df, names)
end

function vec_tuple_to_DataFrame(data; names = string.(keys(first(data))))
    DataFrame(collect(map(idx -> getindex.(data, idx), eachindex(collect(first(data))))), collect(names))
end

