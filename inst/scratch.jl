function op(a, b, ::Type{Mul})
    return a*b
end

using Symbolics
@variables a b c Tair
Esat = a * exp((b * Tair) / (c + Tair))
Symbolics.derivative(Esat, Tair; simplify = true)

# provide example dataset as Parquet and linked local
# look at generated docu

using DataDeps
using RData
import CodecBzip2, CodecXz
register(DataDep(
    "DE_Tha_Jun_2014.rda",
    "downloading exampple dataset DE_Tha_Jun_2014 from bitbucket.org/juergenknauer/bigleaf",
    "https://bitbucket.org/juergenknauer/bigleaf/raw/0ebe11626b4409305951e8add9f6436703c82584/data/DE_Tha_Jun_2014.rda",
    "395f02e1a1a2d175ac7499c200d9d48b1cb58ff4755dfd2d7fe96fd18258d73c"
))
#println(datadep"DE_Tha_Jun_2014.rda")
ENV["DATADEPS_ALWAYS_ACCEPT"]="true" # avoid question to download
DE_Tha_Jun_2014 = first(values(load(joinpath(datadep"DE_Tha_Jun_2014.rda/DE_Tha_Jun_2014.rda"))))
