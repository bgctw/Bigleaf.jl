using DataFrames
Tair = 0:0.25:12
#Tair = [10.0,20.0]
eform_def = Val(:Sonntag_1990)
Esat_def = Esat_from_Tair.(Tair; formula = eform_def)
eforms = (Val(:Sonntag_1990), Val(:Alduchov_1996), Val(:Allen_1998))
eform = eforms[2]
string.(eforms)
df = mapreduce(vcat, eforms) do eform 
    Esat = Esat_from_Tair.(Tair; formula = eform)
    local dff # make sure to not override previous results
    dff = DataFrame(
        formula = eform, Tair = Tair, 
        Esat = Esat,
        dEsat = Esat - Esat_def,
        )
end;
#using Chain
using Pipe
using Plots, StatsPlots
dfw = @pipe df |> select(_, 1,2, :Esat) |> unstack(_, :formula, 3)
dfws = @pipe df |> select(_, 1,2, :dEsat) |> unstack(_, :formula, 3)
@df dfw plot(:Tair, cols(2:4), legend = :topleft, xlab="Tair (degC)", ylab="Esat (kPa)")
savefig("Esat_abs.svg")
@df dfws plot(:Tair, cols(2:4), legend = :topleft, xlab="Tair (degC)", ylab="Esat -ESat_Sonntag_1990 (kPa)")
savefig("Esat_rel.svg")
