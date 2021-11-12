using Bigleaf

using DataFrames, Pipe, Missings
using Dates, TimeZones
using Statistics

using Latexify
using DataDeps, Suppressor
using RData
import CodecBzip2, CodecXz
#@suppress_err # error in github-actions: GitHubActionsLogger has no field stream
register(DataDep(
    "DE_Tha_Jun_2014.rda",
    "downloading exampple dataset DE_Tha_Jun_2014 from bitbucket.org/juergenknauer/bigleaf",
    "https://bitbucket.org/juergenknauer/bigleaf/raw/0ebe11626b4409305951e8add9f6436703c82584/data/DE_Tha_Jun_2014.rda",
    "395f02e1a1a2d175ac7499c200d9d48b1cb58ff4755dfd2d7fe96fd18258d73c"
))
#println(datadep"DE_Tha_Jun_2014.rda")
ENV["DATADEPS_ALWAYS_ACCEPT"]="true" # avoid question to download
DE_Tha_Jun_2014 = first(values(load(joinpath(datadep"DE_Tha_Jun_2014.rda/DE_Tha_Jun_2014.rda"))))
nothing
tha = DE_Tha_Jun_2014
set_datetime_ydh!(tha)

thaf = copy(tha);
setinvalid_qualityflag!(thaf);
setinvalid_range!(thaf, 
     :PPFD => (200, Inf), 
     :ustar => (0.2, Inf), 
     :LE =>(0, Inf), 
     :VPD => (0.01, Inf)
     );
setinvalid_nongrowingseason!(thaf, 0.4);
setinvalid_afterprecip!(thaf; min_precip=0.02, hours_after=24);

thas = subset(thaf, :valid)



function tmpf()
    Dl=0.01
    aerodynamic_conductance!(thas; Gb_model=Val(:Su_2001), 
        Dl, LAI=thal.LAI, zh=thal.zh, zr=thal.zr);
    surface_conductance!(thas, Val(:PenmanMonteith); G=thas.G);
end

# tha48 and thal see runtests.jl


show(thaf.wind[1:48])
show(thaf.VPD[1:48])
show(thaf.LE[1:48])
show(thaf.Rn[1:48])
show(thaf.G[1:48])


dfGPPd = @pipe tha |> 
    transform(_, :datetime => ByRow(yearmonthday) => :ymd, copycols = false) |>
    groupby(_, :ymd) |>
    combine(_, :datetime => (x -> Date(first(x))) => :date, :GPP => mean => :GPPd, ) 
    #x -> x[!,2]

show(dfGPPd.GPPd)


GPPd = [11.288497760271033, 13.013025930772224, 12.851774960756302, 11.996453734696843, 11.635422044472458, 11.155685572574535, 10.774393322790273, 10.181774605065584, 11.257192575993637, 12.9423423493281, 12.352468963712454, 13.402045020057509, 9.53826415212825, 12.071680051895479, 13.692589149111882, 12.845505638824156, 12.378533909407755, 11.672167064607493, 10.401156075240579, 10.705716138705611, 10.207347450816693, 11.052016352768987, 13.54435911634937, 12.060648361220956, 7.758974596237143, 9.869706534050541, 12.998054057980577, 10.627359564105669, 8.685295419767499, 10.874667977293333]
using Plots,StatsPlots
@df dfGPPd plot(:date, [:GPPd], xlab = "Date", ylab="GPP")
