This vignette is a short introduction to the functionalities of the `Bigleaf.jl` package. 
It is directed to first-time package users who are familiar with the basic concepts of Julia. 
After presenting the use of several key functions of the package, 
some useful hints and guidelines are given at the end of the vignette.


# Package scope and important conceptual considerations

`Bigleaf.jl` calculates physical and physiological ecosystem properties from eddy covariance data. Examples for such properties are aerodynamic and surface conductance, surface conditions(e.g. temperature, VPD), wind profile, roughness parameters, vegetation-atmosphere decoupling, potential evapotranspiration, (intrinsic) water-use efficiency, stomatal sensitivity to VPD, or intercellular CO2 concentration.  All calculations in the `Bigleaf.jl` package assume that the ecosystem behaves like a  "big-leaf", i.e. a single, homogenous plane which acts as the only source and sink of the measured fluxes. This assumption comes with the advantages that calculations are simplified considerably and that (in most cases) little ancillary information on the EC sites is required. It is important to keep in mind that these simplifications go hand in hand with critical limitations. All derived variables are bulk ecosystem characteristics and have to be interpreted as such. It is for example not possible to infer within-canopy variations of a certain property.

Please also keep in mind that the `Bigleaf.jl` package does NOT provide formulations for bottom-up modelling. The principle applied here is to use an inversion approach in which ecosystem properties are inferred top-down from the measured fluxes. Such an inversion can, in principle, be also be conducted with more complex models (e.g. sun/shade or canopy/soil models), but keep in mind that these approaches also require that the additional, site-specific parameters are adequately well known. 

The use of more detailed models is not within the scope of the `Bigleaf.jl` package, but it is preferable to use such approaches when important assumptions of the "big-leaf" approach are not met. This is the case in particular when the ecosystem is sparsely covered with vegetation (low LAI, e.g. sparse crops, some savanna systems). 


# Preparing the data

In this tutorial, we will work with a dataset from the eddy covariance site Tharandt (DE-Tha), a spruce forest in Eastern Germany. The DataFrame `DE_Tha_Jun_2014` is downloaded from the `bigleaf` 
[R package](https://bitbucket.org/juergenknauer/Bigleaf/) repository and contains half-hourly data of meteorological and flux measurements made in June 2014. For loading the RData into Julia, see the 
[source](https://github.com/bgctw/Bigleaf.jl/blob/main/docs/src/walkthrough.md?plain=1#L26) of this file. We give the data.frame a shorter name here.

```@example doc
using Bigleaf
using DataFrames
```
```@setup doc
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
```
```@example doc
tha = DE_Tha_Jun_2014
mdtable(select(describe(tha), :variable, :eltype, :min, :max), latex=false) # hide
```


And the first six rows of tha:
```@example doc
mdtable(tha[1:6,:],latex=false) # hide
```

More information on the data (e.g. meaning of column names and units) can be found at the 
[bigleaf R package](https://bitbucket.org/juergenknauer/bigleaf/src/master/man/DE_Tha_Jun_2014.Rd). 
For more information on the site see e.g. Grünwald & Bernhofer 2007.
In addition, we will need some ancillary data for this site throughout this tutorial. To ensure consistency, we define them here at the beginning:

```@example doc
LAI = 7.6   # leaf area index
zh  = 26.5  # average vegetation height (m)
zr  = 42    # sensor height (m)
Dl  = 0.01  # leaf characteristic dimension (m)
nothing # hide
```

# General guidelines on package usage

There are a few general guidelines that are important to consider when using the `Bigleaf.jl` package. 


## Units

It is imperative that variables are provided in the right units, as the plausibility of 
the input units is not checked in most cases. The required units of the input arguments 
can be found in the respective help file of the function. The good news is that units 
do not change across functions. For example, pressure is always required in kPa, 
and temperature always in °c.

## Function arguments

Most functions of `Bigleaf.jl` require a DataFrame, from which the required
variables are extracted. This is usually the first argument of a function. 
Most functions further provide default values for their arguments, 
such that in many cases it is not necessary to provide them explicitly.

The column names in the DataFrame should correspond to the argument names
of the corresponding method that accespts each input individually.

We can demonstrate the usage with a simple example:

```@example doc
# explicit inputs
Tair, pressure, Rn, =  14.81, 97.71, 778.17 
potential_ET(Tair, pressure, Rn, Val(:PriestleyTaylor))
# DataFrame
potential_ET(tha, Val(:PriestleyTaylor))
# DataFrame with a few columns overwritten by user values
potential_ET(transform(tha, :Tair => x -> 25.0; renamecols=false), Val(:PriestleyTaylor))
# varying one input only
Tair_vec =  10.0:1.0:20.0
DataFrame(potential_ET.(Tair_vec, pressure, Rn, Val(:PriestleyTaylor)))
nothing # hide
```

## Ground heat flux and storage fluxes

Many functions require the available energy ($A$), which is defined as ($A = R_n - G - S$, 
all in $\text{W m}^{-2}$), where $R_n$ is the net radiation, $G$ is the ground heat flux, 
and $S$ is the sum of all storage fluxes of the ecosystem 
(see e.g. Leuning et al. 2012 for an overview). For some sites, $G$ is not available, 
and for most sites, only a few components of $S$ are measured. 

In `Bigleaf.jl` it is not a problem if $G$ and/or $S$ are missing (other than the results might be (slightly) biased), but special options exist for the treatment of missing $S$ and $G$ values. 

Note that the default for G and S in the dataframe variant is missing (and assumed zero), 
even if those columns are
present in the DataFrame. You need to explictly pass those columns with the optional
arguments: e.g. `potential_ET(df, Val(:PriestleyTaylor); G = df.G)`

Note that in difference to the bigleaf R package missing entries in a provide
vector are not relaced by zero by default. 
You need to explitly use coalesce when specifying a ground heat flux
for which missings should be replaced by zero: `;G = coalesce(df.G, zero(df.G))`
 

# Function walkthrough #

## Meteorological variables

The `Bigleaf.jl` package provides calculation routines for a number of meteorological variables, which are basic to the calculation of many other variables. A few examples on their usage are given below:

```@example doc
# Saturation vapor pressure (kPa) and slope of the saturation vapor pressure curve (kPa K-1)
Esat_slope(25.0)
```
```@example doc
# psychrometric constant (kPa K-1)
psychrometric_constant(25.0,100.0) # Tair, pressure
```
```@example doc
# air density (kg m-3)
air_density(25.0,100.0) # Tair, pressure
```
```@example doc
# dew point (degC)
dew_point(25.0,1.0) # Tair, VPD
```
```@example doc
# wetbulb temperature (degC)
wetbulb_temp(25.0, 100.0, 1.0) # Tair, pressure, VPD
```
```@example doc
# estimate atmospheric pressure from elevation (hypsometric equation)
pressure_from_elevation(500.0, 25.0) # elev, Tair
```

There are several formulations describing the empirical function `Esat(Tair)`.
The following figure compares them at absole scale and as difference to the 
#default method. The differences are small.

```@setup doc
#using DataFrames
#Tair = 0:0.25:12
##Tair = [10.0,20.0]
#eform_def = Val(:Sonntag_1990)
#Esat_def = Esat_from_Tair.(Tair; formula = eform_def)
#eforms = (Val(:Sonntag_1990), Val(:Alduchov_1996), Val(:Allen_1998))
#eform = eforms[2]
#string.(eforms)
#df = mapreduce(vcat, eforms) do eform 
#    Esat = Esat_from_Tair.(Tair; formula = eform)
#    local dff # make sure to not override previous results
#    dff = DataFrame(
#        formula = eform, Tair = Tair, 
#        Esat = Esat,
#        dEsat = Esat - Esat_def,
#        )
#end;
##using Chain
#using Pipe
#using Plots, StatsPlots
#dfw = @pipe df |> select(_, 1,2, :Esat) |> unstack(_, :formula, 3)
#dfws = @pipe df |> select(_, 1,2, :dEsat) |> unstack(_, :formula, 3)
#@df dfw plot(:Tair, cols(2:4), legend = :topleft, xlab="Tair (degC)", #ylab="Esat (kPa)")
#savefig("Esat_abs.svg")
#@df dfws plot(:Tair, cols(2:4), legend = :topleft, xlab="Tair (degC)", #ylab="Esat -ESat_Sonntag_1990 (kPa)")
#savefig("fig/Esat_rel.svg")
```

![](fig/Esat_abs.svg)

![](fig/Esat_rel.svg)

## Global radiation

Potential radiation for given time and latitude:
```@example doc
doy, hour = 160, 10.5
lat, long = 51.0, 11.5
potrad = potential_radiation(doy, hour, lat, long)
```

Calculation is based on sun's altitude, one of the horizontal coordinates of its position.
```@example doc
using Plots, StatsPlots, DataFrames, Dates, Pipe, Suppressor
hours = 0:24
lat,long = 51.0, 13.6 # Dresden Germany
#deg2second = 24*3600/360
doy = 160
datetimes = DateTime(2021) .+Day(doy-1) .+ Hour.(hours) #.- Second(round(long*deg2second))
res3 = @pipe calc_sun_position_hor.(datetimes, lat, long) |> DataFrame(_)
@df res3 scatter(datetimes, cols([:altitude,:azimuth]), legend = :topleft, xlab="Date and Time", ylab = "rad")
```

The hour-angle at noon represents the difference to
local time. In the following example solar time is
about 55min ahead of local winter time.

```@example doc
summernoon = DateTime(2021) +Day(doy-1) + Hour(12) 
sunpos = calc_sun_position_hor(summernoon, lat, long) 
sunpos.hourangle * 24*60/(2*π) # convert angle to minutes
```

## Unit interconversions

The package further provides a number of useful unit interconversions, which are straightforward to use (please make sure that the input variable is in the right unit, e_g. rH has to be between 0 and 1 and not in percent):

```@example doc
# VPD to vapor pressure (e, kPa)
VPD_to_e(2, 25)
```
```@example doc
# vapor pressure to specific humidity (kg kg-1)
e_to_q(1, 100)
```
```@example doc
# relative humidity to VPD (kPa)
rH_to_VPD(0.6, 25)
```
```@example doc
# conductance from ms-1 to mol m-2 s-1
ms_to_mol(0.01, 25, 100) # mC, Tair, pressure
```
```@example doc
# umol CO2 m-2 s-1 to g C m-2 d-1
umolCO2_to_gC(20)
```

Many functions provide constant empirical parameters. Those can
be changed by overriding the default values with 
[`bigleaf_constants`](@ref) 
and passing this Dictionary to the respective function.


