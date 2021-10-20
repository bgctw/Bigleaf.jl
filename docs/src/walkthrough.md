This vignette is a short introduction to the functionalities of the `bigleaf.jl` package. 
It is directed to first-time package users who are familiar with the basic concepts of Julia. 
After presenting the use of several key functions of the package, 
some useful hints and guidelines are given at the end of the vignette.


# Package scope and important conceptual considerations

`bigleaf.jl` calculates physical and physiological ecosystem properties from eddy covariance data. Examples for such properties are aerodynamic and surface conductance, surface conditions(e.g. temperature, VPD), wind profile, roughness parameters, vegetation-atmosphere decoupling, potential evapotranspiration, (intrinsic) water-use efficiency, stomatal sensitivity to VPD, or intercellular CO2 concentration.  All calculations in the `bigleaf.jl` package assume that the ecosystem behaves like a  "big-leaf", i.e. a single, homogenous plane which acts as the only source and sink of the measured fluxes. This assumption comes with the advantages that calculations are simplified considerably and that (in most cases) little ancillary information on the EC sites is required. It is important to keep in mind that these simplifications go hand in hand with critical limitations. All derived variables are bulk ecosystem characteristics and have to be interpreted as such. It is for example not possible to infer within-canopy variations of a certain property.

Please also keep in mind that the `bigleaf.jl` package does NOT provide formulations for bottom-up modelling. The principle applied here is to use an inversion approach in which ecosystem properties are inferred top-down from the measured fluxes. Such an inversion can, in principle, be also be conducted with more complex models (e.g. sun/shade or canopy/soil models), but keep in mind that these approaches also require that the additional, site-specific parameters are adequately well known. 

The use of more detailed models is not within the scope of the `bigleaf.jl` package, but it is preferable to use such approaches when important assumptions of the "big-leaf" approach are not met. This is the case in particular when the ecosystem is sparsely covered with vegetation (low LAI, e.g. sparse crops, some savanna systems). 


# Installation and Loading

The `bigleaf.jl` R package can be installed with the usual command once:

```julia
using Pkg
Pkg.add(bigleaf)
```

And then importet to the every Julia session by:
```julia
using bigleaf
```


# Preparing the data

In this tutorial, we will work with a dataset from the eddy covariance site Tharandt (DE-Tha), a spruce forest in Eastern Germany. The DataFrame `DE_Tha_Jun_2014` is downloaded from the `bigleaf` 
[R package]() repository and contains half-hourly data of meteorological and flux measurements made in June 2014:

```@setup doc
using Latexify, DataFrames
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
nothing
```
```@example doc
# downloading and caching code not shown
tha = DE_Tha_Jun_2014
mdtable(select(describe(tha), :variable, :eltype, :min, :max), latex=false) # hide
```

And the first six rows of tha:
```@example doc
mdtable(tha[1:6,:],latex=false) # hide
```

We give the data.frame a shorter name here. More information on the data (e.g. meaning of column names and units) can be found at the 
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

There are a few general guidelines that are important to consider when using the `bigleaf.jl` package. 


## Units

It is imperative that variables are provided in the right units, as the plausibility of the input units is not checked in most cases. The required units of the input arguments can be found in the respective help file of the function. The good news is that units do not change across functions. For example, pressure is always required in kPa, and temperature always in °c.

