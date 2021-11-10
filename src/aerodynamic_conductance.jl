"""
    aerodynamic_conductance!(df; 
      Gb_model = Val(:Thom_1972), Ram_model = Val(:wind_zr),
      zr=nothing,zh=nothing, d = isnothing(zh) ? nothing : 0.7*zh,
      ...
      )

Bulk aerodynamic conductance, including options for the boundary layer conductance
formulation and stability correction functions.

# Arguments
- `df`: DataFrame with columns
  - `ustar`     : Friction velocity (m s-1)
  - `wind`      : Wind speed at sensor height (m s-1)
- `Gb_model`  : model for computing boundary layer conductance (see [`compute_Gb!`](@ref))
- `Ram_model` : model for computing aerodynamic resistance (see [`compute_Ram`](@ref))
- `zh`        : canopy height (m)
- `zr`        : Instrument (reference) height (m)

Further required columns of `df` and keyword argument depend on `Gb_model`  
(see [`compute_Gb!`](@ref)) and `Ram_model` (see [`compute_Ram`](@ref)).

If only columns `ustar` and `wind` are available, use default models 
(`Val(:Thom_1972)` and `Val(:wind_zr)`).

# Details

Aerodynamic conductance for heat (Ga_h) is calculated as:

``Ga_h = 1 / (Ra_m + Rb_h)``
  
where ``Ra_m`` is the aerodynamic resistance for momentum and ``Rb_h = 1/Gb_h`` the 
(quasi-laminar) canopy boundary layer resistance ('excess resistance') for heat.

`Ra_m` is computed and described with [`compute_Ram`](@ref) using model `Ram_model`.

`Rb_h` is computed and described with 1/[`compute_Gb!`](@ref) using a given `Gb_model`.


# Value
combined results of [`compute_Gb!`](@ref) and 
- `Ra_m`: Aerodynamic resistance for momentum transfer (s m-1)
- `Ga_m`: Aerodynamic conductance for momentum transfer (m s-1)
- `Ga_h`: Aerodynamic conductance for heat transfer (m s-1)
- `Ra_h`: Aerodynamic resistance for heat transfer (s m-1)
- `Ga_CO2`: Aerodynamic conductance for CO2 transfer (m s-1)
       
# Note
The roughness length for water and heat (z0h) can be computed by [`roughness_z0h`](@ref).

TODO check
Input variables such as LAI, Dl, or zh can be either constants, or
vary with time, i.e. are vectors of the same length as `df`.

Note that boundary layer conductance to water vapor transfer (`Gb_w`) is often 
assumed to equal `Gb_h`. This assumption is also made in `Bigleaf.jl`, for 
example in the function [`surface_conductance`](@ref).

If the roughness length for momentum (`z0m`) is not provided as input, it is estimated 
using [`roughness_parameters`](@ref), which estimates a single `z0m` 
value for the entire time period. If a varying `z0m` value 
(e.g. across seasons or years) is required, `z0m` should be provided as input argument.

# Examples
```jldoctest; output = false
using DataFrames
df = DataFrame(Tair=25,pressure=100,wind=[3,4,5],
  ustar=[0.5,0.6,0.65],H=[200,230,250])   
# simple calculation of Ga  
aerodynamic_conductance!(df;Gb_model=Val(:Thom_1972)) 
# calculation of Ram using a model derived from the logarithmic wind profile
aerodynamic_conductance!(df;Gb_model=Val(:Thom_1972),Ram_model = Val(:wind_profile), 
  zr=40,zh=25,d=17.5,z0m=2) 
# simple calculation of Ga, but a physically based canopy boundary layer model
aerodynamic_conductance!(df,Gb_model=Val(:Su_2001),
  zr=40,zh=25,d=17.5,Dl=0.05,N=2,fc=0.8)
all(isfinite.(df.psi_h))
# output
true
``` 
"""
function aerodynamic_conductance!(df; Gb_model = Val(:Thom_1972), Ram_model = Val(:wind_zr),
  zr=nothing,zh=nothing, d = isnothing(zh) ? nothing : 0.7*zh ,
  z0m=nothing,Dl=nothing,N=2,fc=nothing,LAI=nothing,Cd=0.2,hs=0.01,
  leafwidth=nothing,
  stab_formulation=Val(:Dyer_1970),
  kB_h=nothing,constants=bigleaf_constants()
  )
  # add zeta
  if !isnothing(zr) && !isnothing(d) && !(stab_formulation isa Val{:no_stability_correction})
    stability_parameter!(df::AbstractDataFrame; z=zr, d, constants)
  else
    df[!,:zeta] .= missing
  end
  # adds columns psi_m and psi_h, (:no_stability_correction: just add 0s without using zeta)
  stability_correction!(df; zeta=df.zeta, stab_formulation, constants)
  # pre-estimate z0m to use it in both Gb and Ga
  needs_windprofile = 
    Gb_model isa Union{Val{:Choudhury_1988}, Val{:Su_2001}} || 
    Ram_model isa Val{:wind_profile}
  if needs_windprofile
    if isnothing(z0m) 
      z0m = roughness_parameters(Val(:wind_profile), df; zh, zr, psi_m = df.psi_m).z0m
    end
    wind_zh = wind_profile(zh, df, d, z0m; zh, zr, stab_formulation, constants)
  end
  #
  # calculate canopy boundary layer conductance (Gb)
  Gb_model isa Val{:Thom_1972} && compute_Gb!(df, Gb_model; constants)
  Gb_model isa Val{:constant_kB1} && compute_Gb!(df, Gb_model; kB_h, constants)
  Gb_model isa Val{:Choudhury_1988} && compute_Gb!(
    df, Gb_model; leafwidth, LAI, wind_zh, constants)
  Gb_model isa Val{:Su_2001} && compute_Gb!(
    df, Gb_model; wind_zh, Dl, fc, N, Cd, hs, LAI, constants)
  compute_Gb_quantities!(df)
  #
  # calculate aerodynamic risistance for momentum (Ra_m)
  Ram_model isa Val{:wind_profile} && compute_Ram!(df, Ram_model; zr, d, z0m, constants) 
  Ram_model isa Val{:wind_zr} && compute_Ram!(df, Ram_model) 
  function ft(Ra_m, Gb_h, Gb_CO2) 
    Ga_m = 1/Ra_m
    Ra_h = Ra_m + 1/Gb_h
    Ga_h = 1/Ra_h
    Ga_CO2 =  1/(Ra_m + 1/Gb_CO2)
    (;Ga_m, Ga_h, Ra_h, Ga_CO2)
  end
  transform!(df, [:Ra_m, :Gb_h, :Gb_CO2] => ByRow(ft) => AsTable)
end


"""
    roughness_z0h(z0m, kB_h)

# Arguments
- `z0m` : Roughness length for momentum (m). Can be calculated 
          by [`roughness_parameters`](@ref).
- `kB_h` : kB-1 parameter, Output of [`aerodynamic_conductance!`](@ref)          

# Details
The roughness length for water and heat (z0h) is calculated from the 
relationship (e.g. Verma 1989):
  
``{k_B}_h = ln(z_{0m}/z_{0h})`` 

it follows:

``z_{0h} = z_{0m} / e^{k_{B_h}}``

# References
Verma, S., 1989: Aerodynamic resistances to transfers of heat, mass and momentum.
  In: Estimation of areal evapotranspiration, IAHS Pub, 177, 13-20.
"""
roughness_z0h(z0m, kB_h) = z0m / exp(kB_h)


"""
    compute_Ram(::Val{:wind_profile}, ustar; 
      zr, d, z0m, psi_h, constants=bigleaf_constants())
    compute_Ram!(df, method::Val{:wind_profile};  
      zr, d, z0m, psi_h = df.psi_h, kwargs...)

    compute_Ram(::Val{:wind_zr}, ustar, wind)
    compute_Ram!(df, method::Val{:wind_zr}; kwargs...)

Estimate bulk aerodynamic conductance.

# Arguments
- `ustar`          : Friction velocity (m s-1)
- `wind`           : wind speed at measurement height (m s-1)
- `df`             : DataFrame with above columns
- `zr`             : Instrument (reference) height (m)
- `d`              : Zero-plane displacement height (-), can be estimated using 
                      [`roughness_parameters`](@ref)
- `z0m`            : Roughness length for momentum (m). Can be estimated using 
                     from [`roughness_parameters`](@ref)  
- `psi_h`          : the value of the stability function for heat and water vapor (-)
                     see [`stability_correction`](@ref)

# Details
 
The aerodynamic resistance for momentum ``R_{a_m}`` is given by (`Ram_method = Val(:wind_zr)`):

``R_{a_m} = u/{u^*}^2``

Where u is the horizontal wind velocity. 
Note that this formulation accounts for changes in atmospheric stability, and does not 
require an additional stability correction function. 

An alternative method to calculate ``Ra_m`` is provided (`Ram_method = Val(:wind_profile)`):

``R_{a_m} = (ln((z_r - d)/z_{0m}) - \\psi_h) / (k \\, u^*)``

If the roughness parameters `z0m` and `d` are unknown, they can be estimated using
[`roughness_parameters`](@ref). The argument `stab_formulation` determines the stability 
correction function used to account for the effect of atmospheric stability on `Ra_m` 
(`Ra_m` is lower for unstable and higher for stable stratification). Stratification is based 
on a stability parameter `zeta` ``\\zeta=(z-d/L)``, where `z` is the height, 
`d` the zero-plane displacement height, and `L` the Monin-Obukhov length, calculated with 
[`Monin_Obukhov_length`](@ref)
The stability correction function is chosen by the argument `stab_formulation`. 
Options are `Val(:Dyer_1970)` and `Val(:Businger_1971)` and `Val(:no_stability_correction)`.

# Note
For adding aerodynamic conductance for other species see [`add_Ga!`](@ref).

# Value
Aerodynamic resistance for momentum transfer (s m-1) (``Ra_m``)
       
# References
- Verma, S., 1989: Aerodynamic resistances to transfers of heat, mass and momentum.
  In: Estimation of areal evapotranspiration, IAHS Pub, 177, 13-20.
- Verhoef, A., De Bruin, H., Van Den Hurk, B., 1997: Some practical notes on the parameter kB-1
  for sparse vegetation. Journal of Applied Meteorology, 36, 560-572.
- Hicks, B_B., Baldocchi, D_D., Meyers, T_P., Hosker, J_R., Matt, D_R., 1987:
  A preliminary multiple resistance routine for deriving dry deposition velocities
  from measured quantities. Water, Air, and Soil Pollution 36, 311-330.
- Monteith, J_L., Unsworth, M_H., 2008: Principles of environmental physics.
  Third Edition. Elsevier Academic Press, Burlington, USA. 

# See also
[`aerodynamic_conductance!`](@ref), [`add_Ga!`](@ref)
"""
function compute_Ram(::Val{:wind_profile}, ustar::Union{Missing,Number}; 
  zr, d, z0m, psi_h, constants=bigleaf_constants()
  )
  Ra_m = max((log((zr - d)/z0m) - psi_h),0) / (constants[:k]*ustar)
end
function compute_Ram!(df, method::Val{:wind_profile};  
  zr, d, z0m, psi_h = df.psi_h, kwargs...
  )
  # put keyword arguments to positional arguments for proper broadcast
  ftpos(ustar, zr, d, z0m, psi_h) = 
    compute_Ram(method, ustar; zr, d, z0m, psi_h, kwargs...)
  ft(ustar) = ftpos.(ustar, zr, d, z0m, psi_h)
  transform!(df, :ustar => ft => :Ra_m)
end

function compute_Ram(::Val{:wind_zr}, ustar::Union{Missing,Number}, wind)
  Ra_m = wind / ustar^2
end
function compute_Ram!(df, method::Val{:wind_zr}; kwargs...)
  ft(ustar, wind) = compute_Ram(method, ustar, wind; kwargs...)
  transform!(df, SA[:ustar, :wind] => ByRow(ft) => :Ra_m)
end


"""
    add_Ga(Gb_h, Ga_m, Sc::Vararg{Pair,N}; constants)
    add_Ga!(df::AbstractDataFrame, Sc; Gb_h = df.Gb_h, Ga_m = df.Ga.m, kwargs...) 

compute additional aerodynamic conductance quantities for given Schmidt-numbers

# Arguments  
- `Gb_h`     : Boundary layer conductance for heat transfer (m s-1)
- `Ga_m`     : Aerodynamic conductance for momentum (m s-1)
- `Sc`       : several `Pair{Symbol,Number}` Output name and Schmidt number of 
               additional conductances to be calculated
- `df`       : DataFrame to add output columns               
optional
- `constants=`[`bigleaf_constants`](@ref)`()`: Dictionary with entries 
  - `Pr` - Prandtl number 

# Details
Aerodynamic conductance is calculated as

``G_{a_x} = 1/(1/G_{a_m} + 1/G_{b_x})``

where `Gb_x` is the Boundary layer conductance for other quantities x is calculated 
based on boundary layer for heat transfer, Schmidt-Number, and  Prantl number,
as documented in [`add_Gb!`](@ref).

# Value
a NameTuple or `df` with keys `Ga_x` where `x` are the keys in `Sc` and 
corresponding aerodynamic conductances (m s-1).

# Examples
```jldoctest; output=false
using DataFrames
df = DataFrame(Gb_h=[0.02, missing, 0.055], Ga_m = 1 ./ [0.03, 0.03, 0.03])
add_Ga!(df, :O2 => 0.84, :CH4 => 0.99)
propertynames(df)[3:4] == [:Ga_O2, :Ga_CH4]
# output
true
```
"""
function add_Ga!(df::AbstractDataFrame, Sc::Vararg{Pair,N}; 
  Gb_h = df.Gb_h, Ga_m = df.Ga_m, kwargs...) where N
  N == 0 && return(df)
  Scn, Scv = get_names_and_values("Ga_", Sc...)
  ft() = add_Ga_.(Gb_h, Ga_m, Ref(Scn), Ref(Scv); kwargs...)
  transform!(df, [] => ft => AsTable)
end
function add_Ga(Gb_h::Union{Missing,Number}, Ga_m::Union{Missing,Number}, 
  Sc::Vararg{Pair,N}; kwargs...) where N
  Scn, Scv = get_names_and_values("Ga_", Sc...)
  add_Ga_(Gb_h, Ga_m, Scn, Scv; kwargs...)
end
function add_Ga_(Gb_h::Union{Missing,Number}, Ga_m::Union{Missing,Number}, 
  Scn::NTuple{N,Symbol}, Scv::NTuple{N};  constants=bigleaf_constants()) where N 
  Gbx = add_Gb_(Gb_h, Scn, Scv; constants)
  Gaxv = ntuple(i -> 1/(1/Ga_m + 1/Gbx[i]), length(Gbx))
  Gax = NamedTuple{Scn}(Gaxv)
end




