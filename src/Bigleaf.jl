module Bigleaf

#using DocStringExtensions
using Optim
using FillArrays
using DataFrames
using Tables: columntable
using Dates, TimeZones
using Pipe
using AstroLib
using Suppressor
using Missings
using Statistics, StatsBase # mean, rle
using PaddedViews, StaticArrays
using Infiltrator
using Parameters
using SymbolDispatch

#export ConductanceMethod, Thom1972, Choudhury1988, Su2001, ConstantKB1
export ResistanceMethod, ResistanceWindZr, ResistanceWindProfile
#export RoughnessMethod, Roughness_wind_profile, RoughnessCanopyHeight, RoughnessCanopyHeightLAI
export ETMethod, PriestleyTaylor, PenmanMonteith
export EsatMethod, Sonntag1990, Alduchov1996, Allen1998
export StabilityCorrectionMethod, Dyer1970, Businger1971, NoStabilityCorrection
export SurfaceConductanceMethod, FluxGradient, InversePenmanMonteith

export frac_hour, moving_average, get_nonoverlapping_periods, set_datetime_ydh!
export BigleafConstants
export Esat_slope, Esat_from_Tair, Esat_from_Tair_deriv,
     LE_to_ET, ET_to_LE, ms_to_mol, mol_to_ms, VPD_to_rH, rH_to_VPD,
     e_to_rH, VPD_to_e, e_to_VPD, e_to_q, q_to_e, q_to_VPD, VPD_to_q,
     Rg_to_PPFD, PPFD_to_Rg, kg_to_mol, umolCO2_to_gC, gC_to_umolCO2
export air_density, pressure_from_elevation, psychrometric_constant,
    latent_heat_vaporization, virtual_temp, kinematic_viscosity,
    dew_point_from_e, dew_point,
    wetbulb_temp_from_e_Tair_gamma, wetbulb_temp
export calc_sun_position_MOD, calc_sun_position_hor
export potential_radiation, extraterrestrial_radiation, get_datetime_for_doy_hour
export potential_ET, potential_ET!, equilibrium_imposed_ET, equilibrium_imposed_ET!
export setinvalid_range!, setinvalid_qualityflag!, 
    setinvalid_nongrowingseason!, get_growingseason, setinvalid_afterprecip!
export decoupling, surface_conductance, aerodynamic_conductance
export compute_Gb, compute_Gb!, add_Gb, add_Gb!, 
    Gb_Thom, Gb_Choudhury, Gb_Su, Gb_constant_kB1,
    compute_Gb_quantities, compute_Gb_quantities!
export wind_profile
export Monin_Obukhov_length, Monin_Obukhov_length!, stability_parameter, 
    stability_parameter!, stability_correction, stability_correction!,
    roughness_parameters, Reynolds_Number
export aerodynamic_conductance!, add_Ga, add_Ga!, compute_Ram, compute_Ram!, roughness_length_heat
export surface_conductance, surface_conductance!

# shorthand type: either DataFrame or column table (Tables.columntable)
DFTable = Union{NamedTuple,AbstractDataFrame}

include("util.jl")    
include("bigleaf_constants.jl")
include("unit_conversions.jl")
include("meteorological_variables.jl")
include("sun_position.jl")
include("potential_radiation.jl")
include("evapotranspiration.jl")
include("filter_data.jl")
include("stability_correction.jl")
include("surface_roughness.jl")
include("boundary_layer_conductance.jl")
include("aerodynamic_conductance.jl")
include("surface_conductance.jl")

# testing @symboldispatch
export foo
@symboldispatch foo(::Val{:wind}) = :wind
foo(::Val{:Wutzler2022}) = :Wutzler2022 

end
