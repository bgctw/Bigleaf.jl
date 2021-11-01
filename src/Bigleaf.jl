module Bigleaf

#using DocStringExtensions
using Optim
using FillArrays
using DataFrames
using Dates, TimeZones
using Pipe
using AstroLib
using Suppressor
using Missings
using Statistics, StatsBase # mean, rle
using PaddedViews, StaticArrays
using Infiltrator

export frac_hour, moving_average, get_nonoverlapping_periods, set_datetime_ydh!
export bigleaf_constants
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

include("util.jl")    
include("bigleaf_constants.jl")
include("unit_conversions.jl")
include("meteorological_variables.jl")
include("sun_position.jl")
include("potential_radiation.jl")
include("evapotranspiration.jl")
include("filter_data.jl")

end
