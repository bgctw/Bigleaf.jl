module bigleaf

using DataFrames
using DocStringExtensions

export bigleaf_constants
export Esat_slope, Esat_from_Tair, Esat_from_Tair_deriv,
     LE_to_ET, ET_to_LE, ms_to_mol, mol_to_ms, VPD_to_rH, rH_to_VPD,
     e_to_rH, VPD_to_e, e_to_VPD, e_to_q, q_to_e, q_to_VPD, VPD_to_q,
     Rg_to_PPFD, PPFD_to_Rg, kg_to_mol, umolCO2_to_gC, gC_to_umolCO2
#export air_density, pressure_from_elevation, psychrometric_constant,
#    latent_heat_vaporization, virtual_temp
export latent_heat_vaporization

include("bigleaf_constants.jl")
include("unit_conversions.jl")
include("meteorological_variables.jl")

end
