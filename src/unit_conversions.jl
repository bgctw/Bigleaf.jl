#' Saturation Vapor Pressure (Esat) and Slope of the Esat Curve
#' 
#' Calculates saturation vapor pressure (Esat) over water and the
#'              corresponding slope of the saturation vapor pressure curve_
#' 
#' @param Tair      Air temperature (deg C)
#' @param formula   Formula to be used_ Either \code{:Sonntag_1990} (Default), 
#'                  \code{"Alduchov_1996"}, or \code{"Allen_1998"}_
#' @param constants Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#' 
#' @details Esat (kPa) is calculated using the Magnus equation:
#' 
#'    \deqn{Esat = a * exp((b * Tair) / (c + Tair)) / 1000}
#'  
#'  where the coefficients a, b, c take different values depending on the formula used_
#'  The default values are from Sonntag 1990 (a=611.2, b=17.62, c=243.12)_ This version
#'  of the Magnus equation is recommended by the WMO (WMO 2008; p1.4-29)_ Alternatively,
#'  parameter values determined by Alduchov & Eskridge 1996 or Allen et al_ 1998 can be 
#'  used (see references)_
#'  The slope of the Esat curve (\eqn{\Delta}) is calculated as the first derivative of the function:
#'  
#'    \deqn{\Delta = dEsat / dTair}
#'  
#'  which is solved using \code{\link[stats]{D}}_
#' 
#' @return A dataframe with the following columns: 
#'  \item{Esat}{Saturation vapor pressure (kPa)}
#'  \item{Delta}{Slope of the saturation vapor pressure curve (kPa K-1)}
#'    
#' @references Sonntag D_ 1990: Important new values of the physical constants of 1986, vapor 
#'             pressure formulations based on the ITS-90 and psychrometric formulae_ 
#'             Zeitschrift fuer Meteorologie 70, 340-344_
#'             
#'             World Meteorological Organization 2008: Guide to Meteorological Instruments
#'             and Methods of Observation (WMO-No_8)_ World Meteorological Organization,
#'             Geneva_ 7th Edition_
#'             
#'             Alduchov, O_ A_ & Eskridge, R_ E_, 1996: Improved Magnus form approximation of 
#'             saturation vapor pressure_ Journal of Applied Meteorology, 35, 601-609
#'             
#'             Allen, R_G_, Pereira, L_S_, Raes, D_, Smith, M_, 1998: Crop evapotranspiration -
#'             Guidelines for computing crop water requirements - FAO irrigation and drainage
#'             paper 56, FAO, Rome_
#' 
#' @examples 
#' Esat_slope(seq(0,45,5))[,"Esat"]  # Esat in kPa
#' Esat_slope(seq(0,45,5))[,"Delta"] # the corresponding slope of the Esat curve (Delta) in kPa K-1
#'        
#' @importFrom stats D                  
#' @export
function Esat_slope(Tair::Number; formula=:Sonntag_1990, constants=bigleaf_constants()) 
  Esat = Esat_from_Tair([Tair]; formula, constants)
  Delta = Esat_from_Tair_deriv([Tair]; formula, constants)
  Esat[1], Delta[1]
end

get_EsatCoef(::Val{:Sonntag_1990}) = (a=611.2,b=17.62,c=243.12)
get_EsatCoef(::Val{:Alduchov_1996}) = (a=610.94,b=17.625,c=243.04)
get_EsatCoef(::Val{:Allen_1998}) = (a=610.8,b=17.27,c=237.3)

function Esat_from_Tair(Tair; formula=Val(:Sonntag_1990), constants=bigleaf_constants()) 
  a,b,c = get_EsatCoef(formula)
  Esat = a * exp((b * Tair) / (c + Tair)) * constants[:Pa2kPa]
end

function Esat_from_Tair_deriv(Tair; formula=Val(:Sonntag_1990), constants=bigleaf_constants()) 
  # slope of the saturation vapor pressure curve
  #Delta = eval(D(expression(a * exp((b * Tair) / (c + Tair))),name="Tair"))
  a,b,c = get_EsatCoef(formula)
  #Delta_Pa = @. a*(b / (Tair + c) + (-Tair*b) / ((Tair + c)^2))*exp((Tair*b) / (Tair + c))
  Delta_Pa = @. a * (exp((b * Tair)/(c + Tair)) * (b/(c + Tair) - (b * Tair)/(c + Tair)^2))
  Delta = Delta_Pa .* constants[:Pa2kPa]
end
    

"""
    LE_to_ET(LE,Tair)

    ET_to_LE(ET,Tair)

Converts evaporative water flux from mass (ET=evapotranspiration)
             to energy (LE=latent heat flux) units, or vice versa.

## Arguments
- LE   Latent heat flux (W m-2)
- ET   Evapotranspiration (kg m-2 s-1)
- Tair Air temperature (deg C)

The conversions are given by:
- ``ET = LE/\\lambda``
- ``LE = \\lambda ET``

where ``\\lambda`` is the latent heat of vaporization (J kg-1) as calculated by
[`latent_heat_vaporization`](@ref).

# Examples
```@example
# LE of 200 Wm-2 and air temperature of 25degC
LE_to_ET(200,25)
```

```@example
x = linspace(-π, π) # hide
plot(x, f(x), color = "red")
savefig("f-plot.svg"); 
```
"""
function LE_to_ET(LE,Tair)
  # ![](f-plot.svg)
  lambda = latent_heat_vaporization(Tair)
  ET     = LE/lambda
end,
function ET_to_LE(ET,Tair)
  lambda = latent_heat_vaporization(Tair)
  LE     = ET*lambda
end


# """
# Conversion between Conductance Units

# @description Converts conductances from mass (m s-1)
#              to molar units (mol m-2 s-1), or vice versa_

# @aliases ms_to_mol mol_to_ms

# @param G_ms       Conductance (m s-1)
# @param G_mol      Conductance (mol m-2 s-1)
# @param Tair       Air temperature (deg C)
# @param pressure   Atmospheric pressure (kPa)
# @param constants  Kelvin - conversion degree Celsius to Kelvin \cr
#                   Rgas - universal gas constant (J mol-1 K-1) \cr
#                   kPa2Pa - conversion kilopascal (kPa) to pascal (Pa)

# @details
# The conversions are given by:

# ``G_mol = G_ms * pressure / (Rgas * Tair)``

# ``G_ms = G_mol * (Rgas * Tair) / pressure``

# where Tair is in Kelvin and pressure in Pa (converted from kPa internally)

# @references Jones, H_G_ 1992_ Plants and microclimate: a quantitative approach to environmental plant physiology_
#             2nd Edition_, Cambridge University Press, Cambridge_ 428 p

# @example
# ms_to_mol(0.005,25,100)

# @export
# """
# function ms_to_mol(G_ms,Tair,pressure; constants=bigleaf_constants())
#   Tair     = Tair + constants[:Kelvin] 
#   pressure = pressure * constants[:kPa2Pa] 
#   G_mol  = G_ms * pressure / (constants[:Rgas]  * Tair)
# end
# function mol_to_ms(G_mol,Tair,pressure; constants=bigleaf_constants())
#   Tair     = Tair + constants[:Kelvin] 
#   pressure = pressure * constants[:kPa2Pa] 
#   G_ms  = G_mol * (constants[:Rgas]  * Tair) / (pressure)
# end


# """
# Conversions between Humidity Measures

# @description Conversion between vapor pressure (e), vapor pressure deficit (VPD),
#              specific humidity (q), and relative humidity (rH)_

# @param Tair      Air temperature (deg C)
# @param pressure  Atmospheric pressure (kPa)
# @param e         Vapor pressure (kPa)
# @param q         Specific humidity (kg kg-1)
# @param VPD       Vapor pressure deficit (kPa)
# @param rH        Relative humidity (-)
# @param Esat_formula  Optional: formula to be used for the calculation of esat and the slope of esat_
#                      One of \code{:Sonntag_1990} (Default), \code{"Alduchov_1996"}, or \code{"Allen_1998"}_
#                      See \code{\link{Esat_slope}}_
# @param constants eps - ratio of the molecular weight of water vapor to dry air (-) \cr
#                  Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)

# @family humidity conversion

# @references Foken, T, 2008: Micrometeorology_ Springer, Berlin, Germany_

# @export
# """
function VPD_to_rH(VPD,Tair; Esat_formula=Val(:Sonntag_1990),
                      constants=bigleaf_constants())
  esat = Esat_from_Tair(Tair; formula = Esat_formula,constants)
  rH   = 1 - VPD/esat
end
function rH_to_VPD(rH,Tair; Esat_formula=Val(:Sonntag_1990),
                      constants=bigleaf_constants())
  if !ismissing(rH) && rH > 1 
    @warn("Expected relative humidity (rH) between 0 and 1, but was ", rH)
  end
  esat = Esat_from_Tair(Tair; formula = Esat_formula,constants)
  VPD  = esat - rH*esat
end
function e_to_rH(e,Tair; Esat_formula=Val(:Sonntag_1990),
                    constants=bigleaf_constants())
  esat = Esat_from_Tair(Tair; formula = Esat_formula,constants)
  if !ismissing(e) && (e > esat + sqrt(eps()))
    @warn("Provided vapour pressure that was higher than saturation_
             Returning rH=1 for those cases.")
  end
  rH  = min(1, e/esat)
end
function VPD_to_e(VPD,Tair; Esat_formula=Val(:Sonntag_1990),
                     constants=bigleaf_constants())
  esat = Esat_from_Tair(Tair; formula = Esat_formula,constants)
  e    = esat - VPD
end
function e_to_VPD(e,Tair; Esat_formula=Val(:Sonntag_1990),
                     constants=bigleaf_constants())
  esat = Esat_from_Tair(Tair; formula = Esat_formula,constants)
  VPD  = esat - e
end
function e_to_q(e,pressure; constants=bigleaf_constants())
  q = constants[:eps]  * e / (pressure - (1-constants[:eps] ) * e)
end
function q_to_e(q,pressure; constants=bigleaf_constants())
  e = q * pressure / ((1-constants[:eps] ) * q + constants[:eps] )
end
function q_to_VPD(q,Tair,pressure; Esat_formula=Val(:Sonntag_1990),
                     constants=bigleaf_constants())
  esat = Esat_from_Tair(Tair; formula = Esat_formula,constants)
  e    = q_to_e(q,pressure; constants)
  VPD  = esat - e
end
function VPD_to_q(VPD,Tair,pressure; Esat_formula=Val(:Sonntag_1990),
                     constants=bigleaf_constants())
  esat = Esat_from_Tair(Tair; formula = Esat_formula,constants)
  e    = esat - VPD
  q    = e_to_q(e,pressure;constants)
end

# """
# Conversions between Global Radiation and Photosynthetic Photon Flux Density

# @description Converts radiation from W m-2 to umol m-2 s-1 and vice versa_

# @param Rg       Global radiation = incoming short-wave radiation at the surface (W m-2)
# @param PPFD     Photosynthetic photon flux density (umol m-2 s-1)
# @param J_to_mol Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
# @param frac_PAR Fraction of incoming solar irradiance that is photosynthetically
#                 active radiation (PAR); defaults to 0.5

# @details
# The conversion is given by:

#  ``PPFD = Rg * frac_PAR * J_to_mol}

# by default, the combined conversion factor (\code{frac_PAR * J_to_mol}) is 2.3

# @example
# # convert a measured incoming short-wave radiation of 500 Wm-2 to
# # PPFD in umol m-2 s-1 and backwards
# Rg_to_PPFD(500)
# PPFD_to_Rg(1150)

# @export
# """
function Rg_to_PPFD(Rg,J_to_mol=4.6,frac_PAR=0.5)
  PPFD = Rg * frac_PAR * J_to_mol
end
function PPFD_to_Rg(PPFD,J_to_mol=4.6,frac_PAR=0.5)
  Rg = PPFD / frac_PAR / J_to_mol
end

# """
# Conversion between Mass and Molar Units

# @description Converts mass units of a substance to the corresponding molar units
#              and vice versa_

# @param mass      Numeric vector of mass in kg
# @param molarMass Numeric vector of molar mass of the substance (kg mol-1)
#                  e_g_ as provided by \code{\link{bigleaf_constants}}()$H2Omol
#                  Default is molar mass of Water_

# @return Numeric vector of amount of substance in mol_
# @export
# """
function kg_to_mol(mass, molarMass=bigleaf_constants()[:H2Omol])
  moles = mass / molarMass
end

# """
# Conversion between Mass and Molar Units of Carbon and CO2

# @description Converts CO2 quantities from umol CO2 m-2 s-1 to g C m-2 d-1 and vice versa_

# @param CO2_flux  CO2 flux (umol CO2 m-2 s-1)
# @param C_flux    Carbon (C) flux (gC m-2 d-1)
# @param constants Cmol - molar mass of carbon (kg mol-1) \cr
#                  umol2mol - conversion micromole (umol) to mol (mol) \cr
#                  mol2umol - conversion mole (mol) to micromole (umol)  \cr
#                  kg2g - conversion kilogram (kg) to gram (g) \cr
#                  g2kg - conversion gram (g) to kilogram (kg) \cr
#                  days2seconds - seconds per day

# @example
# umolCO2_to_gC(20)  # gC m-2 d-1

# @export
# """
function umolCO2_to_gC(CO2_flux; constants=bigleaf_constants())
  C_flux = CO2_flux * constants[:umol2mol]  * constants[:Cmol]  * 
  constants[:kg2g]  * constants[:days2seconds] 
end
function gC_to_umolCO2(C_flux; constants=bigleaf_constants())
  CO2_flux = (C_flux * constants[:g2kg]  / constants[:days2seconds] ) / 
  constants[:Cmol]  * constants[:mol2umol] 
end
