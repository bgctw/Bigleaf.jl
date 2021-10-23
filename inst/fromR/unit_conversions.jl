#######################
## Unit Conversions ###
#######################

#' Conversion between Latent Heat Flux and Evapotranspiration
#'
#' converts evaporative water flux from mass (ET=evapotranspiration)
#'              to energy (LE=latent heat flux) units, or vice versa.
#'
#' @aliases LE_to_ET ET_to_LE
#'
#' - LE   Latent heat flux (W m-2)
#' - ET   Evapotranspiration (kg m-2 s-1)
#' - Tair Air temperature (deg C)
#'
#' # Details

#' The conversions are given by:
#'
#' \deqn{ET = LE/\lambda}
#'
#' \deqn{LE = \lambda ET}
#'
#' where \eqn{\lambda} is the latent heat of vaporization (J kg-1) as calculated by
#' `\link{latent_heat_vaporization`}.
#'
#' ```@example; output = false
#' ```
#' # LE of 200 Wm-2 and air temperature of 25degC
#' LE_to_ET(200,25)
#'
"""
"""
function LE_to_ET(LE,Tair)

  lambda = latent_heat_vaporization(Tair)
  ET     = LE/lambda

  return(ET)
end


#' @rdname LE_to_ET
"""
"""
function ET_to_LE(ET,Tair)

  lambda = latent_heat_vaporization(Tair)
  LE     = ET*lambda

  return(LE)
end



#' Conversion between Conductance Units
#'
#' Converts conductances from mass (m s-1)
#'              to molar units (mol m-2 s-1), or vice versa.
#'
#' @aliases ms_to_mol mol_to_ms
#'
#' - G_ms       Conductance (m s-1)
#' - G_mol      Conductance (mol m-2 s-1)
#' - Tair       Air temperature (deg C)
#' - pressure   Atmospheric pressure (kPa)
#' - constants  Kelvin - conversion degree Celsius to Kelvin \cr
#'                   Rgas - universal gas constant (J mol-1 K-1) \cr
#'                   kPa2Pa - conversion kilopascal (kPa) to pascal (Pa)
#'
#' # Details

#' The conversions are given by:
#'
#' \deqn{G_mol = G_ms * pressure / (Rgas * Tair)}
#'
#' \deqn{G_ms = G_mol * (Rgas * Tair) / pressure}
#'
#' where Tair is in Kelvin and pressure in Pa (converted from kPa internally)
#'
#' @references Jones, H_G. 1992. Plants and microclimate: a quantitative approach to environmental plant physiology.
#'             2nd Edition., Cambridge University Press, Cambridge. 428 p
#'
#' ```@example; output = false
#' ```
#' ms_to_mol(0.005,25,100)
#'
"""
"""
function ms_to_mol(G_ms,Tair,pressure,constants=bigleaf_constants())

  Tair     = Tair + constants[:Kelvin]
  pressure = pressure * constants[:kPa2Pa]

  G_mol  = G_ms * pressure / (constants[:Rgas] * Tair)

  return(G_mol)
end


#' @rdname ms_to_mol
"""
"""
function mol_to_ms(G_mol,Tair,pressure,constants=bigleaf_constants())

  Tair     = Tair + constants[:Kelvin]
  pressure = pressure * constants[:kPa2Pa]

  G_ms  = G_mol * (constants[:Rgas] * Tair) / (pressure)

  return(G_ms)
end



#' Conversions between Humidity Measures
#'
#' Conversion between vapor pressure (e), vapor pressure deficit (VPD),
#'              specific humidity (q), and relative humidity (rH).
#'
#' - Tair      Air temperature (deg C)
#' - pressure  Atmospheric pressure (kPa)
#' - e         Vapor pressure (kPa)
#' - q         Specific humidity (kg kg-1)
#' - VPD       Vapor pressure deficit (kPa)
#' - rH        Relative humidity (-)
#' - Esat_formula  Optional: formula to be used for the calculation of esat and the slope of esat.
#'                      One of `"Sonntag_1990"` (Default), `"Alduchov_1996"`, or `"Allen_1998"`.
#'                      See `\link{Esat_slope`}.
#' - constants eps - ratio of the molecular weight of water vapor to dry air (-) \cr
#'                  Pa2kPa - conversion pascal (Pa) to kilopascal (kPa)
#'
#' @family humidity conversion
#'
#' @references Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany.
#'
"""
"""
function VPD_to_rH(VPD,Tair,Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                      constants=bigleaf_constants())

  esat = Esat_slope(Tair,Esat_formula,constants)[,"Esat"]
  rH   = 1 - VPD/esat
  return(rH)
end


#' @rdname VPD_to_rH
#' @family humidity conversion
"""
"""
function rH_to_VPD(rH,Tair,Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                      constants=bigleaf_constants())

  if(any(rH > 1 & !is_na(rH)))
    warning("relative humidity (rH) has to be between 0 and 1.")
end
  
  esat = Esat_slope(Tair,Esat_formula,constants)[,"Esat"]
  VPD  = esat - rH*esat
  return(VPD)
end

#' @rdname VPD_to_rH
#' @family humidity conversion
"""
"""
function e_to_rH(e,Tair,Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                    constants=bigleaf_constants())
  
  esat = Esat_slope(Tair,Esat_formula,constants)[,"Esat"]
  
  if (any(e > esat + .Machine$double_eps^0.5 & !is_na(e)))
    warning("Provided vapour pressure that was higher than saturation.
             Returning rH=1 for those cases.")
end
    
  rH  = pmin(1, e/esat)
  return(rH)
end


#' @rdname VPD_to_rH
#' @family humidity conversion
"""
"""
function VPD_to_e(VPD,Tair,Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                     constants=bigleaf_constants())

  esat = Esat_slope(Tair,Esat_formula,constants)[,"Esat"]
  e    = esat - VPD
  return(e)
end


#' @rdname VPD_to_rH
#' @family humidity conversion
"""
"""
function e_to_VPD(e,Tair,Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                     constants=bigleaf_constants())

  esat = Esat_slope(Tair,Esat_formula,constants)[,"Esat"]
  VPD  = esat - e
  return(VPD)
end


#' @rdname VPD_to_rH
#' @family humidity conversion
"""
"""
function e_to_q(e,pressure,constants=bigleaf_constants())
  q = constants[:eps] * e / (pressure - (1-constants[:eps]) * e)
  return(q)
end


#' @rdname VPD_to_rH
#' @family humidity conversion
"""
"""
function q_to_e(q,pressure,constants=bigleaf_constants())
  e = q * pressure / ((1-constants[:eps]) * q + constants[:eps])
  return(e)
end


#' @rdname VPD_to_rH
#' @family humidity conversion
"""
"""
function q_to_VPD(q,Tair,pressure,Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                     constants=bigleaf_constants())

  esat = Esat_slope(Tair,Esat_formula,constants)[,"Esat"]
  e    = q_to_e(q,pressure,constants)
  VPD  = esat - e
  return(VPD)
end


#' @rdname VPD_to_rH
#' @family humidity conversion
"""
"""
function VPD_to_q(VPD,Tair,pressure,Esat_formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                     constants=bigleaf_constants())

  esat = Esat_slope(Tair,Esat_formula,constants)[,"Esat"]
  e    = esat - VPD
  q    = e_to_q(e,pressure,constants)
  return(q)
end





#' Conversions between Global Radiation and Photosynthetic Photon Flux Density
#'
#' Converts radiation from W m-2 to umol m-2 s-1 and vice versa.
#'
#' - Rg       Global radiation = incoming short-wave radiation at the surface (W m-2)
#' - PPFD     Photosynthetic photon flux density (umol m-2 s-1)
#' - J_to_mol Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
#' - frac_PAR Fraction of incoming solar irradiance that is photosynthetically
#'                 active radiation (PAR); defaults to 0.5
#'
#' # Details

#' The conversion is given by:
#'
#'  \deqn{PPFD = Rg * frac_PAR * J_to_mol}
#'
#' by default, the combined conversion factor (`frac_PAR * J_to_mol`) is 2.3
#'
#' ```@example; output = false
#' ```
#' # convert a measured incoming short-wave radiation of 500 Wm-2 to
#' # PPFD in umol m-2 s-1 and backwards
#' Rg_to_PPFD(500)
#' PPFD_to_Rg(1150)
#'
"""
"""
function Rg_to_PPFD(Rg,J_to_mol=4.6,frac_PAR=0.5)
  PPFD = Rg * frac_PAR * J_to_mol
  return(PPFD)
end


#' @rdname Rg_to_PPFD
"""
"""
function PPFD_to_Rg(PPFD,J_to_mol=4.6,frac_PAR=0.5)
  Rg = PPFD / frac_PAR / J_to_mol
  return(Rg)
end




#' Conversion between Mass and Molar Units
#' 
#' Converts mass units of a substance to the corresponding molar units
#'              and vice versa.
#'
#' - mass      Numeric vector of mass in kg
#' - molarMass Numeric vector of molar mass of the substance (kg mol-1)
#'                  e.g. as provided by `\link{bigleaf_constants`}()$H2Omol
#'                  Default is molar mass of Water.
#'
#' # Value
 Numeric vector of amount of substance in mol.
"""
"""
function kg_to_mol(mass, molarMass=bigleaf_constants()$H2Omol)
  
  moles = mass / molarMass
  
  return(moles)

end




#' Conversion between Mass and Molar Units of Carbon and CO2
#'
#' Converts CO2 quantities from umol CO2 m-2 s-1 to g C m-2 d-1 and vice versa.
#'
#' - CO2_flux  CO2 flux (umol CO2 m-2 s-1)
#' - C_flux    Carbon (C) flux (gC m-2 d-1)
#' - constants Cmol - molar mass of carbon (kg mol-1) \cr
#'                  umol2mol - conversion micromole (umol) to mol (mol) \cr
#'                  mol2umol - conversion mole (mol) to micromole (umol)  \cr
#'                  kg2g - conversion kilogram (kg) to gram (g) \cr
#'                  g2kg - conversion gram (g) to kilogram (kg) \cr
#'                  days2seconds - seconds per day
#'
#' ```@example; output = false
#' ```
#' umolCO2_to_gC(20)  # gC m-2 d-1
#'
"""
"""
function umolCO2_to_gC(CO2_flux,constants=bigleaf_constants())

  C_flux = CO2_flux * constants[:umol2mol] * constants[:Cmol] * constants[:kg2g] * constants[:days2seconds]

  return(C_flux)
end




#' @rdname umolCO2_to_gC
"""
"""
function gC_to_umolCO2(C_flux,constants=bigleaf_constants())

  CO2_flux = (C_flux * constants[:g2kg] / constants[:days2seconds]) / constants[:Cmol] * constants[:mol2umol]

  return(CO2_flux)
end