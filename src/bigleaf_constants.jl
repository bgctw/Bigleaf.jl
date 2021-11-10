"""
  bigleaf_constants(;...)

Constants Used in the Bigleaf.jl Package

Return a Dictionary Dict{Symbol, Real}.
Default values can be overridden by the named arguments:

- cp           Specific heat of air for constant pressure (J K-1 kg-1)
- Rgas         Universal gas constant (J mol-1 K-1)
- Rv           Gas constant of water vapor (J kg-1 K-1) (Stull 1988 p.641)
- Rd           Gas constant of dry air (J kg-1 K-1) (Foken p. 245)
- Md           Molar mass of dry air (kg mol-1)
- Mw           Molar mass of water vapor (kg mol-1)
- eps          Ratio of the molecular weight of water vapor to dry air (=Mw/Md)
- g            Gravitational acceleration (m s-2)
- solar_constant Solar constant (W m-2)
- pressure0    Reference atmospheric pressure at sea level (Pa)
- Tair0        Reference air temperature (K)
- k            von Karman constant
- Cmol         Molar mass of carbon (kg mol-1)
- Omol         Molar mass of oxygen (kg mol-1)
- H2Omol       Molar mass of water (kg mol-1)
- sigma        Stefan-Boltzmann constant (W m-2 K-4)
- Pr           Prandtl number
- Sc_CO2       Schmidt number for CO2
- Kelvin       Conversion degree Celsius to Kelvin
- DwDc         Ratio of the molecular diffusivities for water vapor and CO2
- days2seconds Seconds per day
- kPa2Pa       Conversion kilopascal (kPa) to pascal (Pa)
- Pa2kPa       Conversion pascal (Pa) to kilopascal (kPa)
- umol2mol     Conversion micromole (umol) to mole (mol)
- mol2umol     Conversion mole (mol) to micromole (umol)
- kg2g         Conversion kilogram (kg) to gram (g)
- g2kg         Conversion gram (g) to kilogram (kg)
- kJ2J         Conversion kilojoule (kJ) to joule (J)
- J2kJ         Conversion joule (J) to kilojoule (kJ)
- se_median    Conversion standard error (SE) of the mean to SE of the median
- frac2percent Conversion between fraction and percent
"""
function bigleaf_constants(;
  # make sure to provide Float64 values (rather than Int) so that return is type-stable
  ## Physical constants
  cp         = 1004.834,        # specific heat of air for constant pressure (J K-1 kg-1)
  Rgas       = 8.31451,         # universal gas constant (J mol-1 K-1)
  Rv         = 461.5,           # gas constant of water vapor (J kg-1 K-1) (Stull 1988 p.641)
  Rd         = 287.0586,        # gas constant of dry air (J kg-1 K-1) (Foken 2008 p. 245)
  Md         = 0.0289645,       # molar mass of dry air (kg mol-1)
  Mw         = 0.0180153,       # molar mass of water vapor (kg mol-1)
  eps        = 0.622,           # ratio of the molecular weight of water vapor to dry air (=Mw/Md)
  g          = 9.81,            # gravitational acceleration (m s-2)
  solar_constant = 1366.1,      # solar constant, i.e. solar radation at earth distance from the sun (W m-2)
  pressure0  = 101325.0,        # reference atmospheric pressure at sea level (Pa)
  Tair0      = 273.15,          # reference air temperature (K)
  k          = 0.41,            # von Karman constant
  Cmol       = 0.012011,        # molar mass of carbon (kg mol-1)
  Omol       = 0.0159994,       # molar mass of oxygen (kg mol-1)
  H2Omol     = 0.01801528,      # molar mass of water (kg mol-1)
  sigma      = 5.670367e-08,    # Stefan-Boltzmann constant (W m-2 K-4)
  Pr         = 0.71,            # Prandtl number
  Sc_CO2     = 1.07,            # Schmidt number for CO2 (Hicks et al. 1987)
  ## Conversion constants
  Kelvin       = 273.15,         # conversion degree Celsius to Kelvin
  DwDc         = 1.6,            # Ratio of the molecular diffusivities for water vapor and CO2
  days2seconds = 86400.0,        # seconds per day
  kPa2Pa       = 1000.0,         # conversion kilopascal (kPa) to pascal (Pa)
  Pa2kPa       = 0.001,          # conversion pascal (Pa) to kilopascal (kPa)
  umol2mol     = 1.0e-6,         # conversion micromole (umol) to mole (mol)
  mol2umol     = 1.0e6,          # conversion mole (mol) to micromole (umol)
  kg2g         = 1000.0,         # conversion kilogram (kg) to gram (g)
  g2kg         = 0.001,          # conversion gram (g) to kilogram (kg)
  kJ2J         = 1000.0,         # conversion kilojoule (kJ) to joule (J)
  J2kJ         = 0.001,          # conversion joule (J) to kilojoule (kJ)
  se_median    = 1.253,          # conversion standard error (SE) of the mean to SE of the median (http://influentialpoints.com/Training/standard_error_of_median.htm)
  frac2percent = 100.0           # conversion between fraction and percent
)
  Dict(
:cp => cp, :Rgas => Rgas, :Rv => Rv, :Rd => Rd, :Md => Md, :Mw => Mw, :eps => eps, :g => g,
:solar_constant => solar_constant,
:pressure0 => pressure0, :Tair0 => Tair0, :k => k, :Cmol => Cmol, :Omol => Omol,
:H2Omol => H2Omol,
:sigma => sigma, :Pr => Pr, :Sc_CO2 => Sc_CO2, :Kelvin => Kelvin, :DwDc => DwDc,
:days2seconds => days2seconds, :kPa2Pa => kPa2Pa, :Pa2kPa => Pa2kPa, :umol2mol => umol2mol,
:mol2umol => mol2umol, :kg2g => kg2g, :g2kg => g2kg, :kJ2J => kJ2J, :J2kJ => J2kJ,
:se_median => se_median, :frac2percent => frac2percent
  )
end
