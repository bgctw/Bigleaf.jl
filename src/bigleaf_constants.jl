"""
    BigLeafConstants(;...)

Constants used troughout the BigLeaf.jl Package

Default values can be overridden by the named arguments of the constructor:

- cp           : Specific heat of air for constant pressure (J K-1 kg-1)
- Rgas         : Universal gas constant (J mol-1 K-1)
- Rv           : Gas constant of water vapor (J kg-1 K-1) (Stull 1988 p.641)
- Rd           : Gas constant of dry air (J kg-1 K-1) (Foken p. 245)
- Md           : Molar mass of dry air (kg mol-1)
- Mw           : Molar mass of water vapor (kg mol-1)
- eps          : Ratio of the molecular weight of water vapor to dry air (=Mw/Md)
- g            : Gravitational acceleration (m s-2)
- solar_constant: Solar constant (W m-2)
- pressure0    : Reference atmospheric pressure at sea level (Pa)
- Tair0        : Reference air temperature (K)
- k            : von Karman constant
- Cmol         : Molar mass of carbon (kg mol-1)
- Omol         : Molar mass of oxygen (kg mol-1)
- H2Omol       : Molar mass of water (kg mol-1)
- sigma        : Stefan-Boltzmann constant (W m-2 K-4)
- Pr           : Prandtl number
- Sc_CO2       : Schmidt number for CO2
- Kelvin       : Conversion degree Celsius to Kelvin
- DwDc         : Ratio of the molecular diffusivities for water vapor and CO2
- days2seconds : Seconds per day
- kPa2Pa       : Conversion kilopascal (kPa) to pascal (Pa)
- Pa2kPa       : Conversion pascal (Pa) to kilopascal (kPa)
- umol2mol     : Conversion micromole (umol) to mole (mol)
- mol2umol     : Conversion mole (mol) to micromole (umol)
- kg2g         : Conversion kilogram (kg) to gram (g)
- g2kg         : Conversion gram (g) to kilogram (kg)
- kJ2J         : Conversion kilojoule (kJ) to joule (J)
- J2kJ         : Conversion joule (J) to kilojoule (kJ)
- se_median    : Conversion standard error (SE) of the mean to SE of the median
  (http://influentialpoints.com/Training/standard_error_of_median.htm)
- frac2percent : Conversion between fraction and percent

## Examples
```jldoctest; output=false
BigLeafConstants().g == 9.81
# on the moon change gravity constant to 1/6 that of the earth
BigLeafConstants(g = BigLeafConstants().g/6).g == 9.81/6
# output
true
```
"""
@with_kw struct BigLeafConstants{FT,IT,RT} 
  cp::FT         = 1004.834        # specific heat of air for constant pressure (J K-1 kg-1)
  Rgas::FT       = 8.31451         # universal gas constant (J mol-1 K-1)
  Rv::FT         = 461.5           # gas constant of water vapor (J kg-1 K-1) (Stull 1988 p.641)
  Rd::FT         = 287.0586        # gas constant of dry air (J kg-1 K-1) (Foken 2008 p. 245)
  Md::FT         = 0.0289645       # molar mass of dry air (kg mol-1)
  Mw::FT         = 0.0180153       # molar mass of water vapor (kg mol-1)
  eps::FT        = 0.622           # ratio of the molecular weight of water vapor to dry air
    # (=Mw/Md)
  g::FT          = 9.81            # gravitational acceleration (m s-2)
  solar_constant::FT = 1366.1      # solar constant, i.e. solar radation at earth distance from 
  # the sun (W m-2)
  pressure0::FT  = 101325.0        # reference atmospheric pressure at sea level (Pa)
  Tair0::FT      = 273.15          # reference air temperature (K)
  k::FT          = 0.41            # von Karman constant
  Cmol::FT       = 0.012011        # molar mass of carbon (kg mol-1)
  Omol::FT       = 0.0159994       # molar mass of oxygen (kg mol-1)
  H2Omol::FT     = 0.01801528      # molar mass of water (kg mol-1)
  sigma::FT      = 5.670367e-08    # Stefan-Boltzmann constant (W m-2 K-4)
  Pr::FT         = 0.71            # Prandtl number
  Sc_CO2::FT     = 1.07            # Schmidt number for CO2 (Hicks et al. 1987)
  ## Conversion constants
  Kelvin::FT       = 273.15        # conversion degree Celsius to Kelvin
  DwDc::FT         = 1.6           # Ratio of the molecular diffusivities for water vapor and CO2
  days2seconds::IT = 86400         # seconds per day
  kPa2Pa::IT       = 1000          # conversion kilopascal (kPa) to pascal (Pa)
  Pa2kPa::RT       = 1/1000        # conversion pascal (Pa) to kilopascal (kPa)
  umol2mol::RT     = 1/1e6         # conversion micromole (umol) to mole (mol)
  mol2umol::IT     = 1_000_000     # conversion mole (mol) to micromole (umol)
  kg2g::IT         = 1000          # conversion kilogram (kg) to gram (g)
  g2kg::RT         = 1/1000        # conversion gram (g) to kilogram (kg)
  kJ2J::IT         = 1000          # conversion kilojoule (kJ) to joule (J)
  J2kJ::RT         = 1/1000        # conversion joule (J) to kilojoule (kJ)
  se_median::FT    = 1.253         # conversion standard error (SE) of the mean to SE of the 
  # median (http://influentialpoints.com/Training/standard_error_of_median.htm)
  frac2percent::IT = 100           # conversion between fraction and percent
end

