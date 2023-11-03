"""
    air_density(Tair,pressure; ...)

Air density of moist air from air temperature and pressure_

# Arguments
- Tair:      Air temperature (deg C)
- pressure:  Atmospheric pressure (kPa)
optional 
- `constants=`[`BigleafConstants`](@ref)`()`: Dictionary with entries 
   Kelvin, kPa2Pa

# Details 
Air density ``\\rho`` is calculated as:

  ``\\rho = {pressure \\over Rd * Tair}``

# Value 
air density (kg m-3)

# Examples
```@example doc 
# air density at 25degC and standard pressure (101.325kPa)
air_density(25,101.325)
```

# References  
Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany.
"""
function air_density(Tair,pressure; constants=BigleafConstants())
  Tair_K     = Tair + oftype(Tair,constants.Kelvin)
  pressure_Pa = pressure * oftype(pressure,constants.kPa2Pa)
  rho = pressure_Pa / (oftype(Tair_K,constants.Rd) * Tair_K) 
end


"""
    pressure_from_elevation(elev,Tair,VPD=missing;...)

An estimate of mean pressure at a given elevation as predicted by the
             hypsometric equation.

# Arguments             
- elev:      Elevation a_s_l_ (m)
- Tair:      Air temperature (deg C)
- VPD:       Vapor pressure deficit (kPa); optional
optional 
- `constants=`[`BigleafConstants`](@ref)`()`: Dictionary with entries 
    Kelvin, pressure0, Rd, g, Pa2kPa

# Details 
Atmospheric pressure is approximated by the hypsometric equation:

``pressure = pressure_0 / (exp(g * elevation / (Rd Temp)))``
      
The hypsometric equation gives an estimate of the standard pressure
      at a given altitude.
      If VPD is provided, humidity correction is applied and the
      virtual temperature instead of air temperature is used_ VPD is 
      internally converted to specific humidity.

# Value 
Atmospheric pressure (kPa)
                           
# References  
Stull B_, 1988: An Introduction to Boundary Layer Meteorology.
            Kluwer Academic Publishers, Dordrecht, Netherlands.

# Examples
```@example doc
# mean pressure at 500m altitude at 25 deg C and VPD of 1 kPa
pressure_from_elevation(500,25,1)
```
"""                           
function pressure_from_elevation(elev,Tair,VPD=missing; constants=BigleafConstants())
  Tair_K     = Tair + oftype(Tair,constants.Kelvin)
  if ismissing(VPD)
    pressure = oftype(Tair,constants.pressure0) / exp(oftype(Tair,constants.g) * elev / (oftype(Tair,constants.Rd)*Tair_K))
  else 
    pressure1   = oftype(Tair,constants.pressure0) / exp(oftype(Tair,constants.g) * elev / (oftype(Tair,constants.Rd)*Tair_K))
    Tv          = virtual_temp(
      Tair_K - oftype(Tair,constants.Kelvin), pressure1 * constants.Pa2kPa, VPD;
      Esat_formula=Sonntag1990(),constants) + oftype(Tair,constants.Kelvin)
    pressure    = oftype(Tair,constants.pressure0) / exp(oftype(Tair,constants.g) * elev / (oftype(Tair,constants.Rd)*Tv))
  end
  pressure = pressure * constants.Pa2kPa
end 

"""
    psychrometric_constant(Tair,pressure; ...)

Computes the psychrometric 'constant'.

# Arguments
- Tair:      Air temperature (deg C)
- pressure:  Atmospheric pressure (kPa)
optional 
- `constants=`[`BigleafConstants`](@ref)`()`: Dictionary with entries 
    cp, eps
                 
# Details 
The psychrometric constant (``\\gamma``) is given as:

``\\gamma = cp * pressure / (eps * \\lambda)``
 
 where ``\\lambda`` is the latent heat of vaporization (J kg-1), 
 as calculated from [`latent_heat_vaporization`](@ref).
 
# Value 
the psychrometric constant (kPa K-1)
 
# References  
Monteith J.L., Unsworth M.H., 2008: Principles of Environmental Physics.
            3rd Edition. Academic Press, London.

# Examples            
```@example doc 
psychrometric_constant.(5.0:5.0:25.0, 100)
```
"""
function psychrometric_constant(Tair,pressure; constants=BigleafConstants())
  lambda = latent_heat_vaporization(Tair)
  gamma  = (oftype(pressure,constants.cp) * pressure) / (oftype(lambda,constants.eps) * lambda)
end

"""
    latent_heat_vaporization(Tair) 

Latent heat of vaporization as a function of air temperature
using 

``\\lambda = (2.501 - 0.00237 \\, Tair) 10^6``.

# Arguments:
- Tair:  Air temperature (deg C)

# Value
``\\lambda``: Latent heat of vaporization (J kg-1) 

# References
- Stull, B_, 1988: An Introduction to Boundary Layer Meteorology (p_641)
            Kluwer Academic Publishers, Dordrecht, Netherlands            
- Foken, T, 2008: Micrometeorology_ Springer, Berlin, Germany


```@example
latent_heat_vaporization.(5:5:45)        
```
"""            
function latent_heat_vaporization(Tair) 
  k1 = 2.501
  k2 = 0.00237
  lambda = ( k1 - k2 * Tair ) * 1e+06
end


"""
    wetbulb_temp(Tair, pressure, VPD; ...)

Calculate the wet bulb temperature, i_e_ the temperature
             that the air would have if it was saturated

# Arguments              
- Tair:      Air temperature (deg C)
- pressure:  Atmospheric pressure (kPa)
- VPD:       Vapor pressure deficit (kPa)
optional
- accuracy:  Accuracy of the result (deg C)
- `Esat_formula`: formula used in [`Esat_from_Tair`](@ref)
- `constants=`[`BigleafConstants`](@ref)`()`: Dictionary with entries 
     cp, eps,Pa2kPa

# Details 
Wet-bulb temperature (Tw) is calculated from the following expression:
         
``e = Esat(Tw) - gamma* (Tair - Tw)``
         
The equation is optimized for Tw.
Actual vapor pressure e (kPa) is calculated from VPD using [`VPD_to_e`](@ref).
The psychrometric constant gamma (kPa K-1) is calculated using
[`psychrometric_constant`](@ref).
         
# Value 
wet-bulb temperature (degC)
             
# References  
Monteith J.L., Unsworth M.H., 2008: Principles of Environmental Physics.
            3rd edition. Academic Press, London.
   
# Examples            
```@example doc 
wetbulb_temp.([20,25.0], 100, [1,1.6])             
```
"""
function wetbulb_temp(Tair, pressure, VPD; accuracy=1e-03,
  Esat_formula=Sonntag1990(), constants=BigleafConstants())
  gamma  = psychrometric_constant(Tair,pressure)
  ea     = VPD_to_e(VPD,Tair;Esat_formula)
  wetbulb_temp_from_e_Tair_gamma(ea,Tair,gamma; accuracy,Esat_formula,constants)
end,
function wetbulb_temp_from_e_Tair_gamma(ea, Tair, gamma; accuracy=1e-03,
  Esat_formula=Sonntag1990(), constants=BigleafConstants())
  if accuracy > one(accuracy)
    @warn ("'accuracy' is set to 1 degC")
    accuracy = one(accuracy)
  end
  fopt(Tw) = abs(ea - (Esat_from_Tair(Tw; Esat_formula = Esat_formula,constants) - 
    0.93*gamma*(Tair - Tw)))
  resopt = optimize(fopt, -100, 100, abs_tol = accuracy)
  roundmult(resopt.minimizer, accuracy)
end

"""              
    dew_point(Tair,VPD; ...)
    dew_point_from_e(ea; ...)

Calculate the dew point, the temperature to which air must be 
             cooled to become saturated (i_e_ e = Esat(Td))

# Arguments             
- Tair:     Air temperature (degC)
- VPD:      Vapor pressure deficit (kPa)
- ea:       actual water vapor pressure (kPa)
optional
- accuracy = 1e-03: Accuracy of the result (deg C)
- `Esat_formula`: formula used in [`Esat_from_Tair`](@ref)
- `constants=`[`BigleafConstants`](@ref)`()`: Dictionary with entries 
    Pa2kPa 

# Details 
Dew point temperature (Td) is defined by the temperature for which
saturaed vapour pressure equals current vapour pressure, i.e. below which 
water would start to condensate.

``e = Esat(Td)``

where e is vapor pressure of the air and Esat is the vapor pressure deficit.
         This equation is solved for Td by optimization.
         
# Value 
dew point temperature (degC)

# References  
Monteith J.L., Unsworth M.H., 2008: Principles of Environmental Physics.
            3rd edition. Academic Press, London.

# Examples            
```@example doc
Tair = 20.0
VPD = 1.5
Td = dew_point(Tair, VPD; accuracy = 1e-2)                
(e = VPD_to_e(VPD,Tair), esat_Td = Esat_from_Tair(Td))
```
"""
function dew_point(Tair, VPD; accuracy=1e-03,Esat_formula=Sonntag1990(),
                      constants=BigleafConstants())
  ea = VPD_to_e(VPD,Tair;Esat_formula)
  dew_point_from_e(ea; accuracy, Esat_formula, constants)
end,
function dew_point_from_e(ea;accuracy=1e-03,Esat_formula=Sonntag1990(), constants=BigleafConstants())
  if accuracy > one(accuracy)
    @warn ("'accuracy' is set to 1 degC")
    accuracy = one(accuracy)
  end
  fopt(Td) = abs(ea - Esat_from_Tair(Td;Esat_formula = Esat_formula,constants))
  resopt = optimize(fopt, -100, 100, abs_tol = accuracy)
  roundmult(resopt.minimizer, accuracy)
end
# https://discourse.julialang.org/t/rounding-to-a-multiple-of-float/21295/21
roundmult(val, prec) = (inv_prec = 1 / prec; round(val * inv_prec) / inv_prec)



"""
    virtual_temp(Tair,pressure,VPD; ...)

Virtual temperature, defined as the temperature at which dry air would have the same
             density as moist air at its actual temperature.

# Arguments            
- `Tair`:      Air temperature (deg C)
- `pressure`:  Atmospheric pressure (kPa)
- `VPD`:       Vapor pressure deficit (kPa)
- `Esat_formula`: formula used in [`Esat_from_Tair`](@ref)
optional 
- `constants=`[`BigleafConstants`](@ref)`()`: Dictionary with entries 
  - :Kelvin - conversion degree Celsius to Kelvin
  - :eps - ratio of the molecular weight of water vapor to dry air (-) 

# Details
The virtual temperature is given by:
 
```Tv = Tair / (1 - (1 - eps) e/pressure)```

 where Tair is in Kelvin (converted internally)_ Likewise, VPD is converted 
 to actual vapor pressure (e in kPa) with [`VPD_to_e`](@ref) internally.

# Value
virtual temperature (deg C)

# References
- Monteith J.L., Unsworth M.H., 2008: Principles of Environmental Physics.
            3rd edition. Academic Press, London.
 
# Examples            
```jldoctest; output = false
Tair,pressure,VPD = 25.0,100.0,1.5
vt = virtual_temp(Tair,pressure,VPD)  
≈(vt, 26.9, atol =0.1)
# output
true
```
"""
function virtual_temp(Tair,pressure,VPD;Esat_formula=Sonntag1990(),
                         constants=BigleafConstants())
  e    = VPD_to_e(VPD,Tair;Esat_formula)
  Tair_Kelvin = Tair + oftype(Tair,constants.Kelvin)
  Tv_Kelvin = Tair_Kelvin / (1 - (1 - oftype(pressure,constants.eps)) * e/pressure) 
  Tv = Tv_Kelvin - oftype(Tair,constants.Kelvin)
end


"""         
    kinematic_viscosity(Tair,pressure; ...)

Calculate the kinematic viscosity of air.

# Parameters
- Tair      Air temperature (deg C)
- pressure  Atmospheric pressure (kPa)
optional 
- `constants=`[`BigleafConstants`](@ref)`()`: Dictionary with entries 
    Kelvin, pressure0, Tair0, kPa2Pa

# Details 
Eq where v is the kinematic viscosity of the air (m2 s-1), 
         given by (Massman 1999b):
         
``v = 1.327 * 10^-5(pressure0/pressure)(Tair/Tair0)^{1.81}``
         
# Value 
kinematic viscosity of air (m2 s-1)

# References  
Massman, W.J., 1999b: Molecular diffusivities of Hg vapor in air, 
            O2 and N2 near STP and the kinematic viscosity and thermal diffusivity
            of air near STP. Atmospheric Environment 33, 453-457.      
            
# Examples            
```jldoctest; output = false
Tair,pressure = 25.0,100.0
vis = kinematic_viscosity(Tair,pressure)
≈(vis, 1.58e-5, atol =1e-7)
# output
true
```
"""         
function kinematic_viscosity(Tair,pressure; constants=BigleafConstants())
  (ismissing(Tair) || ismissing(pressure)) && return(missing)
  Tair_Kelvin     = Tair + oftype(Tair, constants.Kelvin)
  pressure_kPa = pressure * oftype(pressure, constants.kPa2Pa)
  v  = oftype(Tair,1.327e-05)*(oftype(pressure,constants.pressure0)/pressure_kPa) *
    (Tair_Kelvin/oftype(Tair,constants.Tair0))^oftype(Tair,1.81)
end
