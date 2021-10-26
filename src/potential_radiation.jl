"""
    extraterrestrial_radiation(doy::Number; ...)
    extraterrestrial_radiation(datetime::TimeType; ...)

Compute the extraterrestrial solar radiation with the
eccentricity correction. Computation follows Lanini, 2010 (Master thesis, Bern University).


# Arguments
- doy: integer vector with day of year (DoY)
optional 
- `constants=`[`bigleaf_constants`](@ref)`()`: Dictionary with entries 
   - solar_constant
- `year`=2030: year to create timestamps. Due to precession results slightly
  change across decades.   

# Value
numeric vector of extraterrestrial radiation (W_m-2)

# Examples
```jldoctest; output = false
ex_rad = extraterrestrial_radiation(1)
≈(ex_rad, 1414, atol = 1)
# output
true
```
"""
function extraterrestrial_radiation(doy::Number;constants = bigleaf_constants(), year = 2030)
  FracYearRad = 2 * π * (doy - 1) / 365.24
  #Eccentricity correction
  ExtRadiation = constants[:solar_constant] * (
    1.00011 + 0.034221 * cos(FracYearRad) + 0.00128 * sin(FracYearRad)
     + 0.000719 * cos(2 * FracYearRad) + 0.000077 * sin(2 * FracYearRad)
     )
end,
function extraterrestrial_radiation(datetime::TimeType; constants = bigleaf_constants())
  # Fractional year in radians
  doy = Dates.dayofyear(datetime)
  extraterrestrial_radiation(doy; constants)
end

"""
    potential_radiation(datetime, lat, long)
    potential_radiation(doy, hour, lat, long; timezone,year)

Compute potential radiation for given geolocation and time.

Because potential radiation does not change across years 
(for a given location and daytime and day of year), time
can be specified alternatively by doy and hour.

# Arguments
- datetime:  UTC timestamp.
- lat:       Latitude (decimal degrees)
- long:      Longitude (decimal degrees)
- doy:       day of year (integer starting from 1)
- hour:      hour within the day (0.0 .. 24.0 float)
optional     
- timezone:  Timezone for doy and hour, defaults to "GMT+x"
   nearest to given longitude.
- year: specific year for doy and hour

# Value
vector of potential radiation (W m-2)

# Examples
```@example
# assume hours in the GMT+x that is closest to given longitude
potrad = potential_radiation(160, 10.5, 51.0, 11.5)
```
"""
function potential_radiation(doy, hour, lat, long; timezone = FixedTimeZone("UTC+"*string(round(Int32, long/15))), year = 2030)
  # the following doctest does not suppress warnings and fails
  # ```jldoctest; output = false
  # # assume hours in the GMT+x that is closest to given longitude
  # potrad = potential_radiation(160, 10.5, 51.0, 11.5)
  # ≈(potrad, 1093, atol=1)
  # # output
  # true
  # ```
  @pipe get_datetime_for_doy_hour(doy,hour; year) |>  
    ZonedDateTime(_, timezone) |> 
    potential_radiation(_, lat, long)
end,
function potential_radiation(datetime::TimeType, lat, long)
  # Calculate potential radiation from solar elevation and extraterrestrial solar radiation
  solElevRad = calc_sun_position_hor(datetime, lat, long).altitude
  extRadiation = extraterrestrial_radiation(datetime)
  potRad = extRadiation * sin(solElevRad)
end

"""
    get_datetime_for_doy_hour(doy, hour=12; year = 2030)

Create DateTime for given day_of_year and hour. Hour defaults
to noon and year to 2030, a near future where earth axis 
precession does not differ too much from year where the function 
is called. Fractional hours can be provided.

# Examples
```jldoctest output = false
get_datetime_for_doy_hour(2; year = 2021)
# output
2021-01-02T12:00:00
```
"""
get_datetime_for_doy_hour(doy, hour::Integer=12; year = 2030) = DateTime(year) + Day(doy-1) + Hour(hour)

get_datetime_for_doy_hour(doy, hour; year = 2030) = DateTime(year) + Day(doy-1) + frac_hour(Dates.Millisecond, hour)





