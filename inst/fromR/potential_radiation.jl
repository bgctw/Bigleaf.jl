
#' Extraterrestrial solar radiation
#'
#' Compute the extraterrestrial solar radiation with the
#' eccentricity correction.
#'
#' - doy integer vector with day of year (DoY)
#' - constants solar_constant - solar constant (W m-2)
#'
#' # Details

#' Computation follows Lanini, 2010 (Master thesis, Bern University)
#'
#' # Value
 numeric vector of extraterrestrial radiation (W_m-2)
#'
#' ```@example; output = false
#' ```
#' plot(1:365, extraterrestrial_radiation(1:365), type = "l"
#'   , ylab = "radiation (W m-2)", xlab = "day of year")
#'   
"""
"""
function extraterrestrial_radiation(doy,constants = bigleaf_constants())
  
  # Fractional year in radians
  FracYearRad = 2 * pi * (doy - 1) / 365.24
  
  #Eccentricity correction
  ExtRadiation = constants[:solar_constant] * (
    1.00011 + 0.034221 * cos(FracYearRad) + 0.00128 * sin(FracYearRad)
     + 0.000719 * cos(2 * FracYearRad) + 0.000077 * sin(2 * FracYearRad)
     )
  
  return(ExtRadiation)
end




#' Potential radiation
#'
#' Compute potential radiation for given geolocation and day of year.
#'
#' - doy          Integer vector with day of year (start at 1), same length as `hour` or length 1.
#' - hour         Numeric vector with daytime as decimal hour of local time zone
#' - latDeg       Latitude (decimal degrees)
#' - longDeg      Longitude (decimal degrees)
#' - timezone     Time zone (hours)
#' - useSolartime by default corrects hour (given in local winter time)
#'                     for latitude to solar time (where noon is exactly at 12:00).
#'                     Set this to `FALSE` to directly use local winter time.
#'
#' # Value
 vector of potential radiation (W m-2)
#' 
#' ```@example; output = false
#' ```
#' hour = seq(5, 18, by = 0.1)
#' potRadApparentLocal = potential_radiation(
#'   160, hour, 39.94, -5.77, timezone = +1)
#' potRadTimezone = potential_radiation(
#'   160, hour, 39.94, -5.77, timezone = +1, useSolartime = FALSE)
#' plot(potRadApparentLocal ~ hour, type = 'l'
#'   , ylab = 'potential radiation (W m-2)')
#' lines(potRadTimezone ~  hour, col = "blue")
#' abline(v = 12, col = "blue", lty = "dotted")
#' legend("bottomright", legend = c("solar time", "local winter time")
#' , col = c("black", "blue"), inset = 0.05, lty = 1)
#' @importFrom solartime computeSunPositionDoyHour
"""
"""
function potential_radiation(doy, hour, latDeg, longDeg, timezone, useSolartime = TRUE)
  
  # Calculate potential radiation from solar elevation and extraterrestrial solar radiation
  solElevRad = computeSunPositionDoyHour(
    doy, hour, latDeg, longDeg, timezone
    , isCorrectSolartime = useSolartime)[,"elevation"]
  
  extRadiation = extraterrestrial_radiation(doy)
  
  potRad = ifelse(
    solElevRad <= 0, 0, extRadiation * sin(solElevRad) )
  
  return(potRad)
end





