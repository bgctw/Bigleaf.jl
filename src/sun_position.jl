"""
    calc_sun_position_hor(datetime, lat, long)

Compute the Sun position at given time and observer coordinates in horizontal coordinates.

# Arguments:
- `datetime`: time: Either a `ZonedDateTime`, or `DateTime` assumed in UTC
- `lat`, `long`: latitude and longitude in degree

# Value
`NamedTuple`: sun position with entries of the same type as `lat` and `long`
- `altitude`: angle above the horizon [rad].
- `azimuth`: angle ange the horizon plane eastwards of north [rad]
- `hourangle`: [rad] as output by AstroLib.eq2hor
   Seems to represent time [day/2pi] after solar noon. 
   Value at local timezone noon provdes (local time - solar time).
"""
function calc_sun_position_hor(datetime::ZonedDateTime, lat, long)
  datetimeUTC = DateTime(datetime,UTC)
  calc_sun_position_hor(datetimeUTC, lat, long)
end
function calc_sun_position_hor(datetime::DateTime, lat::FT, long::FT) where FT
    deg2rad = FT(π/180)
    jd = FT(datetime2julian(datetime))
    pos_eq = calc_sun_position_MOD(jd)
    # precession is already account for in MOD
    pos_hor = eq2hor(pos_eq.α/deg2rad, pos_eq.δ/deg2rad, jd, lat, long; precession = false) .* deg2rad
    (altitude = FT(pos_hor[1]), azimuth = FT(pos_hor[2]), hourangle = FT(pos_hor[3]))
end



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Compute the sun position.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] The Astronomical Almanac for the year 2000 (p. C24).
#
#   [2] http://aa.usno.navy.mil/faq/docs/SunApprox.php
#
#   [3] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorne, CA.
#
#   [4] The Astronomical Almanac for the year 2006.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# modified from https://github.com/JuliaSpace/SatelliteToolbox.jl/blob/c5244c50e55e85f6b3a1702de16f00a9c726ee24/src/sun/sun_position.jl#L24-L30
#
# https://en.wikipedia.org/wiki/Ecliptic_coordinate_system

"""
    calc_sun_position_MOD(JD::Number)

Compute the Sun position at the Julian Day `JD`.

Results are represented in the Mean Equinox of Date (MOD),
i.e. accounting for precession but not for nutation and smaller pertubation
of the polar axes, in spherical ecliptic and equatorial coordinates. 
The algorithm was adapted from [Vallado 2013, p. 277-279].

# Arguments:
- `JD`: time given as Julian Day . 

# Value
`NamedTuple`: sun position where
- Ecliptic coordinates (1:3)
  - λ: Ecliptic longitude of the Sun [rad].
  - β: Ecliptic latitude of the Sun [rad] is assumed 0.
  - r: Distance of the Sun from Earth [m].
- Equatorial coordinate (4:5)  
  - α: ascention [rad]
  - δ: declination [rad]
- ϵ: Obliquity of the ecliptic [rad].

# References
- Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
       Microcosm Press, Hawthorne, CA.
"""
function calc_sun_position_MOD(JD::Number)
    # Constants (apdated from SatelliteToolbox)
    deg2rad = π/180
    au2m = 149597870700.0
    JD_J2000 = 2.451545e6

    # Number of Julian centuries from J2000 epoch.
    T_UT1 = (JD-JD_J2000)/36525.0

    # Here, we will assume that T_TBD ≈ T_UT1. This can be done because we are
    # using a low-precision algorithm [3].
    T_TBD = T_UT1

    # Mean anomaly of the Sun [deg].
    Ms = 357.529_109_2 + 35_999.050_34T_TBD

    # Convert Ms to [rad] and limit to the interval [0,2π].
    Ms = mod2pi(Ms*deg2rad)

    # Compute auxiliary variables.
    sinMs,  cosMs  = sincos(Ms)
    sin2Ms, cos2Ms = sincos(2Ms)

    # Mean longitude of the Sun [deg].
    λ_m = 280.460 + 36_000.771T_UT1

    # Ecliptic latitude of the Sun [deg]. # twutz 2110 longitude?
    λ_e = λ_m + 1.914_666_471sinMs + 0.019_994_643sin2Ms

    # Obliquity of the ecliptic [deg].
    ϵ = 23.439_291 - 0.013_004_2T_TBD

    # Convert λ_e and ϵ to [rad] and limit to the interval [0,2π].
    λ_e = mod2pi(λ_e*deg2rad)
    ϵ   = mod2pi(  ϵ*deg2rad)

    # Auxiliary variables.
    sinϵ   , cosϵ   = sincos(ϵ)
    sinλ_e , cosλ_e = sincos(λ_e)

    # Distance of the Sun from Earth [m].
    r = (1.000_140_612 - 0.016_708_617cosMs - 0.000_139_589cos2Ms )*au2m

    # ascension
    α = atan(cosϵ * sinλ_e, cosλ_e)

    # declination
    δ = asin(sinϵ * sinλ_e)

    S_MOD_rad = (
        λ = λ_e, β = 0.0, r = r, 
        α = α, δ = δ, 
        ϵ = ϵ
        )
end
