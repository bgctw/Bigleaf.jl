@testset "vernal point" begin
    hours = 0:24
    lat,long = 0.0,0.0
    deg2second = 24*3600/360
    datetimes = DateTime(2021,3,20) .+ Hour.(hours) .- Second(round(long*deg2second))
    res2 = @pipe calc_sun_position_MOD.(datetime2julian.(datetimes)) |> toDataFrame(_)
    # latitude over ecliptic of sun is zero
    @test all(res2.β .≈ 0.0) 
    au2m = 149597870700.0
    # distance to sun about one astromical unit (150_000 km)
    @test all(.≈(res2.r,  au2m, rtol = 0.02)) 
    # declination crossing from negative to positive
    @test res2.δ[1] < 0
    @test res2.δ[end] > 0
    # azimuth crossing zero longitude
    @test res2.α[1] < 0
    @test res2.α[end] > 0
    # function tmpf()
    #     using Plots, StatsPlots
    #     @df res2 scatter(hours, cols([:α, :δ]), legend = :topleft, xlab="hours")
    #     @df res3 scatter(hours, cols(1:2), legend = :topleft, xlab="hours")
    # end
    #
    res3 = @pipe calc_sun_position_hor.(datetimes, lat, long) |> toDataFrame(_)
    # altitude (-pi, +pi)
    @test minimum(res3.altitude) ≈ -π/2 atol = 0.1
    @test maximum(res3.altitude) ≈ +π/2 atol = 0.1
    #
    lat,long = 45.0,0.0
    datetimes = DateTime(2021,3,20) .+ Hour.(hours) .- Second(round(long*deg2second))
    res3 = @pipe calc_sun_position_hor.(datetimes, lat, long) |> toDataFrame(_)
    # altitude (-pi, +pi)
    @test minimum(res3.altitude) ≈ -π/4 atol = 0.1
    @test maximum(res3.altitude) ≈ +π/4 atol = 0.1
    # azimuth increasing (aside edge cases)
    @test all(diff(res3.azimuth[2:(end-1)]) .> 0)
end

@testset "Dresden summer" begin
    hours = 0:24
    lat,long = 51.0, 13.6
    deg2second = 24*3600/360
    doy = 160
    datetimes = DateTime(2021) .+Day(doy-1) .+ Hour.(hours) .- Second(round(long*deg2second))
    # function tmpf()
    #     using Plots, StatsPlots
    #     @df res2 scatter(hours, cols([:α, :δ]), legend = :topleft, xlab="hours")
    #     @df res3 scatter(hours, cols(1:2), legend = :topleft, xlab="hours")
    # end
    #
    res3 = @pipe calc_sun_position_hor.(datetimes, lat, long) |> toDataFrame(_)
    # maximum altitude at noon
    @test argmax(res3.altitude) == 1 + (length(hours)-1) ÷ 2
    # azimuth increasing (aside edge cases)
    @test all(diff(res3.azimuth[2:(end-1)]) .> 0)
end
  

