@testset "frac_hour" begin
    @test frac_hour(Minute, 1+1/60) == Hour(1) + Minute(1)
    @test frac_hour(1+1/60) == Hour(1) + Minute(1)
end

@testset "moving_average" begin
    GPPd = @pipe sin.((1:365).*π./365) .+ 0.5 .* rand(365) |> allowmissing(_)
    GPPd[100:120] .= missing
    GPPd_smooth = moving_average(GPPd, 15) 
    #plot(1:365, GPPd)
    #plot!(1:365, GPPd_smooth)
    @test length(GPPd_smooth) == length(GPPd)
    @test GPPd_smooth[1] == mean(GPPd[1:8])
    @test GPPd_smooth[8] == mean(GPPd[1:15])
    @test isfinite(GPPd_smooth[100])
    @test ismissing(GPPd_smooth[101])
end

@testset "get_nonoverlapping_periods" begin
    dfo = DataFrame(
        p_start = DateTime(2021) .+ Hour.(SA[0,4,5,10]),
        p_end = DateTime(2021) .+ Hour.(SA[2,5,6,11]),
        )
    dfno = get_nonoverlapping_periods(dfo)        
    @test dfno.p_start == DateTime(2021) .+ Hour.(SA[0,4,10])
    @test dfno.p_end == DateTime(2021) .+ Hour.(SA[2,6,11])
    #
    # extend last row
    dfo = DataFrame(
        p_start = DateTime(2021) .+ Hour.(SA[0,4,5]),
        p_end = DateTime(2021) .+ Hour.(SA[2,5,6]),
        )
    dfno = get_nonoverlapping_periods(dfo)        
    @test dfno.p_start == DateTime(2021) .+ Hour.(SA[0,4])
    @test dfno.p_end == DateTime(2021) .+ Hour.(SA[2,6])
end

@testset "set_datetime_ydh!" begin
    dfo = DataFrame(year = 2021, doy = repeat(1:3, inner = 48), hour = repeat(0.0:0.5:23.5, outer = 3))
    df = copy(dfo)
    df2 = set_datetime_ydh!(df)
    @test df2 === df
    @test df2.datetime[1] == ZonedDateTime(DateTime(2021), tz"UTC")
    @test all(Dates.year.(df2.datetime) .== dfo.year)
    @test all(Dates.dayofyear.(df2.datetime) .== dfo.doy)
    @test all(Dates.hour.(df2.datetime) .+ Dates.minute.(df2.datetime)/60 .== dfo.hour)
end
