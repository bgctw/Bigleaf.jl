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
