@testset "GPPfilter" begin
    rng = StableRNG(815)
    GPPd = @pipe sin.((1:365).*π./365) .+ 0.5 .* rand(rng,365) |> allowmissing(_)
    GPPd[100:120] .= missing
    #plot(1:365, GPPd)
    growseas = @test_logs (:warn, r"gap in 'GPPd'") filter_growing_season(GPPd, 0.5) 
    #plot(growseas)
    @test rle(growseas) == (Bool[0, 1, 0], [51, 258, 56])
    growseas = filter_growing_season(GPPd, 0.5; min_int = 56, warngap = false) 
    @test rle(growseas) == (Bool[1, 0], [309, 56])
end

