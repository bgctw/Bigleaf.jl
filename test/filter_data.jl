@testset "get_growingseason" begin
    rng = StableRNG(815)
    GPPd = @pipe sin.((1:365).*π./365) .+ 0.5 .* rand(rng,365) |> allowmissing(_)
    GPPd[100:120] .= missing
    #plot(1:365, GPPd)
    growseas = @test_logs (:warn, r"gap in 'GPPd'") get_growingseason(GPPd, 0.5) 
    #plot(growseas)
    @test rle(growseas.is_growingseason) == (Bool[0, 1, 0], [51, 258, 56])
    growseas = get_growingseason(GPPd, 0.5; min_int = 56, warngap = false) 
    @test rle(growseas.is_growingseason) == (Bool[1, 0], [309, 56])
end

@testset "get_growingseason with few data" begin
    rng = StableRNG(815)
    GPPd = @pipe sin.((1:365).*π./365) .+ 0.5 .* rand(rng,365) |> allowmissing(_)
    GPPd[1:200] .= missing
    @test_throws ErrorException get_growingseason(GPPd, 0.5) 
end

@testset "setinvalid_nongrowingseason!" begin
    rng = StableRNG(815)
    nday = 365
    GPPd = @pipe sin.((1:nday).*π./365) .+ 0.5 .* rand(rng,365) |> allowmissing(_)
    GPPd[100:120] .= missing
    df = DataFrame(datetime = DateTime(2021) .+ Hour(1).+ Day.(0:nday-1), GPP = GPPd)
    df2 = copy(df)
    df2a = setinvalid_nongrowingseason!(df2, 0.5; warngap = false)
    @test df2a === df2
    @test propertynames(df2a) == union(propertynames(df), SA[:valid])
    @test rle(df2a.valid) == (Bool[0, 1, 0], [51, 258, 56])
    #
    df2.valid[80:90] .= false
    df2a = setinvalid_nongrowingseason!(df2, 0.5; warngap = false, update_GPPd_smoothed = true)
    @test df2a === df2
    @test propertynames(df2a) == union(propertynames(df), SA[:valid, :GPPd_smoothed])
    @test rle(df2a.valid) == (Bool[0, 1, 0, 1, 0], [51, 28, 11, 219, 56])
end

@testset "setinvalid_qualityflag!" begin
    df = DataFrame(
        NEE = 1:3,
        GPP = 10:10:30,
        NEE_qc = [1,1,2],
        GPP_qc = [1,missing,1],
    )
    df2 = copy(df)
    df2a = setinvalid_qualityflag!(df2; vars = SA["NEE", "GPP"], setvalmissing = false)
    @test df2a === df2
    @test df2a.valid == [true, false, false]
    @test df2a.NEE == 1:3 # not set to missing
    # test with existing :valid column, where already false entries are kept
    df2.valid[1] = false
    df2a = setinvalid_qualityflag!(df2; vars = SA["NEE", "GPP"], setvalmissing = false)
    @test df2a === df2
    @test df2a.valid == [false, false, false]
    @test df2a.NEE == 1:3 # not set to missing
    #
    # test with setting values to missing
    df2 = copy(df)
    df2a = setinvalid_qualityflag!(df2; vars = SA["NEE", "GPP"], setvalmissing = true)
    @test df2a === df2
    @test df2a.valid == [true, false, false]
    @test isequal(df2a.NEE, [1,2,missing])
    @test isequal(df2a.GPP, [10,missing,30])
end

@testset "setinvalid_range!" begin
    df = DataFrame(
        NEE = 1:3,
        GPP = 10:10:30,
    )
    allowmissing!(df, :NEE)
    var_ranges = [:NEE => (-2.0,4.0), :GPP => (8.0,28.0)]
    df2 = copy(df)
    df3 = copy(df2)
    df2a = setinvalid_range!(df2, var_ranges...; setvalmissing = false)
    # @btime setinvalid_range!($df2, $(var_ranges)...; setvalmissing = false)
    @test df2a === df2
    @test df2a.valid == [true, true, false]
    df3a = setinvalid_range!(df3, var_ranges...)
    @test df3a === df3
    @test isequal(df3a.valid, df2a.valid)
    @test ismissing(df3a.GPP[3])
    @test all(.!ismissing.(df3a.GPP[1:2]))
    #
    df2.NEE[1] = missing
    df3 = copy(df2)
    df2b = setinvalid_range!(df2, var_ranges...; setvalmissing = false)
    @test df2b === df2
    @test df2b.valid == [false, true, false]
    df3a = setinvalid_range!(df3, var_ranges...)
    @test isequal(df3a.valid, df2a.valid)
    @test ismissing(df3a.GPP[3])
    @test all(.!ismissing.(df3a.GPP[1:2]))
    @test all(.!ismissing.(df3a.NEE[2:3]))
    #
    df2.valid[2] = false
    df3 = copy(df2)
    df2b = setinvalid_range!(df2, var_ranges...; setvalmissing = false)
    @test df2b === df2
    @test df2b.valid == [false, false, false]
    df3a = setinvalid_range!(df3, var_ranges...)
    @test isequal(df3a.valid, df2a.valid)
end

@testset "setinvalid_afterprecip!" begin
    min_precip = 0.01 
    hours_after = 1.0
    dfo = DataFrame(datetime = DateTime(2021) .+ Hour.(1:13), precip = 0.0)
    dfo.precip[3] = 0.001
    dfo.precip[10:12] .= 2.0
    dfo.precip[6:8] .= 2.0
    allowmissing!(dfo, :precip); dfo.precip[2] = missing
    #
    df = copy(dfo)
    df2 = setinvalid_afterprecip!(df; min_precip, hours_after)
    @test df2 === df 
    @test all(.!df.valid[6:9])
    @test all(.!df.valid[10:13])
    @test all(df.valid[Not(vcat(6:9,10:13))])
end

