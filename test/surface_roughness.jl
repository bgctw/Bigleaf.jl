#@testset "wind_profile" begin
    datetime, ustar, Tair, pressure, H = values(tha48[1,:])
    z = 30
    d=0.7*tha_heights.zh
    z0m=2.65
    u30 = wind_profile(z, ustar, d, z0m)
    @test ≈(u30, 1.93, rtol = 1/100 )
    #
    u30c = wind_profile(Val(:Dyer_1970), z, ustar, Tair,pressure, H, d, z0m)
    @test ≈(u30c, 2.31, rtol = 1/100 )
    #
    df = copy(tha48)
    windz = wind_profile(df, z, d, z0m; stab_formulation = Val(:no_stability_correction))    
    @test length(windz) == 48
    @test windz[1] == u30
    windzc = wind_profile(df, z, d, z0m; stab_formulation = Val(:Dyer_1970))    
    @test windzc[1] == u30c
    #plot(windz)
    #plot!(windz2)
    psi_m = stability_correction!(copy(df, copycols=false), z, d).psi_m
    windzc2 = wind_profile(df, z, d, z0m, psi_m)    
    @test windzc2 == windzc
end