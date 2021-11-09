@testset "Reynolds_Number" begin
    Tair,pressure,ustar,z0m = 25,100,0.5,0.5
    R = Reynolds_Number(Tair,pressure,ustar,z0m)                             
    @test ≈(R, 15870, rtol=1e-3) 
end


@testset "roughness_parameters" begin
    zh = thal.zh
    zr = thal.zr
    LAI = thal.LAI
    keys_exp = (:d, :z0m, :z0m_se)
    rp = roughness_parameters(Val(:canopy_height), zh)
    #round.(values(rp); sigdigits = 4)
    @test keys(rp) == keys_exp
    @test all(isapproxm.(values(rp), (18.55, 2.65, missing), rtol=1e-3))
    #
    rp = roughness_parameters(Val(:canopy_height_LAI), zh, LAI)
    #round.(values(rp); sigdigits = 4)
    @test keys(rp) == keys_exp
    @test all(isapproxm.(values(rp), (21.77, 1.419, missing), rtol=1e-3))
    #
    df = copy(tha48)
    d=0.7*zh
    psi_m = stability_correction!(copy(df, copycols=false), zr, d).psi_m
    rp = roughness_parameters(Val(:wind_profile), df, zh, zr; psi_m)
    #round.(values(rp); sigdigits = 4)
    @test keys(rp) == keys_exp
    #@test all(isapproxm.(values(rp), (18.55, 1.879, 0.3561), rtol=1e-3))
    #from R:
    @test all(isapproxm.(values(rp), (18.55, 1.879402, 0.356108), rtol=1e-3))
    #
    # no stability correction
    rp0 = roughness_parameters(Val(:wind_profile), df, zh, zr; psi_m = 0.0)
    @test keys(rp0) == keys_exp
    # same magnitude as with stability correction
    @test all(isapproxm.(values(rp0), values(rp), rtol=0.5))
    #
    # estimate psi
    rp_psiauto = roughness_parameters(Val(:wind_profile), df, zh, zr)
    @test propertynames(df) == propertynames(tha48) # not changed
    @test rp_psiauto == rp
end

@testset "wind_profile" begin
    datetime, ustar, Tair, pressure, H = values(tha48[1,:])
    z = 30
    d=0.7*thal.zh
    z0m= 2.65
    u30 = wind_profile(z, ustar, d, z0m)
    @test ≈(u30, 1.93, rtol = 1/100 )
    #
    u30c = wind_profile(Val(:Dyer_1970), z, ustar, Tair,pressure, H, d, z0m)
    @test ≈(u30c, 2.31, rtol = 1/100 )
    #
    z0m=1.9 #2.14 #2.65
    u30 = wind_profile(z, ustar, d, z0m) # used below
    u30c = wind_profile(Val(:Dyer_1970), z, ustar, Tair,pressure, H, d, z0m)
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
    #
    # estimate z0m
    # need to give zh and zr in addition to many variables in df
    @test_throws Exception wind_profile(df, z, d)    
    windzc3 = wind_profile(df, z, d; zh=thal.zh, zr=thal.zr)    
    # may have used slightly different estimated z0m
    #windzc3 - windzc
    @test all(isapprox.(windzc3, windzc, atol=0.1))
end

