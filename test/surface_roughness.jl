@testset "Reynolds_Number" begin
    Tair,pressure,ustar,z0m = 25.0,100-0,0.5,0.5
    R = @inferred Reynolds_Number(Tair,pressure,ustar,z0m)                             
    @test ≈(R, 15870, rtol=1e-3) 
end

@testset "roughness_parameters" begin
    zh = thal.zh
    zr = thal.zr
    LAI = thal.LAI
    keys_exp = (:d, :z0m, :z0m_se)
    rp = @inferred roughness_parameters(Val(:canopy_height), zh)
    #round.(values(rp); sigdigits = 4)
    @test keys(rp) == keys_exp
    @test all(isapproxm.(values(rp), (18.55, 2.65, missing), rtol=1e-3))
    #
    rp = @inferred roughness_parameters(Val(:canopy_height_LAI), zh, LAI)
    #round.(values(rp); sigdigits = 4)
    @test keys(rp) == keys_exp
    @test all(isapproxm.(values(rp), (21.77, 1.419, missing), rtol=1e-3))
    #
    df = copy(tha48)
    #df.wind[1] = missing # obscures comparing to R
    d=0.7*zh
    psi_m = stability_correction!(df; z=zr, d).psi_m
    # note: must use columntable for type stability - but needs compilation timede
    # not type-stable if some columns allow missings in Julia 1.6
    #rp = @inferred roughness_parameters(Val(:wind_profile), df.ustar, df.wind, psi_m; zh, zr)
    rp = roughness_parameters(Val(:wind_profile), df.ustar, df.wind, psi_m; zh, zr)
    @test keys(rp) == keys_exp
    #@test all(isapproxm.(values(rp), (18.55, 1.879, 0.3561), rtol=1e-3))
    #from R:
    @test all(isapproxm.(values(rp), (18.55, 1.879402, 0.356108), rtol=1e-3))
    #
    # no stability correction
    # broadcast across Missings is not type-stable 
    # rp0 = @inferred roughness_parameters(Val(:wind_profile), df.ustar, df.wind, 
    #     df.Tair, df.pressure, df.H; zh, zr, 
    #     stab_formulation = Val(:no_stability_correction))
    dfd = disallowmissing(df[!,Not(:LE)])
    rp0 = @inferred roughness_parameters(Val(:wind_profile), dfd.ustar, dfd.wind, 
        dfd.Tair, dfd.pressure, dfd.H; zh, zr, 
        stab_formulation = Val(:no_stability_correction))
    @test keys(rp0) == keys_exp
    # same magnitude as with stability correction
    @test all(isapproxm.(values(rp0), values(rp), rtol=0.5))
    # providing DataFrame is not type stable
    rp0b = roughness_parameters(Val(:wind_profile), df; zh, zr, 
        stab_formulation = Val(:no_stability_correction))
    @test rp0b == rp0         
    #
    # estimate psi
    #@code_warntype stability_correction(columntable(df), zr, 0.7*zh)
    #@code_warntype roughness_parameters(Val(:wind_profile), columntable(df), zh, zr)
    rp_psiauto = roughness_parameters(Val(:wind_profile), df; zh, zr)
    @test rp_psiauto == rp
    #
    # DataFrame with only columns ustar and wind
    rp0c = roughness_parameters(Val(:wind_profile), df[!,Cols(:ustar, :wind)]; zh, zr, 
        stab_formulation = Val(:no_stability_correction))
    @test rp0c == rp0         
end

@testset "wind_profile" begin
    datetime, ustar, Tair, pressure, H = values(tha48[1,:])
    z = 30
    d=0.7*thal.zh
    z0m= 2.65
    u30 = @inferred wind_profile(z, ustar, d, z0m, 0.0)
    @test ≈(u30, 1.93, rtol = 1/100 ) # from R
    #
    u30c = @inferred wind_profile(z, ustar, d, z0m, Tair,pressure, H)
    @test ≈(u30c, 2.31, rtol = 1/100 ) # from R
    #
    z0m=1.9 #2.14 #2.65
    u30 = @inferred wind_profile(z, ustar, d, z0m, 0.0) # used below
    u30c = @inferred wind_profile(z, ustar, d, z0m, Tair,pressure, H)
    df = copy(tha48)
    dfd = disallowmissing(df[!,Not(:LE)])
    windz = @inferred wind_profile(z, columntable(dfd), d, z0m; stab_formulation = Val(:no_stability_correction))    
    windz = wind_profile(z, df, d, z0m; stab_formulation = Val(:no_stability_correction))    
    @test length(windz) == 48
    @test windz[1] == u30
    windzc = @inferred wind_profile(z, columntable(dfd), d, z0m; stab_formulation = Val(:Dyer_1970))    
    @test windzc[1] == u30c
    #plot(windz)
    #plot!(windz2)
    #
    # with providing psi_m
    psi_m = stability_correction(df; z, d).psi_m
    windzc2 = @inferred wind_profile(z, columntable(dfd), d, z0m; psi_m)    
    @test windzc2 == windzc
    #
    # with providing L
    MOL = Monin_Obukhov_length.(df.Tair, df.pressure, df.ustar, df.H)
    windzc3 = @inferred wind_profile(z, columntable(dfd), d, z0m; MOL)    
    @test windzc3 == windzc
end

