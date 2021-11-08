@testset "Monin_Obukhov_length" begin
    datetime, ustar, Tair, pressure, H = values(tha48[24,1:5])
    MOL24 = Monin_Obukhov_length(Tair, pressure, ustar, H)
    @test ≈(MOL24, -104.3, rtol = 1/1000)
    df = copy(tha48)
    Monin_Obukhov_length!(df)
    @test df.MOL[24] == MOL24
    MOL = Monin_Obukhov_length(df)
    @test MOL == df.MOL
end

@testset "stability_parameter" begin
    df = copy(tha48)
    MOL24 = first(Monin_Obukhov_length(df[SA[24],:]))
    zr=40;d=15
    zeta = stability_parameter(zr,d,MOL24)
    @test ≈(zeta , -0.240, rtol = 1/100)
    #
    stability_parameter!(df;zr,d)
    @test df.zeta[24] == zeta
    df = copy(tha48)
    zetas = stability_parameter(df;zr,d)
    @test df == tha48 # did not modify original df
    @test zetas[24] == zeta
end

@testset "stability_correction" begin
    zeta = -2:0.5:0.5
    df2 = DataFrame(stability_correction.(zeta))
    @test all(isapprox.(df2.psi_h, SA[2.431,2.197,1.881,1.386,0,-2.5], rtol=1e-3))
    @test all(isapprox.(df2.psi_m, SA[2.275, 2.061, 1.772, 1.317, 0, -2.5], rtol=1e-3))
    df2 = DataFrame(stability_correction.(zeta; stab_formulation=Val(:Businger_1971)))                         
    @test all(isapprox.(df2.psi_h, SA[2.085, 1.862, 1.564, 1.106,0, -3.9], rtol=1e-3))
    @test all(isapprox.(df2.psi_m, SA[2.418, 2.200, 1.904, 1.435,0, -3], rtol=1e-3))
    #
    resm = stability_correction(missing)
    @test keys(resm) == (:psi_h, :psi_m)
    @test all(ismissing.(values(resm)))
    #
    resm = stability_correction(first(zeta); 
        stab_formulation = Val(:no_stability_correction))
    @test resm == (psi_h = 0.0, psi_m = 0.0)
end

@testset "stability_correction DataFrame variant" begin
    zr=40;d=15
    dfo = DataFrame(Tair=25, pressure=100, ustar=0.2:0.1:1.0, H=40:20:200)
    df = copy(dfo)
    stability_correction!(df, zr, d)
    propertynames(df)[(end-1):end] == SA[:psi_h, :psi_m]
    #
    dfm = copy(dfo)
    stability_correction!(dfm, zr, d; stab_formulation = Val(:no_stability_correction))
    propertynames(dfm)[(end-1):end] == SA[:psi_h, :psi_m]
    @test all(iszero.(dfm.psi_h))
    @test all(iszero.(dfm.psi_m))
    #
    df2 = copy(dfo)
    stability_parameter!(df2; zr, d) # adds zeta
    stability_correction!(df2)
    @test df2.psi_h == df.psi_h
    @test df2.psi_m == df.psi_m
end

@testset "stability_correction from raw" begin
    datetime, ustar, Tair, pressure, H = values(tha48[24,1:5])
    z=40.0;d=15.0
    resm = stability_correction(Tair,pressure,ustar,H, z,d)
    @test resm.psi_h ≈ 0.940 rtol=1/1000
    @test resm.psi_m ≈ 0.902 rtol=1/1000
    #
    resm0 = stability_correction(Tair,pressure,ustar,H, z,d; 
        stab_formulation = Val(:no_stability_correction))
    @test resm0 == (psi_h = 0.0, psi_m = 0.0)
    #
    df = copy(tha48)
    df.ustar[3] = missing
    stability_correction!(df,z,d)
    @test all(ismissing.((df.psi_h[3], df.psi_m[3])))
    @test all(isapprox.((df.psi_h[24], df.psi_m[24]),values(resm)))
end