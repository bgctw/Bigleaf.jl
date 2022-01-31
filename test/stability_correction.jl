@testset "Monin_Obukhov_length" begin
    datetime, ustar, Tair, pressure, H = tha48[24,Cols(:datetime, :ustar, :Tair, :pressure, :H)]
    MOL24 = @inferred Monin_Obukhov_length(Tair, pressure, ustar, H)
    MOL24 = @inferred Monin_Obukhov_length(Tair, pressure, ustar, H; constants =BigleafConstants())
    @test ≈(MOL24, -104.3, rtol = 1/1000)
    #
    df = copy(tha48)
    @inferred Monin_Obukhov_length!(df)
    @test df.MOL[24] == MOL24
end

@testset "stability_parameter" begin
    df = copy(tha48)
    df24 = df[24,:]
    MOL24 = Monin_Obukhov_length(df24.Tair, df24.pressure, df24.ustar, df24.H)
    z=40;d=15
    zeta = @inferred stability_parameter(z,d,MOL24)
    @test ≈(zeta , -0.240, rtol = 1/100)
    zeta2 = @inferred stability_parameter(z,d,df24.Tair, df24.pressure, df24.ustar, df24.H)
    @test zeta2 == zeta
    #
    @inferred stability_parameter!(df;z,d)
    @test df.zeta[24] == zeta
    zeta_scalar = df.zeta
    # 
    # also test zr as a vector
    df = copy(tha48)
    df[!,:zi] .= z
    df.zi[1] = z/2
    zetas = stability_parameter!(df; z=df.zi, d).zeta
    @test zetas[2:24] == zeta_scalar[2:24]
    @test zetas[1] != zeta_scalar[1]
    #
    # specifyngn MOL directly
    df4 = copy(tha48)
    Monin_Obukhov_length!(df4) 
    zetas2 = stability_parameter!(df4; z, d, MOL = df4.MOL./2).zeta
    @test all(zetas2[2:end] .== zetas[2:end] .* 2)
end

@testset "stability_correction zeta" begin
    zeta = -2:0.5:0.5
    @inferred stability_correction(first(zeta))
    df2 = DataFrame(stability_correction.(zeta))
    @test all(isapprox.(df2.psi_h, SA[2.431,2.197,1.881,1.386,0,-2.5], rtol=1e-3))
    @test all(isapprox.(df2.psi_m, SA[2.275, 2.061, 1.772, 1.317, 0, -2.5], rtol=1e-3))
    df2 = DataFrame(stability_correction.(zeta; stab_formulation=Val(:Businger_1971)))                         
    @test all(isapprox.(df2.psi_h, SA[2.085, 1.862, 1.564, 1.106,0, -3.9], rtol=1e-3))
    @test all(isapprox.(df2.psi_m, SA[2.418, 2.200, 1.904, 1.435,0, -3], rtol=1e-3))
    #
    # https://bitbucket.org/juergenknauer/bigleaf/issues/8/logical-inconsistency-possible-bug-with
    resm = @inferred stability_correction(missing)
    @test keys(resm) == (:psi_h, :psi_m)
    @test all(ismissing.(values(resm)))
    #
    resm = @inferred stability_correction(first(zeta); 
        stab_formulation = Val(:no_stability_correction))
    @test resm == (psi_h = 0.0, psi_m = 0.0)
end

@testset "stability_correction metvars" begin
    z=40;d=15
    df = DataFrame(Tair=25.0, pressure=100.0, ustar=0.2:0.1:1.0, H=40:20.0:200)
    df1 = df[1,:]
    zeta1 = stability_parameter(z,d,df1.Tair, df1.pressure, df1.ustar, df1.H)
    res1 = @inferred stability_correction(z,d, df1.Tair, df1.pressure, df1.ustar, df1.H) 
    @test res1 == stability_correction(zeta1)
end

@testset "stability_correction DataFrame variant" begin
    z=40;d=15
    dfo = DataFrame(Tair=25.0, pressure=100.0, ustar=0.2:0.1:1.0, H=40:20.0:200)
    df = copy(dfo)
    res = stability_correction(df; z, d)
    # for type stability, use columntable(df)
    stability_correction!(df; z, d)
    propertynames(df)[(end-1):end] == SA[:psi_h, :psi_m]
    @test DataFrame(res) == df[!, (end-1):end]
    #
    dfm = allowmissing(dfo)
    dfm.ustar[1] = missing
    res = stability_correction(
        dfm; z, d, stab_formulation = Val(:no_stability_correction))
    stability_correction!(
        dfm; z, d, stab_formulation = Val(:no_stability_correction))
    propertynames(dfm)[(end-1):end] == SA[:psi_h, :psi_m]
    @test all(iszero.(dfm.psi_h))
    @test all(iszero.(dfm.psi_m))
    @test DataFrame(res) == dfm[!, (end-1):end]
    #
    # test specifying zeta instead of z and d
    df2 = copy(dfo)
    stability_parameter!(df2; z, d) # adds zeta
    stability_correction!(df2, zeta=df2.zeta)
    @test df2.psi_h == df.psi_h
    @test df2.psi_m == df.psi_m
    #
    # test specifying zr as a vector
    df3 = copy(dfo)
    df3[!,:zi] .= z
    df3.zi[1] = z/2
    res = stability_correction(df3; z=df3.zi, d)
    @inferred stability_correction!(df3; z=df3.zi, d)
    @test df3.psi_h[2:end] == df.psi_h[2:end]
    @test df3.psi_h[1] != df.psi_h[1]
    @test DataFrame(res) == df3[!, (end-1):end]
end

@testset "stability_correction DataFrame variant Float32" begin
    z=40;d=15
    dfo = DataFrame(Tair=25.0, pressure=100.0, ustar=0.2:0.1:1.0, H=40:20.0:200)
    float_cols = names(dfo, All())
    dfo = transform(dfo,  float_cols .=> ByRow(passmissing(Float32)) .=> float_cols)
    
    df = copy(dfo)
    res = stability_correction(df; z, d)
    @test eltype(res.psi_h) == Float32
    # for type stability, use columntable(df)
    stability_correction!(df; z, d)
    propertynames(df)[(end-1):end] == SA[:psi_h, :psi_m]
    @test eltype(df.psi_h) == Float32
    @test DataFrame(res) == df[!, (end-1):end]
    #
    dfm = allowmissing(dfo)
    dfm.ustar[1] = missing
    res = stability_correction(
        dfm; z, d, stab_formulation = Val(:no_stability_correction))
    @test eltype(res.psi_h) == Float32
    stability_correction!(
        dfm; z, d, stab_formulation = Val(:no_stability_correction))
    propertynames(dfm)[(end-1):end] == SA[:psi_h, :psi_m]
    @test all(iszero.(dfm.psi_h))
    @test all(iszero.(dfm.psi_m))
    #@test eltype(dfm.psi_h) == typeof(z) #Float32
    @test DataFrame(res) == dfm[!, (end-1):end]
    #
    # test specifying zeta instead of z and d
    df2 = copy(dfo)
    stability_parameter!(df2; z, d) # adds zeta
    stability_correction!(df2, zeta=df2.zeta)
    @test eltype(df2.psi_h) == Float32
    @test df2.psi_h == df.psi_h
    @test df2.psi_m == df.psi_m
    #
    # test specifying zr as a vector
    df3 = copy(dfo)
    df3[!,:zi] .= z
    df3.zi[1] = z/2
    df3_before = copy(df3)
    res = stability_correction(df3; z=df3.zi, d)
    @test df3 == df3_before
    @test eltype(res.psi_h) == Float32
    @inferred stability_correction!(df3; z=df3.zi, d)
    @test eltype(df3.psi_h) == Float32
    @test df3.psi_h[2:end] == df.psi_h[2:end]
    @test df3.psi_h[1] != df.psi_h[1]
    @test DataFrame(res) == df3[!, (end-1):end]
end



