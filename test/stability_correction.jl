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
end
