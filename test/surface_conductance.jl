@testset "FluxGradient" begin
    ir = 24
    Tair,pressure,VPD,LE = tha48[ir, Cols(:Tair,:pressure,:VPD,:LE)]
    Gs24 = surface_conductance(Val(:FluxGradient), Tair,pressure,VPD,LE)
    @test keys(Gs24) == (:Gs_ms, :Gs_mol)
    # compare with R
    @test all(isapprox.(values(Gs24), (0.00919248, 0.3751482), rtol=1e-3))
    #
    df = copy(tha48)
    surface_conductance!(df, Val(:FluxGradient))
    @test propertynames(df)[(end-1):end] == SA[:Gs_ms, :Gs_mol]
    @test df[ir, Cols(:Gs_ms, :Gs_mol)] == Gs24
end

@testset "PenmanMonteith" begin
    ir = 24
    Tair,pressure,VPD,LE,Rn,Ga_h,G = tha48[ir, Cols(:Tair,:pressure,:VPD,:LE, :Rn,:Ga_h,:G)]
    Gs24 = surface_conductance(Val(:PenmanMonteith), Tair,pressure,VPD,LE,Rn,Ga_h;G)
    @test keys(Gs24) == (:Gs_ms, :Gs_mol)
    # compare with R
    @test all(isapprox.(values(Gs24), (0.006846274, 0.2793988), rtol=1e-3))
    #
    df = copy(tha48)
    surface_conductance!(df, Val(:PenmanMonteith); G=df.G)
    @test propertynames(df)[(end-1):end] == SA[:Gs_ms, :Gs_mol]
    @test df[ir, Cols(:Gs_ms, :Gs_mol)] == Gs24
end