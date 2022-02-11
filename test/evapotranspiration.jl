@testset "potential_ET scalars" begin
    Tair,pressure,Rn = 30.0,100.0,500.0
    ET_pot, LE_pot = potential_ET(PriestleyTaylor(), Tair,pressure,Rn; alpha=1.26)    
    # compare to R
    @test ≈(ET_pot, 0.0002035969; rtol = 1e-5)
    @test ≈(LE_pot, 494.7202; rtol = 1e-5)
    #
    VPD,Ga_h = 2.0, 0.1
    Gs_pot = 0.5 # assume known surface conductance
    ET_pot,LE_pot = potential_ET(PenmanMonteith(), Tair,pressure,Rn,VPD,Ga_h; Gs_pot)    
    Gs_ms, Gs_mol = surface_conductance(InversePenmanMonteith(), Tair,pressure,VPD,LE_pot,Rn,Ga_h)
    @test Gs_mol ≈ Gs_pot
end

@testset "potential_ET dataframe" begin
    dfo = DataFrame(Tair = 20.0:1.0:30.0,pressure = 100.0, Rn = 500.0)
    df = copy(dfo)
    potential_ET!(df, PriestleyTaylor(); alpha=1.26)    
    @test ≈(last(df).ET_pot, 0.0002035969; rtol = 1e-5)
    @test ≈(last(df).LE_pot, 494.7202; rtol = 1e-5)
    #
    df2 = copy(dfo)
    df2[!, :VPD] .= 2.0
    df2[!, :Ga_h] .= 0.1
    potential_ET!(df2, PenmanMonteith(); Gs_pot = 0.4)  
    df2b = rename!(copy(df2, copycols=false),SA[:LE_pot => :LE])  
    surface_conductance!(df2b, InversePenmanMonteith())
    @test all(df2b.Gs_mol .≈ 0.4)
end

@testset "potential_ET dataframe missings in G" begin
    df = @pipe DataFrame(Tair = 20.0,pressure = 100.0, Rn = 500.0, G = 100.0:1:110) |>
      allowmissing(_, Cols(:G))
    df.G[2] = missing
    df_ET = potential_ET!(copy(df), PriestleyTaylor(); G = df.G)    
    @test ncol(df_ET) == ncol(df)+2 # two columns added: ET_pot, LE_pot
    @test nrow(df_ET) == nrow(df)
    @test df_ET.ET_pot[end] < df_ET.ET_pot[1] # smaller: mith higher G less available energy
    @test df_ET.ET_pot[1] < 0.0002035969 # smaller as without G because less energy available
    @test ismissing(df_ET.LE_pot[2])
    #
    df2 = copy(df)
    df2[!, :VPD] .= 2.0
    df2[!, :Ga_h] .= 0.1
    potential_ET!(df2, PenmanMonteith(), S = df.G; Gs_pot = 0.3)    
    @test ncol(df2) == ncol(df)+2+2 # two columns added: ET_pot, LE_pot
    @test nrow(df2) == nrow(df)
    @test ismissing(df2.LE_pot[2])
    df2b = rename!(copy(df2, copycols=false),SA[:LE_pot => :LE])  
    surface_conductance!(df2b, InversePenmanMonteith(); S = df.G)
    @test ismissing(df2b.Gs_mol[2])
    @test all(df2b.Gs_mol[1:end .!= 2] .≈ 0.3)
end

@testset "equilibrium_imposed_ET scalars" begin
    # regression from R package example
    Tair,pressure,Rn, VPD, Gs = 20.0,100.0,50.0, 0.5, 0.01
    ET_eq, ET_imp, LE_eq, LE_imp = equilibrium_imposed_ET(Tair,pressure,VPD,Gs, Rn)    
    @test ≈(ET_eq, 1.399424e-05; rtol = 1e-5)
    @test ≈(ET_imp, 3.695727e-05; rtol = 1e-5)
    @test ≈(LE_eq, 34.33628; rtol = 1e-5)
    @test ≈(LE_imp, 90.67837; rtol = 1e-5)
    #
    df = DataFrame(Tair = 20.0:1.0:22.0, pressure = 100.0, Rn = 50.0,  VPD = 0.5,  Gs = 0.01)
    ncoldf0 = ncol(df)
    # df_ET = @test_logs (:info,r"G is not provided") (:info,r"S is not provided") equilibrium_imposed_ET(df)
    # @test ncol(df) == ncoldf0
    # @test nrow(df_ET) == nrow(df)
    # @test ≈(first(df_ET).ET_eq, ET_eq)
    # @test ≈(first(df_ET).ET_imp, ET_imp)
    #
    dfm = copy(df)
    equilibrium_imposed_ET!(dfm)
    @test ncol(dfm) == ncoldf0 + 4
    @test nrow(dfm) == nrow(df)
    @test ≈(first(dfm).ET_eq, ET_eq)
    @test ≈(first(dfm).ET_imp, ET_imp)
end