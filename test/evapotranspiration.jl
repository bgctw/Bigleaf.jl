@testset "potential_ET scalars" begin
    Tair,pressure,Rn = 30.0,100.0,500.0
    ET_pot, LE_pot = potential_ET(Tair,pressure,Rn, Val(:PriestleyTaylor); alpha=1.26)    
    @test ≈(ET_pot, 0.0002035969; rtol = 1e-5)
    @test ≈(LE_pot, 494.7202; rtol = 1e-5)
    #
    VPD,Ga = 2.0, 0.1
    ET_pot,LE_pot = potential_ET(Tair,pressure,Rn,VPD,Ga, Val(:PenmanMonteith))    
    #
    #TODO surface_conductance(Tair=20,pressure=100,VPD=2,Ga=0.1,Rn=400,LE=LE_pot_PM)
end

@testset "potential_ET dataframe" begin
    df = DataFrame(
    Tair = 20.0:1.0:30.0,pressure = 100.0, Rn = 500.0)
    df_ET = @test_logs (:info,r"G is not provided") (:info,r"S is not provided") potential_ET(df, Val(:PriestleyTaylor); alpha=1.26)    
    # non-mutating
    @test ncol(df) == 3
    @test nrow(df_ET) == nrow(df)
    @test ≈(last(df_ET).ET_pot, 0.0002035969; rtol = 1e-5)
    @test ≈(last(df_ET).LE_pot, 494.7202; rtol = 1e-5)
    #
    df2 = copy(df)
    df2.VPD .= 2.0
    df2.Ga .= 0.1
    # df2 = transform(df,  
    #     [] => ByRow(() -> 2.0) => :VPD,
    #     [] => ByRow(() -> 0.1) => :Ga,
    #     )
    df_ET2 = @test_logs (:info,r"G is not provided") (:info,r"S is not provided") potential_ET(df2, Val(:PenmanMonteith))    
    @test ncol(df2) == 5
    @test nrow(df_ET2) == nrow(df2)
    #TODO surface_conductance(Tair=20,pressure=100,VPD=2,Ga=0.1,Rn=400,LE=LE_pot_PM)
    #
    # mutating
    dfm = copy(df)
    #@test_throws Exception 
    @test_logs (:info,r"G is not provided") (:info,r"S is not provided") potential_ET!(dfm, Val(:PriestleyTaylor); alpha=1.26)   
    #dfm = @test_logs (:info,r"G is not provided") (:info,r"S is not provided") hcat(dfm, Bigleaf.fill_GS_missings(dfm, missing, missing, false, false); copycols = false)        
    #potential_ET!(dfm, Val(:PriestleyTaylor); alpha=1.26)   
    @test ncol(dfm) == 3+2 # two columns added: ET_pot, LE_pot
    dfm = @test_logs (:info,r"G is not provided") (:info,r"S is not provided")  potential_ET!(copy(df2), Val(:PenmanMonteith))   
    @test ncol(dfm) == ncol(df2)+2 
end

@testset "potential_ET dataframe missings in G" begin
    df = @pipe DataFrame(Tair = 20.0:1.0:30.0,pressure = 100.0, Rn = 500.0, G = 105.0) |>
      allowmissing(_, Cols(:G))
    df.G[1] = missing
    df_ET = @test_logs (:info,r"S is not provided") potential_ET(df, Val(:PriestleyTaylor); G = df.G)    
    # non-mutating
    @test ncol(df) == 4
    @test nrow(df_ET) == nrow(df)
    @test last(df_ET).ET_pot < 0.0002035969 # smaller because less energy available
    @test ismissing(first(df_ET).LE_pot)
    #
    df2 = copy(df)
    df2.VPD .= 2.0
    df2.Ga .= 0.1
    df_ET2 = @test_logs (:info,r"G is not provided") potential_ET(df2, Val(:PenmanMonteith), S = df.G)    
    @test ncol(df2) == 6
    @test nrow(df_ET2) == nrow(df2)
    @test ismissing(first(df_ET).LE_pot)
    #TODO surface_conductance(Tair=20,pressure=100,VPD=2,Ga=0.1,Rn=400,LE=LE_pot_PM)
    #
    # mutating
    dfm = copy(df)
    #@test_throws Exception 
    @test_logs (:info,r"G is not provided") potential_ET!(dfm, Val(:PriestleyTaylor); S = dfm.G)   
    #dfm = @test_logs (:info,r"G is not provided") (:info,r"S is not provided") hcat(dfm, Bigleaf.fill_GS_missings(dfm, missing, missing, false, false); copycols = false)        
    #potential_ET!(dfm, Val(:PriestleyTaylor); alpha=1.26)   
    @test ncol(dfm) == ncol(df)+2 # two columns added: ET_pot, LE_pot
    dfm = @test_logs (:info,r"S is not provided")  potential_ET!(copy(df2), Val(:PenmanMonteith); G = dfm.G)   
    @test ncol(dfm) == ncol(df2)+2 
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
    df_ET = @test_logs (:info,r"G is not provided") (:info,r"S is not provided") equilibrium_imposed_ET(df)
    @test ncol(df) == ncoldf0
    @test nrow(df_ET) == nrow(df)
    @test ≈(first(df_ET).ET_eq, ET_eq)
    @test ≈(first(df_ET).ET_imp, ET_imp)
    #
    dfm = copy(df)
    dfm_ET = equilibrium_imposed_ET!(dfm; infoGS = false)
    @test ncol(dfm) == ncoldf0 + 4
    @test nrow(dfm) == nrow(df)
    @test ≈(first(dfm_ET).ET_eq, ET_eq)
    @test ≈(first(dfm_ET).ET_imp, ET_imp)
    #TODO surface_conductance(Tair=20,pressure=100,VPD=2,Ga=0.1,Rn=400,LE=LE_pot_PM)
end