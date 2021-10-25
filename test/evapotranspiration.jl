@testset "potential_ET scalars" begin
    Tair,pressure,Rn = 30.0,100.0,500.0
    ET_pot, LE_pot = @test_logs (:info,r"G is not provided") (:info,r"S is not provided") potential_ET(Tair,pressure,Rn, Val(:PriestleyTaylor); alpha=1.26)    
    @test ≈(ET_pot, 0.0002035969; rtol = 1e-5)
    @test ≈(LE_pot, 494.7202; rtol = 1e-5)
    #
    VPD,Ga = 2.0, 0.1
    ET_pot,LE_pot = @test_logs (:info,r"G is not provided") (:info,r"S is not provided") potential_ET(Tair,pressure,Rn,VPD,Ga, Val(:PenmanMonteith))    
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
    @test_throws Exception potential_ET!(
        dfm, Val(:PriestleyTaylor); alpha=1.26)   
    dfm = @test_logs (:info,r"G is not provided") (:info,r"S is not provided") hcat(dfm, Bigleaf.fill_GS_missings(dfm, missing, missing, false, false); copycols = false)        
    potential_ET!(dfm, Val(:PriestleyTaylor); alpha=1.26)   
    @test ncol(dfm) == 3+4 # four columns added: G,S, ET_pot, LE_pot
    dfm = @test_logs (:info,r"G is not provided") (:info,r"S is not provided")  fill_GS_missings!(copy(df2), missing, missing)
    potential_ET!(dfm, Val(:PenmanMonteith))   
    @test ncol(dfm) == ncol(df2)+4 # four columns added: G,S, ET_pot, LE_pot
end

