@testset "add_Ga" begin
    Gb_h = 0.0347
    Ga_m = 0.3
    @test @inferred add_Ga(Gb_h, Ga_m) == NamedTuple()
    #
    # here cannot @inferred the names in the return
    # TODO
    Ga = add_Ga(Gb_h, Ga_m, :N20 => 2, :CH4 => 4)
    @test keys(Ga) == (:Ga_N20, :Ga_CH4)
    @test ≈(Ga.Ga_N20, 0.0164, rtol = 1/100)
    @test ≈(Ga.Ga_CH4, 0.0105, rtol = 1/100)
    #
    Gam = add_Ga(Gb_h, missing, :N20 => 2, :CH4 => 4)
    @test keys(Gam) == (:Ga_N20, :Ga_CH4)
    @test all(ismissing.(values(Gam)))
    # 
    # DataFrame variant
    dfo = DataFrame(Gb_h=SA[Gb_h, missing, 0.055], Ga_m = SA[0.3,0.3,0.3])
    df = copy(dfo)
    df2 = @inferred add_Gb!(df)
    @test isequal(df2, dfo)
    df2 = @inferred add_Ga!(df, :N20 => 2, :CH4 => 4)
    @test df2 === df # mutating
    @test propertynames(df2)[end-1:end] == [:Ga_N20, :Ga_CH4]
    @test df2.Ga_N20[1] == Ga.Ga_N20
end

@testset "compute_Ram" begin
    zr, zh, d = thal.zr, thal.zh, 0.7*thal.zh
    ustar, wind = tha48[24, [:ustar, :wind]]
    df = copy(tha48)
    stability_correction!(df; z=zr, d)
    z0m = roughness_parameters(Roughness_wind_profile(), df; zh, zr, psi_m = df.psi_m).z0m
    Ram_r = @inferred compute_Ram(ResistanceWindZr(), ustar, wind)
    @test Ram_r ≈ 6.31 rtol = 1/100 # TODO check with R
    Ram_p =  @inferred compute_Ram(
        ResistanceWindProfile(), ustar; zr, d, z0m, psi_h = df.psi_h[24])
    @test Ram_p ≈ 5.41 rtol = 1/100 # TODO check with R
    #
    #Gb_h, Gb_CO2 = compute_Gb!(df, Thom1972())[24, [:Gb_h, :Gb_CO2]]
    @inferred compute_Ram!(df, ResistanceWindZr())
    #@test all(propertynames(df)[(end+1-length(Ram)):end] .== keys(Ram))
    @test propertynames(df)[end] == :Ra_m
    @test df.Ra_m[24] == Ram_r
    @inferred compute_Ram!(df, ResistanceWindProfile(); zr, d, z0m)
    @test df.Ra_m[24] == Ram_p
    Ram_skalar = df.Ra_m
    #
    # test zr as a vecttor
    df[!,:zri] .= zr
    df.zri[1] = zr/2
    @inferred compute_Ram!(df, ResistanceWindProfile(); zr=df.zri, d, z0m)
    @test df.Ra_m[2:end] == Ram_skalar[2:end]
    @test df.Ra_m[1] != Ram_skalar[1]
end

@testset "roughness_length_heat" begin
    z0m, kB_h = 1.2, 3.4
    z0h = @inferred roughness_length_heat(z0m, kB_h)
    @test z0h ≈ 0.0400 rtol = 1e-2
end

@testset "aerodynamic_conductance! only ustar and wind" begin
    # using only ustar and wind
    dfo = DataFrame(ustar = [0.1,missing,0.3], wind = [3.4, 2.8, 3.3])
    df = copy(dfo)
    @inferred aerodynamic_conductance!(df; Gb_model=Thom1972())
    # missing due to missing zr
    @test all(ismissing.(df.psi_h))        
    @test all(ismissing.(df.psi_m))        
    @test all(ismissing.(df.z0h))        
    @inferred aerodynamic_conductance!(
        df; Gb_model=Thom1972(), zr = thal.zr, zh = thal.zh,
        stab_formulation = NoStabilityCorrection())
    @test all(iszero.(df.psi_h))        
    @test all(iszero.(df.psi_m))        
    @test all(ismissing.(df.z0h))        
    @test propertynames(df)[(end-9):end] == 
        SA[:Gb_h, :Rb_h, :kB_h, :Gb_CO2, :Ra_m, :Ga_m, :Ga_h, :Ra_h, :Ga_CO2, :z0h]
    # compare with R        
    @test all(isapproxm.(df.Gb_h, (0.0347, missing, 0.0723), rtol = 1e-3))
    @test all(isapproxm.(df.Ga_m, (0.00294, missing, 0.0273), rtol = 1e-3))
    @test all(isapproxm.(df.Ga_h, (0.00271, missing, 0.0198), rtol = 1e-3))
end

@testset "aerodynamic_conductance! Gb_Thom" begin
    df = tha48[:,Not(Cols(:Gb_h, :Ga_h))]
    # test Ga_m by wind_profile
    @inferred aerodynamic_conductance!(
        df; Gb_model=Thom1972(), Ram_model=ResistanceWindProfile(), 
        zh = thal.zh, zr = thal.zr)
    @test propertynames(df)[(end-9):end] == 
        SA[:Gb_h, :Rb_h, :kB_h, :Gb_CO2, :Ra_m, :Ga_m, :Ga_h, :Ra_h, :Ga_CO2, :z0h]
end

@testset "aerodynamic_conductance! constant_kB1" begin
    # specify kB1 as a vecttor
    df = tha48[:,Not(Cols(:Gb_h, :Ga_h))]
    df[!,:kB_hi] .= 1.18
    df.kB_hi[1] = 1.18/2
    @inferred aerodynamic_conductance!(df; Gb_model=ConstantKB1(), kB_h = df.kB_hi, 
        zh = thal.zh, zr = thal.zr)
    @test propertynames(df)[(end-9):end] == 
        SA[:Gb_h, :Rb_h, :kB_h, :Gb_CO2, :Ra_m, :Ga_m, :Ga_h, :Ra_h, :Ga_CO2, :z0h]
    @test all(isapproxm.(df.Gb_h[1:3], (0.375, 0.17, 0.167), rtol = 1e-2))
end

@testset "aerodynamic_conductance! Gb_Choudhury" begin
    leafwidth=0.1
    df = tha48[:,Not(Cols(:Gb_h, :Ga_h))]
    # test Ga_m by wind_profile
    @inferred aerodynamic_conductance!(
        df; Gb_model=Choudhury1988(), Ram_model=ResistanceWindProfile(), 
        leafwidth, LAI=thal.LAI, zh = thal.zh, zr = thal.zr)
    @test propertynames(df)[(end-9):end] == 
        SA[:Gb_h, :Rb_h, :kB_h, :Gb_CO2, :Ra_m, :Ga_m, :Ga_h, :Ra_h, :Ga_CO2, :z0h]
    # compare with R        
    @test all(isapproxm.(df.Gb_h[1:3], (0.157, 0.15, 0.151), rtol = 1e-2))
    @test all(isapproxm.(df.Ga_m[1:3], (0.0709, 0.0649, 0.0603), rtol = 1e-2))
end

@testset "aerodynamic_conductance! Gb_Su" begin
    Dl=0.01
    df = tha48[:,Not(Cols(:Gb_h, :Ga_h))]
    df[!,:LAI] .= thal.LAI
    df.LAI[1] = thal.LAI/2
    # test providing changing LAI
    @inferred aerodynamic_conductance!(df; Gb_model=Su2001(), 
        Dl, LAI=df.LAI, zh=thal.zh, zr=thal.zr);
    @test propertynames(df)[(end-9):end] == 
        SA[:Gb_h, :Rb_h, :kB_h, :Gb_CO2, :Ra_m, :Ga_m, :Ga_h, :Ra_h, :Ga_CO2, :z0h]
    # compare with R
    @test all(isapproxm.(df.Gb_h[1:3], (0.202, 0.177, 0.167), rtol = 1e-2))
    @test all(isapproxm.(df.Ga_m[1:3], (0.0693, 0.0538, 0.0507), rtol = 1e-2))
end

