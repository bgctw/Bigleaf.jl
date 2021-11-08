@testset "add_Ga" begin
    Gb_h = 0.0347
    Ga_m = 0.3
    @test add_Ga(Gb_h, Ga_m) == NamedTuple()
    #
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
    df2 = add_Gb!(df)
    @test isequal(df2, dfo)
    df2 = add_Ga!(df, :N20 => 2, :CH4 => 4)
    @test df2 === df # mutating
    @test propertynames(df2)[end-1:end] == [:Ga_N20, :Ga_CH4]
    @test df2.Ga_N20[1] == Ga.Ga_N20
end

@testset "compute_Ram" begin
    zr, zh, d = thal.zr, thal.zh, 0.7*thal.zh
    ustar, wind = tha48[24, [:ustar, :wind]]
    df = copy(tha48)
    stability_parameter!(df; zr, d)
    stability_correction!(df)
    z0m = roughness_parameters(Val(:wind_profile), df, zh, zr; psi_m = df.psi_m).z0m
    Ram_r = compute_Ram(Val(:wind_zr), ustar, wind)
    @test Ram_r ≈ 6.31 rtol = 1/100 # TODO check with R
    Ram_p = compute_Ram(Val(:wind_profile), ustar; zr, d, z0m, psi_h = df.psi_h[24])
    @test Ram_p ≈ 5.41 rtol = 1/100 # TODO check with R

    #Gb_h, Gb_CO2 = compute_Gb!(df, Val(:Thom_1972))[24, [:Gb_h, :Gb_CO2]]
    compute_Ram!(df, Val(:wind_zr))
    #@test all(propertynames(df)[(end+1-length(Ram)):end] .== keys(Ram))
    @test propertynames(df)[end] == :Ra_m
    @test df.Ra_m[24] == Ram_r
    compute_Ram!(df, Val(:wind_profile); zr, d, z0m)
    @test df.Ra_m[24] == Ram_p
end

@testset "aerodynamic_conductance! only ustar and wind" begin
    # using only ustar and wind
    dfo = DataFrame(ustar = SA[0.1,missing,0.3], wind = [3.4, 2.8, 3.3])
    df = copy(dfo)
    aerodynamic_conductance!(df; Gb_model=Val(:Thom_1972))
    # missing due to missing zr
    @test all(ismissing.(df.psi_h))        
    @test all(ismissing.(df.psi_m))        
    aerodynamic_conductance!(df; Gb_model=Val(:Thom_1972), zr = thal.zr, zh = thal.zh
        ,stab_formulation = Val(:no_stability_correction))
    @test all(iszero.(df.psi_h))        
    @test all(iszero.(df.psi_m))        
    @test propertynames(df)[(end-8):end] == 
        SA[:Rb_h, :Gb_h, :kB_h, :Gb_CO2, :Ra_m, :Ga_m, :Ga_h, :Ra_h, :Ga_CO2]
    # TODO compare with R        
end

@testset "compute_Gb Gb_Thom" begin
    df = copy(tha48)
    aerodynamic_conductance!(df; Gb_model=Val(:Thom_1972), zh = thal.zh, zr = thal.zr)
    @test propertynames(df)[(end-8):end] == 
        SA[:Rb_h, :Gb_h, :kB_h, :Gb_CO2, :Ra_m, :Ga_m, :Ga_h, :Ra_h, :Ga_CO2]
end

@testset "constant_kB1" begin
    kB_h = 1.18
    # DataFrame variant
    df = copy(tha48)
    aerodynamic_conductance!(df; Gb_model=Val(:constant_kB1), kB_h, zh = thal.zh, zr = thal.zr)
    @test propertynames(df)[(end-8):end] == 
        SA[:Rb_h, :Gb_h, :kB_h, :Gb_CO2, :Ra_m, :Ga_m, :Ga_h, :Ra_h, :Ga_CO2]
end

@testset "compute_Gb Gb_Choudhury" begin
    leafwidth=0.1
    df = copy(tha48)
    aerodynamic_conductance!(df; Gb_model=Val(:Choudhury_1988),
        leafwidth, LAI=thal.LAI, zh = thal.zh, zr = thal.zr)
    @test propertynames(df)[(end-8):end] == 
        SA[:Rb_h, :Gb_h, :kB_h, :Gb_CO2, :Ra_m, :Ga_m, :Ga_h, :Ra_h, :Ga_CO2]
end

@testset "compute_Gb Gb_Su" begin
    Dl=0.01
    df = copy(tha48)
    aerodynamic_conductance!(df; Gb_model=Val(:Su_2001), Dl, LAI=thal.LAI, zh=thal.zh, zr=thal.zr)
    @test propertynames(df)[(end-8):end] == 
        SA[:Rb_h, :Gb_h, :kB_h, :Gb_CO2, :Ra_m, :Ga_m, :Ga_h, :Ra_h, :Ga_CO2]
end

