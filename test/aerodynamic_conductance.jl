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




@testset "compute_Gb :no_stability_correction" begin
    dfo = DataFrame(ustar = SA[0.1,missing,0.3], wind = [3.4, 2.8, 3.3])
    df = copy(dfo)
    aerodynamic_conductance!(df, Val(:Thom_1972); zh = thal.zh, zr = thal.zr,
        stab_formulation = Val(:no_stability_correction))
    @test all(iszero.(df.psi_h))        
    @test all(iszero.(df.psi_m))        
    @test propertynames(df)[(end-3):end] == SA[:Rb_h, :Gb_h, :kB_h, :Gb_CO2]
end


@testset "compute_Gb Gb_Thom" begin
    df = copy(tha48)
    aerodynamic_conductance!(df, Val(:Thom_1972); zh = thal.zh, zr = thal.zr)
    @test propertynames(df)[(end-3):end] == SA[:Rb_h, :Gb_h, :kB_h, :Gb_CO2]
end

@testset "constant_kB1" begin
    kB_h = 1.18
    # DataFrame variant
    df = copy(tha48)
    aerodynamic_conductance!(df, Val(:constant_kB1); kB_h, zh = thal.zh, zr = thal.zr)
    @test propertynames(df)[(end-3):end] == SA[:Rb_h, :Gb_h, :kB_h, :Gb_CO2]
end

@testset "compute_Gb Gb_Choudhury" begin
    leafwidth=0.1
    df = copy(tha48)
    aerodynamic_conductance!(df, Val(:Choudhury_1988); 
        leafwidth, LAI=thal.LAI, zh = thal.zh, zr = thal.zr)
    @test propertynames(df)[(end-3):end] == SA[:Rb_h, :Gb_h, :kB_h, :Gb_CO2]
end

@testset "compute_Gb Gb_Su" begin
    Dl=0.01
    df = copy(tha48)
    aerodynamic_conductance!(df, Val(:Su_2001); Dl, LAI=thal.LAI, zh=thal.zh, zr=thal.zr)
    @test propertynames(df)[(end-3):end] == SA[:Rb_h, :Gb_h, :kB_h, :Gb_CO2]
end

