@testset "add_Gb" begin
    Gb_h = 0.0347
    @test add_Gb(Gb_h) == NamedTuple()
    #
    Gb = add_Gb(Gb_h, :Gb_N20 => 2, :Gb_CH4 => 4)
    @test keys(Gb) == (:Gb_N20, :Gb_CH4)
    @test ≈(Gb.Gb_N20, 0.0173, rtol = 1/100)
    @test ≈(Gb.Gb_CH4, 0.0109, rtol = 1/100)
    #
    Gb = add_Gb(missing, :Gb_N20 => 2, :Gb_CH4 => 4)
    @test keys(Gb) == (:Gb_N20, :Gb_CH4)
    @test all(ismissing.(values(Gb)))
end

@testset "Gb_Thom" begin
    Gb = Gb_Thom(0.1)
    @test keys(Gb) == (:Rb_h, :Gb_h, :kB_h, :Gb_CO2)
    @test all(isapprox.(values(Gb), values((Rb_h = 28.8, Gb_h = 0.0347, kB_h = 1.18, Gb_CO2 = 0.0264)), rtol = 1e-2))
    Gb = Gb_Thom(missing)
    @test keys(Gb) == (:Rb_h, :Gb_h, :kB_h, :Gb_CO2)
    @test all(ismissing.(values(Gb)))
    #
    dfo = DataFrame(ustar = SA[0.1,missing,0.3])
    df = copy(dfo)
    compute_Gb!(df, Val(:Thom_1972))
    @test propertynames(df) == [:ustar, :Rb_h, :Gb_h, :kB_h, :Gb_CO2]
    @test all(isapprox.(values(df[1,2:end]), values((Rb_h = 28.8, Gb_h = 0.0347, kB_h = 1.18, Gb_CO2 = 0.0264)), rtol = 1e-2))
    compute_Gb!(df, Val(:Thom_1972); Sc = SA[:Gb_N20 => 2, :Gb_CH4 => 4])
    @test propertynames(df) == [:ustar, :Rb_h, :Gb_h, :kB_h, :Gb_CO2, 
    :Gb_N20, :Gb_CH4]
    @test all(isapprox.(values(df[1,(end-1):end]), values((Gb_N20 = 0.0173, Gb_CH4 = 0.0109)), rtol = 1e-2))
end

