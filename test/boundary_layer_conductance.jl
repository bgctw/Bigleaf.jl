@testset "add_Gb" begin
    Gb_h = 0.0347
    @test add_Gb(Gb_h) == NamedTuple()
    #
    Gb = add_Gb(Gb_h, :N20 => 2, :CH4 => 4)
    @test keys(Gb) == (:Gb_N20, :Gb_CH4)
    @test ≈(Gb.Gb_N20, 0.0173, rtol = 1/100)
    @test ≈(Gb.Gb_CH4, 0.0109, rtol = 1/100)
    #
    Gbm = add_Gb(missing, :N20 => 2, :CH4 => 4)
    @test keys(Gbm) == (:Gb_N20, :Gb_CH4)
    @test all(ismissing.(values(Gbm)))
    # 
    # DataFrame variant
    dfo = DataFrame(Gb_h=SA[Gb_h, missing, 0.055])
    df = copy(dfo)
    df2 = add_Gb!(df)
    @test isequal(df2, dfo)
    df2 = add_Gb!(df, :N20 => 2, :CH4 => 4)
    @test df2 === df # mutating
    @test propertynames(df2)[end-1:end] == [:Gb_N20, :Gb_CH4]
    @test df2.Gb_N20[1] == Gb.Gb_N20
end

@testset "compute_Gb Gb_constant_kB1" begin
    kB_h = 1.18
    Gb = Gb_constant_kB1(0.1, kB_h)
    @test keys(Gb) == (:Rb_h, :Gb_h, :kB_h, :Gb_CO2)
    @test all(isapprox.(values(Gb), values((Rb_h = 28.8, Gb_h = 0.0347, kB_h = 1.18, Gb_CO2 = 0.0264)), rtol = 1e-2))
    Gb = Gb_constant_kB1(missing, kB_h)
    @test keys(Gb) == (:Rb_h, :Gb_h, :kB_h, :Gb_CO2)
    @test Gb.kB_h == kB_h
    @test all(ismissing.(values(getindex.(Ref(Gb), SA[:Rb_h, :Gb_h, :Gb_CO2]))))
    #
    # DataFrame variant
    dfo = DataFrame(ustar = SA[0.1,missing,0.3])
    df = copy(dfo)
    compute_Gb!(df, Val(:constant_kB1); kB_h)
    @test propertynames(df) == [:ustar, :Rb_h, :Gb_h, :kB_h, :Gb_CO2]
    @test all(isapprox.(values(df[1,2:end]), values((Rb_h = 28.8, Gb_h = 0.0347, kB_h = 1.18, Gb_CO2 = 0.0264)), rtol = 1e-2))
end

@testset "compute_Gb Gb_Thom" begin
    Gb = Gb_Thom(0.1)
    @test keys(Gb) == (:Rb_h, :Gb_h, :kB_h, :Gb_CO2)
    @test all(isapprox.(values(Gb), values((Rb_h = 28.8, Gb_h = 0.0347, kB_h = 1.18, Gb_CO2 = 0.0264)), rtol = 1e-2))
    Gb = Gb_Thom(missing)
    @test keys(Gb) == (:Rb_h, :Gb_h, :kB_h, :Gb_CO2)
    @test all(ismissing.(values(Gb)))
    #
    # DataFrame variant
    dfo = DataFrame(ustar = SA[0.1,missing,0.3])
    df = copy(dfo)
    compute_Gb!(df, Val(:Thom_1972))
    @test propertynames(df) == [:ustar, :Rb_h, :Gb_h, :kB_h, :Gb_CO2]
    @test all(isapprox.(values(df[1,2:end]), values((Rb_h = 28.8, Gb_h = 0.0347, kB_h = 1.18, Gb_CO2 = 0.0264)), rtol = 1e-2))
end

@testset "compute_Gb Gb_Choudhury" begin
    zh, zr, LAI = thal.zh, thal.zr, thal.LAI
    leafwidth=0.1
    wind_zh = 1.2
    Gb = Gb_Choudhury(tha48.ustar[24]; leafwidth, LAI, wind_zh)
    @test keys(Gb) == (:Rb_h, :Gb_h, :kB_h, :Gb_CO2)
    @test Gb.Rb_h ≈ 8.533 rtol=1e-3 # regression to first implementation
    Gbm = Gb_Choudhury(missing; leafwidth, LAI, wind_zh)
    @test keys(Gb) == (:Rb_h, :Gb_h, :kB_h, :Gb_CO2)
    @test Gbm.Rb_h == Gb.Rb_h
    #
    # DataFrame variant
    df = copy(tha48)
    compute_Gb!(df, Val(:Choudhury_1988); leafwidth, LAI, zh, zr)
    @test propertynames(df)[(end-3):end] == SA[:Rb_h, :Gb_h, :kB_h, :Gb_CO2]
    all(isapprox.(values(df[24,SA[:Rb_h, :Gb_h, :kB_h, :Gb_CO2]]), 
        (7.534, 0.1327, 2.255, 0.1008), rtol=1e-3))
end

@testset "compute_Gb Gb_Su" begin
    zh, zr, LAI = thal.zh, thal.zr, thal.LAI
    Dl=0.01; leafwidth=0.1; fc = (1-exp(-LAI/2)) 
    wind_zh = 1.2
    Tair, pressure, ustar = values(tha48[24, SA[:Tair, :pressure, :ustar]])
    Gb = Gb_Su(Tair, pressure, ustar; wind_zh, Dl, fc)
    @test keys(Gb) == (:Rb_h, :Gb_h, :kB_h, :Gb_CO2)
    @test Gb.Rb_h ≈ 1.221 rtol=1e-3 # regression to first implementation
    Gb = Gb_Su(missing, pressure, ustar; wind_zh, Dl, fc)
    @test keys(Gb) == (:Rb_h, :Gb_h, :kB_h, :Gb_CO2)
    @test all(ismissing.(values(Gb)))
    #
    # DataFrame variant
    df = copy(tha48)
    compute_Gb!(df, Val(:Su_2001); Dl, LAI, zh, zr)
    @test propertynames(df)[(end-3):end] == SA[:Rb_h, :Gb_h, :kB_h, :Gb_CO2]
    all(isapprox.(values(df[24,SA[:Rb_h, :Gb_h, :kB_h, :Gb_CO2]]), 
        (1.767, 0.5659, 0.5289, 0.4299), rtol=1e-3)) # regression to first implementation
end

