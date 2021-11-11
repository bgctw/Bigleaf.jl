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
    df2 = @inferred add_Gb!(df)
    @test isequal(df2, dfo)
    df2 = @inferred add_Gb!(df, :N20 => 2, :CH4 => 4)
    @test df2 === df # mutating
    @test propertynames(df2)[end-1:end] == [:Gb_N20, :Gb_CH4]
    @test df2.Gb_N20[1] == Gb.Gb_N20
end

@testset "compute_Gb Gb_constant_kB1 and compute_Gb_quantities" begin
    kB_h = 1.18
    ustar = 0.1
    Gb = @inferred Gb_constant_kB1(ustar, kB_h)
    @test Gb ≈ 0.0347 rtol = 1e-2
    Gbq = compute_Gb_quantities(Gb, ustar)
    @test keys(Gbq) == (:Rb_h, :Gb_h, :kB_h, :Gb_CO2)
    @test all(isapprox.(values(Gbq), values(
        (Rb_h = 28.8, Gb_h = 0.0347, kB_h = 1.18, Gb_CO2 = 0.0264)), rtol = 1e-2))
    Gb = @inferred Gb_constant_kB1(missing, kB_h)
    @test ismissing(Gb)
    Gbq = compute_Gb_quantities(Gb, ustar)
    @test keys(Gbq) == (:Rb_h, :Gb_h, :kB_h, :Gb_CO2)
    @test all(ismissing.(values(getindex.(Ref(Gbq), SA[:Rb_h, :Gb_h, :Gb_CO2]))))
    #
    # DataFrame variant
    dfo = DataFrame(ustar = SA[0.1,missing,0.3])
    df = copy(dfo)
    @inferred compute_Gb!(df, Val(:constant_kB1); kB_h)
    @test propertynames(df) == [:ustar, :Gb_h]
    compute_Gb_quantities!(df)
    @test propertynames(df) == [:ustar, :Gb_h, :Rb_h, :kB_h, :Gb_CO2]
    @test all(isapprox.(values(df[1,2:end]), values(
        (Gb_h = 0.0347, Rb_h = 28.8, kB_h = 1.18, Gb_CO2 = 0.0264)), rtol = 1e-2))
end

@testset "compute_Gb Gb_Thom" begin
    ustar = 0.1
    Gb1 = @inferred Gb_Thom(ustar)
    @test Gb1 ≈ 0.0347 rtol = 1e-2
    Gbm = Gb_Thom(missing)
    @test ismissing(Gbm)
    #
    # DataFrame variant
    dfo = DataFrame(ustar = SA[ustar,missing,0.3])
    df = copy(dfo)
    @inferred compute_Gb!(df, Val(:Thom_1972))
    @test propertynames(df) == [:ustar, :Gb_h]
    @test df.Gb_h[1] == Gb1
end

@testset "compute_Gb Gb_Choudhury and Gb_Su" begin
    zh, zr, LAI = thal.zh, thal.zr, thal.LAI
    leafwidth=0.1
    wind_zh = 2.1662688
    Gb_Choud = @inferred Gb_Choudhury(; leafwidth, LAI, wind_zh)
    @test Gb_Choud ≈ 0.157 rtol=1e-2 # from R
    Gbm = @inferred Gb_Choudhury(; leafwidth, LAI, wind_zh=missing)
    @test ismissing(Gbm)
    #
    Tair, pressure, ustar = values(tha48[1, SA[:Tair, :pressure, :ustar]])
    Dl=0.01; fc = (1-exp(-LAI/2)) 
    Gb_S = @inferred Gb_Su(Tair, pressure, ustar; wind_zh, Dl, fc)
    @test Gb_S ≈ 0.185 rtol=1e-2 # from R
    Gbm = @inferred Gb_Su(missing, pressure, ustar; wind_zh, Dl, fc)
    @test ismissing(Gbm)
    #
    # DataFrame variant
    df = tha48[:,Not(:Gb_h)]
    # sum(skipmissing(...)) not inferrable in Julia 1.6
    # @descend_code_warntype roughness_parameters(
    #     Val(:wind_profile), df.ustar, df.wind, df.Tair, df.pressure, df.H; 
    #     zh, zr)
    # @code_warntype roughness_parameters(
    #     Val(:wind_profile), df.ustar, df.wind, df.Tair, df.pressure, df.H; 
    #     zh, zr)
    # z0m = (@inferred roughness_parameters(
    #     Val(:wind_profile), df.ustar, df.wind, df.Tair, df.pressure, df.H; 
    #     zh, zr)).z0m
    z0m = roughness_parameters(
        Val(:wind_profile), df.ustar, df.wind, df.Tair, df.pressure, df.H; 
        zh, zr).z0m
    # function f1(zh, ustar, z0m, Tair, pressure, H) 
    #     wind_zh = wind_profile.(zh, ustar, 0.7*zh, z0m, Tair, pressure, H)
    #     wind_zh * 2
    # end
    # @code_warntype f(zh, df.ustar, z0m,df.Tair, df.pressure, df.H)    
    wind_zh = wind_profile.(zh, df.ustar, 0.7*zh, z0m, df.Tair, df.pressure, df.H)
    @inferred compute_Gb!(df, Val(:Choudhury_1988); leafwidth, LAI, wind_zh)
    @test last(propertynames(df)) == :Gb_h
    @test df.Gb_h[1] ≈ Gb_Choud rtol=1e-6
    #
    df = tha48[:,Not(:Gb_h)]
    @inferred compute_Gb!(df, Val(:Su_2001); wind_zh, Dl, LAI)
    @test last(propertynames(df)) == :Gb_h
    @test df.Gb_h[1] ≈ Gb_S rtol=1e-6
end

