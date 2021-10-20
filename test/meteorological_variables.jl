@testset "air_density" begin
  ad = air_density(25.0,100.0) # Tair, pressure
  # regression test
  @test ≈(ad, 1.168, atol =0.001)
end

@testset "pressure_from_elevation" begin
  pressure = pressure_from_elevation(500.0, 25.0) # elev, Tair
  # regression test
  @test ≈(pressure, 95.681, atol =0.001)
  #
  pressure = pressure_from_elevation(500.0, 25.0, 1.5) # elev, Tair, VPD
  # regression test
  @test ≈(pressure, 95.717, atol =0.001)
end

@testset "virtual_temp" begin
  vt = virtual_temp(25,100,1.5)  
  # regression test
  @test ≈(vt, 26.883, atol =0.001)
end

@testset "kinematic_viscosity" begin
  Tair,pressure = 25,100
  vis = kinematic_viscosity(Tair,pressure)
  # regression test
  @test ≈(vis, 1.58e-5, atol =1e-7)
end

@testset "air_density" begin
  ad = air_density(25.0,100.0) # Tair, pressure
  # regression test
  @test ≈(ad, 1.168, atol =0.001)
end




@testset "dew_point" begin
  VPD = 1.5
  Tair = 25.0
  accuracy = 1e-2
  ea = VPD_to_e(VPD,Tair)
  Td = dew_point_from_e(ea,accuracy = accuracy)                
  @test ≈(ea, Esat_from_Tair(Td), atol = accuracy)
  @test ≈(dew_point(Tair, VPD; accuracy = accuracy), Td, atol = accuracy)
  Td1 = @test_logs (:warn,r"set to 1")  dew_point(Tair, VPD; accuracy = 1.2)
  @test ≈(Td1, Td, atol = 1.0)
end

@testset "wetbulb temperature" begin
  Tair = 25.0
  VPD = 1.0
  pressure = 100
  accuracy = 1e-2
  #
  gamma  = psychrometric_constant(Tair,pressure)
  ea     = VPD_to_e(VPD,Tair)
  Tw = wetbulb_temp_from_e_Tair_gamma(ea, Tair, gamma; accuracy = accuracy)
  #@test ≈(ea, Esat_from_Tair(Tw) - gamma* (Tair - Tw), atol = accuracy)
  @test ≈(ea, Esat_from_Tair(Tw) - 0.93*gamma* (Tair - Tw), atol = accuracy)
  @test ≈(wetbulb_temp(Tair, pressure, VPD; accuracy = accuracy), Tw, atol = accuracy)
  Tw1 = @test_logs (:warn,r"set to 1") wetbulb_temp(Tair, pressure, VPD; accuracy = 1.2)
  @test ≈(Tw1, Tw, atol = 1.0)
end
