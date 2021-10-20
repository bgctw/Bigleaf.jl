@testset "dew_point" begin
  VPD = 1.5
  Tair = 25.0
  accuracy = 1e-2
  ea = VPD_to_e(VPD,Tair)
  Td = dew_point_from_e(ea,accuracy = accuracy)                
  @test ≈(ea, Esat_from_Tair(Td), atol = accuracy)
  @test ≈(dew_point(Tair, VPD; accuracy = accuracy), Td, atol = accuracy)
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
end
