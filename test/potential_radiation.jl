@testset "frac_hour" begin
    p = frac_hour(1+1/60)
    @test p == Hour(1) + Minute(1)
end

@testset "get_datetime_for_doy_hour summer" begin
    hours = [0,π,24]
    dts = get_datetime_for_doy_hour.(1, hours; year = 2021)
    @test dts[1] == DateTime(2021,1,1)
    @test dts[3] == DateTime(2021,1,2)
    @test dts[1] < dts[2]  < dts[3]
end

@testset "Dresden summer" begin
    hours = 8:16
    potRadSolar = potential_radiation.(160, hours, 39.94, -5.77; timezone = tz"UTC+1")
    expSolar = [
      484.152670743821, 717.876981534078, 925.130678985721,
      1091.78976612035, 1206.4967015669, 1261.43439723686, 1252.85893995917,
      1181.35473297567, 1051.79466982602]
    @test all(.≈(potRadSolar, expSolar, rtol = 0.02))
end





