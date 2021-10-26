@testset "frac_hour" begin
    @test frac_hour(Minute, 1+1/60) == Hour(1) + Minute(1)
    @test frac_hour(1+1/60) == Hour(1) + Minute(1)
end
