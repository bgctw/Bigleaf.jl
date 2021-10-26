using Bigleaf
using Test
using Pipe, DataFrames, Dates, TimeZones

@testset "Bigleaf" begin
    @testset "util" begin
        include("util.jl")
    end
    @testset "unit_conversions" begin
        include("unit_conversions.jl")
    end
    @testset "meteorological_variables" begin
        include("meteorological_variables.jl")
    end
    @testset "sun_position" begin
        include("sun_position.jl")
    end
    @testset "potential_radiation" begin
        include("potential_radiation.jl")
    end
    @testset "evapotranspiration" begin
        include("evapotranspiration.jl")
    end
end
