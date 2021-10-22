using Bigleaf
using Test
using Pipe, DataFrames, Dates

@testset "Bigleaf" begin
    @testset "unit_conversions" begin
        include("unit_conversions.jl")
    end
    @testset "meteorological_variables" begin
        include("meteorological_variables.jl")
    end
    @testset "sun_position" begin
        include("sun_position.jl")
    end
end
