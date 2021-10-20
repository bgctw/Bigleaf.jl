using bigleaf
using Test

@testset "bigleaf" begin
    @testset "unit_conversions" begin
        include("unit_conversions.jl")
    end
    @testset "meteorological_variables" begin
        include("meteorological_variables.jl")
    end
end
