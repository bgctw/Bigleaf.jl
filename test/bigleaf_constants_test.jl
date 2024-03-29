using Bigleaf, Test

@testset "bigleaf_constants" begin
    cst = BigleafConstants()
    cst_moon = BigleafConstants(g = 9.81/6)
    f1(;cst = BigleafConstants()) = 2*cst.g
    @test f1() == 2*cst.g
    @test f1(;cst=cst_moon) == 2*cst_moon.g
    #@code_llvm f1()
    #@code_llvm f1(;cst=BigleafConstants(g = 9.81/6))
    #const cst_moon2 = BigleafConstants(g = 9.81/6)
    #@code_llvm f1(;cst=cst_moon2)
    #@code_llvm f1(;cst=BigleafConstants()) # not optimized
    f2(Tair;cst = BigleafConstants()) = Tair*oftype(Tair,cst.g)
    @test f2(2.0) == 2*cst.g
    #@code_llvm f2(Float32(2.0)) 
    Tair32 = Float32(2.0)
    @test f2(Tair32; cst=cst_moon) == Float32(2*cst_moon.g)
    #@code_llvm f2(Tair32;cst=cst_moon)
end

