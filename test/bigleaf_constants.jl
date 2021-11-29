#@testset "bigleaf_constants" begin
    cst = BigleafConstants()
    twog(cst::SymbolValueCall) = 2*cst(Val(:g))
    @test twog(cst) == 2*Bigleaf.BigleafConstantsDef.g
    f1() = cst(Val(:g))
    f1()
    #@code_llvm f1()
end

