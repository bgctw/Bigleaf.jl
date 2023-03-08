module mymod
    const g = 9.81
end

module mymod2
    const g = 9.81/6
end


struct Constants1 end
(c::Constants1)(::Val{s}) where s = getproperty(mymod, s)
cst1 = Constants1()
cst1(Val(:g))
@code_llvm (() -> 2*Constants1()(Val(:g)))()

f1(::Val{cst}) where cst = 2*cst(Val(:g))
f1(Val(cst1))
@code_llvm f1(Val(cst1))

struct Constants2 end
(c::Constants2)(::Val{s}) where s = getproperty(mymod2, s)
cst3 = Constants2()
f1(Val(cst3))
@code_llvm f1(Val(cst3))

struct ConstantsB end
(c::ConstantsB)(::Val{s}) where s = getproperty(BigLeaf.BigLeafConstantsDef, s)
cstB = ConstantsB()
cstB(Val(:g))
f1(Val(cstB))
@code_llvm f1(Val(cstB))

using BigLeaf
f3(;cst=BigLeafConstants()) = 2 * cst.g
f3()
@code_llvm f3()
@code_llvm f3(cst=BigLeafConstants())
@code_warntype f3()
@code_warntype f3(cst=BigLeafConstants())
@code_warntype bigleaf_constants.g

