function op(a, b, ::Type{Mul})
    return a*b
end

using Symbolics
@variables a b c Tair
Esat = a * exp((b * Tair) / (c + Tair))
Symbolics.derivative(Esat, Tair; simplify = true)