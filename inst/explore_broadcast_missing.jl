x = [1.0, 2.0, missing, missing]
y = [3.0, 3.0, 3.0, 3.0]
f(x,y) = x+y
function g(x,y)
    z = f.(x,y)
    if z isa Union
        return typeof(z)
        convert(typeof(z).b.a, z)
    else
        z
    end
end
z = g(x,y)
@code_warntype g(x,y)


typeof(x), typeof(x[1:2]), typeof(x[3:4]) # preserves eltype
typeof(x.*1.0), typeof(x[1:2].*1.0), typeof(x[3:4].*1.0) # does not
typeof(x*1.0), typeof(x[1:2]*1.0), typeof(x[3:4]*1.0) # does not


x = [1.0, 2.0, missing, missing]
y = 1.0
typeof(x.*y), typeof(x[1:2].*y), typeof(x[3:4].*y) 
# (Vector{Union{Missing, Float64}}, Vector{Float64}, Vector{Missing})
res = x[1:2].*y
res_eltype = Union{eltype(x),typeof(y)}
convert(Vector{res_eltype}, res)::Vector{res_eltype}

f1(x,y) = x[1:2].*y
@code_warntype f1(x,y)

function f2(x,y)
    res_eltype = Union{eltype(x),typeof(y)}
    convert(Vector{res_eltype}, x[1:2].*y)::Vector{res_eltype}
end
@inferred f2(x,y)
