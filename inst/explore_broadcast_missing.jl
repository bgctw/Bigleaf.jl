x = [1.0, 2.0, missing]
y = [3.0, 3.0, 3.0]
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