#Firzandi Ahsan Dwi Styana
#21537144016
function quad_midpoint(f,a,b,N) 
    h = (b-a)/N
    int = 0.0
    for k=1:N
        xk_mid = (b-a) * (2k-1)/(2N) + a
        int = int + h*f(xk_mid)
    end
    return int
end

function quad_trap(f,a,b,N) 
    h = (b-a)/N
    int = h * ( f(a) + f(b) ) / 2
    for k=1:N-1
        xk = (b-a) * k/N + a
        int = int + h*f(xk)
    end
    return int
end

using Printf
f(x) = exp(cos(x)^3)
ref = quad_trap(f,0,2π,10000);

for n=2:5
    error = abs(quad_trap(f,0,2π,2^n) - ref)
    @printf("n=%3i, %4.2e\n", 2^n,error)
end

x_true = LinRange(0,2π,1600)
x_trap = LinRange(0,2π,16)
plot(x_trap,f.(x_trap),fillrange=0,label="Trapezoidal rule")
plot!(x_true,f.(x_true),line=4,legend=:bottomright,label="Function values")

using QuadGK
f(x) = sqrt(x)
a=0
b=1
I,est = quadgk(f, a, b, rtol=1e-8)


using FastGaussQuadrature
@time x, w = gausslegendre( 10_000_000 );

x, w = gausslegendre( 7 );
exact = exp(1)-exp(-1);
numer = w'exp.(x);
abs(numer-exact)

f(x,y) = exp(-(x^2+y^2))

function Quad2D(f,n)
    x, w = gausslegendre( n );
    y = x;
    sum = 0
    for i=1:n, j=1:n
        sum = sum + f(x[i],y[j]) * w[i] * w[j]
    end
    return sum
end

function simps(f::Function, a::Number, b::Number, n::Number)
    n % 2 == 0 || error("`n` must be even")
    h = (b-a)/n
    s = f(a) + f(b)
    s += 4sum(f(a + collect(1:2:n) * h))
    s += 2sum(f(a + collect(2:2:n-1) * h))
    return h/3 * s
end

function simps(y::Vector, h::Number)
    n = length(y)-1
    n % 2 == 0 || error("`y` length (number of intervals) must be odd")
    s = sum(slice(y,1:2:n) + 4slice(y,2:2:n) + slice(y,3:2:n+1))
    return h/3 * s
end

function simps(y::Vector, x::Union{Vector,Range})
    n = length(y)-1
    n % 2 == 0 || error("`y` length (number of intervals) must be odd")
    length(x)-1 == n || error("`x` and `y` length must be equal")
    h = (x[end]-x[1])/n
    s = sum(slice(y,1:2:n) + 4slice(y,2:2:n) + slice(y,3:2:n+1))
    return h/3 * s
end

function simps(y::Matrix, h::Number, axis::Integer=2)
    n = size(y)[axis]-1
    n % 2 == 0 || error("`y` length (number of intervals) must be odd")
    axis in [1,2] || error("`axis` must be 1 or 2")
    if axis == 2
        inds = [(:,1:2:n),(:,2:2:n),(:,3:2:n+1)]
    else
        inds = [(1:2:n,:),(2:2:n,:),(3:2:n+1,:)]
    end
    s = sum(slice(y,inds[1]...) + 4slice(y,inds[2]...) + slice(y,inds[3]...),axis)
    return h/3 * s
end

function simps(y::Matrix, x::Union{Matrix,Vector,Range}, axis::Integer=2)
    n = size(y)[axis]-1
    n % 2 == 0 || error("`y` length (number of intervals) must be odd")
    axis in [1,2] || error("`axis` must be 1 or 2")
    if length(size(x)) == 1
        length(x)-1 == n || error("`x` and `y` length must be equal in axis `axis`")
        h = (x[end]-x[1])/n
    else
        size(x)[axis]-1 == n || error("`x` and `y` length must be equal in axis `axis`")
        h = nothing
    end
    if axis == 2
        h != nothing || (h = (x[1,end]-x[1,1])/n)
        inds = [(:,1:2:n),(:,2:2:n),(:,3:2:n+1)]
    else
        h != nothing || (h = (x[end,1]-x[1,1])/n)
        inds = [(1:2:n,:),(2:2:n,:),(3:2:n+1,:)]
    end
    s = sum(slice(y,inds[1]...) + 4slice(y,inds[2]...) + slice(y,inds[3]...),axis)
    return h/3 * s
end


