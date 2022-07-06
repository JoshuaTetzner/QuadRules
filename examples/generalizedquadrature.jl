using SpecialPolynomials
using FastGaussQuadrature
using QuadGK

function logfct(n, x)
    if iseven(n)
        return basis.(SpecialPolynomials.Chebyshev, Int(n/2))(x)
    else
        return basis.(SpecialPolynomials.Chebyshev, Int(floor(n/2)))(x)*log(x) 
    end
end

function chebyshevpolynomials(n, x)
    return basis.(SpecialPolynomials.Chebyshev, n)(x)
end

##
N = 10
order = 2*N-1 
println("nestedapprox")
@time sys = nestedsystem(order, 2, 30, chebyshevpolynomials, precision=1e-64)
@time sgsys = nestedsystem(order, 50, 50, logfct, precision=1e-64)
println("gramschmidt")
osys = gramschmidt(sys)
osgsys = gramschmidt(sgsys)
println("quadrature")
x, w = gausslegendre(N) 
#x = big.(x)
x = x .* big(0.5) .+ 0.5
#sys.segments

@time x, w, e = nestedquadrature(osgsys, osys, x, tol=1e-64)

##

convert(Vector{Matrix{Float64}}, osys.systems[:])

sys.intpl
sgsys.intpl