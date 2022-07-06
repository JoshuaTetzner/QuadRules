using SpecialPolynomials
using FastGaussQuadrature

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
N = 5
order = 2*N-1 
println("nestedapprox")
@time sys = nestedsystem(order, 2, 50, chebyshevpolynomials)
@time sgsys = nestedsystem(order, 50, 50, logfct)
println("gramschmidt")
osys = gramschmidt(sys)
osgsys = gramschmidt(sgsys)
println("quadrature")
x, w = gausslegendre(N) 
#x = big.(x)
x = x .* big(0.5) .+ 0.5
#sys.segments
@time x, w, e = nestedquadrature(osgsys, osys, x)

