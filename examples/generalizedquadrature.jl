using SpecialPolynomials
using FastGaussQuadrature

function logfct(n, x)
    if iseven(n)
        return basis.(SpecialPolynomials.Chebyshev, Int(n/2))(x)
    else
        return basis.(SpecialPolynomials.Chebyshev, Int(floor(n/2)))(x)*log(x+1) 
    end
end

function chebyshevpolynomials(n, x)
    return basis.(SpecialPolynomials.Chebyshev, n)(x)#x^n#
end

#
N = 15
order = 2*N-1 
println("nestedapprox")
@time sys = nestedsystem(order, 2, 50, chebyshevpolynomials)
@time sgsys = nestedsystem(order, 50, 50, logfct)
println("gramschmidt")
osys = gramschmidt(sys)
osgsys = gramschmidt(sgsys)
println("quadrature")
x, w = gausslegendre(N) 
#x = x .* big(0.5) .+ 0.5
#x = big.(x)
@time x, w, e = nestedquadrature(osgsys, sys, x, tol=1e-12)
#osys.intpl

