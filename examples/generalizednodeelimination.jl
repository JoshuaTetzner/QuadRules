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
    return basis.(SpecialPolynomials.Chebyshev, n)(x)
end

N1 = 6
order1 = 2*N1-1 
@time sysa = nestedsystem(order1, 50, 50, logfct)
N2 = 4
order2 = 2*N2-1 
@time sysb = nestedsystem(order2, 2, 22, chebyshevpolynomials)

sysa = gramschmidt(sysa)
sysb = gramschmidt(sysb)

xa = Float64.(logquadx[N1-2].*2 .- 1)
wa = Float64.(logquadw[N1-2].*2)
#xa, wa = gausslegendre(N1)
xb, wb = gausslegendre(N2) 

x, w = nonsymmetricquad2(sysa, sysb, xa, xb, wa, wb, order1-1, order2-1)

println(x)
println(length(w))
##
function chebyshevpolynomials(n, x)
    return basis.(SpecialPolynomials.Chebyshev, n)(x)
end

n = 6
order = 2*n-1 
@time sys = nestedsystem(order, 2, 30, chebyshevpolynomials)

sys = gramschmidt(sys)
xa, wa = gausslegendre(n)

x, w = nonsymmetricquad2(sys, sys, xa, xa, wa, wa, order, order)
printnln(length(weights))