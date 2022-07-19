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

N1 = 8
order1 = 2*N1-1 
@time sysa = nestedsystem(order1, 2, 22, chebyshevpolynomials)
N2 = 4
order2 = 2*N2-1 
@time sysb = nestedsystem(order2, 2, 22, chebyshevpolynomials)

sysa = gramschmidt(sysa)
sysb = gramschmidt(sysb)

xa = Float64.(logquadx[N1-2].*2 .- 1)
wa = Float64.(logquadw[N1-2].*2)
# = gausslegendre(N1)
xb, wb = gausslegendre(N2) 

x, w = nonsymmetricquad2(sysa, sysb, xa, xb, wa, wb, 14, 6)

println(x)
println((w))
##

N = 5
order = 2*N-1
xa, wa = gausslegendre(N1) 
nodes, weights = tensorrule(xa, wa, xa, wa, 2)
nodes, weights = nonsymmetricquad(nodes, weights, order)
