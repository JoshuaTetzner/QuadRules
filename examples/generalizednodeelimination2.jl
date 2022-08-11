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

##

# example rectangel generalized cubature degree 4
N1 = 10
order1 = 2*N1-1 
sysa = nestedsystem(order1, 50, 50, logfct)
sysa = gramschmidt(sysa)
N2 = 5
order2 = 2*N2-1 
sysb = nestedsystem(order2, 2, 30, chebyshevpolynomials)
sysb = gramschmidt(sysb)

xa = Float64.(logquadx[N1-2].*2 .- 1)
wa = Float64.(logquadw[N1-2].*2)
xb, wb = gausslegendre(N2) 

x, w = nonsymmetricquad3(sysa, sysb, xa, xb, wa, wb, order1-2, order2-1)
print("Degree p = 4")
println(Float64.(x))
println(Float64.(w))

##

# reconsruction of rectangel generalized cubature rules
# significance index has to be adjusted in generalizednodeelimination2
for N = 2:5
    N1 = 2*N
    order1 = 2*N1-1 
    @time sysa = nestedsystem(order1, 50, 50, logfct)
    sysa = gramschmidt(sysa)
    N2 = N
    order2 = 2*N2-1 
    @time sysb = nestedsystem(order2, 2, 30, chebyshevpolynomials)
    sysb = gramschmidt(sysb)

    xa = Float64.(logquadx[N1-2].*2 .- 1)
    wa = Float64.(logquadw[N1-2].*2)
    xb, wb = gausslegendre(N2) 

    x, w = nonsymmetricquad3(sysa, sysb, xa, xb, wa, wb, order1-2, order2-1)
    print("Degree p = ")
    print(2*N-2)
    println(x)
    println(w)

    x, w = nonsymmetricquad3(sysa, sysb, xa, xb, wa, wb, order1, order2)
    print("Degree p = ")
    print(2*N-1)
    println(x)
    println(w)
end