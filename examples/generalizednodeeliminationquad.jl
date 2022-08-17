using SpecialPolynomials
using FastGaussQuadrature
using QuadRules

# generalized cubature for a complete quadratic systems of degree d
d = 5
x, w = generalizedcubature(d, type = :logquad)

##

# definition of systems of functions
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

# example quadratic generalized cubature of polynomial degree 4
p = 3
order = 2*p-1 
sysa = nestedsystem(order, 50, 50, logfct)
sysa = gramschmidt(sysa)
sysb = nestedsystem(order, 2, 30, chebyshevpolynomials)
sysb = gramschmidt(sysb)

xa, wa = generalizedquadrature(p)
xa = Float64.(xa) .* 2 .- 1
wa = Float64.(wa) .* 2
xb, wb = gausslegendre(p) 

x, w = nonsymmetricquadquad(sysa, sysb, xa, xb, wa, wb, order)
print("Degree p = ")
println(2*p-1)
println(Float64.(x))
println(Float64.(w))

##

# reconsruction of quadratic generalized cubature rules
# significance index has to be adjusted in generalizednodeelimination2
for p = 3:5
    p = 3
    order = 2*p-1 
    sysa = nestedsystem(order, 50, 50, logfct)
    sysa = gramschmidt(sysa)
    sysb = nestedsystem(order, 2, 30, chebyshevpolynomials)
    sysb = gramschmidt(sysb)
    
    xa, wa = generalizedquadrature(p)
    xa = Float64.(xa) .* 2 .- 1
    wa = Float64.(wa) .* 2
    xb, wb = gausslegendre(p) 
    
    x, w = nonsymmetricquadquad(sysa, sysb, xa, xb, wa, wb, order-1)
    print("Degree p = ")
    println(2*p-2)
    println(Float64.(x))
    println(Float64.(w))

    x, w = nonsymmetricquadquad(sysa, sysb, xa, xb, wa, wb, order)
    print("Degree p = ")
    println(2*p-1)
    println(Float64.(x))
    println(Float64.(w))
end