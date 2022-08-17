using SpecialPolynomials
using FastGaussQuadrature
using QuadRules

# generalized cubature for a complete rectangle systems of degree p
p = 5
x, w = generalizedcubature(p, type = :logrect)

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

# example rectangle generalized cubature of polynomial degree 4
p1 = 6
order1 = 2*p1-1 
sysa = nestedsystem(order1, 50, 50, logfct)
sysa = gramschmidt(sysa)
p2 = 3
order2 = 2*p2-1 
sysb = nestedsystem(order2, 2, 30, chebyshevpolynomials)
sysb = gramschmidt(sysb)

xa, wa = generalizedquadrature(p1)
xb, wb = gausslegendre(p2) 

x, w = nonsymmetricquadrect(sysa, sysb, Float64.(xa), xb, Float64.(wa), wb, order1, order2)
print("Degree p = ")
print(2*p2-1)
println(x)
println(w)

##

# reconsruction of rectangle generalized cubature rules
# significance index has to be adjusted in generalizednodeelimination2
for p = 3:5
    p1 = 2*p
    order1 = 2*p1-1 
    @time sysa = nestedsystem(order1, 50, 50, logfct)
    sysa = gramschmidt(sysa)
    p2 = p
    order2 = 2*p2-1 
    @time sysb = nestedsystem(order2, 2, 30, chebyshevpolynomials)
    sysb = gramschmidt(sysb)

    xa, wa = generalizedquadrature(p1)
    xb, wb = gausslegendre(p2) 

    x, w = nonsymmetricquadquad(
        sysa,
        sysb,
        Float64.(xa),
        xb,
        Float64.(wa),
        wb,
        order1-2,
        order2-1
    )
    print("Degree p = ")
    print(2*p-2)
    println(x)
    println(w)

    x, w = nonsymmetricquadquad(sysa, sysb, Float64.(xa), xb, Float64.(wa), wb, order1, order2)
    print("Degree p = ")
    print(2*p-1)
    println(x)
    println(w)
end