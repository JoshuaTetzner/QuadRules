using SpecialPolynomials
using FastGaussQuadrature

# generalized quadrature for logarithmic system with n = 5 points
n = 5
x, w = generalizedquadrature(n)

##

# generalized quadrature for log(x+1) with N = 5 points
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

N = 5
order = 2*N-1 
println("nestedapprox")
sys = nestedsystem(order, 2, 50, chebyshevpolynomials, T=BigFloat)
sgsys = nestedsystem(order, 50, 50, logfct, T = BigFloat)
println("gramschmidt")
sys = gramschmidt(sys)
sgsys = gramschmidt(sgsys)
println("quadrature")
x, w = gausslegendre(N) 
x, w, e = nestedquadrature(sgsys, sys, big.(x), tol=1e-64)

##

# generalized quadrature with second continuation with N = 10 points

function logfct(n, x)
    if iseven(n)
        return basis.(SpecialPolynomials.Chebyshev, Int(n/2))(x)
    else
        return basis.(SpecialPolynomials.Chebyshev, Int(floor(n/2)))(x)*log(x+4) 
    end
end

function chebyshevpolynomials(n, x)
    return basis.(SpecialPolynomials.Chebyshev, n)(x)
end

N = 10
order = 2*N-1 
println("nestedapprox")
@time sys = nestedsystem(order, 2, 50, chebyshevpolynomials)
@time sgsys = nestedsystem(order, 50, 50, logfct)
println("gramschmidt")
osys = gramschmidt(sys)
osgsys = gramschmidt(sgsys)
println("quadrature")
x, w = gausslegendre(N) 
x = big.(x)
@time x, w, e = nestedquadrature(osgsys, sys, x)
step = 1
c1 = 4
c2 = c1-step
while c2 >= 1
    print("C1: ")
    println(c1)
    print("C2: ")
    println(c2)
    function logfct1(n, x)
        if iseven(n)
            return basis.(SpecialPolynomials.Chebyshev, Int(n/2))(x)
        else
            return basis.(SpecialPolynomials.Chebyshev, Int(floor(n/2)))(x)*log(x+c1) 
        end
    end
    function logfct2(n, x)
        if iseven(n)
            return basis.(SpecialPolynomials.Chebyshev, Int(n/2))(x)
        else
            return basis.(SpecialPolynomials.Chebyshev, Int(floor(n/2)))(x)*log(x+c2) 
        end
    end

    println("nestedapprox")
    @time sys = nestedsystem(order, 50, 50, logfct1)
    @time sgsys = nestedsystem(order, 50, 50, logfct2)
    println("gramschmidt")
    osys = gramschmidt(sys)
    osgsys = gramschmidt(sgsys)
    println("quadrature")
    @time xnew, wnew, e = nestedquadrature(osgsys, osys, x)

    if e < 1e-16
        x = xnew
        w = wnew
        c1 = c2
        if step < 0.5
            step .* 2
        end
        if (c2 + step) < 1
            c2 = 1
        else
            c2-=step
        end
    else
        c2 += step
        step = step/2
        c2-=step
    end
end
