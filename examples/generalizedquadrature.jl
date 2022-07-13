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
N = 3
order = 2*N-1 
println("nestedapprox")
@time sys = nestedsystem(order, 2, 50, chebyshevpolynomials, T=BigFloat)
@time sgsys = nestedsystem(order, 70, 70, logfct, T = BigFloat)
println("gramschmidt")
osys = gramschmidt(sys)
osgsys = gramschmidt(sgsys)
println("quadrature")
x, w = gausslegendre(N) 
x = big.(x)
@time x, w, e = nestedquadrature(osgsys, sys, x, tol=1e-64)
#osys.intpl
##

##
function logfct(n, x)
    if iseven(n)
        return basis.(SpecialPolynomials.Chebyshev, Int(n/2))(x)
    else
        return basis.(SpecialPolynomials.Chebyshev, Int(floor(n/2)))(x)*log(x+2) 
    end
end

function chebyshevpolynomials(n, x)
    return basis.(SpecialPolynomials.Chebyshev, n)(x)#x^n#
end

#
N = 8
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
c1 = 2
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
        #step *= 2
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

##
0.5 .*x .+ 0.5
w./2


##
for i=3:10
    x, w = correctlog(logquadx[i-2])
    println(x)
    println(w)
end