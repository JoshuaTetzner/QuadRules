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

N = 7
order = 2*N-1 

@time sysa = nestedsystem(order, 2, 50, chebyshevpolynomials)
@time sysb = nestedsystem(order, 50, 50, logfct)

sysa = gramschmidt(sysa)
sysb = gramschmidt(sysb)

xa, wa = gausslegendre(N) 
xa = xa# .* big(0.5) .+ 0.5
wa = wa# ./ 2

#xb = [
#    0.005652228205080097035552221875020730091975742103370603085205,
#    0.07343037174265227307828329122173937037355904523379546523871, 
#    0.28495740446255815329163031667161033867641141702377923505920, 
#    0.61948226408477838111959593386326987077397680519302133362163, 
#    0.91575808300469833371375577163284827270461436931693579356679
#]
#wb = [
#    0.021046945791854628902309427287856552365627336619228036320496, 
#    0.130705540744446697393661892987517663559602966486760316450186, 
#    0.289702301671314156867388901514800418962601391914838606731876, 
#    0.350220370120398710502902411482726506001064257050041742596754, 
#    0.208324841671985806333737366727098859111104047929131297900686
#]
##

x, w = nonsymmetricquad2(sysa, sysb, (xa), Float64.(logquadx[5] .*2 .-1), (wa), Float64.(logquadw[5].*2), 13)

println(x)
println((w))
##
N = 5
order = 2*N-1
xa, wa = gausslegendre(N) 
nodes, weights = tensorrule(xa, wa, xa, wa, 2)
nodes, weights = nonsymmetricquad(nodes, weights, order)