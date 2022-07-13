using FastGaussQuadrature
using Plots

polynomial = []
sin2 = []
logf = []
sqrtx = []

#functions observed are x^10, x^10 * log(x+1), x^10 * log(x+0.1), x^10 * log(x)
for i = 1:5:50
    x, w = gausslegendre(i)
    push!(polynomial, abs(sum([w[j]*x[j]^10 for j = 1:i]) - 2/11))
    push!(sin2, abs(sum([w[j]*(sin(x[j])^2) for j = 1:i]) + (sin(2)-2)/2))
    push!(sqrtx, abs(sum([w[j]*sqrt(x[j]+1) for j = 1:i]) - 
        2^(5/2)/3))
    push!(logf, abs(sum([w[j]*log(x[j]+1) for j = 1:i]) - (2*log(2)-2)))
end

plot(1:5:50, polynomial .+ 1e-17, yaxis=:log)
plot!(1:5:50, sin2, yaxis=:log)
plot!(1:5:50, sqrtx, yaxis=:log)
plot!(1:5:50, logf, yaxis=:log)
##
x = [i for i = 1:5:50]
polynomial .+ 1e-18

polynomial.+ 1e-18
sin2.+ 1e-18
sqrtx.+ 1e-18
logf.+ 1e-18

##

function tester(a, tol)
    if tol < 1e-16
        T = BigFloat
    else
        T = Float64
    end

    return T(a)
end

eltype(tester(1, 1e-25))