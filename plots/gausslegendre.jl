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

using FastGaussQuadrature
using Plots

n = []
polynomial = []
sin2 = []
logf = []
sqrtx = []

#functions observed are x^10, x^10 * log(x+1), x^10 * log(x+0.1), x^10 * log(x)
for i = 1:length(logquadw)
    x = Float64.(logquadx[i]).*2 .-1
    w = Float64.(logquadw[i]).*2
    push!(n , length(w))
    push!(polynomial, abs(sum([w[j]*x[j]^10 for j = 1:length(w)]) - 2/11))
    push!(sin2, abs(sum([w[j]*(sin(x[j])^2) for j = 1:length(w)]) + (sin(2)-2)/2))
    push!(sqrtx, abs(sum([w[j]*sqrt(x[j]+1) for j = 1:length(w)]) - 
        2^(5/2)/3))
    push!(logf, abs(sum([w[j]*log(x[j]+1) for j = 1:length(w)]) - (2*log(2)-2)))
end

plot(n, polynomial .+ 1e-17, yaxis=:log)
plot!(n, sin2 .+ 1e-17, yaxis=:log)
plot!(n, sqrtx .+ 1e-17, yaxis=:log)
plot!(n, logf .+ 1e-17, yaxis=:log)
##
n
polynomial.+ eps(Float64)
sin2.+ eps(Float64)
sqrtx.+ eps(Float64)
logf.+ eps(Float64)