using FastGaussQuadrature
using Plots

polynomial = []
log1polynomial = []
log01polynomial = []
log0polynomial = []

#functions observed are x^10, x^10 * log(x+1), x^10 * log(x+0.1), x^10 * log(x)
for i = 1:5:50
    x, w = gausslegendre(i)
    push!(polynomial, abs(sum([w[j]*x[j]^10 for j = 1:i]) - 2/11))
    push!(log1polynomial, abs(sum([w[j]*(x[j]^5 * log(x[j]+2)) for j = 1:i]) + 
        (945*log(3) - 1052)/90))
    push!(log01polynomial, abs(sum([w[j]*(x[j]^5 * log(x[j]+1.1)) for j = 1:i]) - 
        0.3665530927935019))
    push!(log0polynomial, abs(sum([w[j]* x[j]^5 * log(x[j]+1) for j = 1:i]) - 23/45))
end

plot(1:5:50, polynomial .+ 1e-17, yaxis=:log)
plot!(1:5:50, log1polynomial, yaxis=:log)
plot!(1:5:50, log01polynomial, yaxis=:log)
plot!(1:5:50, log0polynomial, yaxis=:log)

x = [i for i = 1:5:50]
polynomial .+ 1e-18

log1polynomial
log1polynomial
log01polynomial
log0polynomial