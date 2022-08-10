using FastGaussQuadrature
using Plots

polynomial = []
sin2 = []
logf = []
sqrtx = []
npoints = []

#functions observed are x^10, x^10 * log(x+1), x^10 * log(x+0.1), x^10 * log(x)
for i = 3:9
    x, w = gausslegendre(i)
    push!(polynomial, sum([w[j]*x[j]^10 for j in eachindex(w)]))
    push!(sin2, sum([w[j]*(sin(x[j])^2) for j in eachindex(w)]))
    push!(sqrtx, sum([w[j]*sqrt(x[j]+1) for j in eachindex(w)]))
    push!(logf, sum([w[j]*log(x[j]+1) for j in eachindex(w)]))
    push!(npoints, length(w))
end

for i = 10:5:30
    x, w = gausslegendre(i)
    push!(polynomial, sum([w[j]*x[j]^10 for j in eachindex(w)]))
    push!(sin2, sum([w[j]*(sin(x[j])^2) for j in eachindex(w)]))
    push!(sqrtx, sum([w[j]*sqrt(x[j]+1) for j in eachindex(w)]))
    push!(logf, sum([w[j]*log(x[j]+1) for j in eachindex(w)]))
    push!(npoints, length(w))
end

errpol = 10 .* log10.(abs.((polynomial .- 2/11)/(2/11) .+ eps(Float64)))
errsin2 = 10 .* log10.(abs.((sin2 .+ (sin(2)-2)/2)/((sin(2)-2)/2) .+ eps(Float64)))
errsqrtx = 10 .* log10.(abs.((sqrtx .- 2^(5/2)/3)/(2^(5/2)/3) .+ eps(Float64)))
errlogf = 10 .* log10.(abs.((logf .- (2*log(2)-2))/(2*log(2)-2) .+ eps(Float64)))

plot(npoints, errpol)
plot!(npoints, errsin2)
plot!(npoints, errsqrtx)
plot!(npoints, errlogf)
