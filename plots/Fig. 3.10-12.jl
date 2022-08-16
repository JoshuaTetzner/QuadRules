using Plots
using FastGaussQuadrature

# Fig. 3.10
#f(x, y) = x^4*y^4*log(x+1) + x^2*y^2
# Fig. 3.11
#f(x, y) = x^2*y^2*log(x+1) + x^4 + y^4
# Fig. 3.12
f(x, y) = x^2*y^4*log(x+1) + x^4*log(x+1) + y^8

#quadratic nodeelimination
asgg = []
nasgg = []
for d = 4:11
    x, w = generalizedcubature(d, type = :logquad)
    val = sum([f(x[i,1], x[i,2])*w[i] for i in eachindex(w)])
    push!(asgg, val)
    push!(nasgg, length(w))
end

#product rule
pgg = []
npgg = []

for i = 1:8
    xa = Float64.(gqlogx[i]).*2 .-1
    wa =  Float64.(gqlogw[i]).*2
    if iseven(length(xa))
        xb, wb = gausslegendre(Int(length(xa)/2))
        nodes, weights = tensorrule(xa, wa, xb, wb, 2)
        val = sum([f(nodes[i,1], nodes[i,2])*weights[i] for i in eachindex(weights)])
        push!(pgg, val)
        push!(npgg, length(weights))
    end
end

#rectangel nodeelimination
asrgg = []
nasrgg = []
for d = 4:9
    x, w = generalizedcubature(d, type = :logrect)
    val = sum([f(x[i,1], x[i,2])*w[i] for i in eachindex(w)])
    push!(asrgg, val)
    push!(nasrgg, length(w))
end

# Fig. 3.10
#trueval = 4*(56 + 45* log(2))/1125
# Fig. 3.11
#trueval = 4/135 * (34 + 5*log(8))
# Fig. 3.12
trueval = 16/225 * (15*log(2)-16)

#Plots 
pggval = abs.((pgg .- trueval)) ./ trueval .+ eps(Float64)
asggval = abs.((asgg .- trueval)) ./ trueval .+ eps(Float64)
asrggval = abs.((asrgg .- trueval)) ./ trueval .+ eps(Float64)

pggval = 10 .* log10.(abs.(pggval))
asggval = 10 .* log10.(abs.(asggval))
asrggval = 10 .* log10.(abs.(asrggval))

p2 = plot(npgg, pggval, label = "tensorproduct",  marker = :a)
plot!(nasgg, asggval, label = "eliminated",  marker = :a)
plot!(nasrgg, asrggval, label = "rectangle eliminated",  marker = :a)
##

pggval
npgg

asggval
nasgg

asrggval
nasrgg