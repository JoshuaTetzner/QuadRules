using Plots
using PyCall
using JLD2
using FastGaussQuadrature

# Fig. 3.3
f(x, y) = x^6*y^6
# Fig. 3.4
#f(x, y) = x^6*y^6 + x^12 + y^12
# Fig. 3.5
#f(x, y) = x^55*y^55
# Fig. 3.6
#f(x, y) = sin(x)^2 + sin(y)^2
# Fig. 3.7
#f(x, y) = x^2*y^2*log(abs(x+2))
# Fig. 3.8
#f(x, y) = x^2*y^2*log(abs(x+1))

#Tensorrule
tenquadval = []
nnodesten = []
for i=2:10
    xa, wa = gausslegendre(i) 
    nodes, weights = tensorrule(xa, wa, xa, wa, 2)
    val = sum([f(nodes[i,1], nodes[i,2])*weights[i] for i in eachindex(weights)])
    push!(tenquadval, val)
    push!(nnodesten, length(weights))
end

#Dunavant
dunavantquadval = []
nnodesdunavant = []

for i = 4:10
    quadpy = pyimport("quadpy")
    num = ""
    if i != 10
        num = string(0) * string(i)
    else
        num = string(i)
    end
    a = quadpy.c2.schemes["dunavant_" * num]
    nodes = a().points
    weights = 4 .* a().weights
    val = sum([f((nodes[1, j]), (nodes[2, j]))*weights[j] for j in eachindex(weights)])
    push!(dunavantquadval, val)
    push!(nnodesdunavant, length(weights))
end

#nonsymmetric(Tetzner)
nonsymquadval = []
nnonsym = []

for d = 3:15
    x, w = asymmetriccubature(d)
    val = sum([f(x[j, 1], x[j, 2])*w[j] for j in eachindex(w)])
    push!(nonsymquadval, val)
    push!(nnonsym, length(w))
end

#symmetric(Tetzner)
symquadval = []
nsym = []

for d = 5:2:21
    x, w = symmetriccubature(d)
    val = sum([f(x[j, 1], x[j, 2])*w[j] for j in eachindex(w)])
    push!(symquadval, val)
    push!(nsym, length(w))
end

# Fig. 3.3
trueval = 4/49
# Fig. 3.4
#trueval = 444/637
# Fig. 3.5
#trueval = 0
# Fig. 3.6
#trueval = 4-2*sin(2)
# Fig. 3.7
#trueval = log(9) - 52/27
# Fig. 3.8
#trueval = 4/27*(log(8)-4)

#Plots 
if trueval != 0
    tenval = abs.((tenquadval .- trueval)) / trueval .+ eps(Float64)
    dunavantval = abs.((dunavantquadval .- trueval)) / trueval .+ eps(Float64)
    nonsymval = abs.((nonsymquadval .- trueval)) / trueval .+ eps(Float64)
    symval = abs.((symquadval .- trueval)) / trueval .+ eps(Float64)
else
    tenval = abs.((tenquadval .- trueval)) .+ eps(Float64)
    dunavantval = abs.((dunavantquadval .- trueval)) .+ eps(Float64)
    nonsymval = abs.((nonsymquadval .- trueval)) .+ eps(Float64)
    symval = abs.((symquadval .- trueval)) .+ eps(Float64)
end

tenval = 10 .* log10.(abs.(tenval))
dunavantval = 10 .* log10.(abs.(dunavantval))
nonsymquadval = 10 .* log10.(abs.(nonsymval))
symval = 10 .* log10.(abs.(symval))

p2 = plot(nnodesten, tenval, label = "tensorproduct",  marker = :a)
plot!(nnodesdunavant, dunavantval, label = "quadpy(Dunavant)",  marker = :a)
plot!(nnonsym, nonsymquadval, label = "nonsymmertric(Tetzner)",  marker = :a)
plot!(nsym, symval, label = "symmertric(Tetzner)",  marker = :a)
