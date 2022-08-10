using Plots
using PyCall
using HCubature
using JLD2
using FastGaussQuadrature

x_val = 1
y_val = 1

f(x, y) = x^10*y^10
#f(x, y) = x^10*y^10 + x^20 + y^20
#f(x, y) = x^55*y^55
#f(x, y) = sin(x)^2 + sin(y)^2
#f(x, y) = x^2*y^2*log(abs(x+2))
#f(x, y) = x^2*y^2*log(abs(x+1))

#Tensorrule
tenquadval = []
nnodesten = []
for i=2:15
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

for i = 1:length(cplw)
    val = sum([f(cplx[i][j, 1], cplx[i][j, 2])*cplw[i][j] for j=1:length(cplw[i])])
    push!(nonsymquadval, val)
    push!(nnonsym, length(cplw[i]))
end

#symmetric(Tetzner)
symquadval = []
nsym = []

for i = 1:length(csplw)
    val = sum([f(csplx[i][j][1], csplx[i][j][2])*csplw[i][j] for j in eachindex(csplw[i])])
    push!(symquadval, val)
    push!(nsym, length(csplw[i]))
end

trueval = 4/121
#trueval = 1052/2541
#trueval = 0
#trueval = 4-2*sin(2)
#trueval = log(9) - 52/27
#trueval = 4/27*(log(8)-4)

##

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
##

(tenval)
(nnodesten)
(dunavantval)
(nnodesdunavant)
(nonsymquadval)
(nnonsym)
(symval)
(nsym)
#
