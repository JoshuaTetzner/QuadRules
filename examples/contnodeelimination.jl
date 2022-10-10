using JLD2
using FastGaussQuadrature
using HCubature
##

# asymmetric cubature for a complete polynomial system of degree n.
#n = 3
#order = 2*n-1 
#xa, wa = gausslegendre(n) 
#nodes, weights = tensorrule(xa, wa, xa, wa, 2)
#nodes, weights = asymmetriccubature(12)
nodes = c23pls1x
weights = c23pls1w
order = 23
print("Knoten: ")
println(length(weights))
nodes2, weights2 = contnonsymmetricquad(nodes, weights, order)
print("Knoten danach: ")
println(length(weights2))
##
f(x, y) = x^11*y^4
println(length(weights2))
println(sum([f(nodes2[i, 1], nodes2[i,2])* weights2[i] for i = 1:length(weights2)]))
hcubature(x->f(x[1], x[2]), [-1,-1], [1, 1])

using Plots
scatter(nodes[:, 1], nodes[:,2])

sum(nodes .* weights)
sum(nodes)
