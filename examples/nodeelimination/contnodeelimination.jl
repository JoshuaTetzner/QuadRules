using JLD2
using FastGaussQuadrature
using HCubature
##

# asymmetric cubature for a complete polynomial system of degree n.
n = 8
order = 2*n-1 
#xa, wa = gausslegendre(n) 
#nodes, weights = tensorrule(xa, wa, xa, wa, 2)
nodes, weights = asymmetriccubature(order)
#nodes = c6plSFx
#weights = c6plSFw
#order = 6
print("Knoten: ")
println(length(weights))
@time nodes2, weights2 = contnonsymmetricquad2(nodes, weights, order-1)
print("Knoten danach: ")
println(length(weights2))
##

for n = 3:11
    nodes = cplxb[n-2]
    weights = cplwb[n-2]
    order = n
    println("\n")
    print("Order: ")
    println(order)
    print("Knoten: ")
    println(length(weights))
    nodes2, weights2 = contnonsymmetricquad2(nodes, weights, order)
    print("Knoten danach: ")
    println(length(weights2))
end

##

