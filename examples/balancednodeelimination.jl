using JLD2
using FastGaussQuadrature

# asymmetric cubature for a complete polynomial system of degree n.
n = 9
order = 2*n-1 
xa, wa = gausslegendre(n) 
nodes, weights = tensorrule(xa, wa, xa, wa, 2)
nodes, weights = balancednonsymmetricquad(nodes, weights, order-1)
println(nodes)
println(weights)
## 


