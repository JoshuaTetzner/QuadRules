using FileIO
using JLD2
using FastGaussQuadrature

## 
# asymmetric cubature for a complete polynomial system of degree n.
n = 5
order = 2*n-1 
xa, wa = gausslegendre(n) 
nodes, weights = tensorrule(xa, wa, xa, wa, 2)
nodes, weights = nonsymmetricquad(nodes, weights, order)

## 
# reconstruction of asymmetric cubatures for complete polynomials from order "a" upto order "b"
# the significance indices can be selected in src/nodelimination/nodelimination.jl
function asymcub(a::Int, b::Int)
    order = a
    points = Int(ceil((order+1)/2))
    xa, wa = gausslegendre(points) 
    nodes, weights = tensorrule(xa, wa, xa, wa, 2)
    nodes, weights = nonsymmetricquad(nodes, weights, order)
    dict = Dict{String, Any}(string(order) => Dict("weights" => weights, "nodes" => nodes))
    save("nonsymmetric" * string(order) * ".jld2", dict)

    for order = (a+1):b
        points = Int(ceil((order+1)/2))
        xa, wa = gausslegendre(points) 
        nodes, weights = tensorrule(xa, wa, xa, wa, 2)
        nodes, weights = nonsymmetricquad(nodes, weights, order)
        dict = Dict{String, Any}(string(order) => Dict("weights" => weights, "nodes" => nodes))
        save("nonsymmetric" * string(orderS) * ".jld2", dict)
    end
end

asymcub(3,15)
