using FileIO
using JLD2
using FastGaussQuadrature

## Single cubature for a complete polynomial system up to given order.
order = 2
points = Int(ceil((order+1)/2))
xa, wa = gausslegendre(points) 
nodes, weights = tensorrule(xa, wa, xa, wa, 2)
nodes, weights = nonsymmetricquad(nodes, weights, order)

## Multible cubatures for complete polynomials from order "a" upto order "b"
function start(a::Int, b::Int)
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

start(2,10)
