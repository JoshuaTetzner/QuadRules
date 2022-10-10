using JLD2
using FastGaussQuadrature

# asymmetric cubature for a complete polynomial system of degree n.
n = 6
order = 2*n-1 
xa, wa = gausslegendre(n) 
nodes, weights = tensorrule(xa, wa, xa, wa, 2)
@time nodes, weights = balancednonsymmetricquad(nodes, weights, order-1)
#println(nodes)
#println(weights)
## 

function asymcub(a::Int, b::Int)
    order = a
    points = Int(ceil((order+1)/2))
    xa, wa = gausslegendre(points) 
    nodes, weights = tensorrule(xa, wa, xa, wa, 2)
    nodes, weights = balancednonsymmetricquad(nodes, weights, order)
    dict = Dict{String, Any}(string(order) => Dict("weights" => weights, "nodes" => nodes))
    save("nonsymmetric" * string(order) * ".jld2", dict)

    for order = (a+1):b
        points = Int(ceil((order+1)/2))
        xa, wa = gausslegendre(points) 
        nodes, weights = tensorrule(xa, wa, xa, wa, 2)
        nodes, weights = balancednonsymmetricquad(nodes, weights, order)
        dict = Dict{String, Any}(string(order) => Dict("weights" => weights, "nodes" => nodes))
        save("nonsymmetric" * string(order) * ".jld2", dict)
    end
end

asymcub(15,23)

