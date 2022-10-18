using JLD2
using FastGaussQuadrature
using QuadRules 

# asymmetric cubature for a complete system of polynomial depgree p = 5
p = 5
x, w = asymmetriccubature(5)

##

# asymmetric cubature for a complete polynomial system of degree n.
n = 4
order = 2*n-1 
xa, wa = gausslegendre(n) 
nodes, weights = tensorrule(xa, wa, xa, wa, 2)
nodes, weights = nonsymmetricquad(nodes, weights, order-1)
println(nodes)
println(weights)
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

##

# alternative way of generating cubature cubatures
for order = 3:15
    print("\n Order: ")
    println(order) 
    cmin = order^2
    min = Int(round(2.5*order))
    max = 4*order 
    for ord = min:max
        nodes, weights = initialquad(order, ord)
        println(ord)
        nodes, weights = nonsymmetricquad(nodes, weights, order)
        print("Points: ")
        println(length(weights))
        if length(weights) < cmin
            cmin = length(weights)
            dict = Dict{String, Any}(string(order) => Dict("weights" => weights, "nodes" => nodes))
            save("nonsymmetric" * string(order) * ".jld2", dict)
        end
    end
end

##

# improving existing cubature rules by continuation elimination
for i = 3:15
    quad = load("nonsymmetric" * string(i) * ".jld2")
    nodes = quad[string(i)]["nodes"]
    weights = quad[string(i)]["weights"]
    print("\n Order: ")
    println(i) 
    print("Points: ")
    println(length(weights))
    nodes, weights = contnonsymmetricquad2(nodes, weights, i)
    print("cont.-Points: ")
    println(length(weights))
end

