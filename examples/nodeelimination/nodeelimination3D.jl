using JLD2
using FastGaussQuadrature

# asymmetric cubature for a complete polynomial system of degree n.
order = 5
xa, wa = gausslegendre(3) 
nodes, weights = tensorrule3DX(c5plSFx, c5plSFw, xa, wa, 3)
#nodes, weights = tensorrule3D(xa, wa, xa, wa, xa, wa, 3)
weights
nodes, weights = nonsymmetricquad3D(nodes, weights, order)
#nodes2, weights2 = contnonsymmetricquad3D(nodes, weights, order)

##
order = 4
ord = 10
nodes, weights = initialquad3D(order, ord)
println(ord)
nodes, weights = nonsymmetricquad3D(nodes, weights, order)
println(length(weights))

##
for order = 3:7
    print("\n Order: ")
    println(order) 
    cmin = order^2
    min = Int(round(2*order))
    max = 4*order 
    for ord = min:max
        nodes, weights = initialquad3D(order, ord)
        println(ord)
        nodes, weights = nonsymmetricquad3D(nodes, weights, order)
        print("Points: ")
        println(length(weights))
        if length(weights) < cmin
            cmin = length(weights)
            dict = Dict{String, Any}(string(order) => Dict("weights" => weights, "nodes" => nodes))
            save("3Dnonsymmetric" * string(order) * ".jld2", dict)
        end
    end
end


