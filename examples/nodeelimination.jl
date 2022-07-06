using FileIO
using JLD2
using Base.Threads
using FastGaussQuadrature

##
function start()
    d = 3
    points = Int(ceil((d+1)/2))
    xa, wa = gausslegendre(points) 
    nodes, weights = tensorrule(xa, wa, xa, wa, 2)
    nodes, weights = nonsymmetricquad(nodes, weights, d)
    dict = Dict{String, Any}(string(d) => Dict("weights" => weights, "nodes" => nodes))
    save("nonsymmetric" * string(d) * ".jld2", dict)

    for d = 4:25
        print("Order: ")
        println(d)
        points = Int(ceil((d+1)/2))
        xa, wa = gausslegendre(points) 
        nodes, weights = tensorrule(xa, wa, xa, wa, 2)
        nodes, weights = nonsymmetricquad(nodes, weights, d)
        dict = Dict{String, Any}(string(d) => Dict("weights" => weights, "nodes" => nodes))
        save("nonsymmetric" * string(d) * ".jld2", dict)
    end
end

##
start()
##
d = 5
points = Int(ceil((N+1)/2))
xa, wa = gausslegendre(points) 
nodes, weights = tensorrule(xa, wa, xa, wa, 2)
nodes, weights = nonsymmetricquad(nodes, weights, d)