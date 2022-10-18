using QuadRules

function start()
    for order = 3:6
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
                save("nonsymmetricnew" * string(order) * ".jld2", dict)
            end
        end
    end
end

start()