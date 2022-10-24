using QuadRules
using ArgParse
using JLD2

function start(a, b)
    for order = a:b
        print("\n Order: ")
        println(order) 
        cmin = order^2
        min = Int(round(2.5*order))
        max = 4*order 
        for ord = min:max
            nodes, weights = initialquad3D(order, ord)
            println(ord)
            nodes, weights = nonsymmetricquad3D(nodes, weights, order)
            print("Points: ")
            println(length(weights))
            nodes, weights = contnonsymmetricquad3D(nodes, weights, order)
            print("cont.-Points: ")
            println(length(weights))
            if length(weights) < cmin
                cmin = length(weights)
                dict = Dict{String, Any}(string(order) => Dict("weights" => weights, "nodes" => nodes))
                save("nonsymmetric3D" * string(order) * ".jld2", dict)
            end
        end
    end
end

start(3,10)
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--start"
            help = "another option with an argument"
            arg_type = Int
            default = 3
        "--end"
            help = "another option with an argument"
            arg_type = Int
            default = 3

    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    start(parsed_args["start"], parsed_args["end"])
end
main()