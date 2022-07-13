using JLD2

symquad = load("symquad5.jld2")
for i = 7:2:21
    dict = load("symquad" * string(i) * ".jld2")
    merge!(symquad, dict)
end
save("symmetricquad.jld2", symquad)

##
quadrature = load("symmetricquad.jld2")
for i = 5:2:21
    nodes = quadrature[string(i)]["nodes"]
    weights = quadrature[string(i)]["weights"]
    println(i)
    println(nodes)
    println(weights)
    #println(length(weights))
end

