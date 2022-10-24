using JLD2

symquad = load("3Dnonsymmetric3.jld2")
for i = 4:6
    dict = load("3Dnonsymmetric" * string(i) * ".jld2")
    merge!(symquad, dict)
end
save("3Dnonsymmetric.jld2", symquad)

#
quadrature = load("3Dnonsymmetric.jld2")
for i = 3:6
    nodes = quadrature[string(i)]["nodes"]
    weights = quadrature[string(i)]["weights"]
    println(i)
    #println(nodes)
    #println(weights)
    println(length(weights))
end

##
quadrature = load("symmetricquad.jld2")
for i = 5:2:21
    nodes = quadrature[string(i)]["nodes"]
    weights = quadrature[string(i)]["weights"]
    print((i))
    for j in eachindex(weights)
        print("&")
        print(Float64(nodes[j][1]))
        print("&")
        print(Float64(nodes[j][2]))
        print("&")
        print(Float64(weights[j]))
        println("\\\\")
    end
end

##
for i = 3:9
    print(i)
    x, w = generalizedcubature(i)
    for j in eachindex(w)
        print("&")
        print(Float64(w[j]))
        print("&")
        print(Float64(x[j,1]))
        print("&")
        print(Float64(x[j,2]))
        println("\\\\")
    end
end
##
for i = 4:11
    print(i)
    x, w = generalizedcubature(i, type = :logquad)
    for j in eachindex(w)
        print("&")
        print(Float64(w[j]))
        print("&")
        print(Float64(x[j,1]))
        print("&")
        print(Float64(x[j,2]))
        println("\\\\")
    end
end