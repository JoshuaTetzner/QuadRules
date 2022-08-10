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
for i in eachindex(clogw)
    print(i+3)
    for j in eachindex(clogw[i])
        print("&")
        print(Float64(clogw[i][j]))
        print("&")
        print(Float64(clogx[i][j,2]))
        print("&")
        print(Float64(clogx[i][j,1]))
        println("\\\\")
    end
end
##
for i in eachindex(gclogw)
    print(i+2)
    for j in eachindex(gclogw[i])
        print("&")
        print(Float64(gclogw[i][j]))
        print("&")
        print(Float64(gclogx[i][j,1]))
        print("&")
        print(Float64(gclogx[i][j,2]))
        println("\\\\")
    end
end