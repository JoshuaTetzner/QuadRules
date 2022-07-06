order = 5

n = length([[j, k] for k = 0:2:order for j = k:2:order-k])
nxvec = getcombinations(order, 0)

for i = 1:length(nxvec)
    n = nxvec[i][1] + 2*(nxvec[i][2] + nxvec[i][3]) + 3*nxvec[i][4] 
    X = rand(Float64, n)
    println(nxvec[i][:])
    X, nodes, weights = symquadratur(nxvec[i][:], X, order)
    
    f(x,y) = x^2*y^2
    println(sum([f(nodes[j][1],nodes[j][2])*weights[j] for j = 1:length(weights)]))
    #println(norm(getsystem(nxvec[i,:], X, order)))
    #println((weights))
    #println(X)
end

##

function symmetricparallel(maxorder::Int)
    f(x,y) = x^2*y^2
    found = zeros(Bool, maxorder)
    for order = 3:2:maxorder
        println(order)
        add = 0
        while add <= 2 && !found[order]
            nxvec = getcombinations(order, add)
            @threads for i = 1:length(nxvec)
                println(nxvec[i])
                count = 0
                while !found[order] && count < 20
                    count += 1
                    n = nxvec[i][1] + 2*(nxvec[i][2] + nxvec[i][3]) + 3*nxvec[i][4] 
                    X = rand(Float64, n)
                    X, nodes, weights = symquadratur(nxvec[i][:], X, order, 100)
                    intval = sum([f(nodes[i][1],nodes[i][2])*weights[i] for i = 1:length(weights)])
                    if norm(getsystem(nxvec[i][:], X, order)) < 1e-14 && 
                        isapprox(intval, 0.44444444444444, atol=1e-10)
                        if !found[order]
                            found[order] = true
                            dict = Dict{String, Any}(string(order) => Dict("weights" => weights, "nodes" => nodes))
                            save("symquad" * string(order)  * ".jld2", dict)
                        else
                            break
                        end
                    end
                end
            end
            add += 2
        end
    end
end


##
@time symmetricparallel(23)
##
using JLD2
symquad = load("symquad5.jld2")
for i = 7:2:23
    dict = load("symquad" * string(i) * ".jld2")
    merge!(symquad, dict)
end
save("symmetricquad.jld2", symquad)