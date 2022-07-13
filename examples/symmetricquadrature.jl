using FileIO
using JLD2
using Base.Threads

## Singel symmetric cubature rule for a complete set of polynomials with given order. 
# Symmetric Cubatures are not always found for random start values. 
# Function might have to be called more than once.
order = 4
nxvecs = getcombinations(order, 0)

for nxvec in nxvecs
    n = nxvec[1] + 2*(nxvec[2] + nxvec[3]) + 3*nxvec[4] 
    X = rand(Float64, n)
    println(nxvec[:])
    X, nodes, weights = symquadratur(nxvec[:], X, order)
    println(nodes)
    println(weights)
    #Check cubature for f(x,y) = x^2 y^2. True value = 0.4444444444444444. 
    f(x,y) = x^2*y^2
    println(sum([f(nodes[j][1],nodes[j][2])*weights[j] for j = 1:length(weights)]))
end

## Mulitble symmetric cubature rules for complete sets of polynomials up to maxorder. 
# Funciton is called multibel times until cubature is found. 
# Function stops solving combinations of n_0,...,n_3, if cubature is found. 
function symmetricparallel(maxorder::Int)
    f(x,y) = x^2*y^2
    found = zeros(Bool, maxorder)
    npoints = zeros(maxorder)
    for order = 3:2:maxorder
        println(order)
        add = 0
        while add <= 2 && !found[order]
            nxvec = getcombinations(order, add)
            for i = 1:length(nxvec)
                println(nxvec[i])
                count = 0
                while count < 20
                    count += 1
                    n = nxvec[i][1] + 2*(nxvec[i][2] + nxvec[i][3]) + 3*nxvec[i][4] 
                    X = rand(Float64, n)
                    X, nodes, weights = symquadratur(nxvec[i][:], X, order, iterations=200)
                    intval = sum([f(nodes[i][1],nodes[i][2])*weights[i] for i = 1:length(weights)])
                    if norm(getsystem(nxvec[i][:], X, order)) < 1e-14 && 
                        isapprox(intval, 0.44444444444444, atol=1e-10)
                        if !found[order]
                            found[order] = true
                            npoints[order] = length(weights)
                            dict = Dict{String, Any}(string(order) => Dict("weights" => weights, "nodes" => nodes))
                            save("symquad" * string(order)  * ".jld2", dict)
                        else
                            if length(weights) < npoints[order]
                                npoints[order] = length(weights)
                                dict = Dict{String, Any}(string(order) => Dict("weights" => weights, "nodes" => nodes))
                                save("symquad" * string(order)  * ".jld2", dict)
                            end
                        end
                    end
                end
            end
            if found[order]
                break
            else
                add += 2
            end
        end
    end
end

@time symmetricparallel(21)
