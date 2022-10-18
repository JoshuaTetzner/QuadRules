using LinearAlgebra
using Base.Threads

function nonsymmetricquad(
    nodes::Matrix{T},
    weights::Vector{T},
    order::Int
) where {T <: AbstractFloat}

    #print("Order: ")
    #println(order)
    Φ = getpolynomes(order)
    dΦ_x = getpolynomes_dx(order)
    dΦ_y = getpolynomes_dy(order)
    int_f = zeros(T, length(Φ))
    int_f[1] = 2


    n = length(weights)

    x = zeros(T, 3*length(weights))
    for i = 1:length(weights)
        x[(i-1)*3+1] = nodes[i,1]
        x[(i-1)*3+2] = nodes[i,2]
        x[(i-1)*3+3] = weights[i]
    end
    #println(norm(int_f - getA(x, Φ) * x[3:3:end]))
    for k = (n-1):-1:1        
        delnode = x[(3*k+1):(3*k+3)]
        pop!(x)
        pop!(x)
        pop!(x)
        savex = x
        
        # descending order
        #sindex = [i for i = k:-1:1]
        # sum(phi(x))
        #sindex = [sum([Φ[j](x[3*i+1], x[3*i+2])^2 for j in eachindex(Φ)]) for i = 0:(k-1)] 
        # sum(phi(x)) * w
        sindex = [sum([Φ[j](x[3*i+1], x[3*i+2])^2 for j in eachindex(Φ)]) * x[3*i+3] for i = 0:(k-1)] 

        counter = 1
        while counter <= k
            counter += 1
            currentindex = argmin(sindex)-1
            sindex[currentindex+1] = maximum(sindex)+1
            iter = 0
            ϵ=1
            while iter < 10 && !isapprox(ϵ, 0, atol=1e-14)
                iter += 1
                J = jacobian(dΦ_x, dΦ_y, Φ, x) 
                fct = fmat(Φ, x)
                x -= pinv(Float64.(J))*(fct-int_f)
                # check interior points and positive weights
                @threads for i = 0:(k-1)
                    if abs(x[i*3+1]) > 1
                        x[i*3+1] = sign(x[i*3+1])*0.9
                    end
                end

                @threads for i = 0:(k-1)
                    if abs(x[i*3+2]) > 1
                        x[i*3+2] = sign(x[i*3+2])*0.9
                    end
                end

                @threads for i = 0:(k-1)
                    if x[i*3+3] < 0
                        x[i*3+3] = eps(Float64)
                    end
                end
                
                ϵ = norm(int_f - getA(x, Φ) * x[3:3:(3*k)]) 
            end

            if !isapprox(ϵ, 0, atol=1e-14) && iter == 10 
                x=savex
                x[(3*currentindex+1):(3*currentindex+3)], delnode = 
                    delnode, x[(3*currentindex+1):(3*currentindex+3)] 
            else
                print("n-Points: ")
                println(k)
                nodes = [x[1:3:3*(k)] x[2:3:3*(k)]]
                weights = x[3:3:3*(k)]
                break
            end
        end
        if length(weights) == k+1
            return nodes, weights
        end
    end
end