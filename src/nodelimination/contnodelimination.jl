using LinearAlgebra
using Base.Threads

function contnonsymmetricquad(
    nodes::Matrix{T},
    weights::Vector{T},
    order::Int
) where {T <: AbstractFloat}

    Φ = getpolynomes(order)
    dΦ_x = getpolynomes_dx(order)
    dΦ_y = getpolynomes_dy(order)
    int_f = zeros(T, length(Φ))
    int_f[1] = 2

    x = zeros(T, 3*length(weights))
    for i = 1:length(weights)
        x[(i-1)*3+1] = nodes[i,1]
        x[(i-1)*3+2] = nodes[i,2]
        x[(i-1)*3+3] = weights[i]
    end

    totalfailed = false
    while !totalfailed
        totalfailed = true
        
        # descending order
        #sindices = [i for i = 0:(Int(length(x)/3)-1)]
        # sum(phi(x))
        #sindices = [sum([Φ[j](x[3*i+1], x[3*i+2])^2 for j in eachindex(Φ)]) for i = 0:(Int(length(x)/3)-1)] 
        # sum(phi(x)) * w
        sindices = [sum([Φ[j](x[3*i+1], x[3*i+2])^2 for j in eachindex(Φ)]) * x[3*i+3] for i = 0:(Int(length(x)/3)-1)] 

        elimind = argmin(sindices)
        sindices[elimind] = maximum(sindices)+1
        
        delnode = x[(length(x)-2):end]
        pop!(x)
        pop!(x)
        pop!(x)
        
        #print("n-Points: ")
        #println(length(x)/3 + 1)
        Ient = 1*Matrix(I, length(x), length(x)) 
        J = zeros(Float64, length(dΦ_x), length(x))
        H = zeros(Float64, length(x), length(x))
        F = zeros(Float64, length(Φ))
        firstind = elimind
        if elimind != Int(length(x)/3 + 1)
            x[(3*elimind-2):(3*elimind)], delnode = 
                delnode, x[(3*elimind-2):(3*elimind)]
        end
        saverx = zeros(length(x)) 
        saver = zeros(3)
        saver .= delnode
        for ind = 1:Int(length(x)/3)
            
            #print("Elim.-index: ")
            #println(elimind)
            failed = false
            savex = x
            factor = 0.2
            while !isapprox(delnode[3], 0) && !failed 
                if delnode[3] <= 1e-4
                    delnode[3] = 0
                end
                saverx .= x
                iter = 0
                ϵ = norm(int_f - (getA(x, Φ) * x[3:3:end] + fmat(Φ, delnode)))
               
                λ = 0.01
                while iter != 100 && !isapprox(ϵ, 0, atol=1e-14)

                    F .= fmat(Φ, x) + fmat(Φ, delnode) -int_f
                    J .= jacobian(dΦ_x, dΦ_y, Φ, x)
                   
                    #Levenberg-Marquart
                    H .= transpose(J)*J 
                    xdiff = pinv((H + λ .* (diag(H).* Ient))) * transpose(J) * F
                    ϵ1 = norm(int_f - (getA(x-xdiff, Φ) * x[3:3:end] + fmat(Φ, delnode)))
                
                    if ϵ1 < ϵ
                        ϵ = ϵ1
                        x -= xdiff
                        λ = λ / 3
                    else
                        λ = λ * 2
                        if λ > 10 
                            break
                        end
                    end

                    x = checkinterior(x)

                    iter += 1
                end

                if isapprox(ϵ, 0, atol=1e-14) && iter != 50
                    if isapprox(delnode[3], 0)
                        println("eliminated")
                        break
                    else
                        if factor < 0.2
                            factor = factor * 2
                        end
                        delnode[3] = delnode[3] / (1+factor)
                    end
                else
                    if factor < 0.05
                        failed = true
                    else
                        x = saverx
                        delnode[3] = delnode[3] * (1+factor)
                        factor = factor / 2
                        delnode[3] = delnode[3] / (1+factor)
                    end
                end
            end
            
            if !failed
                totalfailed = false
                break
            else
                x = savex
                elimind = argmin(sindices)
                sindices[elimind] = maximum(sindices)+1

                if elimind != Int(length(x)/3 + 1)
                    x[(3*elimind-2):(3*elimind)], delnode = 
                        saver, x[(3*elimind-2):(3*elimind)]
                    saver .= delnode
                else
                    x[(3*firstind-2):(3*firstind)], delnode = 
                        saver, x[(3*firstind-2):(3*firstind)]
                    saver .= delnode
                end
            end
        end
        if totalfailed
            push!(x, saver[1])
            push!(x, saver[2])
            push!(x, saver[3])
        end
    end
    nodes = [x[1:3:end] x[2:3:end]]
    weights = x[3:3:end]
    return nodes, weights
end


