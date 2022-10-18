using LinearAlgebra
using Base.Threads

function contnonsymmetricquad3D(
    nodes::Matrix{T},
    weights::Vector{T},
    order::Int
) where {T <: AbstractFloat}

    #print("Order: ")
    #println(order)
    Φ = getpolynomes3D(order)
    dΦ_x = getpolynomes_dx3D(order)
    dΦ_y = getpolynomes_dy3D(order)
    dΦ_z = getpolynomes_dz3D(order)
    int_f = zeros(T, length(Φ))
    int_f[1] = 8

    x = zeros(T, 4*length(weights))
    for i = 1:length(weights)
        x[(i-1)*4+1] = nodes[i,1]
        x[(i-1)*4+2] = nodes[i,2]
        x[(i-1)*4+3] = nodes[i,3]
        x[(i-1)*4+4] = weights[i]
    end
    println(norm(int_f - (getA3D(x, Φ) * x[4:4:end])))
    totalfailed = false
    while !totalfailed
        totalfailed = true
        
        sindices = [sum([Φ[j](x[4*i+1], x[4*i+2], x[4*i+3])^2 for j in eachindex(Φ)]) * 
            x[4*i+4] for i = 0:(Int(length(x)/4)-1)] 
 

        #print("NewPoints ")
        #println((Int(length(x)/3)-1))
        elimind = argmin(sindices)
        sindices[elimind] = maximum(sindices)+1
        
        delnode = x[(length(x)-3):end]
        pop!(x)
        pop!(x)
        pop!(x)
        pop!(x)
        
        print("n-Points: ")
        println(length(x)/4 + 1)
        Ient = 1*Matrix(I, length(x), length(x)) 
        J = zeros(Float64, length(dΦ_x), length(x))
        H = zeros(Float64, length(x), length(x))
        F = zeros(Float64, length(Φ))
        firstind = elimind
        if elimind != Int(length(x)/4 + 1)
            x[(4*elimind-3):(4*elimind)], delnode = 
                delnode, x[(4*elimind-3):(4*elimind)]
        end
        saverx = zeros(length(x)) 
        saver = zeros(4)
        saver .= delnode
        for ind = 1:Int(length(x)/4)
            
            print("Elim.-index: ")
            println(elimind)
            failed = false
            savex = x
            factor = 0.1
            while !isapprox(delnode[4], 0) && !failed 
                if delnode[4] <= 1e-4
                    delnode[4] = 0
                end
                saverx .= x
                iter = 0
                ϵ = norm(int_f - (getA3D(x, Φ) * x[4:4:end] + fmat3D(Φ, delnode)))
                λ = 0.01
                while iter != 200 && !isapprox(ϵ, 0, atol=1e-14)
                    iter += 1
                    F .= fmat3D(Φ, x) + fmat3D(Φ, delnode) - int_f
                    J .= jacobian3D(dΦ_x, dΦ_y, dΦ_z, Φ, x)

                    H .= transpose(J)*J 
                    xdiff = pinv((H + λ .* (diag(H).* Ient))) * transpose(J) * F
                    ϵ1 = norm(int_f - (getA3D(x-xdiff, Φ) * x[4:4:end] + fmat3D(Φ, delnode)))
                
                    if ϵ1 < ϵ
                        ϵ = ϵ1
                        x -= xdiff
                        if λ > 1e-8
                            λ = λ / 5
                        end
                    else
                        if λ < 1000
                            λ = λ * 4
                        else
                            iter = 200
                            break
                        end
                    end

                    x = checkinterior3D(x)

                   
                    
                end
                println(delnode[4])
                if isapprox(ϵ, 0, atol=1e-14)
                    if isapprox(delnode[4], 0)
                        println("eliminated")
                        break
                    else
                        if factor < 0.5
                            factor = factor * 2
                        end
                        delnode[4] = delnode[4] / (1+factor)
                    end
                else
                    if factor < 1e-2
                        failed = true
                    else
                        x = saverx
                        delnode[4] = delnode[4] * (1+factor)
                        factor = factor / 2
                        delnode[4] = delnode[4] / (1+factor)
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

                if elimind != Int(length(x)/4 + 1)
                    x[(4*elimind-3):(4*elimind)], delnode = 
                        saver, x[(4*elimind-3):(4*elimind)]
                    saver .= delnode
                else
                    x[(4*firstind-3):(4*firstind)], delnode = 
                        saver, x[(4*firstind-3):(4*firstind)]
                    saver .= delnode
                end
            end
        end
        if totalfailed
            push!(x, saver[1])
            push!(x, saver[2])
            push!(x, saver[3])
            push!(x, saver[4])
        end
    end
    nodes = [x[1:4:end] x[2:4:end] x[3:4:end]]
    weights = x[4:4:end]
    return nodes, weights
end