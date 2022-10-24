using LinearAlgebra
using Base.Threads

function nonsymmetricquad3D(
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


    n = length(weights)

    x = zeros(T, 4*length(weights))
    for i = 1:length(weights)
        x[(i-1)*4+1] = nodes[i,1]
        x[(i-1)*4+2] = nodes[i,2]
        x[(i-1)*4+3] = nodes[i,3]
        x[(i-1)*4+4] = weights[i]
    end

    ϵ = norm(int_f - getA3D(x, Φ) * x[4:4:end]) 
    #println(ϵ)
    for k = (n-1):-1:1        
        delnode = x[(4*k+1):(4*k+4)]
        pop!(x)
        pop!(x)
        pop!(x)
        pop!(x)
        savex = x
        
        Ient = 1*Matrix(I, length(x), length(x)) 
        # allocation
        J = zeros(Float64, length(dΦ_x), length(x))
        H = zeros(Float64, length(x), length(x))
        F = zeros(Float64, length(Φ))

        sindex = [sum([Φ[j](x[4*i+1], x[4*i+2], x[4*i+3])^2 for j in eachindex(Φ)]) * x[4*i+4] for i = 0:(k-1)] 

        counter = 1
        while counter <= k
            counter += 1
            currentindex = argmin(sindex)-1
            sindex[currentindex+1] = maximum(sindex)+1
            iter = 0
            ϵ = norm(int_f - getA3D(x, Φ) * x[4:4:(4*k)])
            λ = 0.1
            #println(currentindex)
            while iter < 500 && !isapprox(ϵ, 0, atol=1e-14)
                iter += 1
                #J = jacobian(dΦ_x, dΦ_y, dΦ_z, Φ, x) 
                #fct = fmat(Φ, x)
                #x -= pinv(Float64.(J))*(fct-int_f)

                #Levenberg-Marquart
                F .= fmat3D(Φ, x) - int_f
                J .= jacobian3D(dΦ_x, dΦ_y, dΦ_z, Φ, x)
                H .= transpose(J)*J 

                xdiff = pinv((H + λ .* (diag(H).* Ient))) * transpose(J) * F
                ϵ1 = norm(int_f - getA3D(x-xdiff, Φ) * x[4:4:(4*k)])
            
                if ϵ1 < ϵ
                    ϵ = ϵ1
                    x -= xdiff
                    if λ > 1e-3
                        λ = λ / 5
                    end
                else
                    if λ < 100
                        λ = λ * 4
                    else
                        iter = 500
                        break
                    end
                end
            
                x = checkinterior3D(x)
                ϵ = norm(int_f - getA3D(x, Φ) * x[4:4:(4*k)]) 
            end
            if !isapprox(ϵ, 0, atol=1e-14) && iter == 500 
                x=savex
                x[(4*currentindex+1):(4*currentindex+4)], delnode = 
                    delnode, x[(4*currentindex+1):(4*currentindex+4)] 
            else
                #print("n-Points: ")
                println(k)
                nodes = [x[1:4:4*(k)] x[2:4:4*(k)] x[3:4:4*(k)]]
                weights = x[4:4:4*(k)]
                break
            end
        end
        if length(weights) == k+1
            return nodes, weights
        end
    end
end
