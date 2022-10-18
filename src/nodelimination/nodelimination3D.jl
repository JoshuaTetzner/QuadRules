using LinearAlgebra
using Base.Threads

function fmat(fct, x::Vector{T}) where {T <: AbstractFloat}
    F = zeros(T, length(fct))
    for i in eachindex(fct)
        for j = 1:4:length(x)
            F[i] += fct[i](x[j], x[j+1], x[j+2])*x[j+3]
        end
    end

    return F
end

function getpolynomes(p::Int)
    Φ=[]
    for i = 0:p
        for j = 0:p-i
            for k = 0:(p-(i+j))
                ϕ(x,y,z) = legendre(i,x) * legendre(j,y) * legendre(k,z) 
                    sqrt((2*i + 1)/2) * sqrt((2*j + 1)/2) * sqrt((2*k + 1)/2)  
                push!(Φ, ϕ)   
            end 
        end
    end
    return Φ
end

function getpolynomes_dx(p::Int)
    Φ=[]
    for i = 0:p
        for j = 0:p-i
            for k = 0:(p-(i+j))
                ϕ(x,y,z) = (i + 1)/(1 - x^2) * (x * legendre(i, x) - legendre(i+1, x)) * 
                    legendre(j,y) * legendre(k,z) *
                    sqrt((2*i + 1)/2) * sqrt((2*j + 1)/2) * sqrt((2*k + 1)/2)   
                push!(Φ, ϕ)
            end    
        end
    end
    return Φ
end

function getpolynomes_dy(p::Int)
    Φ=[]
    for i = 0:p
        for j = 0:p-i
            for k = 0:(p-(i+j))
                ϕ(x,y,z) = (j + 1)/(1 - y^2) * (y * legendre(j, y) - legendre(j+1, y)) * 
                    legendre(i,x) * legendre(k,z) *
                    sqrt((2*i + 1)/2) * sqrt((2*j + 1)/2) * sqrt((2*k + 1)/2)
                    push!(Φ, ϕ) 
            end   
        end
    end
    return Φ
end

function getpolynomes_dz(p::Int)
    Φ=[]
    for i = 0:p
        for j = 0:p-i
            for k = 0:(p-(i+j))
                ϕ(x,y,z) = (k + 1)/(1 - z^2) * (z * legendre(k, y) - legendre(k+1, z)) * 
                    legendre(i,x) * legendre(j,y) *
                    sqrt((2*i + 1)/2) * sqrt((2*j + 1)/2) * sqrt((2*k + 1)/2)
                    push!(Φ, ϕ) 
            end   
        end
    end
    return Φ
end

function jacobian(fdx, fdy, fdz, fdw, x::Vector{T}) where {T <: AbstractFloat}
    J = zeros(T, length(fdx), length(x))
    @threads for i in eachindex(fdx)
        @threads for j = 1:4:length(x)
            J[i, j] = fdx[i](x[j], x[j+1], x[j+2]) * x[j+3]
            J[i, j+1] = fdy[i](x[j], x[j+1], x[j+2]) * x[j+3]
            J[i, j+2] = fdz[i](x[j], x[j+1], x[j+2]) * x[j+3]
            J[i, j+3] = fdw[i](x[j], x[j+1], x[j+2])
        end
    end
    return J
end

function getA(x::Vector{T}, Φ) where {T <: AbstractFloat}
    A = zeros(T, length(Φ), Int(length(x)/4))
    @threads for j = 1:Int(length(x)/4)
        @threads for i in eachindex(Φ)
            A[i, j] = Φ[i](x[(j-1)*4+1], x[(j-1)*4+2], x[(j-1)*4+3])
        end
    end
    return A
end

function checkinterior(x::Vector)
    k = Int(length(x)/4)
    @threads for i = 0:(k-1)
        if abs(x[i*4+1]) >= 1
            x[i*4+1] = sign(x[i*4+1])*0.99
        end
    end

    @threads for i = 0:(k-1)
        if abs(x[i*4+2]) >= 1
            x[i*4+2] = sign(x[i*4+2])*0.99
        end
    end

    @threads for i = 0:(k-1)
        if abs(x[i*4+3]) >= 1
            x[i*4+3] = sign(x[i*4+3])*0.99
        end
    end

    @threads for i = 0:(k-1)
        if x[i*4+4] < 0
            x[i*4+4] = eps(Float64)
        end
    end

    return x
end

function nonsymmetricquad3D(
    nodes::Matrix{T},
    weights::Vector{T},
    order::Int
) where {T <: AbstractFloat}

    print("Order: ")
    println(order)
    Φ = getpolynomes(order)
    dΦ_x = getpolynomes_dx(order)
    dΦ_y = getpolynomes_dy(order)
    dΦ_z = getpolynomes_dz(order)
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

    ϵ = norm(int_f - getA(x, Φ) * x[4:4:end]) 
    println(ϵ)
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
            ϵ = norm(int_f - getA(x, Φ) * x[4:4:(4*k)])
            λ = 0.1
            println(currentindex)
            while iter < 500 && !isapprox(ϵ, 0, atol=1e-14)
                iter += 1
                #J = jacobian(dΦ_x, dΦ_y, dΦ_z, Φ, x) 
                #fct = fmat(Φ, x)
                #x -= pinv(Float64.(J))*(fct-int_f)

                #Levenberg-Marquart
                F .= fmat(Φ, x) - int_f
                J .= jacobian(dΦ_x, dΦ_y, dΦ_z, Φ, x)
                H .= transpose(J)*J 

                xdiff = pinv((H + λ .* (diag(H).* Ient))) * transpose(J) * F
                ϵ1 = norm(int_f - getA(x-xdiff, Φ) * x[4:4:(4*k)])
            
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
            
                x = checkinterior(x)
                ϵ = norm(int_f - getA(x, Φ) * x[4:4:(4*k)]) 
            end
            if !isapprox(ϵ, 0, atol=1e-14) && iter == 500 
                x=savex
                x[(4*currentindex+1):(4*currentindex+4)], delnode = 
                    delnode, x[(4*currentindex+1):(4*currentindex+4)] 
            else
                print("n-Points: ")
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
