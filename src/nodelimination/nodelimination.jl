using LinearAlgebra

function fmat(fct, x::Vector{T}) where {T <: AbstractFloat}
    F = zeros(T, length(fct))
    for i = 1:length(fct)
        for j = 1:3:length(x)
            F[i] += fct[i](x[j], x[j+1])*x[j+2]
        end
    end

    return F
end

function getpolynomes(p::Int)
    Φ=[]
    for i = 0:p
        for j = 0:p-i
            ϕ(x,y) = legendre(i,x)*legendre(j,y) * sqrt((2*i + 1)/2)* sqrt((2*j + 1)/2)  
            push!(Φ, ϕ)    
        end
    end
    return Φ
end

function getpolynomes_dx(p::Int)
    Φ=[]
    for i = 0:p
        for j = 0:p-i
            ϕ(x,y) = (i + 1)/(1 - x^2) * (x * legendre(i, x) - legendre(i+1, x)) * 
                legendre(j,y) * sqrt((2*i + 1)/2) * sqrt((2*j + 1)/2)   
            push!(Φ, ϕ)    
        end
    end
    return Φ
end

function getpolynomes_dy(p::Int)
    Φ=[]
    for i = 0:p
        for j = 0:p-i
            ϕ(x,y) = (j + 1)/(1 - y^2) * (y * legendre(j, y) - legendre(j+1, y)) * 
                legendre(i,x) * sqrt((2*i + 1)/2) * sqrt((2*j + 1)/2)  
            push!(Φ, ϕ)    
        end
    end
    return Φ
end

function jacobian(fdx, fdy, fdw, x::Vector{T}) where {T <: AbstractFloat}
    J = zeros(T, length(fdx), length(x))
    for i = 1:length(fdx)
        for j = 1:3:length(x)
            J[i, j] = fdx[i](x[j], x[j+1]) * x[j+2]
            J[i, j+1] = fdy[i](x[j], x[j+1]) * x[j+2]
            J[i, j+2] = fdw[i](x[j], x[j+1])
        end
    end
    return J
end

function getA(x::Vector{T}, Φ) where {T <: AbstractFloat}
    A = zeros(T, length(Φ), Int(length(x)/3))
    for j = 1:Int(length(x)/3)
        for i = 1:length(Φ)
            A[i, j] = Φ[i](x[(j-1)*3+1], x[(j-1)*3+2])
        end
    end
    return A
end


function nonsymmetricquad(
    nodes::Matrix{T},
    weights::Vector{T},
    order::Int
) where {T <: AbstractFloat}

    print("Order: ")
    println(order)
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

    for k = (n-1):-1:1        
        delnode = x[(3*k+1):(3*k+3)]
        pop!(x)
        pop!(x)
        pop!(x)
        savex = x
        currentindex = k
        
        while currentindex != 1
            currentindex -= 1
            iter = 0
            ϵ=1
            while iter < 10 && !isapprox(ϵ, 0, atol=1e-14)
                iter += 1
                J = jacobian(dΦ_x, dΦ_y, Φ, x) 
                fct = fmat(Φ, x)
                x -= pinv(Float64.(J))*(fct-int_f)
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