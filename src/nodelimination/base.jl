using Base.Threads
using LinearAlgebra

function fmat(fct, x::Vector{T}) where {T <: AbstractFloat}
    F = zeros(T, length(fct))
    for i in eachindex(fct)
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
    @threads for i in eachindex(fdx)
        @threads for j = 1:3:length(x)
            J[i, j] = fdx[i](x[j], x[j+1]) * x[j+2]
            J[i, j+1] = fdy[i](x[j], x[j+1]) * x[j+2]
            J[i, j+2] = fdw[i](x[j], x[j+1])
        end
    end
    return J
end

function getA(x::Vector{T}, Φ) where {T <: AbstractFloat}
    A = zeros(T, length(Φ), Int(length(x)/3))
    @threads for j = 1:Int(length(x)/3)
        @threads for i in eachindex(Φ)
            A[i, j] = Φ[i](x[(j-1)*3+1], x[(j-1)*3+2])
        end
    end
    return A
end

function checkinterior(x::Vector)
    k = Int(length(x)/3)
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

    return x
end

# 3D
function fmat3D(fct, x::Vector{T}) where {T <: AbstractFloat}
    F = zeros(T, length(fct))
    for i in eachindex(fct)
        for j = 1:4:length(x)
            F[i] += fct[i](x[j], x[j+1], x[j+2])*x[j+3]
        end
    end

    return F
end

function getpolynomes3D(p::Int)
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

function getpolynomes_dx3D(p::Int)
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

function getpolynomes_dy3D(p::Int)
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

function getpolynomes_dz3D(p::Int)
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

function jacobian3D(fdx, fdy, fdz, fdw, x::Vector{T}) where {T <: AbstractFloat}
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

function getA3D(x::Vector{T}, Φ) where {T <: AbstractFloat}
    A = zeros(T, length(Φ), Int(length(x)/4))
    @threads for j = 1:Int(length(x)/4)
        @threads for i in eachindex(Φ)
            A[i, j] = Φ[i](x[(j-1)*4+1], x[(j-1)*4+2], x[(j-1)*4+3])
        end
    end
    return A
end

function checkinterior3D(x::Vector)
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