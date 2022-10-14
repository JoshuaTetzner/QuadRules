using LinearAlgebra


function getpolynomes(p::Int)
    Φ=[]
    for i = 0:p
        for j = 0:p-i
            ϕ(x,y) = legendre(i,x)*legendre(j,y)* sqrt((2*i + 1)/2)* sqrt((2*j + 1)/2)  
            push!(Φ, ϕ)    
        end
    end
    return Φ
end

function getA(x::Vector{T}, Φ) where {T <: AbstractFloat}
    A = zeros(T, Int(length(x)/3), length(Φ))
    @threads for j = 1:Int(length(x)/3)
        @threads for i in eachindex(Φ)
            A[j, i] = Φ[i](x[(j-1)*3+1], x[(j-1)*3+2]) * sqrt(x[(j-1)*3+3]) 
        end
    end
    return A
end

function pivotedGS(A::Matrix{T}) where {T <: AbstractFloat}
    ncol = size(A)[1]
    k = ncol
    p = collect(1:ncol)
    for i = 1:ncol
        μ = argmax([norm(A[j,:]) for j = i:ncol])
        if !isapprox(norm(A[μ, :]), 0, atol=1e-15)
            println(norm(A[μ, :]))
            A[i,:], A[μ,:] = A[μ,:], A[i,:]
            p[i], p[μ] = p[μ], p[i]
            for j = (i+1):ncol
                A[j, :] = A[j, :] - (dot(A[j, :], A[i, :]) / dot(A[i, :], A[i, :])) * A[i, :]
            end
        else
            k = i-1
        end
    end 

    return A, k, p
end


##
n = 8
order = 2*n-1 
xa, wa = gausslegendre(n) 
nodes, weights = tensorrule(xa, wa, xa, wa, 2)
x = zeros(Float64, 3*length(weights))
for (i, w) in enumerate(weights)
    x[(i-1)*3+1] = nodes[i,1]
    x[(i-1)*3+2] = nodes[i,2]
    x[(i-1)*3+3] = w
end
Φ = getpolynomes(order)
A = getA(x, Φ)

A, k, p = pivotedGS(A)

dot(A[64, :], A[57,:])

p
