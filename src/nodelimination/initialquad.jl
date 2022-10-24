using LinearAlgebra
using Base.Threads

function getwA(nodes::Matrix, weights::Vector{T}, Φ) where {T <: AbstractFloat}
    
    x = zeros(Float64, 3*length(weights))
    for (i, w) in enumerate(weights)
        x[(i-1)*3+1] = nodes[i,1]
        x[(i-1)*3+2] = nodes[i,2]
        x[(i-1)*3+3] = w
    end

    A = zeros(T, length(Φ), Int(length(x)/3))
    @threads for j = 1:Int(length(x)/3)
        @threads for i in eachindex(Φ)
            A[i, j] = Φ[i](x[(j-1)*3+1], x[(j-1)*3+2]) * sqrt(x[(j-1)*3+3]) 
        end
    end
    return A
end

function getA(nodes::Matrix{T}, Φ::Vector) where {T <: AbstractFloat}
    
    A = zeros(T, length(Φ), Int(length(nodes)/2))
    @threads for j = 1:Int(length(nodes)/2)
        @threads for i in eachindex(Φ)
            A[i, j] = Φ[i](nodes[j, 1], nodes[j, 2]) 
        end
    end
    return A
end

function initialquad(plorder, quadorder)
    xa, wa = gausslegendre(quadorder) 
    nodes, weights = tensorrule(xa, wa, xa, wa, 2)

    x = zeros(Float64, 3*length(weights))
    for (i, w) in enumerate(weights)
        x[(i-1)*3+1] = nodes[i,1]
        x[(i-1)*3+2] = nodes[i,2]
        x[(i-1)*3+3] = w
    end
    Φ = getpolynomes(plorder)
    wA = getwA(nodes, weights, Φ)
    
    _, _, p = qr(wA, Val(true))#ColumnNorm())
    newnodes = zeros(Float64, length(Φ), 2)
    for (i, val) in enumerate(p[1:length(Φ)])
        newnodes[i,1] = x[(val-1)*3+1] 
        newnodes[i,2] = x[(val-1)*3+2]
    end
    
    A = getA(newnodes, Φ)
    int = zeros(length(Φ))
    int[1] = 2
    newweights = inv(A) * int

    return newnodes, newweights
end

function getwA3D(nodes::Matrix, weights::Vector{T}, Φ) where {T <: AbstractFloat}
    
    x = zeros(Float64, 4*length(weights))
    for (i, w) in enumerate(weights)
        x[(i-1)*4+1] = nodes[i,1]
        x[(i-1)*4+2] = nodes[i,2]
        x[(i-1)*4+3] = nodes[i,3]
        x[(i-1)*4+4] = w
    end

    A = zeros(T, length(Φ), Int(length(x)/4))
    @threads for j = 1:Int(length(x)/4)
        @threads for i in eachindex(Φ)
            A[i, j] = Φ[i](x[(j-1)*4+1], x[(j-1)*4+2], x[(j-1)*4+3]) * sqrt(x[(j-1)*4+4]) 
        end
    end
    return A
end

function getA3D(nodes::Matrix{T}, Φ::Vector) where {T <: AbstractFloat}
    
    A = zeros(T, length(Φ), Int(length(nodes)/3))
    @threads for j = 1:Int(length(nodes)/3)
        @threads for i in eachindex(Φ)
            A[i, j] = Φ[i](nodes[j, 1], nodes[j, 2], nodes[j, 3]) 
        end
    end
    return A
end

function initialquad3D(plorder, quadorder)
    xa, wa = gausslegendre(quadorder) 
    nodes, weights = tensorrule3D(xa, wa, xa, wa, xa, wa)

    x = zeros(Float64, 4*length(weights))
    for (i, w) in enumerate(weights)
        x[(i-1)*4+1] = nodes[i,1]
        x[(i-1)*4+2] = nodes[i,2]
        x[(i-1)*4+3] = nodes[i,3]
        x[(i-1)*4+4] = w
    end
    Φ = getpolynomes3D(plorder)
    wA = getwA3D(nodes, weights, Φ)
    
    _, _, p = qr(wA, Val(true))
    newnodes = zeros(Float64, length(Φ), 3)
    for (i, val) in enumerate(p[1:length(Φ)])
        newnodes[i,1] = x[(val-1)*4+1] 
        newnodes[i,2] = x[(val-1)*4+2]
        newnodes[i,3] = x[(val-1)*4+3]
    end
    
    A = getA3D(newnodes, Φ)
    int = zeros(length(Φ))
    int[1] = 8
    newweights = inv(A) * int

    return newnodes, newweights
end
