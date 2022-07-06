using LinearAlgebra

function pl(n, x::T) where {T <: AbstractFloat}
    return legendre(n, x)
end
function dpl(n, x::T)  where {T <: AbstractFloat}
    return dlegendre(n, x)
end

function ln0(w::T, j, k)  where {T <: AbstractFloat}
    return w * pl(j, 0.0) * pl(k, 0.0)
end

function ln0_dw(j, k)
    return pl(j, 0.0) * pl(k, 0.0)
end

function ln1(w::T, ξ::T, j , k) where {T <: AbstractFloat}
    return 2 * w * (pl(j, ξ) * pl(k, 0.0)  + pl(j, 0.0)  * pl(k, ξ))
end

function ln1_dw(ξ::T, j, k) where {T <: AbstractFloat}
    return 2 * (pl(j, ξ) * pl(k, 0.0)  + pl(j, 0.0)  * pl(k, ξ))
end

function ln1_dξ(w::T, ξ::T, j , k)  where {T <: AbstractFloat}
    return 2 * w * (dpl(j, ξ) * pl(k, 0.0) + pl(j, 0.0) * dpl(k, ξ))
end

function ln2(w::T, ξ::T, j , k) where {T <: AbstractFloat}
    return 4 * w * pl(j, ξ)*pl(k, ξ)
end

function ln2_dw(ξ::T, j , k) where {T <: AbstractFloat}
    return 4 * pl(j, ξ)*pl(k, ξ)
end

function ln2_dξ(w::T, ξ::T, j , k) where {T <: AbstractFloat}
    return 4 * w *(dpl(j, ξ)*pl(k, ξ) + dpl(k, ξ)*pl(j, ξ))
end

function ln3(w::T, ξ::T, η::T, j, k) where {T <: AbstractFloat}
    return 4 * w * (pl(j, ξ) * pl(k, η) + pl(k, ξ) * pl(j, η))
end

function ln3_dw(ξ::T, η::T, j, k) where {T <: AbstractFloat}
    return 4 * (pl(j, ξ) * pl(k, η) + pl(k, ξ) * pl(j, η))
end

function ln3_dξ(w::T, ξ::T, η::T, j, k) where {T <: AbstractFloat}
    return 4 * w * (dpl(j, ξ) * pl(k, η) + dpl(k, ξ) * pl(j, η))
end

function ln3_dη(w::T, ξ::T, η::T, j, k) where {T <: AbstractFloat}
    return 4 * w * (pl(j, ξ) * dpl(k, η) + pl(k, ξ) * dpl(j, η))
end

function getcombinations(p, add)
    n = length([[j, k] for k = 0:2:p for j = k:2:p-k])
    n += add
    nxvec = []
    #calculate n3min
    if p > 7
        num = Int(ceil((p - 7) / 2))
        n3diff = 0
        if iseven(num)
            #upper
            num = Int(ceil(num/2))
            n3diff = num^2+num
        else
            #lower
            num = Int(ceil(num/2))
            n3diff = (num-1)^2+(num-1) +num
        end
        println(n3diff)
        n3 = Int(ceil(n3diff/3))
    else
        n3 = 0
    end
    
    for n3 = n3:Int(floor((n-4)/3))
        nleft = n - n3*3
        for n0 = 0:1
            nleft -= n0 
            for n1 = 1:Int(ceil((nleft-2)/2))
                n2 = Int(ceil((nleft-2*n1)/2))
                push!(nxvec, [n0, n1, n2, n3])     
            end
        end
    end
    
    return nxvec
end

function getsystem(nx, X::Vector{T}, order) where {T <: AbstractFloat}
    jk = [[j, k] for k = 0:2:order for j = k:2:order-k]
    F = zeros(T, length(jk))
    for (x, xjk) in enumerate(jk)
        xind = 0
        if nx[1] != 0
            F[x] += ln0(X[xind+1], xjk[1], xjk[2])
            xind += 1
        end
        for i = 1:nx[2]
            F[x] += ln1(X[xind+1], X[xind+2], xjk[1], xjk[2])
            xind += 2
        end
        for i = (nx[2] + 1):(nx[2] + nx[3])
            F[x] += ln2(X[xind+1], X[xind+2], xjk[1], xjk[2])
            xind += 2
        end
        for i = (nx[2] + nx[3] + 1):(nx[2] + nx[3] + nx[4])
            F[x] += ln3(X[xind+1], X[xind+2], X[xind+3], xjk[1], xjk[2])
            xind += 3
        end
        F[1] -= 4 / length(jk)
    end

    return F
end

function jacobian(nx, X::Vector{T}, order) where {T <: AbstractFloat}
    jk = [[j, k] for k = 0:2:order for j = k:2:order-k]
    J = zeros(T, length(jk), length(X))
    xind = 0
    for nxind = 1:4
        for nind = 1:nx[nxind]
            for y = 1:length(jk)
                if nxind == 1
                    J[y, xind+1] = ln0_dw(jk[y][1], jk[y][2])
                elseif nxind == 2
                    J[y, xind+1] = ln1_dw(X[xind+2], jk[y][1], jk[y][2])
                    J[y, xind+2] = ln1_dξ(X[xind+1],X[xind+2], jk[y][1], jk[y][2])
                elseif nxind == 3
                    J[y, xind+1] = ln2_dw(X[xind+2], jk[y][1], jk[y][2])
                    J[y, xind+2] = ln2_dξ(X[xind+1],X[xind+2], jk[y][1], jk[y][2])
                elseif nxind == 4
                    J[y, xind+1] = ln3_dw(X[xind+2], X[xind+3], jk[y][1], jk[y][2])
                    J[y, xind+2] = ln3_dξ(X[xind+1], X[xind+2], X[xind+3], jk[y][1], jk[y][2])
                    J[y, xind+3] = ln3_dη(X[xind+1], X[xind+2], X[xind+3], jk[y][1], jk[y][2])
                end
            end
            if nxind == 1
                xind += 1
            elseif nxind == 2
                xind += 2
            elseif nxind == 3
                xind += 2
            elseif nxind == 4
                xind += 3
            end
        end
    end

    return J
end

function expand(X::Vector{T}, nx) where {T <: AbstractFloat}
    nodes = []
    weights = []
    xind = 0
    for i = 1:4
        for j =1:nx[i]
            if i == 1
                push!(weights, X[xind+1])
                push!(nodes, [0, 0])
                xind += 1
            elseif i == 2
                push!(weights, X[xind+1])
                push!(nodes, [X[xind+2], 0])
                push!(weights, X[xind+1])
                push!(nodes, [-X[xind+2], 0])
                push!(weights, X[xind+1])
                push!(nodes, [0, X[xind+2]])
                push!(weights, X[xind+1])
                push!(nodes, [0, -X[xind+2]])
                xind += 2    
            elseif i == 3
                push!(weights, X[xind+1])
                push!(nodes, [X[xind+2], X[xind+2]])
                push!(weights, X[xind+1])
                push!(nodes, [-X[xind+2], X[xind+2]])
                push!(weights, X[xind+1])
                push!(nodes, [X[xind+2], -X[xind+2]])
                push!(weights, X[xind+1])
                push!(nodes, [-X[xind+2], -X[xind+2]])
                xind += 2
            elseif i == 4
                push!(weights, X[xind+1])
                push!(nodes, [X[xind+2], X[xind+3]])
                push!(weights, X[xind+1])
                push!(nodes, [X[xind+2], -X[xind+3]])
                push!(weights, X[xind+1])
                push!(nodes, [-X[xind+2], X[xind+3]])
                push!(weights, X[xind+1])
                push!(nodes, [-X[xind+2], -X[xind+3]])
                push!(weights, X[xind+1])
                push!(nodes, [X[xind+3], X[xind+2]])
                push!(weights, X[xind+1])
                push!(nodes, [X[xind+3], -X[xind+2]])
                push!(weights, X[xind+1])
                push!(nodes, [-X[xind+3], X[xind+2]])
                push!(weights, X[xind+1])
                push!(nodes, [-X[xind+3], -X[xind+2]])
                xind += 3
            end
        end
    end
    return nodes, weights
end

function symquadratur(
    nx,
    X::Vector{T},
    order;
    iterations=1000,
    tol=1e-16
) where {T <: AbstractFloat}

    Ient = 1*Matrix(I, length(X), length(X)) 
    λ = 0.01
    ϵ = 1
    iter = 0

    while iter < iterations && !isapprox(ϵ, 0, atol=tol)    
        F = getsystem(nx, X, order)
        J = jacobian(nx, X, order)
        
        #Levenberg-Marquart
        H = transpose(J)*J 
        Xdiff = pinv((H + λ .* (diag(H).* Ient))) * transpose(J) * F
        ϵ1 = norm(getsystem(nx, X-Xdiff, order))
    
        if ϵ1 < 2*ϵ
            ϵ = ϵ1
            X -= Xdiff
            λ = λ / 3
        else
            λ = λ * 2
        end
        iter+=1
    end

    nodes, weights = expand(X, nx)
    return X, nodes, weights
end