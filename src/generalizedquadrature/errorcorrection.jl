using Base.Threads
using LinearAlgebra

mutable struct TrueSystem
    pl::Function
    dpl::Function
    intpl::Vector
end

function correctlog(x::Vector{BigFloat})
    function logfct(n, x)
        if iseven(n)
            return x^(floor(n/2))
        else
            return x^(floor(n/2))*log(x) 
        end
    end
    
    function intlogfct(n, x)
        n = big(n)
        if iseven(n)
            1/big(floor(n/2)+1) * x^Int(floor(n/2)+1)
        else
            if x == 0
                return 0
            else
                newn = floor(n/2)+1
                x^newn*(newn*log(big(x))-1)/big(newn)^2
            end
        end
    end
    
    function dlogfct(n, x)
        n = big(n)
        newn = floor(n/2)
        if iseven(n)
            if newn == 0
                return 0
            else
                return newn*x^(newn-1)
            end
        else
            if newn > 0
                return newn*x^(newn-1)*log(big(x)) + x^(newn-1)
            else
                return 1/big(x)
            end
        end
    end
    
    function opf(
        sgsys,
        sys,
        x::T,
        n::Int,
        t
    ) where {T <: AbstractFloat}
    
        return (1 - t) * sys.pl(n, x) + 
            t * sgsys.pl(n, x)
    end
    
    function dopf(
        sgsys,
        sys,
        x::T,
        n::Int,
        t
    ) where {T <: AbstractFloat}
    
        return (1 - t) * sys.dpl(n, x) + 
            t * sgsys.dpl(n, x)
    end
    
    function intopf(sgsys, sys, n::Int, t)
        return (1 - t) * sys.intpl[n] + t * sgsys.intpl[n]
    end
    
    function getαβ(
        sgsys,
        sys,
        x::Vector{T},
        t
    ) where {T <: AbstractFloat}
    
        n = length(x)
        len = length(sys.intpl)
        α = zeros(T, n, 2*n)
        β = zeros(T, n, 2*n)
        ση = zeros(T, 2*n, 2*n)
            
        @threads for i = 1:n
            @threads for j = 1:len
                ση[i, j] = opf(sgsys, sys, x[i], j-1, t)
                ση[i+n, j] = dopf(sgsys, sys, x[i], j-1, t)
            end
        end
    
        @threads for index = 1:n
            fα = zeros(T, 2*n)
            fα[n+index] = T(1.0)
            fβ = zeros(T, 2*n)
            fβ[index] = T(1.0)
            α[index, :] = ση \ fα
            β[index, :] = ση \ fβ
        end
    
        return α, β
    end
    
    function intση(
        sgsys,
        sys,
        α::Matrix{T},
        β::Matrix{T},
        t
    ) where {T<:AbstractFloat}
      
        len = size(α)[1]
        lensys = length(sys.intpl)
        integralσ = zeros(T, len)
        integralη = zeros(T, len)
    
        for j = 1:lensys
            integral = intopf(sgsys, sys, j, t)
            @threads for i = 1:len
                integralσ[i] += α[i, j] * integral
                integralη[i] += β[i, j] * integral
            end
        end
    
        return integralσ, integralη
    end

    ints = []
    for i = 0:(2*length(x)-1)
        push!(ints, intlogfct(i, big(1))-intlogfct(i,big(0)))
    end
    t=0
    sys = TrueSystem(logfct, dlogfct, ints)
    sgsys = sys
    i=0
    α, β = getαβ(sgsys, sys, x, t)
    intσ, intη = intση(sgsys, sys, α, β, t)
    while i < 10 && norm(intσ) > 1e-64
        i += 1
        α, β = getαβ(sgsys, sys, x, t)
        intσ, intη = intση(sgsys, sys, α, β, t)
        x += intσ ./ (intη)  
    end
    α, β = getαβ(sgsys, sys, x, t)
    intσ, intη = intση(sgsys, sys, α, β, t)

    return x, intη
end

