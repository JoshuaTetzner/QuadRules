using LinearAlgebra

function opf(
    sgsys::NestedSystem,
    sys::NestedSystem,
    x::T,
    n::Int,
    t
) where {T <: AbstractFloat}
    
    return (1 - t) * sys.pl(sys.systems[n], sys.segments, x) + 
        t * sgsys.pl(sgsys.systems[n], sgsys.segments, x)
end

function dopf(
    sgsys::NestedSystem,
    sys::NestedSystem,
    x::T,
    n::Int,
    t
) where {T <: AbstractFloat}

    return (1 - t) * sys.dpl(sys.systems[n], sys.segments, x) + 
        t * sgsys.dpl(sgsys.systems[n], sgsys.segments, x)
end

function intopf(sgsys::NestedSystem, sys::NestedSystem, n::Int, t)
    return (1 - t) * sys.intpl[n] + t * sgsys.intpl[n]
end

function getαβ(
    sgsys::NestedSystem,
    sys::NestedSystem,
    x::Vector{T},
    t
) where {T <: AbstractFloat}

    n = length(x)
    len = length(sys.systems)
    α = zeros(T, n, 2*n)
    β = zeros(T, n, 2*n)
    ση = zeros(T, 2*n, 2*n)
        
    @threads for i = 1:n
        @threads for j = 1:len
            ση[i, j] = opf(sgsys, sys, x[i], j, t)
            ση[i+n, j] = dopf(sgsys, sys, x[i], j, t)
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
    sgsys::NestedSystem,
    sys::NestedSystem,
    α::Matrix{T},
    β::Matrix{T},
    t
) where {T <: AbstractFloat}

    len = size(α)[1]
    lensys = length(sys.systems)
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

function nestedquadrature(
    sgsys::NestedSystem,
    sys::NestedSystem,
    x::Vector{T};
    tol=1e-16
) where {T <: AbstractFloat}
    toltype = Float64
    if tol < eps(Float64)
        println("Switched to BigFloat")
        toltype = BigFloat
    end
    x = toltype.(x)
    sgsys.systems = convert(Vector{Matrix{toltype}}, sgsys.systems)
    sys.systems = convert(Vector{Matrix{toltype}}, sys.systems)
    sgsys.segments = convert(Vector{toltype}, sgsys.segments)
    sys.segments = convert(Vector{toltype}, sys.segments)
    sgsys.intpl = convert(Vector{toltype}, sgsys.intpl)
    sys.intpl = convert(Vector{toltype}, sys.intpl)


    t = 0.0
    α, β = getαβ(sgsys, sys, x, t)
    intσ, intη = intση(sgsys, sys, α, β, t)

    println("Progress:")
    print(".")
    step = 0.01
    ϵ = 1
    while t <= 1.1
        print(".")
        xnew = x
        iter = 0
        ϵ = 1
        ϵ2 = 1
        mult = 1
        
        while iter < 20 && !isapprox(0, ϵ, atol = tol) && 
            ϵ <= 10 && ϵ2 <= 10 && maximum(abs.(xnew)) < 1.0 

            iter+=1
            α, β = getαβ(sgsys, sys, xnew, t)
            intσ, intη = intση(sgsys, sys, α, β, t)
            diff = mult .* intσ ./ intη 
            ϵ = norm(diff)
            ϵ2 = norm(intσ)
            xnew += diff 
        end

        if isapprox(ϵ, 0, atol=tol)
            x = xnew
            if t >= 1            
                break
            end
            if step < 0.01
                step = step * 10 
            end
            t += step
        else
            if isapprox(step, 0, atol = 1e-10)
                break
            else
                t -= step             
                step = step/10
                t += step
            end
        end
        
        if t > 1
            step -= t-1
            t = 1 
        end
        println(ϵ)
        println(x)
        println(t)    
    end

    α, β = getαβ(sgsys, sys, x, t)
    intσ, w = intση(sgsys, sys, α, β , t)
    return x, w, ϵ
end
