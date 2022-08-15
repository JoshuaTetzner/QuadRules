using SpecialPolynomials
using LinearAlgebra
using Base.Threads

function fmat(fct, x::Vector)
    F = zeros(Float64, length(fct))
    for i in eachindex(fct)
        for j = 1:3:length(x)
            F[i] += fct[i](x[j], x[j+1])*x[j+2]
        end
    end

    return F
end

function getpolynomes(sysa, sysb, order)
    Φ=[]
    intf=[]
    for i = 0:order
        for j = 0:(order-i)
            ϕ(x,y) = sysa.pl(
                sysa.systems[i+1], 
                sysa.segments, 
                big.(x)
            ) * sysb.pl(sysb.systems[j+1], sysb.segments, big.(y))
            push!(Φ, ϕ)
            push!(intf, sysa.intpl[i+1]*sysb.intpl[j+1])    
        end
    end

    return Φ, intf
end

function getpolynomes_dx(sysa, sysb, order)
    Φ=[]
    for i = 0:order
        for j = 0:(order-i)
            ϕ(x,y) = sysa.dpl(
                sysa.systems[i+1], 
                sysa.segments, 
                big.(x)
            ) * sysb.pl(sysb.systems[j+1], sysb.segments, big.(y))
            push!(Φ, ϕ)    
        end
    end

    return Φ
end

function getpolynomes_dy(sysa, sysb, order)
    Φ=[]
    for i = 0:order
        for j = 0:(order-i)
            ϕ(x,y) = sysa.pl(
                sysa.systems[i+1], 
                sysa.segments, 
                big(x)
            ) * sysb.dpl(sysb.systems[j+1], sysb.segments, big(y))
            push!(Φ, ϕ)    
        end
    end

    return Φ
end

function jacobian(fdx, fdy, fdw, x::Vector)
    J = zeros(Float64, length(fdx), length(x))
    @threads for i in eachindex(fdx)
        @threads for j = 1:3:length(x)
            J[i, j] = fdx[i](x[j], x[j+1]) * x[j+2]
            J[i, j+1] = fdy[i](x[j], x[j+1]) * x[j+2]
            J[i, j+2] = fdw[i](x[j], x[j+1])
        end
    end
    return J
end

function getA(x::Vector, Φ)
    A = zeros(Float64, length(Φ), Int(length(x)/3))
    @threads for j = 1:Int(length(x)/3)
        @threads for i in eachindex(Φ)
            A[i, j] = Φ[i](x[(j-1)*3+1], x[(j-1)*3+2])
        end
    end
    return A
end

function nonsymmetricquad2(
    sysa,
    sysb, 
    xa::Vector{T},
    xb::Vector{T},
    wa::Vector{T},
    wb::Vector{T},
    order
) where {T <: AbstractFloat}
    if length(sysa.intpl) != length(sysb.intpl)
        return 0, 0
    else
        println(order)
        Φ, int_f = getpolynomes(sysa, sysb, order)
        dΦ_x = getpolynomes_dx(sysa, sysb, order)
        dΦ_y = getpolynomes_dy(sysa, sysb, order)

        nodes, weights = tensorrule(xa, wa, xb, wb, 2)
        n = length(weights)
        x = zeros(T, 3*length(weights))
        for i in eachindex(weights)
            x[(i-1)*3+1] = nodes[i,1]
            x[(i-1)*3+2] = nodes[i,2]
            x[(i-1)*3+3] = weights[i]
        end
        #println(norm(int_f - getA(x, Φ) * x[3:3:(3*n)])) 
        
        for k = (n-1):-1:1        
            delnode = x[(3*k+1):(3*k+3)]
            pop!(x)
            pop!(x)
            pop!(x)
            savex = x

            println(k)
            # descending order
            #sindex = [i for i = k:-1:1]
            # sum(phi(x))
            #sindex = [sum([Φ[j](x[3*i+1], x[3*i+2])^2 for j in eachindex(Φ)]) for i = 0:(k-1)] 
            # sum(phi(x)) * w
            sindex = [sum([Φ[j](x[3*i+1], x[3*i+2])^2 for j in eachindex(Φ)]) * x[3*i+3] for i = 0:(k-1)] 
            counter = 1
            while counter <= k
                counter += 1
                currentindex = argmin(sindex)-1
                sindex[currentindex+1] = maximum(sindex)+1

                iter = 0
                ϵ=1
                print("cindex: ")
                println(currentindex+1)
                while iter < 25 && !isapprox(ϵ, 0, atol=1e-14) && ϵ < 50
                    iter += 1
                    J = jacobian(dΦ_x, dΦ_y, Φ, x) 
                    fct = fmat(Φ, x)
                    x -= pinv(Float64.(J))*(fct-int_f)
                    for i = 1:3:3*k
                        if abs(x[i]) > 1
                            x[i] = sign(x[i])*1
                        end
                    end
                    for i = 2:3:3*k
                        if abs(x[i]) > 1
                            x[i] = sign(x[i])*1
                        end
                    end
                    for i = 3:3:3*k
                        if x[i] < 0
                            x[i] = eps(Float64)
                        end
                    end
                    ϵ = norm(int_f - getA(x, Φ) * x[3:3:(3*k)]) 
                end
                if !isapprox(ϵ, 0, atol=1e-14)
                    x=savex
                    x[(3*currentindex+1):(3*currentindex+3)], delnode = 
                        delnode, x[(3*currentindex+1):(3*currentindex+3)] 
                else
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
end

