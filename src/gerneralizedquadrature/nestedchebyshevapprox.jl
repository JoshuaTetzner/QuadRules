using Base.Threads
using QuadGK

mutable struct NestedSystem{T <: AbstractFloat}
    systems::Vector{Matrix}
    segments::Vector
    pl::Function
    dpl::Function
    intpl::Vector{T}
end

#function NestedSystem(system, segments, pl, dpl, intpl)
#    return NestedSystem(system, segments, pl, dpl, intpl)
#end


function segments(n::Int)
    s = zeros(Float64, n)
    s[1] = 0
    s[n] = 1
    for i = 2:n-1
        s[i] = s[1] + 3/big(2)^(n-i)
    end
    return s
end

function nestednodes(nsegments::Int, nnodes::Int, fct; T=BigFloat)
    function segments(n)
        s = zeros(n)
        s[1] = -1
        s[n] = 1
        for i = 2:n-1
            s[i] = s[1] + 3/2^(n-i)
        end
        return s
    end

    δ(a) = ==(a,0)
    s = segments(nsegments)
    segmentnodes = zeros(T, nnodes, nsegments-1)
    
    for (nc, col) in enumerate(eachcol(segmentnodes))
        xn = [
            1/2 * T(s[nc+1] + s[nc]) + 
            1/2 * T(s[nc+1] - s[nc]) * cos(pi * T(k - 1 / 2) / T(nnodes))
            for k = 1:nnodes
        ]
        for (nn, node) in enumerate(col)
            segmentnodes[nn, nc] = sum([
                (2-δ(nn-1)) / T(nnodes) * fct(xn[i]) * schebychev(nn-1, xn[i], s[nc], s[nc+1]) 
                for i = 1:nnodes
            ])
        end
    end

    return segmentnodes, s
end

function nestedsystem(order::Int, nsegments::Int, nnodes::Int, fct; T=BigFloat)

    systems = Vector{Matrix}(undef, order+1)
    integrals = zeros(T, order+1)
    seg = []

    @threads for nf = 0:order 
        f(x) = fct(nf, x)
        systems[nf+1], seg = nestednodes(nsegments, nnodes, f, T)
    end
    function approxpl(segmentnodes, seg, x)
        nnodes = size(segmentnodes)[1]
        nseg = size(segmentnodes)[2]
        for n = 1:nseg
            if x <= seg[n+1]
                return sum([segmentnodes[k, n]*schebychev(k-1, x, seg[n], seg[n+1]) for k = 1:nnodes])
                break
            end
        end
    end

    function approxdpl(segmentnodes, seg, x)
        nnodes = size(segmentnodes)[1]
        nseg = size(segmentnodes)[2]
        for n = 1:nseg
            if x <= seg[n+1]
                return sum([
                    segmentnodes[k, n] * 
                    dschebychev(k-1, x, seg[n], seg[n+1]) for k = 1:nnodes
                ])
                break
            end
        end
    end

    return NestedSystem(systems, seg, approxpl, approxdpl, integrals)
end
