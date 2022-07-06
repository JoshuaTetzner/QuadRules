using Base.Threads

mutable struct NestedSystem
    systems::Vector{Matrix}
    segments::Vector{BigFloat}
    pl::Function
    dpl::Function
    intpl::Vector
end

function nestednodes(nsegments::Int, nnodes::Int, fct)
    function segments(n)
        s = zeros(BigFloat, n)
        s[1] = big(0)
        s[n] =big(1)
        for i = 2:n-1
            s[i] = s[1] + 1/2^(n-i)
        end
        return s
    end

    δ(a) = ==(a,0)
    s = (segments(nsegments))
    segmentnodes = zeros(BigFloat, nnodes, nsegments-1)
    for (nc, col) in enumerate(eachcol(segmentnodes))
        xn = [
            1/2 * (s[nc+1] + s[nc]) + 
            1/2 * (s[nc+1] - s[nc]) * cos(pi * big(k - 1 / 2) / big(nnodes))
            for k = 1:nnodes
        ]
        for (nn, node) in enumerate(col)
            segmentnodes[nn, nc] = sum([
                (2-δ(nn-1)) / big(nnodes) * fct(xn[i]) * schebychev(nn-1, xn[i], s[nc], s[nc+1]) 
                for i = 1:nnodes
            ])
        end
    end

    return segmentnodes, s
end

function nestedsystem(order, nsegments, nnodes, fct)

    systems = Vector{Matrix}(undef, order+1)
    integrals = zeros(BigFloat, order+1)
    seg = []

    @threads for nf = 0:order 
        f(x) = fct(nf, x)
        systems[nf+1], seg = nestednodes(nsegments, nnodes, f)
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
##
