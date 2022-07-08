using Test

natlog(x) = log(x+1)
segmentnodes, seg = nestednodes(50, 50, natlog, T=BigFloat)

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

function checkapprox(x)
    return log(big(x)+1) - approxpl(segmentnodes, seg, big(x))
end

@test isapprox(0, checkapprox(1e-6), atol=1e-16)  
@test isapprox(0, checkapprox(1e-3), atol=1e-16)  
@test isapprox(0, checkapprox(0.01), atol=1e-16)  
@test isapprox(0, checkapprox(0.1), atol=1e-16)
@test isapprox(0, checkapprox(0.5), atol=1e-16)
@test isapprox(0, checkapprox(1.0), atol=1e-16)


