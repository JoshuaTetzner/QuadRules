using Test

using Base.Threads

f(x) = log(big(x+1))
function approxpl(segmentnodes, seg, x)
    incheb(x, a, b) = 2*(x-a)/big(b-a) - 1
    nnodes = size(segmentnodes)[1]
    nseg = size(segmentnodes)[2]
    for n = 1:nseg
        if x <= seg[n+1]
            return  sum([segmentnodes[k, n]*chebychev(k-1, incheb(x, seg[n], seg[n+1])) for k = 1:nnodes])
            break
        end
    end
end


@time system, seg = nestednodes(50, 50, f, T=BigFloat)
elements= 25    
error = []
seg[1] += 1e-64
for ns = 1:length(seg)-1
    step = (seg[ns+1]-seg[ns])/elements
    er(x) = abs(approxpl(system, seg, x) - log(big(x+1)))^2
    l2 = 0
    for j = 1:elements
        l2 += (er(seg[ns]+(j-1)*step) + er(seg[ns]+(j)*step))/2 * step
    end
    push!(error, (l2))
end

(Float64.(error))
##
for i = 2:length(seg)
    println(seg[i-1]-seg[i])
end
##
function logfct(n, x)
    if iseven(n)
        return x^Int(n/2)
    else
        return x^Int(floor(n/2))*log(x+1) 
    end
end
