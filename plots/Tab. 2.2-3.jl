

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

elements= 25
segm= [10, 20, 30, 50]
polyn = [10, 20, 50]
for nseg in segm
    for npol in polyn
        system, seg = nestednodes(nseg, npol, f, T=BigFloat)    
        error = []
        seg[1] += 1e-64
        for ns = 1:2
            step = (seg[ns+1]-seg[ns])/elements
            er(x) = abs(approxpl(system, seg, x) - log(big(x+1)))^2
            l2 = 0
            for j = 1:elements
                l2 += (er(seg[ns]+(j-1)*step) + er(seg[ns]+(j)*step))/2 * step
            end
            push!(error, (l2))
        end
        print("Segments: ")
        print(nseg)
        print(", Polynomials: ")
        println(npol)
        print("Error 1. segment: ")
        print(Float64.(error[1]))
        print(", Size 1. segment: ")
        print(Float64(seg[2]-seg[1]))
        println(" ")
        print("Error 2. segment: ")
        print(Float64.(error[2]))
        print(", Size 2. segment: ")
        print(Float64(seg[3]-seg[2]))
        println("\n")
    end
end
