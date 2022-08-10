function logfct(n, x)
    if iseven(n)
        return x^Int(n/2)
    else
        return x^Int(floor(n/2))*log(x+1) 
    end
end

segm = [10, 20 , 30 , 50]
trueintegrals = [2, 2*log(big(2))-2, 0, 1, 2/big(3), (6*log(big(2))-8)/big(9)]
for nseg in segm 
    print("Segments: ")
    println(nseg)
    N = 3
    order = 2*N-1 
    sgsys = nestedsystem(order, nseg, 50, logfct, T = BigFloat)
    for (i, val) in enumerate(trueintegrals)
        if val != 0
            println(Float64(abs((sgsys.intpl[i]-val)/val)))
        else
            println(Float64(abs(sgsys.intpl[i]-val)))
        end
    end
end
