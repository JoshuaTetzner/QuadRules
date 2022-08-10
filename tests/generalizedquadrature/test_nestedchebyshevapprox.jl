function logfct(n, x)
    if iseven(n)
        return x^Int(n/2)
    else
        return x^Int(floor(n/2))*log(x+1) 
    end
end

N = 3
order = 2*N-1 
sys = nestedsystem(order, 50, 50, logfct, T=BigFloat)

sys.intpl .= 0
for (i, eachsys) in  enumerate(sys.systems)
    for (ns, seg) in enumerate(eachcol(eachsys))
        sys.intpl[i] += sum([
            coeff * 
            intschebychev(nc-1, sys.segments[ns], sys.segments[ns+1]) 
            for (nc, coeff) in enumerate(seg)
        ]) 
    end
end

trueints = [2, (2*log(big(2))-2), 0, 1, 2/big(3), ((6*log(big(2))-8)/9)]
for (i, int) in enumerate(trueints)
    @show @test sys.intpl[i] - int ≈ 0 atol=eps(Float64)
end
