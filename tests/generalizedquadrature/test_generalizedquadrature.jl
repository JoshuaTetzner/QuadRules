using Test

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

#Test for quadratures up to 10 points
for n = 3:10
    println(n)
    for tn = 0:(2*n - 1)
        solution = sum([logquadw[n-2][j]*logfct(tn, logquadx[n-2][j]) for j = 1:length(logquadw[n-2])])
        @test solution ≈ intlogfct(tn, 1)-intlogfct(tn, 0) atol=eps(Float64)
        println(solution - (intlogfct(tn, 1)-intlogfct(tn, 0)))
    end
end


##
logfct(x)=log(x+1)
n=3
solution = sum([w[j]*logfct(x[j]) for j = 1:length(w)]) - (2*log(big(2))-2)