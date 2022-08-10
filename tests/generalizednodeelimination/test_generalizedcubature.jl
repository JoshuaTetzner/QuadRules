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
            newn = big(floor(n/2)+1)
            x^newn*(newn*log(x)-1)/newn^2
        end
    end
end

function monomial(n, x)
    return x^n
end

function intmonom(n, x)
    return 1/(n+1) * x^(n+1)
end

#Test for quadratures up to 10 points
for n = 4:10
    x = clogx[n-3] .* 0.5 .+0.5
    w = clogw[n-3] ./ 4
    for tn1 = 0:n
        for tn2 = 0:(n - tn1)
            
            solution = sum([
                logfct(tn1, x[j,2])*monomial(tn2,x[j,1]) * w[j]
                for j = 1:length(clogw[n-3])
                ])
            @test solution ≈ ((intlogfct(tn1, 1)-intlogfct(tn1, 0))*(intmonom(tn2, 1)-intmonom(tn2, 0))) atol=eps(Float64)
        end
    end
end