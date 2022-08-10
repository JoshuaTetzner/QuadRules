using Test

function monomial(n, x)
    return x^n
end

function intmonom(n, x)
    return 1/(n+1) * x^(n+1)
end

for (i , cub) in enumerate(csplx)
    order = 3 + i*2
    for xn = 0:order
        for yn = 0:(order-xn)
            integral =  sum([monomial(xn, cub[j][1]) * monomial(yn, cub[j][2]) * csplw[i][j] for j in eachindex(cub)])
            trueint = (intmonom(xn, 1) - intmonom(xn, -1))*(intmonom(yn, 1) - intmonom(yn, -1)) 
            @test integral ≈ trueint atol=3e-15
        end
    end
end

##