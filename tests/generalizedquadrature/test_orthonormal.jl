using Test

function chebyshevpolynomials(n, x)
    return x^n
end

N = 5
order = 2*N-1 
sys = nestedsystem(order, 2, 20, chebyshevpolynomials)
osys = gramschmidt(sys)

for i = 2:size(osys.systems)[1]
    @show @test dot(osys.systems[1,1], osys.systems[i,1]) ≈ 0 atol=1e-64
    @show @test dot(osys.systems[i,1], osys.systems[i,1]) ≈ 0.5 atol=1e-64
end