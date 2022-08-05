nXiao = [6, 7, 10, 12, 16, 18, 22, 27]
print("Oder \t")
print("nmin \t")
print("n \t")
print("nXiao \t")
println("prod")
for i = 1:length(clogx)
    Nd = 0
    if iseven(i+3)
        k=(i+3)/2
        Nd = (k+1) * (k+2) / 2
    else 
        k=(i+2)/2
        Nd = (k+1) * (k+2) / 2 + floor((k+1)/2)
    end
    print(i+3)
    print("\t")
    print(Nd)
    print("\t")
    print(length(clogw[i]))
    print("\t")
    print(nXiao[i])
    print("\t")
    println((floor((i+5)/2)^2))
end
##
using Plots
using FastGaussQuadrature

x_val = 1
y_val = 1

f(x, y) = x^2*y^4*log(x+1) + x^4*log(x+1) + y^8

#Tensorrule
asgg = []
nasgg = []
for i in eachindex(clogw)
    nodes = Float64.(clogx[i])
    weights = Float64.(clogw[i])
    val = sum([f(nodes[i,2], nodes[i,1])*weights[i] for i in eachindex(weights)])
    println(val)
    push!(asgg, val)
    push!(nasgg, length(weights))
end

#nonsymmetric(Tetzner)
pgg = []
npgg = []

for i = 1:8
    xa = Float64.(logquadx[i]).*2 .-1
    wa =  Float64.(logquadw[i]).*2
    if iseven(length(xa))
        xb, wb = gausslegendre(Int(length(xa)/2))
        nodes, weights = tensorrule(xa, wa, xb, wb, 2)
        val = sum([f(nodes[i,1], nodes[i,2])*weights[i] for i in eachindex(weights)])
        println(val)
        push!(pgg, val)
        push!(npgg, length(weights))
    end
end


trueval = 16/225 * (15*log(2)-16)#4/135 * (34 + 5*log(8)) #4*(56 + 45* log(2))/1125

##

#Plots 
pggval = abs.((pgg .- trueval)) ./ trueval .+ eps(Float64)
asggval = abs.((asgg .- trueval)) ./ trueval .+ eps(Float64)

println(pggval)
println(asggval)


pggval = 10 .* log10.(abs.(pggval))
asggval = 10 .* log10.(abs.(asggval))


p2 = plot(npgg, pggval, label = "tensorproduct",  marker = :a)
plot!(nasgg, asggval, label = "eliminated",  marker = :a)
##
(pggval)
(asggval)
npgg
nasgg