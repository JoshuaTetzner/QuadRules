#Comparisson number of points symmetric and nonsymmetric
println("Order\tnonsym\tsym")
for i = 1:length(cplw)
    print(i+2)
    print("\t")
    print(length(cplw[i]))
    if (2*length(csplw)) >= (i-1) && i >= 2
        print("\t")
        println(length(csplw[Int(floor(((i)/2)))]))
    else
        print("\t")
        println("-")
    end
end

##
for i in eachindex(logquadx)
    print(length(logquadx[i]))
    for j in eachindex(logquadx[i])
        print("&")
        print(Float64(logquadx[i][j]))
        print("&")
        print(Float64(logquadw[i][j]))
        println("\\\\")
    end
end