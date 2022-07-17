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

