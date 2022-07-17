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