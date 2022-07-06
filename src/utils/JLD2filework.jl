using JLD2

symquad = load("symquad5.jld2")
for i = 7:2:23
    dict = load("symquad" * string(i) * ".jld2")
    merge!(symquad, dict)
end
save("symmetricquad.jld2", symquad)