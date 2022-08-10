function tensorrule(
    x::Vector,
    w::Vector,
    x1::Vector,
    w1::Vector, 
    dim::Int
) 
    len = length(x)*length(x1)
    weights = ones(Float64, len)
    nodes = zeros(Float64, len, dim)
    index = 1
    for (i, ex) in enumerate(x)
        for (i1, ex1) in enumerate(x1)
            weights[index] = w1[i1]*w[i]
            nodes[index, 1] = ex
            nodes[index, 2] = ex1
            index += 1
        end
    end
   
    return nodes, weights
end
