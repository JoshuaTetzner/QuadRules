function tensorrule(
    x::Vector,
    w::Vector,
    x1::Vector,
    w1::Vector, 
    dim::Int
) 
    n = length(x)
    weights = ones(Float64, n^dim)
    nodes = zeros(Float64, n^dim, dim)
    count = 1
    count1 = 1
    for i = 1:n^dim
        nodes[i, 1] = x[count]
        nodes[i, 2] = x1[count1]
        weights[i] = w[count]*w1[count1]
        if count == n
            count = 0
            count1 += 1
        end

        count +=1
    end
   
    return nodes, weights
end
