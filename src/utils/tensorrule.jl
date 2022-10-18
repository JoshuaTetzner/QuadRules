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

function tensorrule3D(
    x::Vector,
    w::Vector,
    x1::Vector,
    w1::Vector, 
    x2::Vector,
    w2::Vector, 
    dim::Int
) 
    len = length(x) * length(x1) * length(x2)
    weights = ones(Float64, len)
    nodes = zeros(Float64, len, dim)
    index = 1
    for (i, ex) in enumerate(x)
        for (i1, ex1) in enumerate(x1)
            for (i2, ex2) in enumerate(x2)
                weights[index] = w[i] * w1[i1] * w2[i2]
                nodes[index, 1] = ex
                nodes[index, 2] = ex1
                nodes[index, 3] = ex2
                index += 1
            end
        end
    end
   
    return nodes, weights
end

function tensorrule3DX(
    x::Matrix,
    w::Vector,
    x1::Vector,
    w1::Vector, 
    dim::Int
) 
    len = length(w) * length(w1)
    weights = ones(Float64, len)
    nodes = zeros(Float64, len, dim)
    index = 1
    for i in eachindex(w)
        for (i1, ex1) in enumerate(x1)
            weights[index] = w[i] * w1[i1]
            nodes[index, 1] = x[i, 1]
            nodes[index, 2] = x[i, 2]
            nodes[index, 3] = ex1
            index += 1
            
        end
    end
   
    return nodes, weights
end