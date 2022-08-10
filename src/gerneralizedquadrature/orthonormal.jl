using LinearAlgebra
using Base.Threads

function dotproduct(veca, vecb)
    lena = length(veca)
    lenb = length(vecb)
    if lena != lenb
        a = zeros(max(lena, lenb))
        b = zeros(max(lena, lenb))
        a[1:lena] = veca
        b[1:lenb] = vecb
        return dot(a, b)
    else 
        return dot(veca, vecb)
    end
end

function substract(veca, vecb)
    lena = length(veca)
    lenb = length(vecb)
    if lena != lenb
        a = zeros(max(lena, lenb))
        b = zeros(max(lena, lenb))
        a[1:lena] = veca
        b[1:lenb] = vecb
        return a .- b
    else 
        return veca .- vecb
    end
end

function gramschmidt(sys::NestedSystem)

    for (i, fs) in enumerate(sys.systems)
        normval = sum([dotproduct(seg, seg) * 
            (sys.segments[nseg+1] - sys.segments[nseg])  for (nseg, seg) in enumerate(eachcol(fs))])

        @threads for ns = 1:size(fs)[2]
            sys.systems[i][:, ns] = sys.systems[i][:, ns] ./ big(sqrt(normval))
        end

        @threads for j = (i+1):length(sys.systems)
            dotval = sum([dotproduct(
                sys.systems[i][:, l],
                sys.systems[j][:, l]
            ) * (sys.segments[l+1] - sys.segments[l])  for l = 1:size(fs)[2]])

            @threads for ns = 1:size(fs)[2]
                sys.systems[j][:, ns] = substract(
                    sys.systems[j][:, ns], 
                    sys.systems[i][:, ns] .* dotval
                )
            end  
        end
    end
    sys.intpl .= 0
    for (i, eachsys) in  enumerate(sys.systems)
        for (ns, seg) in enumerate(eachcol(eachsys))
            sys.intpl[i] += sum([
                coeff * 
                intschebychev(nc-1, sys.segments[ns], sys.segments[ns+1]) 
                for (nc, coeff) in enumerate(seg)
            ]) 
        end
    end
    
    return sys
end
