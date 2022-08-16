
#Error messages
mutable struct GCExeptionDegree <: Exception
end

mutable struct GCRExeptionDegree <: Exception
end

mutable struct GCExeptionType <: Exception
end

mutable struct CExeptionDegree <: Exception
end

mutable struct ASCExeptionDegree <: Exception
end

mutable struct GQExeptionDegree <: Exception
end


Base.showerror(io::IO, e::GCExeptionDegree) = print(
    io,
    "Only cubatures of degree 4..11 available."
)

Base.showerror(io::IO, e::GCRExeptionDegree) = print(
    io,
    "Only cubatures of degree 3..9 available."
)

Base.showerror(io::IO, e::GCExeptionType) = print(
    io,
    "Possible types :logquad, :logrect and :pol."
)

Base.showerror(io::IO, e::CExeptionDegree) = print(
    io,
    "Only cubatures of odd degree 5..21 available."
)

Base.showerror(io::IO, e::ASCExeptionDegree) = print(
    io,
    "Only cubatures of degree 3..15 available."
)

Base.showerror(io::IO, e::GQExeptionDegree) = print(
    io,
    "Only cubatures of degree 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25 and 30 available."
)

function generalizedcubature(degree::Int; type = :logrect)
    if type == :logrect
        if degree < 3 || degree > 9
            throw(GCRExeptionDegree())
        else
            return gclogrectxb[degree-2], gclogrectwb[degree-2] 
        end
    else
        if degree < 4 || degree > 11
            throw(GCExeptionDegree())
        else
            if type == :logquad
                return gclogquadxb[degree-3], gclogquadwb[degree-3] 
            else
                if type ==:pol
                    return gcpolxb[degree-3], gcpolwb[degree-3] 
                else
                    throw(GCExeptionType())
                end    
            end
        end
    end
end

function symmetriccubature(degree::Int)
    if degree < 5 || degree > 21 || iseven(degree)
        throw(CExeptionDegree())
    else
        xn = zeros(Float64, length(csplx[Int(ceil((degree - 4)/2))]), 2)
        for (i, ex) in enumerate(csplx[Int(ceil((degree - 4)/2))])
            xn[i, 1] = ex[1]
            xn[i, 2] = ex[2]
        end 
        return xn, csplw[Int(ceil((degree - 4)/2))]
    end
end

function asymmetriccubature(degree::Int)
    if degree < 3 || degree > 15
        throw(ASCExeptionDegree())
    else
        return cplxb[degree-2], cplwb[degree-2]
    end
end

function generalizedquadrature(n::Int)
    rules = [15, 20, 25, 30]
    if n >= 3 && n <= 10
        return gqlogx[n-2], gqlogw[n-2] 
    else
        if n in rules
            gqlogx[8+Int((n-10)/5)], gqlogw[8+Int((n-10)/5)] 
        else
            throw(GQxeptionDegree())
        end
    end
end
