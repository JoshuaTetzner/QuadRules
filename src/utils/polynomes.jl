using SpecialPolynomials

# Chebyshev polynomials
function chebychev(n, x::T) where {T <: AbstractFloat}
    return basis.(SpecialPolynomials.Chebyshev, n)(x)
end

function dchebychev(n, x::T) where {T <: AbstractFloat}
    if n > 0
        return n * basis.(SpecialPolynomials.ChebyshevU, n-1)(x)
    else
        return 0
    end
end

# Shifted chebychev polynomials
function schebychev(n, x::T, a::T, b::T) where {T <: AbstractFloat}
    xs = 2 * (x - a) / (b - a) - 1
    return basis.(SpecialPolynomials.Chebyshev, n)(xs)
end

function dschebychev(n, x::T, a::T, b::T) where {T <: AbstractFloat}
    xs = 2 * (x - a) / (b - a) - 1
    din = 2 / (b - a)
    if n > 0
        return din * n * basis.(SpecialPolynomials.ChebyshevU, n-1)(xs)
    else
        return 0
    end
end

function intschebychev(n, a::T, b::T) where {T <: AbstractFloat}
    dincheb = 2 / (b-a)
    if n == 1
        return 0
    else
        return 1/dincheb * ((-1)^n + 1) / T(1 - n^2)
    end
end

# Legendere polynomials
function legendre(n, x::T) where {T <: AbstractFloat}
    return basis.(Legendre, n)(x)
end

function dlegendre(n, x::T) where {T <: AbstractFloat}
    if n > 1
        return n / (x^2 - 1) * (x * legendre(n, x) - legendre(n-1,x))
    elseif n > 0
        return T(1.0)
    else
        return 0.0
    end
end
