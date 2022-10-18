module QuadRules

using Base.Threads 
using FastGaussQuadrature
using LinearAlgebra
using SpecialPolynomials

include("generalizedquadrature/nestedchebyshevapprox.jl")
include("generalizedquadrature/generalizedquadrature.jl")
include("generalizedquadrature/orthonormal.jl")
include("generalizedquadrature/errorcorrection.jl")

include("generalizednodeelimination/generalizednodeeliminationquad.jl")
include("generalizednodeelimination/generalizednodeeliminationrect.jl")

include("nodelimination/nodelimination.jl")
include("nodelimination/contnodelimination.jl")

include("symmetricquadrature/symquad.jl")

include("quadratures/quadrature/gqlog.jl")
include("quadratures/cubature/cpl.jl")
include("quadratures/cubature/cspl.jl")
include("quadratures/cubature/gclogquad.jl")
include("quadratures/cubature/gclogrect.jl")
include("quadratures/cubature/gcpol.jl")
include("quadratures/getrules.jl")

include("utils/polynomes.jl")
include("utils/tensorrule.jl")

export nonsymmetricquadquad
export nonsymmetricquadrect

export correctlog
export nestedquadrature
export nestedsystem
export gramschmidt

export nonsymmetricquad

export symquadratur
export getcombinations

export generalizedcubature
export symmetriccubature
export asymmetriccubature
export generalizedquadrature

export chebychev
export dchebychev
export schebychev
export dschebychev
export intschebychev
export legendre
export dlegendre

export tensorrule

end
