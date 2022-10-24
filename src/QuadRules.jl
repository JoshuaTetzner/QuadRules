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
include("nodelimination/nodelimination3D.jl")
include("nodelimination/contnodelimination3D.jl")
include("nodelimination/initialquad.jl")
include("nodelimination/base.jl")

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
export nonsymmetricquad3D
export contnonsymmetricquad
export contnonsymmetricquad3D
export initialquad
export initialquad3D

export fmat
export fmat3D
export getpolynomes, getpolynomes_dx, getpolynomes_dy
export getpolynomes3D, getpolynomes_dx3D, getpolynomes_dy3D, getpolynomes_dz3D
export getA
export getA3D
export jacobian
export jacobian3D
export checkinterior
export checkinterior3D

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
export tensorrule3D

end
