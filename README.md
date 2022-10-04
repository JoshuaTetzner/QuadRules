# Generalized Quadrature and Cubature Rules
This code presents in the knowledge of the author new generalized cubature rules for different systems of functions on the unit square.

The following types of Quadrature and cubature rules can be calculated and are saved in this code:

1. Generalized Gaussian Quadrature Rules for the system $1, \ln(x), x, x \ln(x), x^2, x^2 \ln(x), \dots$ 
2. Cubature Rules calculated by the Gaussian Product Rule of two arbitrary quadratures
3. Symmetric Cubature Rules for complete systems of polynomials
4. Asymmetric Cubature Rules for complete systems of polynomials
5. Generalized Cubature Rules for Rectangular Systems of polynomials with logarithmic singularity
6. Generalized Cubature Rules for Quadratic Systems of polynomials with logarithmic singularity

## Gaussian Quadrature Rules
Gaussian quadrature is a numerical integration technique integrating polynomials up to order $2n-1$, where $n$ is the order of the quadrature rule, exactly. Generalized gaussian quadrature extends this technique to integrate not only systems of polynomials exactly, including square-root functions and functions with endpoint-singularities.

## Generalized Gaussian Quadrature Rules
The Generalized Gaussian Quadrature Rules of degree n for the system $1, \ln(x), x, x \ln(x), x^2, x^2 \ln(x), \dots, x^n$ can be obtained by 
```
x, w = generalizedquadrature(n::Int)
```

## Cubature Rules
The term cubature is used for quadrature rules in two dimensions. An one dimensional quadrature rule can be expanded in two dimensions by integrating each direction by itself. If this technique is used $n^2$ evaluations are needed. 
This number can significantly be reduced if we only consider a complete set of functions.

## Completeness of a System of Functions
A set of monomials up to order 5 contains the functions $\{1, x, x^2, x^3, x^4, x^5\}$. In two dimensions we differentiate between complete to order $p$ or incomplete to order $p$.

1. A complete system of polynomials up to order $p$ contains all polynomials $f(x,y) = x^i y^j$, where $i+j \leq 5$
2. An incomplete system of polynomials up to order $p$ contains all polynomials $f(x, y) = x^i y^j$, where $ i, j \leq 5$

A Gaussian product rule of two one dimensional quadrature rules of order $p$ integrates an incomplete set of functions up to order $2p$ in two dimensions. If only the complete system of polynomials should be integrated, then the number of therms can significantly be reduced. This is done by either using symmetries of the integration region or simply randomly neglecting points and adjusting the remaining ones. 

## Cubature Rules
Cubature rules calculated by the Gaussian product rule of two arbitrary quadratures can be obtained by 
```
nodes, weights = tensorrule(x::Vector,w::Vector,x1::Vector,w1::Vector, 2) 
```
The variables $x, w$ are the points and weights of the first quadrature rule and $x1, w1$ the points and weights of the second rule.

## Symmetric Cubature
The symmetric cubature rules for complete systems of polynomials of degree $p$ can be obtained by
```
nodes, weights = symmetriccubature(p::Int)
```

## Asymmetric Cubature
The asymmetric cubature rules for complete systems of polynomials of degree $p$ can be obtained by
```
nodes, weights = asymmetriccubature(p::Int)
```

## Generalized Cubature Rules for Rectangular Systems
Rectangular systems are complete systems of polynomials with a logarithmic singularity in x-direction fp to degree $p$. Such a cubature of degree 3 integrates all linear combinations of $1,\ln(x),  x, x\ln(x),  x^2, x^2\ln(x),  x^3, x^3\ln(x), y,\ln(x) y,  x y, x\ln(x) y,  x^2 y, x^2\ln(x) y,\ln(x) y^2,  x y^2, x\ln(x) y^2, y, y^2, y^3$. In short all polynomials with or without singularity in x-directions, where the sum of the exponents is smaller or equal to $3$. 
These rules can be obtained by
```
nodes, weights = generalizedcubature(degree::Int; type = :logrect)
```

## Generalized Cubature Rules of Quadratic Systems
A quadratic system is a system with the same amount of terms in x and y-direction. Since every second term in x direction contains the logarithmic singularity, we integrate polynomials with twice the order in y than in x direction. A generalized Cubature rule of degree 3 of a quadratic system integrates all linear combinations of the functions $1,\ln(x),  x, x\ln(x), y,\ln(x) y,  x y, x\ln(x) y,  x^2 y, x^2\ln(x) y, y^2,\ln(x) y^2, y, y^2, y^3$.
These rules can be obtained by
```
nodes, weights = generalizedcubature(degree::Int; type = :logquad)
```