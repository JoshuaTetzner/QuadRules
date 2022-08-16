using BasisFunctions, DomainSets, GeneralizedGauss

# algorithm of Huybrechs
N = 5
cheb = ChebyshevT(N) → big(0)..big(1)        # Chebyshev polynomials scaled to [0,1]
basis = cheb ⊕ log*cheb                # space of the form p(x) + log(x)*q(x) with p and q polynomials
wk, xk = compute_gauss_rules(basis, verbose=true)