module OscarHomotopyContinuation

import Oscar: Oscar, gens, parent, coeff, exponent_vector, monomials, MPolyIdeal
import HomotopyContinuation: HomotopyContinuation, Variable
import LinearAlgebra: LinearAlgebra

export poly_to_expr, nsolve, ndim, ndegree, nid

"""
    poly_to_expr(f)

Takes a polynomial `f` in an `Oscar.MPolyRing` and converts it into a
`HomotopyContinuation.Expression`. The `HomotopyContinuation.Variable`
objects used in the expression have the same names as the generators
in the `Oscar.MPolyRing` (including brackets as in `x[1]`).
"""
function poly_to_expr(f)
    # Get a list of HomotopyContinuation variables whose names are the
    # same as the ones in the Oscar polynomial f.
    v = map(x -> Variable(string(x)), gens(parent(f)))
    # Get an array of representations of the monomials in f in the form
    # [C, A] where C is the coefficient and A is the exponent vector with
    # respect to gens(parent(f)) which is exactly how v is ordered, too.
    poly = map(m -> [coeff(f, m), exponent_vector(m, 1)], monomials(f))
    # Make the HomotopyContinuation expression
    +([*([Rational(c), [v[i]^e for (i,e) in enumerate(a)]...]...) for (c,a) in poly ]...)
end

"""
    System(I::MPolyIdeal; args...)

Takes an `Oscar.MPolyIdeal` and turns it into a `HomotopyContinuation.System`
containing the ideal generators. It forwards all `args` to `HomotopyContinuation.System`.
"""
function System(I::MPolyIdeal; args...)
    HomotopyContinuation.System(map(f -> poly_to_expr(f), gens(I)); args...)
end

"""
    nsolve(I::MPolyIdeal; args...)

Call `HomotopyContinuation.solve` on the `HomotopyContinuation.System` derived
from `I` forwarding all `args`.
"""
function nsolve(I::MPolyIdeal; show_progress=false, args...)
    HomotopyContinuation.solve(System(I); show_progress, args...)
end

function witness_set(I::MPolyIdeal; show_progress=false, args...)
    HomotopyContinuation.witness_set(System(I); show_progress, args...)
end

"""
    ndim(I::MPolyIdeal)

Compute the dimension of `I` numerically.
"""
function ndim(I::MPolyIdeal)
    # This is provided by HomotopyContinuation.jl and computes the
    # dimension based on the Jacobian rank at a random point.
    F = System(I)
    HomotopyContinuation.nvariables(F) - LinearAlgebra.rank(HomotopyContinuation.fixed(F))
end

"""
    ndegree(I::MPolyIdeal; args...)

Compute the degree of `I` numerically, forwarding all `args` to
`HomotopyContinuation.witness_set` if the dimension is positive
and to `HomotopyContinuation.solve` if the dimension is zero.
"""
function ndegree(I::MPolyIdeal; show_progress=false, args...)
    if ndim(I) == 0
        F = System(I)
        sols = HomotopyContinuation.solve(F; show_progress, args...)
	cert = HomotopyContinuation.certify(F, sols; show_progress)
        HomotopyContinuation.ncertified(cert)
    else
        HomotopyContinuation.degree(witness_set(I; show_progress, args...))
    end
end

"""
    nid(I::MPolyIdeal; args...)

Call `HomotopyContinuation.nid` on the `HomotopyContinuation.System` derived
from `I` forwarding all `args`.
"""
function nid(I::MPolyIdeal; show_progress=false, args...)
    HomotopyContinuation.nid(System(I); show_progress, args...)
end

# Also add methods to HomotopyContinuation in case the user has it loaded at global scope.
#HomotopyContinuation.System(I::MPolyIdeal; args...) = System(I; args...)
#HomotopyContinuation.solve(I::MPolyIdeal; args...)  = nsolve(I; args...)
#HomotopyContinuation.nid(I::MPolyIdeal; args...)    = nid(I; args...)

# Also add some methods to Oscar.
#Oscar.nsolve(I::MPolyIdeal; args...)  = nsolve(I; args...)
#Oscar.ndim(I::MPolyIdeal; args...)    = ndim(I; args...)
#Oscar.ndegree(I::MPolyIdeal; args...) = ndegree(I; args...)

end
