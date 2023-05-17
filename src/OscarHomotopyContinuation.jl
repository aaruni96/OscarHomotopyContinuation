module OscarHomotopyContinuation

import Oscar: Oscar, gens, parent, coeff, exponent_vector, monomials, MPolyIdeal
import HomotopyContinuation: HomotopyContinuation, Variable

export poly_to_expr, numerical_solve, nid

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
    numerical_solve(I::MPolyIdeal; args...)

Call `HomotopyContinuation.solve` on the `HomotopyContinuation.System` derived
from `I` forwarding all `args`.
"""
function numerical_solve(I::MPolyIdeal; args...)
    HomotopyContinuation.solve(System(I); args...)
end

"""
    nid(I::MPolyIdeal; args...)

Call `HomotopyContinuation.nid` on the `HomotopyContinuation.System` derived
from `I` forwarding all `args`.
"""
function nid(I::MPolyIdeal; args...)
    HomotopyContinuation.nid(System(I); args...)
end

# Also add methods to HomotopyContinuation in case the user has it
# loaded at global scope.
HomotopyContinuation.System(I::MPolyIdeal; args...) = System(I; args...)
HomotopyContinuation.solve(I::MPolyIdeal; args...) = numerical_solve(I; args...)
HomotopyContinuation.nid(I::MPolyIdeal; args...) = nid(I; args...)

# Also provide a method in Oscar.
Oscar.solve(I::MPolyIdeal; args...) = numerical_solve(I; args...)

end
