# The KZG polynomial commitment scheme (PCS).
from functools import reduce
from operator import add
from bls12 import bls12_381_symm
from field import Field

# A curated selection of a finite field and subgroup generator G.
FieldElement = Field(bls12_381_symm.order)
G = bls12_381_symm.G
e = bls12_381_symm.pair  # the pairing operator `e`

# Trusted setup
# - input `n`: the maximum number of polynomial coefficients allowed in the scheme.
# - The return value `pp` is a list of public parameters.
# - The return value `alpha` must be kept secret or discarded.
def setup( n ):
    alpha = FieldElement.random_element()
    pp = [G * alpha ** i for i in range(n)]
    return (pp, alpha)

# `f`: a polynomial to commit
# returns `f(alpha) * G`
def commit( pp, f ) -> bls12_381_symm:
    assert len(pp) >= len(f.coeff)
    # Compute f0 * H0 + f1 * H1 + ... + fd * Hd, where
    # - fn is the n-th coefficient of f (`f.coeff[i]`).
    # - Hn is the n-th public parameter (`pp[i]`).
    # - Using `reduce`, instead of `sum` to add up without the start value.
    return reduce(add, [f.coeff[i] * pp[i] for i in range(len(f.coeff))])

# Prove `f(u) = v`.
def prove( pp, f, u, v ) -> (bls12_381_symm, bls12_381_symm):
    # the quotient polynomial `q`
    # - `q` should be a polynomial, if f(u) = v.
    X = FieldElement.X
    q = (f - v) / (X - u)

    # commit `f` and `q` as the proof
    com_f = commit( pp, f )
    com_q = commit( pp, q )
    return (com_f, com_q)

# Verify `(alpha - u) * com_q == com_f - v * G`.
def verify( pp, com_f, com_q, u, v ) -> bool:
    return e(pp[1], com_q) == e(com_f + -(v * G) + u * com_q, G)


# ============================================================================
# Unit test

def unit_test():
    print( "Setting up..." )
    size = 32  # the size of polynomial
    pp, alpha = setup( size )
    print( "alpha:", alpha )

    # Generate a random polynomial `f` and a random point `u`.
    coeffs = [FieldElement.random_element() for _ in range(size)]
    f = FieldElement.Polynomial( coeffs )
    u = FieldElement.random_element()
    v = f(u)

    print( "Proving..." )
    com_f, com_q = prove( pp, f, u, v )
    print( "com_f:", com_f )
    print( "com_q:", com_q )

    print( "Verifying..." )
    assert verify( pp, com_f, com_q, u, v )

    wrong_u = u + 1
    assert not verify( pp, com_f, com_q, wrong_u, v )

    wrong_v = v - 1
    assert not verify( pp, com_f, com_q, u, wrong_v )

    u2 = FieldElement.random_element()
    v2 = f(u2)
    assert not verify( pp, com_f, com_q, u2, v2 )

    print( "Success!" )

if __name__ == "__main__":
    unit_test()
