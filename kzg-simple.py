# The KZG polynomial commitment scheme (PCS).
from bls12 import bls12_381_symm as bls12_381
from field import Field

# A curated selection of a finite field and subgroup generator G.
FieldElement = Field(bls12_381.order)
G = bls12_381.G
e = bls12_381.pair  # the pairing operator `e`

# Trusted setup
# - input `n`: the maximum number of polynomial coefficients allowed in the scheme.
# - The return value `pp` is a list of public parameters.
# - The return value `alpha` must be kept secret or discarded.
def setup( n ):
    alpha = FieldElement.random_element()
    pp = [G * alpha ** i for i in range(n)]
    return (pp, alpha)

# Commit polynomial `f`.
# - `H`: the sequence from `pp`.
# - returns `f0 * H0 + f1 * H1 + ... + fd * Hd` (`d` is the degree of `f`).
def commit( H, f ) -> bls12_381:
    assert len(H) >= len(f.coeff)
    return sum([f.coeff[i] * H[i] for i in range(len(f.coeff))], G*0)

# Prove `f(u) = v`.
# returns the commitment of `(f(x) - v) / (x - u)`.
def prove_eval( pp, f, u, v ) -> bls12_381:
    # the quotient polynomial `q`
    # - `q` should be a polynomial, if f(u) = v.
    X = FieldElement.X
    q = (f - v) / (X - u)
    return commit( pp, q )

# Verify the proof `com_q` for `f(u) = v`.
def verify_eval( pp, com_f, u, v, com_q ) -> bool:
    # Checks `(alpha - u) * com_q == com_f - v * G` using `e`.
    G_alpha = pp[1]
    return e(G_alpha, com_q) == e(com_f + -(v * G) + u * com_q, G)


# ============================================================================
# Unit test

def unit_test():
    print( "Setting up..." )
    size = 32  # the size of polynomial
    pp, alpha = setup( size )
    print( "alpha:", alpha )
    alpha = None # discard alpha

    # Assume Prover and Verifier aggreed on some function `f`.
    # - Use a random function `f` as an example.
    coeffs = [FieldElement.random_element() for _ in range(size)]
    f = FieldElement.Polynomial( coeffs )

    # 1. Prover sends the commitment of its `f` to Verifier.
    com_f = commit( pp, f )
    print( "com_f:", com_f )

    # 2. Interactive Oracle Proof
    for _ in range(2):
        # 2.1. Verifier sends requests to prove `f(u) = v` at random `u`.
        u = FieldElement.random_element()
        v = f(u)  # Verifier evaluates `f(u)` with its own `f`.
        print( "u:", u )
        print( "v:", v )

        # 2.2. Prover sends the proof `pi` for `f(u) = v`.
        print( "Proving..." )
        pi = prove_eval( pp, f, u, v )
        print( "pi:", pi )

        # 2.3. Verifier checks the proof `pi` against the commitment.
        print( "Verifying..." )
        assert verify_eval( pp, com_f, u, v, pi )

    # Following are extra negative tests.
    print( "Running extra tests..." )

    wrong_u = u + 1
    assert not verify_eval( pp, com_f, wrong_u, v, pi )

    wrong_v = v - 1
    assert not verify_eval( pp, com_f, u, wrong_v, pi )

    u2 = FieldElement.random_element()
    v2 = f(u2)
    assert not verify_eval( pp, com_f, u2, v2, pi )

    coeffs2 = [FieldElement.random_element() for _ in range(size)]
    wrong_pi = commit( pp, FieldElement.Polynomial( coeffs2 ) )
    assert not verify_eval( pp, com_f, u, v, wrong_pi )

    print( "Success!" )

if __name__ == "__main__":
    unit_test()
