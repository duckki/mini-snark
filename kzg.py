# The KZG polynomial commitment scheme (PCS) - complete version.
# Added the following features to the basic KZG version:
# - Multiple evaluation points using the target polynomial `t`.
# - Polynomial restriction by adding a shifted calculation as a "check-sum".
from bls12 import bls12_381_symm as bls12_381
from field import Field

# A curated selection of a finite field and subgroup generator G.
FieldElement = Field(bls12_381.order)
G = bls12_381.G
e = bls12_381.pair  # the pairing operator `e`

# Trusted setup
# - input `n`: the maximum number of polynomial coefficients allowed in the scheme.
# - input `t`: the target polynomial
# - `s` is the secret key and `alpha` is the shift value.
# - The return value `pk` is the proving key.
# - The return value `vk` is the verification key.
# - The return value `s` and `alpha` must be kept secret or discarded.
def setup( n, t ):
    s = FieldElement.random_element()
    alpha = FieldElement.random_element()
    pk = ( [G * s ** i for i in range(n)] \
         , [G * s ** i * alpha for i in range(n)] )
    vk = ( G * s, G * t(s), G * alpha )
    return (pk, vk, s, alpha)

# Commit polynomial `f`.
# - `H`: one of the sequences from `pk`.
# - returns `f0 * H0 + f1 * H1 + ... + fd * Hd` (`d` is the degree of `f`).
def commit( H, f ) -> bls12_381:
    assert len(H) >= len(f.coeff)
    return sum([f.coeff[i] * H[i] for i in range(len(f.coeff))], G*0)

# Prove `f(u) = v`.
# returns the commitment of `(f(x) - v) / (x - u)`.
def prove_eval( pk, f, u, v ) -> bls12_381:
    # the quotient polynomial `q`
    # - `q` should be a polynomial, if f(u) = v.
    X = FieldElement.X
    q = (f - v) / (X - u)
    return commit( pk[0], q )

# Verify the proof `com_h` for `f(u) = v`.
def verify_eval( vk, com_f, u, v, com_h ) -> bool:
    # Checks `(s - u) * com_h == com_f - v * G` using `e`.
    G_s = vk[0]
    return e(G_s, com_h) == e(com_f + -(v * G) + u * com_h, G)

# Prove `f(r_i) = 0` for all roots `r_i` of `t`.
def prove_roots( pk, f, t ) -> bls12_381:
    # the quotient polynomial `h`
    # - `h` should be a polynomial, if `f` share the same roots of `t`.
    h = f / t
    return commit( pk[0], h )

# Verify the proof `com_h` for `f(r_i) = 0` for all roots `r_i` of `t`.
def verify_roots( vk, com_f, com_h ) -> bool:
    # Checks: `com_f == com_h * t(s)` using `e`.
    G_ts = vk[1]
    return e( com_f, G ) == e( com_h, G_ts )

# Verify shifted value: fs(X) == f(X) * alpha
def verify_shift( vk, com_f, com_fs ) -> bool:
    # Checks: `com_fs == com_f * alpha` using `e`.
    G_alpha = vk[2]
    return e( com_fs, G ) == e( com_f, G_alpha )


# ============================================================================
# Unit test

def unit_test():
    print( "Setting up..." )
    size = 32  # the maximum size of polynomial
    k = 30  # the # of roots
    roots = [FieldElement(7) ** i for i in range(k)]
    t = FieldElement.vanishing_polynomial(roots)
    pk, vk, s, alpha = setup( size, t )
    print( "s:", s )
    print( "alpha:", alpha )
    s, alpha = None, None # discard s and alpha

    # Assume Prover and Verifier aggreed on some function `f`.
    # - Use a random polynomial `f` that shares the same roots with `t`.
    # - The quotient polynomial `h` is defined as `f / t`.
    f = (FieldElement.X - FieldElement.random_element()) * t
    h = f / t

    # 1. Prover sends commitments of f, h, fs to Verifier.
    print( "Proving..." )
    com_f = commit( pk[0], f )
    com_h = prove_roots( pk, f, t )
    com_fs = commit( pk[1], f ) # f's shift
    print( "com_f:", com_f )
    print( "com_h:", com_h )
    print( "com_fs:", com_fs )

    # 2. Verifier checks the relationship between commitments.
    assert verify_roots( vk, com_f, com_h ) \
       and verify_shift( vk, com_f, com_fs )

    # 3. Interactive Oracle Proof
    for _ in range(2):
        # 3.1. Verifier sends requests to prove `h(r) = v` at random `r`.
        r = FieldElement.random_element()
        v = h(r)  # Verifier evaluates `h(u)` with its own `h`.
        print( "r:", r )
        print( "v:", v )

        # 3.2. Prover sends the proof `pi` for `h(r) = v`.
        pi_h = prove_eval( pk, h, r, v )
        print( "pi_h:", pi_h )

        # 3.3 Verifier checks the proof `pi_h` against the commitment.
        assert verify_eval( vk, com_h, r, v, pi_h )

    # Following are extra negative tests.
    print( "Running extra tests..." )

    def verify( vk, pi ) -> bool:
        com_f, com_h, com_fs = pi
        return  verify_roots( vk, com_f, com_h ) \
            and verify_shift( vk, com_f, com_fs )

    wrong_com_f = -com_f
    assert not verify( vk, (wrong_com_f, com_h, com_fs) )

    wrong_com_h = -com_h
    assert not verify( vk, (com_f, wrong_com_h, com_fs) )

    wrong_com_f2 = -com_fs
    assert not verify( vk, (com_f, com_h, wrong_com_f2) )

    print( "Success!" )

if __name__ == "__main__":
    unit_test()
