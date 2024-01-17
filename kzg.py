# The KZG polynomial commitment scheme (PCS) - complete version.
# Added the following features to the basic KZG version:
# - Multiple evaluation points using the target polynomial `t`.
# - Polynomial restriction by adding a shifted calculation as a "check sum".
# - Improved Zero-Knowledge property by randomly shifting the commitments.
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

# input `f`: a polynomial to commit
# returns `f0 * H0 + f1 * H1 + ... + fd * Hd` (`d` is the degree of `f`).
def commit( H, f ) -> bls12_381:
    assert len(H) >= len(f.coeff)
    return sum([f.coeff[i] * H[i] for i in range(len(f.coeff))], G*0)

# Prove `f(r_i) = 0` for all roots `r_i` of `t`.
def prove( pk, f, t ) -> (bls12_381, bls12_381, bls12_381):
    # the quotient polynomial `h`
    # - `h` should be a polynomial, if `f` share the same roots of `t`.
    h = f / t

    # commit `f` and `h` as the proof
    com_f = commit( pk[0], f )
    com_h = commit( pk[0], h )
    com_f2 = commit( pk[1], f ) # f's shift

    # shift all commitments by a random `delta`
    delta = FieldElement.random_element()
    return (com_f * delta, com_h * delta, com_f2 * delta)

# Verify at one point: `f(u) = v` via h s.t h(X) = (f(X) - v) / (X - u).
def verify_point( E_s, com_f, com_h, u, v ) -> bool:
    # Goal: `(s - u) * com_h == com_f - v * G`
    return e(E_s, com_h) == e(com_f + -(v * G) + u * com_h, G)

# Verify at roots of t: `f(r) = 0` for all roots `r` of `t` via `h` where
# - h(X) = f(X) / t(X)
def verify_roots( E_ts, com_f, com_h ) -> bool:
    return e( com_f, G ) == e( E_ts, com_h )

# Verify shifted value: f2(X) == f(X) * alpha
def verify_shift( E_alpha, com_f, com_f2 ) -> bool:
    return e( com_f2, G ) == e( com_f, E_alpha )

# Verify both roots and shift at the same time.
def verify( vk, pi ) -> bool:
    com_f, com_h, com_f2 = pi
    return  verify_roots( vk[1], com_f, com_h ) \
        and verify_shift( vk[2], com_f, com_f2 )


# ============================================================================
# Unit test

def unit_test():
    print( "Setting up..." )
    size = 32  # the maximum size of polynomial
    num_roots = 30  # the # of roots
    roots = [FieldElement(3) ** i for i in range(num_roots)]
    t = FieldElement.vanishing_polynomial(roots)
    pk, vk, s, alpha = setup( size, t )
    print( "s:", s )
    print( "alpha:", alpha )
    s, alpha = None, None # discard s and alpha

    # Generate a random polynomial `f` with the same roots as `t`.
    f = (FieldElement.X - FieldElement.random_element()) * t

    print( "Proving..." )
    com_f, com_q, com_f2 = prove( pk, f, t )
    print( "com_f:", com_f )
    print( "com_q:", com_q )
    print( "com_f2:", com_f2 )

    print( "Verifying..." )
    assert verify( vk, (com_f, com_q, com_f2) )

    wrong_com_f = -com_f
    assert not verify( vk, (wrong_com_f, com_q, com_f2) )

    wrong_com_q = -com_q
    assert not verify( vk, (com_f, wrong_com_q, com_f2) )

    wrong_com_f2 = -com_f2
    assert not verify( vk, (com_f, com_q, wrong_com_f2) )

    print( "Success!" )

if __name__ == "__main__":
    unit_test()
