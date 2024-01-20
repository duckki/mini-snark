# A single operation zk-SNARK example optimized by removing unnecessary
# commitment and verification.
# The enforcing of the first input value and the output value as
# "statement" was added to make the example more interesting.
import kzg

FieldElement = kzg.FieldElement
G = kzg.G
e = kzg.e

def setup( r ):
    max_size = 2 # max size of polynomials
    t = (FieldElement.X - r) # target polynomial
    pk, vk, _, _ = kzg.setup( max_size, t )
    return (pk, vk)


def prover( pk, r, a ):
    # compute the witness and statement values
    wit = FieldElement.random_element() # witness
    stmt = a * wit # statement
    print( "stmt:", stmt )
    print( "wit:", wit )

    # polynomials representing the operation `l * r = o`
    X = FieldElement.X
    t = (X - r)
    p_l = X * (a / r)
    p_r = X * (wit / r)
    p_o = X * (stmt / r)
    h_op = (p_l * p_r - p_o) / t    # operational relationship

    # commitments
    pi_r = kzg.commit( pk[0], p_r )
    pi_r2 = kzg.commit( pk[1], p_r ) # a shift of p_r
    pi_op = kzg.commit( pk[0], h_op )
    return (stmt, (pi_r, pi_r2, pi_op))
# end def prover


def verifier( vk, r, a, stmt, pi ):
    pi_r, pi_r2, pi_op = pi

    # polynomial restriction check
    assert kzg.verify_shift( vk, pi_r, pi_r2 )

    # operation check along with the known values
    # Checks `p_l * p_r - p_o == h_op * t(s)` using `e`.
    G_s = vk[0]
    G_ts = vk[1]
    assert e(G_s, pi_r) ** (a/r).val == e(pi_op, G_ts) * e(G_s, G) ** (stmt/r).val
# end def verifier


def main():
    print( "Field order:", FieldElement.order )

    # Fixed public configuration
    # - r: A random generator for the polynomial roots
    # - a: One of the input operands
    r = FieldElement(7)
    a = FieldElement(11748457154244067814715434074073542859880617394130524462)

    # Setup
    pk, vk = setup( r )

    # Prover
    stmt, pi = prover( pk, r, a )

    # Verifier
    verifier( vk, r, a, stmt, pi )

    print( "Success!" )


if __name__ == "__main__":
    main()
