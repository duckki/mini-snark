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
    # the operation `l * r = o`
    wit = FieldElement.random_element() # witness
    stmt = a * wit # statement
    print( "stmt:", stmt )
    print( "wit:", wit )

    # polynomials
    X = FieldElement.X
    t = (X - r)

    p_l = X * (a / r)
    p_r = X * (wit / r)
    p_o = X * (stmt / r)
    h_a = (p_l - a) / t             # first input value
    h_stmt = (p_o - stmt) / t       # output statement
    h_op = (p_l * p_r - p_o) / t    # operational relationship

    # commitments
    pi_l = kzg.commit( pk[0], p_l )
    pi_r = kzg.commit( pk[0], p_r )
    pi_o = kzg.commit( pk[0], p_o )
    pi_l2 = kzg.commit( pk[1], p_l ) # a shift of p_l
    pi_r2 = kzg.commit( pk[1], p_r ) # a shift of p_r
    pi_o2 = kzg.commit( pk[1], p_o ) # a shift of p_o
    pi_a = kzg.commit( pk[0], h_a )
    pi_stmt = kzg.commit( pk[0], h_stmt )
    pi_op = kzg.commit( pk[0], h_op )
    return (stmt, (pi_l, pi_r, pi_o, pi_l2, pi_r2, pi_o2, pi_a, pi_stmt, pi_op))
# end def prover


def verifier( vk, r, a, stmt, pi ):
    pi_l, pi_r, pi_o, pi_l2, pi_r2, pi_o2, pi_a, pi_stmt, pi_op = pi

    # polynomial restriction check
    assert kzg.verify_shift( vk[2], pi_l, pi_l2 )
    assert kzg.verify_shift( vk[2], pi_r, pi_r2 )
    assert kzg.verify_shift( vk[2], pi_o, pi_o2 )

    # input/output checks
    # assert: l'(r) == a
    assert kzg.verify_point( vk[0], pi_l, pi_a, r, a )
    # assert: o'(r) == stmt
    assert kzg.verify_point( vk[0], pi_o, pi_stmt, r, stmt )

    # operation check
    # assert: l'(X) * r'(X) - o'(X) == op'(X)
    assert e(pi_l, pi_r) == e(vk[1], pi_op) * e(pi_o, G)
# end def verifier


def main():
    print( "Field order:", FieldElement.order )

    # Fixed public configuration
    r = FieldElement(7) # generator of polynomial roots
    a = FieldElement(1) # fixed input value

    # Setup
    pk, vk = setup( r )

    # Prover
    stmt, pi = prover( pk, r, a )

    # Verifier
    verifier( vk, r, a, stmt, pi )

    print( "Success!" )


if __name__ == "__main__":
    main()
