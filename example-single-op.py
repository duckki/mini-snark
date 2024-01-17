# A single operation zk-SNARK example, albeit incomplete.
# The design of computation/operation is based on the chapters 4.1 through 4.4
# from "Why and how zk-SNARK works" (https://arxiv.org/abs/1906.07221).
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
    pi_l = kzg.commit( pk[0], p_l )
    pi_r = kzg.commit( pk[0], p_r )
    pi_o = kzg.commit( pk[0], p_o )
    pi_l2 = kzg.commit( pk[1], p_l ) # a shift of p_l
    pi_r2 = kzg.commit( pk[1], p_r ) # a shift of p_r
    pi_o2 = kzg.commit( pk[1], p_o ) # a shift of p_o
    pi_op = kzg.commit( pk[0], h_op )
    return (stmt, (pi_l, pi_r, pi_o, pi_l2, pi_r2, pi_o2, pi_op))
# end def prover


def verifier( vk, pi ):
    pi_l, pi_r, pi_o, pi_l2, pi_r2, pi_o2, pi_op = pi

    # polynomial restriction check
    assert kzg.verify_shift( vk[2], pi_l, pi_l2 )
    assert kzg.verify_shift( vk[2], pi_r, pi_r2 )
    assert kzg.verify_shift( vk[2], pi_o, pi_o2 )

    # operation check
    assert e(pi_l, pi_r) == e(vk[1], pi_op) * e(pi_o, G)
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
    verifier( vk, pi )

    print( "Success!" )


if __name__ == "__main__":
    main()
