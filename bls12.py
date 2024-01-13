# A symmetric adaptor for the BLS12-381 curve.
# - Adopted from the babysnark implementation.

from py_ecc import optimized_bls12_381 as bls12_381
from field import FieldElement

class bls12_381_symm:
    order = bls12_381.curve_order

    def __init__(self, m1, m2):
        self.m1 = m1
        self.m2 = m2

    def in_group(self):
        return bls12_381.pairing(self.m2, bls12_381.G1) == bls12_381.pairing(
            bls12_381.G2, self.m1
        )

    def __neg__(self):
        return bls12_381_symm(bls12_381.neg(self.m1), bls12_381.neg(self.m2))

    def __add__(self, other):
        assert type(other) is bls12_381_symm
        return bls12_381_symm(bls12_381.add(self.m1, other.m1), bls12_381.add(self.m2, other.m2))

    def __mul__(self, x):
        assert type(x) in (int, FieldElement)
        if type(x) is FieldElement:
            x = x.val
        return bls12_381_symm(bls12_381.multiply(self.m1, x), bls12_381.multiply(self.m2, x))

    def __rmul__(self, x):
        return self.__mul__(x)

    def __eq__(self, other):
        return bls12_381.eq(self.m1, other.m1) and bls12_381.eq(self.m2, other.m2)

    # `pair` is the pairing operator `e`.
    @staticmethod
    def pair(x, y):
        t1 = bls12_381.pairing(x.m2, y.m1)
        # t2 = bls12_381.pairing(y.m2, x.m1)
        # assert t1 == t2 # symmetry
        return t1

    def __str__(self):
        return "(" + str(self.m1) + "," + str(self.m2) + ")"


# The generator G
bls12_381_symm.G = bls12_381_symm(bls12_381.G1, bls12_381.G2)
