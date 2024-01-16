# A simple finite field library.
# - Field class defines a finite field of an arbitrary order.
# - FieldElement class defines an element of a finite field.
# - Polynomial class defines a polynomial over a finite field.
# - Adopted from the stark101 implementation.

# ============================================================================
# imports

from random import randint
try:
     from itertools import zip_longest
except ImportError:
     from itertools import izip_longest as zip_longest


# ============================================================================
# class Field: a finite field of an arbitrary order

class Field:
    # members
    # order: int
    # Polynomial: a shorthand constructor for Polynomial

    def __init__(self, order):
        self.order = order
        self.Polynomial = self.make_Polynomial_ctor()

    def zero(self):
        """
        Obtains the zero element of the field.
        """
        return FieldElement(self, 0)

    def one(self):
        """
        Obtains the unit element of the field.
        """
        return FieldElement(self, 1)

    def random_element(self):
        fe = FieldElement(self, randint(0, self.order - 1))
        return fe

    # A shorthand constructor to allow `<field-value>(val)`,
    # instead of calling `FieldElement(<field-value>, val)`.
    def __call__(self, val):
        if isinstance(val, int):
            return FieldElement(self, val)
        assert isinstance(val, FieldElement), f'Type mismatch: FieldElement expected, but got {type(val)}.'
        return val

    # Internal utility method used to initialize the Polynomial member.
    def make_Polynomial_ctor(self):
        """
        Returns a shorthand constructor for Polynomial
        """
        field = self

        def ctor( *args ):
            return Polynomial(field, *args )

        return ctor

    # Constructs a monomial.
    def monomial(self, degree: int, coefficient: 'FieldElement'):
        """
        Constructs the monomial `coefficient * x**degree`.
        """
        return Polynomial(self, [self.zero()] * degree + [coefficient])

    # Constructs a specific polynomial `x`.
    @property
    def X(self):
        """
        Constructs the polynomial `x`.
        """
        return self.monomial(1, self.one())

    def vanishing_polynomial( self, roots ):
        """
        returns: p(X) = (X-r1)*(X-r2)*...*(X-rn)
        """
        X = self.X
        p = self.monomial(0, self.one()) # p(x) = 1
        for r in roots:
            p *= (X - r)
        return p

# end class Field


# ============================================================================
# class FieldElement

class FieldElement:
    def __init__(self, field: Field, val: int):
        self.field = field
        self.val = val % field.order

    @property
    def order(self):
        return self.field.order

    def __repr__(self):
        return repr(self.val)

    def __hash__(self):
        return hash(self.val)

    def __eq__(self, other):
        if isinstance(other, int):
            other = FieldElement(self.field, other)
        return isinstance(other, FieldElement) and self.val == other.val

    def __neg__(self):
        return self.field.zero() - self

    def __add__(self, other):
        other = self.field(other)
        return FieldElement(self.field, (self.val + other.val) % self.order)

    __radd__ = __add__

    def __sub__(self, other):
        other = self.field(other)
        return FieldElement(self.field, (self.val - other.val) % self.order)

    def __rsub__(self, other):
        return -(self - other)

    def __mul__(self, other):
        try:
            other = self.field(other)
        except AssertionError:
            return NotImplemented
        return FieldElement(self.field, (self.val * other.val) % self.order)

    __rmul__ = __mul__

    def inverse(self):
        t, new_t = 0, 1
        r, new_r = self.order, self.val
        while new_r != 0:
            quotient = r // new_r
            t, new_t = new_t, (t - (quotient * new_t))
            r, new_r = new_r, r - quotient * new_r
        assert r == 1
        return FieldElement(self.field, t)

    def __truediv__(self, other):
        other = self.field(other)
        return self * other.inverse()

    def __pow__(self, n):
        assert n >= 0
        cur_pow = self
        res = FieldElement(self.field, 1)
        while n > 0:
            if n % 2 != 0:
                res *= cur_pow
            n = n // 2
            cur_pow *= cur_pow
        return res

# end class FieldElement


# ============================================================================
# class Polynomial

def trim_trailing_zeros(l: list):
    """
    Removes zeros from the end of a list.
    """
    if len(l) == 0: return l
    i = len(l) - 1
    while i >= 0 and l[i] == 0:
        i -= 1
    return l[:i+1]

class Polynomial:
    # coefficients: a list of FieldElement values
    def __init__(self, field, coefficients, var='x'):
        # Internally storing the coefficients in self.coeff, least-significant (i.e. free term)
        # first, so $9 - 3x^2 + 19x^5$ is represented internally by the list  [9, 0, -3, 0, 0, 19].
        # Note that coefficients is copied, so the caller may freely modify the given argument.
        self.field = field
        self.coeff = trim_trailing_zeros(coefficients)
        self.var = var

    def __repr__(self):
        result = ""
        for i, coeff in enumerate(self.coeff[::-1]):
            if coeff == 0:
                continue
            if len(result) > 0:
                result += " + "
            deg = len(self.coeff) - i - 1
            if deg == 0 or coeff != 1: # if coffecient is not elided
                result += str(coeff)
            if deg == 1:
                result += self.var
            elif deg > 1:
                result += self.var + "^" + str(deg)
            # otherwise, the variable is elided

        if len(result) == 0:
            result = "0" # zero polynomial
        return result

    def typecast(self, other):
        """
        Constructs a Polynomial from `FieldElement` or `int`.
        """
        if isinstance(other, int):
            other = FieldElement(self.field, other)
        if isinstance(other, FieldElement):
            other = Polynomial(self.field, [other])
        assert isinstance(other, Polynomial), f'Type mismatch: Polynomial and {type(other)}.'
        return other

    def __eq__(self, other):
        other = self.typecast(other)
        return self.coeff == other.coeff

    def __neg__(self):
        return Polynomial(self.field, [-c for c in self.coeff])

    def __add__(self, other):
        other = self.typecast(other)
        zip_poly = zip_longest(self.coeff, other.coeff, fillvalue=self.field.zero())
        return Polynomial(self.field, [l + r for (l, r) in zip_poly])

    __radd__ = __add__  # To support <int> + <Polynomial> (as in `1 + x + x**2`).

    def __sub__(self, other):
        other = self.typecast(other)
        zip_poly = zip_longest(self.coeff, other.coeff, fillvalue=self.field.zero())
        return Polynomial(self.field, [l - r for (l, r) in zip_poly])

    def __rsub__(self, other):  # To support <int> - <Polynomial> (as in `1 - x + x**2`).
        return -(self - other)

    def scalar_mul(self, scalar):
        """
        Multiplies polynomial by a scalar factor
        """
        return Polynomial(self.field, [c * scalar for c in self.coeff])

    def __mul__(self, other):
        other = self.typecast(other)
        res = [0] * (self.degree() + other.degree() + 1)
        # perform integer operations for performance
        pol1, pol2 = ([x.val for x in p.coeff] for p in (self, other))
        for i, c1 in enumerate(pol1):
            for j, c2 in enumerate(pol2):
                res[i + j] += c1 * c2
        res = [FieldElement(self.field, x) for x in res]
        return Polynomial(self.field, res)

    def __pow__(self, other):
        assert other >= 0
        # Calculates self**other using repeated squaring.
        res = Polynomial(self.field, [FieldElement(self.field,1)])
        cur = self
        while True:
            if other % 2 != 0:
                res *= cur
            other >>= 1
            if other == 0:
                break
            cur = cur * cur
        return res

    def __divmod__(self, other):
        """
        Returns q, r the quotient and remainder polynomials respectively, such that
        f = q * g + r, where deg(r) < deg(g).
        * Assert that g is not the zero polynomial.
        """
        other = self.typecast(other)
        pol2 = trim_trailing_zeros(other.coeff)
        assert pol2, 'Dividing by zero polynomial.'
        pol1 = trim_trailing_zeros(self.coeff)
        if not pol1:
            return [], []
        rem = pol1
        deg_dif = len(rem) - len(pol2)
        quotient = [self.field.zero()] * (deg_dif + 1)
        g_msc_inv = pol2[-1].inverse()
        while deg_dif >= 0:
            tmp = rem[-1] * g_msc_inv
            quotient[deg_dif] = quotient[deg_dif] + tmp
            last_non_zero = deg_dif - 1
            for i, coef in enumerate(pol2, deg_dif):
                rem[i] = rem[i] - (tmp * coef)
                if rem[i] != self.field.zero():
                    last_non_zero = i
            # Eliminate trailing zeroes (i.e. make r end with its last non-zero coefficient).
            rem = rem[:last_non_zero + 1]
            deg_dif = len(rem) - len(pol2)
        return Polynomial(self.field, trim_trailing_zeros(quotient)) \
             , Polynomial(self.field, rem)

    def __truediv__(self, other):
        div, mod = divmod(self, other)
        assert mod == 0, 'Polynomials are not divisible.'
        return div

    def __mod__(self, other):
        return divmod(self, other)[1]

    def degree(self):
        """
        The polynomials are represented by a list so the degree is the length of the list minus the
        number of trailing zeros (if they exist) minus 1.
        This implies that the degree of the zero polynomial will be -1.
        """
        return len(trim_trailing_zeros(self.coeff)) - 1

    def get_nth_degree_coefficient(self, n):
        """
        Returns the coefficient of x**n
        """
        if n > self.degree():
            return self.field.zero()
        else:
            return self.coeff[n]

    def compose(self, other):
        """
        Composes this polynomial with `other`.
        Example:
        >>> f = X**2 + X
        >>> g = X + 1
        >>> f.compose(g) == (2 + 3*X + X**2)
        True
        """
        other = self.typecast(other)
        res = Polynomial(self.field, [])
        for coef in self.coeff[::-1]:
            res = (res * other) + Polynomial(self.field, [coef])
        return res

    def eval(self, point):
        """
        Evaluates the polynomial at the given point using Horner evaluation.
        """
        # Using int operations within the loop to speed up.
        point = self.field(point).val
        val = 0
        for coef in self.coeff[::-1]:
            val = (val * point + coef.val) % self.field.order
        return FieldElement(self.field, val)

    def __call__(self, other):
        """
        If `other` is an int or a FieldElement, evaluates the polynomial on `other` (in the field).
        If `other` is a polynomial, composes self with `other` as self(other(x)).
        """
        if isinstance(other, (int)):
            other = FieldElement(self.field, other)
        if isinstance(other, FieldElement):
            return self.eval(other)
        if isinstance(other, Polynomial):
            return self.compose(other)
        raise NotImplementedError(str(type(other)))

# end class Polynomial


# ============================================================================
# Lagrange Polynomial Interpolation

def prod(values):
    """
    Computes a product.
    """
    len_values = len(values)
    if len_values == 0:
        return 1
    if len_values == 1:
        return values[0]
    return prod(values[:len_values // 2]) * prod(values[len_values // 2:])

def calculate_lagrange_polynomials(x_values):
    """
    Given the x_values for evaluating some polynomials, it computes part of the lagrange polynomials
    required to interpolate a polynomial over this domain.
    """
    assert len(x_values) > 0
    field = x_values[0].field
    lagrange_polynomials = []
    monomials = [field.monomial(1, field.one()) -
                 field.monomial(0, x) for x in x_values]
    numerator = prod(monomials)
    for j in range(len(x_values)):
        # In the denominator, we have:
        # (x_j-x_0)(x_j-x_1)...(x_j-x_{j-1})(x_j-x_{j+1})...(x_j-x_{len(X)-1})
        denominator = prod([x_values[j] - x for i, x in enumerate(x_values) if i != j])
        # Numerator is a bit more complicated, since we need to compute a poly multiplication here.
        # Similarly to the denominator, we have:
        # (x-x_0)(x-x_1)...(x-x_{j-1})(x-x_{j+1})...(x-x_{len(X)-1})
        cur_poly, _ = divmod(numerator, monomials[j].scalar_mul(denominator))
        lagrange_polynomials.append(cur_poly)
    return lagrange_polynomials

def interpolate_poly_lagrange(y_values, lagrange_polynomials):
    """
    :param y_values: y coordinates of the points.
    :param lagrange_polynomials: the polynomials obtained from calculate_lagrange_polynomials.
    :return: the interpolated poly/
    """
    assert len(y_values) > 0
    field = y_values[0].field
    poly = field.Polynomial([])
    for j, y_value in enumerate(y_values):
        poly += lagrange_polynomials[j].scalar_mul(y_value)
    return poly


def interpolate_poly(x_values, y_values):
    """
    Returns a polynomial of degree < len(x_values) that evaluates to y_values[i] on x_values[i] for
    all i.
    """
    assert len(x_values) == len(y_values)
    assert all(isinstance(val, FieldElement) for val in x_values),\
        'Not all x_values are FieldElement'
    lp = calculate_lagrange_polynomials(x_values)
    assert all(isinstance(val, FieldElement) for val in y_values),\
        'Not all y_values are FieldElement'
    return interpolate_poly_lagrange(y_values, lp)


# ============================================================================
# Unit test

def unit_test():
    FF = Field(7)

    print( FF.zero(), FF.one() )
    print( FF(6) + FF(3) )
    print( FF(6) + 3 )
    print( FF(6) - 1 )
    print( FF(6) * FF(2) )

    print( FF.X )
    print( FF.monomial(2, 2) )
    print( -FF.monomial(2, 2) )
    print( FF.monomial(2, 2) + FF.X + 1 )
    print( FF.monomial(2, 2) - FF.X - 1 )
    print( 1 + FF.X )
    print( 1 - FF.X )

    print( FF.X.scalar_mul(3) )
    print( FF.X * FF.X )
    print( divmod(FF.X + 2, FF.X) )

    print( FF.X ** 2 - 1 )
    print( divmod(FF.X ** 2 - 1, FF.X - 1) )
    print( (FF.X ** 2 - 1) / (FF.X - 1) )
    print( FF.monomial(2, 1)(FF.X + 1) )

    X = FF.X
    f = Polynomial(FF, [0, 0, 1])
    v = FF(1)
    u = FF(1)
    print( (f - v) / (X - u) )

if __name__ == '__main__':
    unit_test()
