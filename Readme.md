# mini-SNARK

mini-SNARK is an implementation of KZG Polynomial Commitment Scheme (aka Kate PCS), closely following Dan Boneh's presentation at [ZK Whiteboard Sessions](https://zkhack.dev/whiteboard/) recorded back in August 2022. Among those sessions, the module two - "Building a SNARK (Part I)" is the most relevant to this implementation.

It's implemented in Python and the only external dependency is `py-ecc`, which implements the BLS12 elliptic curve group and its pairing operator. I wanted something simpler and easier to understand than the [babySNARK](https://github.com/initc3/babySNARK). And the goal is to keep the code as simple as possible.

I think I got the simplest implementation of KZG PCS, complementing Dan's presentation. But, it's still just the PCS part. The IOP part of SNARK is not there, yet. I hope to implement the rest of SNARK and then expand it to zk-SNARK in the future.

## Directory layout

* [`bls12.py`](bls12.py): A symmetric adaptation of the BLS12-381 curve.
* [`field.py`](field.py): Implementation of finite field arithmetic and polynomial.
* [`kzg.py`](kzg.py): KZG PCS implementation.

### Credit

The BLS part is beyond my head and I took the code from [the babySNARK repo](https://github.com/initc3/babySNARK). Also, I took the `field.py` partially from [the stark101 tutorial repo](https://github.com/starkware-industries/stark101) with some modifications.

## Explanation of KZG implementation

The `kzg.py` is the main KZG PCS implementation. It's only 22 lines of code, excluding comments and empty lines.

### Background
Assuming you watched Dan's presentation linked above, most of the code should be easy to follow, except for the pairing part. KZG uses elliptic curve cryptography with pairing, which is the magic behind it. The polynomials we are going to prove/verify will be over a finite field `F` of order `p` and we will use a subgroup `G` as described in Dan's talk. But, one thing to note is that those `p`, which defines `F`, and the generator of `G` are not just any field or generator. They are specific ones handpicked by cryptographers so the `G` has a superpower of being on an elliptic curve (thus, one can't derive `x` from `x * G`) and having a pairing property. One well-known such curve is BLS12 and it's implemented in the `py-ecc` package from the Ethereum foundation.

While Dan's presentation skips over on the pairing part, this [YouTube video](https://youtu.be/9TFEBuANioo?si=3CFyCDRLioTnpVGG) by OpenZeppelin does a good job of explaining the concept of elliptic curve pairing. In summary, the magical `G` of BLS12 has a special pairing function `e`. The `e` takes in two elements of group `G` and returns an element of another group `GT`. It's not critical to understand what `GT` is. The important thing to understand is that `e` has a special algebraic property called "bilinearity".

For any `a` and `b`,
```python
    e((a + b)*G, G) == e(a*G, G) * e(b*G, G)  # in Python
    e(G, (a + b)*G) == e(G, a*G) * e(G, b*G)
```

As corollary, the following are also true:

```python
    e(a*G, G) == e(G, G) ** a
    e(G, a*G) == e(G, G) ** a
    e(a*G, b*G) == e(G, G) ** (a*b)
    e(a*b*G, G) == e(b*G, G) ** a == e(G, G) ** (a*b)
```

That's all we need to know about `e` to understand KZG.

### The Setup Process

```python
def setup( n ):
    alpha = FieldElement.random_element()
    pp = [G * alpha ** i for i in range(n)]
    return (pp, alpha)
```

`G` is a generator of an additive group and `pp` holds a bunch of `G`'s group elements. That's exactly as explained in Dan's presentation.

### Prover

The prover has two functions: `commit` and `prove`. Their implementation follows the description of Dan's presentation.

```python
# Warning: slightly simplified for clarity.
def commit( pp, f ):
    # Compute f0 * H0 + f1 * H1 + ... + fd * Hd, where
    # - fn is the n-th coefficient of f (`f.coeff[i]`).
    # - Hn is the n-th public parameter (`pp[i]`).
    return sum([f.coeff[i] * pp[i] for i in range(len(f.coeff))])
```

```python
def prove( pp, f, u, v ):
    # the quotient polynomial `q`
    # - `q` should be a polynomial, if f(u) = v.
    X = FieldElement.X
    q = (f - v) / (X - u)

    # commit `f` and `q` as the proof
    com_f = commit( pp, f )
    com_q = commit( pp, q )
    return (com_f, com_q)
```

In `prove`, `FieldElement.X` stands for a simple polynomial `x` (aka the identify function like `id(x) = x`).

### Verifier

The verifier is really interesting. Its function body is just one line of code.

```python
def verify( pp, com_f, com_q, u, v ) -> bool:
    return e(pp[1], com_q) == e(com_f + -(v * G) + u * com_q, G)
```

The actual verification condition is `(alpha - u) * com_q == com_f - v * G`, which is described in Dan's talk. How did I end up with that implementation? Here are the derivation steps:

1. `(alpha - u) * com_q == com_f - v * G`   (premise)
2. `alpha * com_q - u * com_q == com_f - v * G`  (distributivity of `*` over `-`)
3. `alpha * com_q == com_f - v * G + u * com_q`  (add `u * com_q` on both sides of equation)
4. `e(alpha * com_q, G) == e(com_f - v * G + u * com_q, G)` (pseudo-injectivity of `e`)
5. `e(alpha * q(alpha) * G, G) == e(com_f - v * G + u * com_q, G)` (definition of `com_q`)
6. `e(alpha * G, q(alpha) * G) == e(com_f - v * G + u * com_q, G)` (bilinearity of `e`)
7. `e(H1, com_q) == e(com_f - v * G + u * com_q, G)` (definition of `setup` and `com_q`)

The step (5) uses the fact that `com_q` equals `q(alpha) * G` from the definition of `commit`.

At the end, `H1` is stored in `pp[1]` and the other variables are also directly available. Thus, the equality on step (7) can be directly computed, which is implemented in the `verify` function.

Notice that the decision procedure for (7) must be sound and complete against (1), which is the ultimate goal. That means (1) must imply (7) and (7) must imply (1) as well.

The chain of implication from (1) to (7) is fine. However, the other direction of the chain is a bit problematic. Look at the step (4). The (3) does imply (4), but the other way around is unclear. It relies on the injectivity of `e`. What we want is that `e(a, G)` = `e(c, G)` implies `a` = `c`. But, I could not find a proof of that. The `e` is required to have "non-degeneracy", which sounds somewhat similar to injectivity, but it doesn't necessarily mean injectivity. I'm guessing that, if `e(a, G)` = `e(c, G)` is true, then `a` = `c` is likely to be true with a high probability, while it may not be 100%.
