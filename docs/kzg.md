# Explanation of KZG implementations

## The simplified version of KZG

The `kzg-simple.py` is a simplified version of KZG PCS implementation. It's only 22 lines of code, excluding comments and empty lines. It's implemented closely following [Dan Boneh's presentation](https://zkhack.dev/whiteboard/). Among those sessions, the module two - "Building a SNARK (Part I)" is the most relevant to this implementation. I think I got the simplest implementation of KZG PCS, complementing Dan's presentation.

### Background
Assuming you watched the Dan's presentation, most of the code should be easy to follow, except for the pairing part. KZG uses elliptic curve cryptography with pairing, which is the magic behind it. The polynomials we are going to prove/verify will be over a finite field `F` of order `p` and we will use a subgroup `G` as described in Dan's talk. But, one thing to note is that those `p`, which defines `F`, and the generator of `G` are not just any field or generator. They are specific ones handpicked by cryptographers so the `G` has a superpower of being on an elliptic curve (thus, one can't derive `x` from `x * G`) and having a pairing property. One well-known such curve is BLS12 and it's implemented in the `py-ecc` package.

While Dan's presentation skips over on the pairing part, this [YouTube video](https://youtu.be/9TFEBuANioo?si=3CFyCDRLioTnpVGG) by OpenZeppelin does a good job of explaining the concept of elliptic curve pairing. In summary, the magical `G` of BLS12 has a special pairing function `e`. The `e` takes in two elements of group `G` and returns an element of another group `GT`. It's not critical to understand what `GT` is. The important thing to understand is that `e` has a special algebraic property called "bilinearity".

For any `a` and `b`,
```python
    e((a + b)*G, G) == e(a*G, G) * e(b*G, G)  # expressed in Python
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

The prover has two functions: `commit` and `prove_eval`. Their implementation follows the description of Dan's presentation.

```python
# Commit polynomial `f`.
# - `H`: the sequence from `pp`.
# - returns `f0 * H0 + f1 * H1 + ... + fd * Hd` (`d` is the degree of `f`).
def commit( H, f ):
    return sum([f.coeff[i] * H[i] for i in range(len(f.coeff))])
```

```python
# Prove `f(u) = v`.
# returns the commitment of `(f(x) - v) / (x - u)`.
def prove_eval( pp, f, u, v ):
    # the quotient polynomial `q`
    # - `q` should be a polynomial, if f(u) = v.
    X = FieldElement.X
    q = (f - v) / (X - u)
    return commit( pp, q )
```

In `prove_eval`, `FieldElement.X` stands for a simple polynomial `x` (aka the identify function like `id(x) = x`).

### Verifier

The verifier is really interesting. Its function body is just one line of code.

```python
# Verify the proof `com_q` for `f(u) = v`.
def verify_eval( pp, com_f, com_q, u, v ) -> bool:
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

At the end, `H1` is stored in `pp[1]` and the other variables are also directly available. Thus, the equality on step (7) can be directly computed, which is implemented in the `verify_eval` function.

Notice that the decision procedure for (7) must be sound and complete against (1), which is the ultimate goal. That means (1) must imply (7) and (7) must imply (1) as well.

The chain of implication from (1) to (7) is fine. However, the other direction of the chain is a bit problematic. Look at the step (4). The (3) does imply (4), but the other way around is unclear. It relies on the injectivity of `e`. What we want is that `e(a, G)` = `e(c, G)` implies `a` = `c`. But, I could not find a proof of that. The `e` is required to have "non-degeneracy", which sounds somewhat similar to injectivity, but it doesn't necessarily mean injectivity. I'm guessing that, if `e(a, G)` = `e(c, G)` is true, then `a` = `c` is likely to be true with a high probability, while it may not be 100%.


## Complete Version of KZG

The `kzg.py` is the reusable version of KZG PCS. It's following a paper by Maksym Petkus, ["Why and how zk-SNARK works" (2019)](https://arxiv.org/abs/1906.07221). To follow the paper's notation, I renamed the secret key to `s`, from `alpha` and the new `alpha` represents the random shift value.

Following features were added on top of the simplified version of KZG:

- Multiple evaluation points using the target polynomial `t`.
- Polynomial restriction by adding a shifted calculation as a "check sum".

### Trusted Setup

```python
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
    vk = ( G * t(s), G * alpha )
    return (pk, vk, s, alpha)
```

The setup now returns proving key and verification key separately.
The proving key contains two series of encrypted keys: one for secret group elements and one for the same elements shifted by `alpha`. The verification key contains the values of `G * t(s)` and `G * alpha`.

### Prover

The `commit` and `prove_eval` functions are the same as before.

```python
# Commit polynomial `f`.
# - `H`: one of the sequences from `pk`.
# - returns `f0 * H0 + f1 * H1 + ... + fd * Hd` (`d` is the degree of `f`).
def commit( H, f ):
    return sum([f.coeff[i] * H[i] for i in range(len(f.coeff))])

# Prove `f(u) = v`.
# returns the commitment of `(f(x) - v) / (x - u)`.
def prove_eval( pk, f, u, v ):
    # the quotient polynomial `q`
    # - `q` should be a polynomial, if f(u) = v.
    X = FieldElement.X
    q = (f - v) / (X - u)
    return commit( pk[0], q )
```

The `prove_roots` function takes the target polynomial `t`. Instead of checking `f(u) = v`, the function now proves that `f` shares the same roots of `t`'s.

```python
# Prove `f(r_i) = 0` for all roots `r_i` of `t`.
def prove_roots( pk, f, t ):
    # the quotient polynomial `h`
    # - `h` should be a polynomial, if `f` share the same roots of `t`.
    h = f / t
    return commit( pk[0], h )
```


### Verifier

The `verify_eval` function is the same as before.

```python
# Verify the proof `com_h` for `f(u) = v`.
def verify_eval( vk, com_f, u, v, com_h ) -> bool:
    # Checks `(s - u) * com_h == com_f - v * G` using `e`.
    G_s = vk[0]
    return e(G_s, com_h) == e(com_f + -(v * G) + u * com_h, G)
```

This verifer does not check the `f`'s evaluation. Instead, it checks `h`'s evaluations. The benefit is that the polynomial `h` is much simpler than `f`. Thus, the verifier can be simpler and more efficient.

The new functions `verify_roots` checks the relationship between the commitments of `f` and `h`.

```python
# Verify the proof `com_h` for `f(r_i) = 0` for all roots `r_i` of `t`.
def verify_roots( vk, com_f, com_h ) -> bool:
    # Checks: `com_f == com_h * t(s)` using `e`.
    G_ts = vk[1]
    return e( com_f, G ) == e( com_h, G_ts )
```

The verifier also checks that prover's evaluation was restricted to the encrypted `s` values by verifying the cryptographic check sum: `com_fs == com_f * alpha`.

```python
# Verify shifted value: fs(X) == f(X) * alpha
def verify_shift( vk, com_f, com_fs ) -> bool:
    # Checks: `com_fs == com_f * alpha` using `e`.
    G_alpha = vk[2]
    return e( com_fs, G ) == e( com_f, G_alpha )
```

Both `verify_roots` and `verify_shift` can only be determined using the pairing function `e`.

Deriving the `verify_roots`'s computation:

1. `com_f == com_h * t(s)`
2. `e(com_f, G) == e(t(s) * com_h, G)`   -- injectivity of `e`
3. `e(com_f, G) == e(t(s) * h(s) * G, G)`  -- definition of `com_h`
4. `e(com_f, G) == e(h(s) * G, t(s) * G)`  -- bilinearity of `e`
4. `e(com_f, G) == e(com_h, t(s) * G)`  -- definition of `com_h`
5. `e(com_f, G) == e(com_h, vk[1])`  -- definition of `vk[1]`

Deriving the `verify_shift`'s computation:

1. `com_fs == com_f * alpha`
2. `e(com_fs, G) == e(com_f * alpha, G)`   -- injectivity of `e`
3. `e(com_fs, G) == e(com_f, G * alpha)`  -- bilinearity of `e`
4. `e(com_fs, G) == e(com_f, vk[2])`  -- definition of `vk[2]`
