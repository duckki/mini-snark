# Example Applications

## Single Operation Example

This example is based on the chapters 4.1 through 4.4 from ["Why and how zk-SNARK works" (2019)](https://arxiv.org/abs/1906.07221). It is designed to reuse the `kzg.py` module.

### Setup

This example is so simple that the target polynomial has only one root. Otherwise, it uses KZG's setup directly.

```python
def setup( r ):
    max_size = 2 # max size of polynomials
    t = (FieldElement.X - r) # target polynomial
    pk, vk, _, _ = kzg.setup( max_size, t )
    return (pk, vk)
```

### Prover

The prover is implemented the same as the reference paper. To make the example a bit more interesting, I fixed the left operand `l`'s value (the parameter `a`). The right operand `r` is considered the secret witness and the result `o` is considered the public statement produced from the prover.

```python
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
```

### Verifier

```python
def verifier( vk, pi ):
    pi_l, pi_r, pi_o, pi_l2, pi_r2, pi_o2, pi_op = pi

    # polynomial restriction check
    assert kzg.verify_shift( vk[2], pi_l, pi_l2 )
    assert kzg.verify_shift( vk[2], pi_r, pi_r2 )
    assert kzg.verify_shift( vk[2], pi_o, pi_o2 )

    # operation check
    assert e(pi_l, pi_r) == e(vk[1], pi_op) * e(pi_o, G)
```

### Review

It's nice that the `kzg` module is reusable at least partially. Those `kzg.commit` and `kzg.verify_shift` abstract away some math/complexity out of the prover and verifier implementation.

Notice that the `kzg.prove` and `kzg.verify` functions are not used here. It turned out that the operation check part is specific to each application. It means the setup, prover and verifier code need to be re-built for every application. Thus, it makes sense to use a compiler to generate them from a circuit definition.

By the way, this implementation of verifier actually does not check the `pi_l` and `pi_o` are correct. This will be fixed in the next optimized version below.


## Optimized Single Operation Example

Since the `a` and `stmt` values are public and we know what the polynomials for `l` and `o` should be, there is no need to check the commitments for `l` and `o`. Thus, we can drop `pi_l`, `pi_o`, `pi_l2` and `pi_o2` from the proof `pi`. Then, the verifier can be simplified as following.

```python
def verifier( vk, r, a, stmt, pi ):
    pi_r, pi_r2, pi_op = pi

    # polynomial restriction check
    assert kzg.verify_shift( vk[2], pi_r, pi_r2 )

    # operation check along with the known values
    assert e(vk[0], pi_r) ** (a/r).val == e(vk[1], pi_op) * e(vk[0], G) ** (stmt/r).val
```

The operation check needs some explanation. The original assertion was this:

```python
    assert e(pi_l, pi_r) == e(vk[1], pi_op) * e(pi_o, G)
```

Since we don't have `pi_l` and `pi_o` any longer, we need to compute them in the verifier. Let's recall how they were computed in the unoptimized prover:

```python
    p_l = X * (a / r)
    p_o = X * (stmt / r)

    # commitments
    pi_l = kzg.commit( pk[0], p_l )
    pi_o = kzg.commit( pk[0], p_o )
```

The new verifier could do the same computation:

```python
    pi_l = kzg.commit( pk[0], X * (a / r) )
    pi_o = kzg.commit( pk[0], X * (stmt / r) )
    assert e(pi_l, pi_r) == e(vk[1], pi_op) * e(pi_o, G)
```

The catch is that the verifier would now need `pk[0]` (the value of `G * s`) as its parameter unless it is not! The cool thing about pairing is that the verifier can check the assertion without `pk[0]`. Note that the commitment `pi_l` equals `G * s * (a/r)` (by the definition of `pk[0]` and `kzg.commit`). Thus, that assertion above can be re-written as following:

1. `assert e(pi_l, pi_r) == e(vk[1], pi_op) * e(pi_o, G)` (original assertion)
2. `assert e(G * s * (a/r), pi_r) == e(vk[1], pi_op) * e(pi_o, G)` (definition of `pi_l`)
3. `assert e(G * s, pi_r) ** (a/r) == e(vk[1], pi_op) * e(pi_o, G)` (bilinearity of `e`)
4. `assert e(vk[0], pi_r) ** (a/r) == e(vk[1], pi_op) * e(pi_o, G)` (definition of `vk[0]`)
5. `assert e(vk[0], pi_r) ** (a/r) == e(vk[1], pi_op) * e(vk[0], G) ** (stmt/r)` (similar steps 2-4 for `pi_o`)

Therefore, the assertion in (5) can be computed directly by the verifier without `pi_l` nor `pk[0]`.


## Missing elements

The examples lack some Zero-Knowledge aspect of `kzg.prove` (shifting the proof by a random delta), which needs to be implemented.

Also, the examples are still not quite SNARK, since the simulation of interactive proof is not implemented, yet.
