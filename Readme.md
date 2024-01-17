# mini-SNARK

mini-SNARK is a simple implementation of zk-SNARK protocol. It's implemented in Python and the only external dependency is `py-ecc` from the Ethereum foundation. I wanted something simpler and easier to understand than the [babySNARK](https://github.com/initc3/babySNARK). And the goal is to keep the code as simple as possible.

## Repo Layout

### Math libraries
* [`field.py`](field.py): Implementation of finite field arithmetic and polynomial.
* [`bls12.py`](bls12.py): A symmetric adaptation of the BLS12-381 curve.

#### Credit
The BLS12 code is from [the babySNARK repo](https://github.com/initc3/babySNARK). Also, I took the `field.py` partially from [the stark101 tutorial repo](https://github.com/starkware-industries/stark101) with some modifications.

### Polynomial Commitment Scheme (PCS)

KZG is one PCS that is commonly used. I've implemented two versions.

* [`kzg-simple.py`](kzg-simple.py): A simplified KZG PCS implementation.
* [`kzg.py`](kzg.py): A reusable implementation of KZG PCS.

The detailed explanation of their implementation is in [`docs/kzg.md`](docs/kzg.md).

### Example Applications

Complete examples combining PCS with a simple circuit:

* [`example-single-op.py`](example-single-op.py): A single operation circuit example.
* [`example-single-op-optimized.py`](example-single-op-optimized.py): An optimized version of the example above.

The detailed explanation of their implementation is in [`docs/example.md`](docs/kzg.md).


## References

* [ZK Whiteboard lectures](https://zkhack.dev/whiteboard/) by Dan Boneh (2022 Aug)
* ["Why and how zk-SNARK works" (2019)](https://arxiv.org/abs/1906.07221) by Maksym Petkus (2019)
