# Exact finite-field reference implementation

This standard-library-only Python package executes selected finite constructions
from the all-fields proof. It is a correctness-oriented prototype, with no speed
claim and **no generator for the decompositions promised by the 9/4 theorem**.
The existing Lean proof and its protected specification files are unchanged.
Python 3.10 or newer is required. Run from the repository root:

```sh
python3 -m reference.demo
python3 -m unittest reference.test_reference -v
```

## What is implemented

| Executable component | Proof source / correspondence |
| --- | --- |
| `BilinearScheme.verify()` | Exact coordinate identities in `Tensor.RankAtMost`; checks every coefficient rather than only input samples. |
| `multiply()` and `multiply_recursive()` | The scalar bilinear identity and square recursive block construction in [RecursiveBlockPrograms.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/Arithmetic/RecursiveBlockPrograms.lean). |
| `tensor_product()` and `tensor_power()` | Tensor powers with product coordinates reindexed as larger matrix coordinates. |
| `separation_period()`, `primitive_root()`, `fourier_filter()` | The characteristic-aware period and normalized root-of-unity filter in [Fourier.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Separation/Fourier.lean). |
| `recover_leading_coefficient()` | Nonzero-node Lagrange interpolation after removing a leading power, from [Interpolation.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Polynomial/Interpolation.lean). |
| `descend()` | The explicit basis/projection formula in [FieldDescent.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/Arithmetic/FieldDescent.lean). |

This correspondence is an implementation guide, not a formal proof that the
Python code refines the Lean definitions. Tests independently check tensor
coefficients, compare with ordinary matrix multiplication, exhaust all 256 pairs
of 2x2 matrices over F_2, exercise several characteristics, and reject malformed
certificates and invalid interpolation/Fourier hypotheses.

## Coefficients and usage

`FiniteField(p, modulus)` represents F_p[X]/(modulus). The characteristic must
be prime and the supplied monic polynomial must be irreducible; both are checked.
The default modulus `(0, 1)` represents F_p. Coefficients of the modulus are in
ascending order. Elements are canonical integer encodings of power-basis
coordinates, **not ordinary integer scalars**:

```python
from reference import FiniteField, descend, strassen_scheme

extension = FiniteField(2, (1, 1, 1))  # F_4, alpha^2 = alpha + 1
assert extension.coordinates(2) == (0, 1)  # encoding 2 means alpha
assert extension.embed(2) == 0             # the integer scalar 2 is zero

# Known seven-term Strassen decomposition, used only as a test seed.
# Rescaling introduces extension coefficients without changing its tensor.
seed = strassen_scheme(extension).rescale_terms(2, 1)
scheme = descend(seed.tensor_power(2))
assert scheme.verify()
assert (scheme.n, scheme.terms) == (4, 196)

left = [[1, 0, 1, 1], [0, 1, 1, 0], [1, 1, 0, 1], [1, 0, 0, 1]]
right = [[0, 1, 1, 0], [1, 1, 0, 1], [1, 0, 1, 1], [0, 1, 0, 1]]
print(scheme.multiply(left, right))  # exact F_2 arithmetic
```

`BilinearScheme(field, n, m, k, a, b, c)` also accepts externally supplied
coefficients for `(n x m) * (m x k)`. All coefficient vectors are row-major;
`c[q][i*k+j]` contributes to output `(i,j)`. Lean writes that coordinate as
`(j,i)`. Input dimensions and canonical element encodings are checked.
Execution first requires the full coefficient certificate to pass. Verification
is cached on an immutable scheme.

`multiply_recursive()` accepts square schemes with block size at least two and
input sizes that are powers of that block size, including size one. Other sizes
need explicit padding, which is not implemented. Direct evaluation can use a
rectangular scheme of exactly the supplied dimensions.

## The fixed overhead matters

For E=F_p[X]/(f) of degree d, the basis is `1,X,...,X^(d-1)` and the
F_p-linear projection pi takes the constant coordinate, so pi(1)=1. Each term
is replaced by d^2 terms with coefficients

```text
a'[q,i,j](x) = coord_i(a[q](x))
b'[q,i,j](y) = coord_j(b[q](y))
c'[q,i,j](z) = pi(X^i * X^j * c[q](z))
```

Fix E first, take the tensor power in E, then project. The decomposition length
is `d^2 * r^t`. Descending once and recursively tensoring that descended scheme
would instead give `(d^2 * r)^t`. Both compute the correct tensor, but the latter
repeats the overhead. The demo verifies the projected schemes at powers 0, 1,
and 2, and compares both construction orders at power 2: **196 vs 784 terms**.
Terms with zero coefficients are retained, so these are rank upper bounds, not
minimal ranks. The whole supplied extension is used; no minimal coefficient
algebra is found automatically. Scalar-operation counts are not instrumented.

## Scope and remaining work

- Only explicitly represented finite fields are supported. Arithmetic is exact;
  floating-point numbers, arbitrary fields, and algebraic closures are absent.
- The user supplies an extension large enough for the roots and interpolation
  nodes. Missing roots or too few nodes raise errors rather than inventing a
  root or reusing colliding nodes. The demo uses F_16 for period 5 and F_25 for
  the adjusted period 6 in characteristic 5.
- Interpolation recovers a polynomial coefficient. General polynomial tensor
  degenerations and their local restriction maps are not implemented.
- The full auxiliary-separation construction, determinant/sector argument,
  entropy limits, and spectral existence argument are not translated.
- Strassen and naive decompositions are input fixtures. Neither is a 9/4
  decomposition extracted from this proof. The example scaling is deliberately
  artificial; descent is demonstrated, not proposed as an improvement to them.
- No running-time, memory-efficiency, bit-complexity, or 9/4 arithmetic-operation
  bound is claimed. Dense tensor powers and exhaustive coefficient checking
  become large quickly. Use small dimensions and powers.

A future coefficient generator can supply an exact scheme to the existing
verification/execution interface. Turning the current spectral existence proof
into such a generator remains a separate mathematical and implementation task.
