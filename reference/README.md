# Exact finite-field reference implementation

This standard-library-only Python package executes selected finite constructions
from the all-fields proof. It is a correctness-oriented prototype, with no speed
claim and **no generator for the decompositions promised by the 9/4 theorem**.
The existing Lean proof and its protected specification files are unchanged.
Python 3.10 or newer is required. Run from the repository root:

```sh
python3 -m reference.demo
python3 -m unittest discover -s reference -t . -v
```

## What is implemented

| Executable component | Proof source / correspondence |
| --- | --- |
| `BilinearScheme.verify()` | Exact coordinate identities in `Tensor.RankAtMost`; checks every coefficient rather than only input samples. |
| `multiply()` and `multiply_recursive()` | The scalar bilinear identity and square recursive block construction in [RecursiveBlockPrograms.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/Arithmetic/RecursiveBlockPrograms.lean). |
| `tensor_product()` and `tensor_power()` | Tensor powers with product coordinates reindexed as larger matrix coordinates. |
| `separation_period()`, `primitive_root()`, `fourier_filter()` | The characteristic-aware period and normalized root-of-unity filter in [Fourier.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Separation/Fourier.lean). |
| `recover_leading_coefficient()` | Nonzero-node Lagrange interpolation after removing a leading power, from [Interpolation.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Polynomial/Interpolation.lean). |
| `PolynomialDegeneration.verify()` and `tensor_power()` | Formal polynomial coefficients, low-degree vanishing, the matrix-tensor leading coefficient, and powers of [PolynomialApproximation](../lean/OAI/LinearAlgebra/MatrixMultiplication/Polynomial/ComplexPolynomialApproximation.lean). |
| `PolynomialDegeneration.recover_scheme()` | Generate an exact `BilinearScheme` by the normalized nonzero-node recovery identity in [Interpolation.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Polynomial/Interpolation.lean). |
| `descend()` | The explicit basis/projection formula in [FieldDescent.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/Arithmetic/FieldDescent.lean). |

This correspondence is an implementation guide, not a formal proof that the
Python code refines the Lean definitions. Tests independently check tensor
coefficients, compare with ordinary matrix multiplication, exhaust all 256 pairs
of 2x2 matrices over F_2, exercise several characteristics, and reject malformed
certificates and invalid interpolation/Fourier hypotheses.
The polynomial tests also check cancellation between terms, tensor-valued
higher-order noise, and false certificates whose evaluations agree over F_2.

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

## Polynomial degeneration to generated coefficients

`Polynomial(field, coefficients)` stores a formal polynomial with coefficients
in ascending degree. The encodings are the same as for scalar field elements.
It supports exact addition, subtraction, multiplication, nonnegative powers,
coefficient inspection, and evaluation. Trailing zeros are removed.

`PolynomialDegeneration(field, n, m, k, leading, a, b, c)` accepts polynomial
vectors instead of the scalar vectors of `BilinearScheme`. These specify

```text
P(t)[x,y,z] = sum_q a[q][x](t) * b[q][y](t) * c[q][z](t)
```

`verify()` sums the terms in the **formal polynomial ring** and checks every
tensor coordinate: coefficients below `leading` must be zero, and the
coefficient at `leading` must equal the matrix multiplication tensor.
Higher coefficients are allowed to contain other tensors. Vanishing is checked
after summation; individual terms may have low coefficients that cancel.
This is stronger than checking a few parameter values, including all values of
a finite field: over F_2, `t` and `t^2` evaluate identically but have different
linear coefficients.

After the certificate passes, `recover_scheme(nodes=None)` evaluates the
polynomial families at distinct nonzero nodes. For each node t it appends the
evaluated terms and scales their output vectors by `w_t / t^leading`, where
`w_t` is the Lagrange constant-recovery weight. The resulting scalar scheme
is then independently checked by `BilinearScheme.require_exact()` before it
is returned. `evaluate_at(t)` alone returns coefficient families for P(t);
they usually fail the exact matrix-tensor certificate and cannot be used to
multiply matrices. They are intermediate tensors, not the recovered target.

The demo supplies a synthetic fixture built from a known Strassen scheme.
For each left vector, it replaces coordinate j by
`t*a[q][j] + t^2*a[q][(j+1) mod (n*m)]`. This constructs
`P(t)=t*T+t^2*S`, where S cyclically permutes the left tensor coordinates.
For this matrix example S is not proportional to T: it includes nonzero
coefficients where T has zero coefficients. This is a test input for recovery,
not a proof-specific degeneration or a new matrix multiplication algorithm.

```python
from reference import FiniteField, strassen_scheme, descend
from reference.degenerations import left_perturbation_fixture

field = FiniteField(2, (1, 1, 1))
seed = strassen_scheme(field).rescale_terms(2, 1)
degeneration = left_perturbation_fixture(seed)
assert degeneration.verify()

# Use the SAME parameter t in all tensor factors; interpolate after powering.
powered = degeneration.tensor_power(2)
generated = powered.recover_scheme()
assert (powered.leading, powered.degree_bound) == (2, 4)
assert (powered.recovery_node_count, generated.terms) == (3, 147)
assert generated.verify()  # an executable 4x4 matrix scheme over F_4

base = descend(generated)
assert (base.field, base.terms) == (FiniteField(2), 588)
assert base.verify()       # an executable 4x4 matrix scheme over F_2
```

External polynomial coefficient families can be passed to the same interface;
they do not need to come from this fixture. `from_scheme(scheme, leading=d)`
provides a constant/monomial lift as a control. There is no search for a
low-rank polynomial degeneration.

### Interpolate after taking powers

Let D be the sum of the maximum degrees of the three polynomial families and
d the leading degree. For a valid certificate, each normalized coordinate
`P(t)/t^d` has degree at most `h=D-d`. Recovery therefore uses h+1 nodes by
default. This uses the normalized degree bound; the Lean power-rank theorem
states the looser `(D*m+1)*r^m` bound, while the normalized recovery identity
also justifies `(h*m+1)*r^m`.

With r terms, recovering after the m-th power gives `(h*m+1)*r^m` terms.
Recovering first and then powering gives `((h+1)*r)^m` terms. The fixture has
r=7 and h=1, so at power 2 the two paths produce **147 vs 196 terms**. Both
schemes pass exact coefficient checks and compute the same matrix product.
Following recovery with descent through a fixed degree-e extension gives the
upper bound `e^2*(h*m+1)*r^m`; the demo obtains **588 terms** at e=2 and m=2.
None of these term counts measure minimum rank or execution time.

The supplied finite field must have at least h*m+1 distinct nonzero elements.
Thus a **fixed finite field supports only a finite range of powers** when h>0.
Too few nodes raise an error. These examples do not replace the infinite
algebraic closure or establish an asymptotic bound. No extension field is
selected automatically, and degree bounds are conservative even when further
cancellation could reduce the actual tensor degree.

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
- Polynomial matrix-tensor degenerations supplied as coefficient families are
  verified and recovered. General source tensors and polynomial local
  restriction maps from those source tensors are not implemented.
- The full auxiliary-separation construction, determinant/sector argument,
  entropy limits, and spectral existence argument are not translated.
- Strassen and naive decompositions are input fixtures. Neither is a 9/4
  decomposition extracted from this proof. The example scaling is deliberately
  artificial; descent is demonstrated, not proposed as an improvement to them.
- No running-time, memory-efficiency, bit-complexity, or 9/4 arithmetic-operation
  bound is claimed. Dense tensor powers and exhaustive coefficient checking
  become large quickly. Use small dimensions and powers.

Proof-specific degeneration generators can supply polynomial families to the
verification/recovery interface. Translating the finite determinant/sector
constructions and turning the current spectral existence proof into a generator
for the 9/4 schemes remain separate mathematical and implementation tasks.
