# Exact finite-field reference implementation

This standard-library-only Python package executes selected finite constructions
from the all-fields proof. It is a correctness-oriented prototype, with no speed
claim and **no generator for the decompositions promised by the 9/4 theorem**.
The existing Lean proof and its protected specification files are unchanged.
Python 3.10 or newer is required. Run from the repository root:

```sh
python3 -m reference.demo
python3 -m reference.proof_demo  # just the proof-derived construction pipeline
python3 -m reference.extraction_demo  # branch extraction to a small matrix scheme
python3 -m reference.catalyst_demo  # actual maps from a supplied catalyst
python3 -m reference.witnessed_demo  # finite dual -> catalyst -> matrix coefficients
python3 -m unittest discover -s reference -t . -v
```

## What is implemented

| Executable component | Proof source / correspondence |
| --- | --- |
| `FiniteTensor` and `TensorScheme` | Finite three-leg coefficient tensors, exact `RankAtMost` decompositions, and bilinear contraction for arbitrary targets in [ComplexTensor.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/Tensor/ComplexTensor.lean). |
| `BilinearScheme.verify()` | Exact coordinate identities in `Tensor.RankAtMost`; checks every coefficient rather than only input samples. |
| `multiply()` and `multiply_recursive()` | The scalar bilinear identity and square recursive block construction in [RecursiveBlockPrograms.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/Arithmetic/RecursiveBlockPrograms.lean). |
| `tensor_product()` and `tensor_power()` | Tensor powers with product coordinates reindexed as larger matrix coordinates. |
| `separation_period()`, `primitive_root()`, `fourier_filter()` | The characteristic-aware period and normalized root-of-unity filter in [Fourier.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Separation/Fourier.lean). |
| `recover_leading_coefficient()` | Nonzero-node Lagrange interpolation after removing a leading power, from [Interpolation.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Polynomial/Interpolation.lean). |
| `PolynomialDegeneration.verify()` and `tensor_power()` | Formal polynomial coefficients, low-degree vanishing, the matrix-tensor leading coefficient, and powers of [PolynomialApproximation](../lean/OAI/LinearAlgebra/MatrixMultiplication/Polynomial/ComplexPolynomialApproximation.lean). |
| `PolynomialDegeneration.recover_scheme()` | Generate an exact `BilinearScheme` by the normalized nonzero-node recovery identity in [Interpolation.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Polynomial/Interpolation.lean). |
| `TensorDegeneration` | The same formal certificate and recovery operations for any supplied finite tensor target, via `PolynomialApproximation`. The matrix API delegates certificate/recovery to it. |
| `ThreeSectorConstruction` | The concrete source `C(a,3*h+a-1)`, interval weights, three ambient branches, and diagonal polynomial restriction from [Sector/Weights.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Sector/Weights.lean) and [Sector/Degeneration.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Sector/Degeneration.lean). |
| `convolution_scheme()` and `lagrange_basis()` | Evaluation of both input polynomials and full coefficient interpolation from [Convolution/Rank.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Convolution/Rank.lean). |
| `LinearMap`, `restrict()`, `direct_sum()`, and `permute_axes()` | Independent local substitutions, complete independent blocks, and coordinate relabelings from the finite tensor operations. |
| `shared_first_tensor()` and `tag_shared_scheme()` | Shared X with matching Y/Z branch labels from [Separation/BranchTagging.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Separation/BranchTagging.lean). |
| `FiniteSeparation` | The three explicit Fourier maps, integer square weights, formal leading coefficient, and full branch/dot-product direct sum from [Separation/Basic.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Separation/Basic.lean) and [SquareWeights.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Separation/SquareWeights.lean). |
| `extract_branch()`, `factor_singleton_leg()`, `square_from_boundary()` | Exact branch selections, boundary coefficient-matrix factorization, and cyclic tensor products reindexed as a matrix multiplication tensor. The small example yields eight terms, not a 9/4 witness. |
| `CatalyticCertificate` and conditional planners | Supplied local maps for a finite positive-gain catalytic restriction, full coefficient verification, gain removal, and integer cost planning. Planning returns numerical conditions; a separate compiler generates the arrays. |
| `compile_powered_catalyst()` and `compile_absorption()` | Explicit powered local maps with one fixed auxiliary, repeated use of one catalyst, and independently verified matrix coefficient arrays. |
| `compile_gain_powered_catalyst()`, `compile_gain_absorption()`, and `compile_gain_target()` | Retain independent scalar gains, remove them by exact substitution, select integer parameters, generate matrix coefficients, and descend once. Array limits report incomplete construction. |
| `WitnessedStateSystem`, `find_nonnegative_dual()`, `assemble_catalyst()` | Exact finite-cone search, integer dual checking, and constructive assembly from retained restriction witnesses and tensor-block provenance. This is a derived finite construction, not a new Lean theorem. |
| `export_sector_order()` and `lift_witnessed_order()` | Actual finite interpolation restrictions and tensor-context local maps usable as witnessed LP rows. |
| `export_exact_type_order()` and `export_finite_separation_order()` | Coordinate pulls selecting exact branch-count words, and finite Fourier/square restrictions with all source copies retained. These yield ordinary witnessed rows without assuming state multiplicativity. |
| `DeterminantFiltration` | Explicit quotient/kernel bases and their output duals, the formal determinant degeneration, generated coefficient families, and an ordinary finite interpolation row. |
| `build_geometric_catalyst()` | Convert a supplied exact M_(d^b) scheme with r<k^b into a D=0 certificate over the same field, using a specified direct sum of matrix powers. It does not discover the required input scheme. |
| `projective_convolution_scheme()` | Finite-node interpolation plus a leading-coefficient/infinity term, verified as an exact convolution scheme. |
| `descend()` | The explicit basis/projection formula in [FieldDescent.lean](../lean/OAI/LinearAlgebra/MatrixMultiplication/Arithmetic/FieldDescent.lean). |

This correspondence is an implementation guide, not a formal proof that the
Python code refines the Lean definitions. Tests independently check tensor
coefficients, compare with ordinary matrix multiplication, exhaust all 256 pairs
of 2x2 matrices over F_2, exercise several characteristics, and reject malformed
certificates and invalid interpolation/Fourier hypotheses.
The polynomial tests also check cancellation between terms, tensor-valued
higher-order noise, and false certificates whose evaluations agree over F_2.
Generic tensor tests cover arbitrary coefficients, non-matrix bilinear maps,
empty axes, tensor products, projection, and round trips to the matrix API.
Sector tests check the formal identity on every coordinate, independently
contract the full local map matrices, and reject missing shifts/reversals.
Map and separation tests also contract actual source copies, reject hidden
overlaps, preserve genuine extension coefficients, detect missing Fourier
normalization and incorrect weights, and exercise adjusted periods in
characteristics three and five.

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

## Arbitrary finite tensor targets

`FiniteTensor(field, shape, coefficients)` represents a three-leg tensor with
independent axis sizes `(dim_X, dim_Y, dim_Z)`. Each axis enumerates a finite
coordinate set as `0,...,dim-1`; a zero-length axis is allowed. These sizes are
vector-space dimensions, not the `(n,m,k)` dimensions of matrix multiplication.
Coefficients are a flat immutable tuple in `(x,y,z)` lexicographic order, with
z varying fastest. `from_function()` can construct this representation from a
coefficient function. Any canonical field coefficients are allowed, including
extension elements; the target is not limited to zero-one support tensors.

`contract(left, right)` directly computes the bilinear map

```text
output[z] = sum_{x,y} T[x,y,z] * left[x] * right[y]
```

`TensorScheme(target, a, b, c)` declares an explicit target and supplies rank-one
coefficient families. `verify()` checks every tensor coordinate against this
target, including its coefficient values, rather than merely checking matching
dimensions. `apply(left, right)` requires that certificate and evaluates the
same bilinear map through the decomposition. `from_tensor(target)` constructs
a trivial exact decomposition with at most `dim_X*dim_Y` terms, omitting pairs
whose entire output vector is zero. It does not search for low rank. The zero
tensor has a valid zero-term decomposition.

For example, ordinary multiplication of polynomials with two and three
coefficients is a finite tensor of shape `(2,3,4)`:

```python
from reference import FiniteField, FiniteTensor, TensorScheme, TensorDegeneration

field = FiniteField(5)
target = FiniteTensor.from_function(field, (2, 3, 4),
                                   lambda i, j, k: int(i + j == k))
scheme = TensorScheme.from_tensor(target)
assert scheme.verify()
assert scheme.terms == 6
assert scheme.apply((1, 2), (3, 4, 1)) == (3, 0, 4, 2)
assert target.contract((1, 2), (3, 4, 1)) == (3, 0, 4, 2)

# Monomial lift and recovery for a target that is NOT a matrix tensor.
degeneration = TensorDegeneration.from_scheme(scheme, leading=1)
recovered = degeneration.recover_scheme()
assert recovered.target == target
assert recovered.apply((1, 2), (3, 4, 1)) == (3, 0, 4, 2)
```

`TensorDegeneration(target, leading, a, b, c)` takes the same polynomial
families as the matrix-specific API but compares the leading coefficient to
the supplied target. Low coefficients must still vanish formally. Recovery
returns a `TensorScheme`; tensor products and powers preserve explicit targets.
Arbitrary polynomial coefficient families remain valid inputs. The proof's
specific three-sector diagonal maps can now generate such families directly,
as described below; general polynomial local restrictions are still absent.

### General local maps, independent blocks, and shared inputs

`LinearMap(field, input_size, rows)` stores an immutable matrix in
`[output][input]` order. `identity()` and `selection()` make coordinate maps;
selection permits repetition and reordering. `second.compose(first)` applies
first, then second. Maps may have empty input or output axes.

`tensor.restrict(A,B,C)` independently substitutes on the three legs:

```text
restricted[x,y,z] = sum_{i,j,k} A[x,i] * B[y,j] * C[z,k] * tensor[i,j,k]
```

`scheme.restrict(A,B,C)` applies the same maps to its exact decomposition and
checks all resulting coefficients. `degeneration.restrict(A,B,C)` applies
constant maps to the target and polynomial families, preserving the leading
order. These are general scalar linear maps; arbitrary polynomial matrices
are not accepted by `LinearMap`.

`FiniteTensor.direct_sum(tensors)` and `TensorScheme.direct_sum(schemes)`
concatenate independent blocks on **all three axes** and set cross-block
coefficients to zero. Blocks may have different shapes. The family must be
nonempty and use one field representation. `permute_axes(order)` exchanges
leg roles; coordinate permutations within a leg use `LinearMap.selection()`.

By contrast, `shared_first_tensor(branches)` leaves X unchanged and tags only
Y and Z. It requires a nonempty family with one common ambient shape.
`tag_shared_scheme(scheme, branches, y_labels, z_labels)` derives those tags
by independent maps when both side coordinates determine every supported
branch's label. It checks each branch before field addition, preventing
characteristic-two cancellation from hiding invalid overlaps. The scheme must
target the coefficientwise sum of the branches. For the three-sector
construction, `side_labels` provides the two label lists and `tag_scheme()`
inserts them into a supplied exact decomposition of R without increasing its
term count or duplicating its first input.

### Convolution by evaluation and full interpolation

`convolution_scheme(field,a,b,nodes=None)` now generates the decomposition
used in `Convolution/Rank.lean`. With n=a+b-1 and distinct nodes v_r,

```text
A[r,i] = v_r^i
B[r,j] = v_r^j
C[r,k] = coefficient k of prod_{s != r} (X-v_s)/(v_r-v_s)
```

The default uses n canonical field encodings, including zero, and gives n
terms when both input axes are positive. Empty input axes give zero terms.
Extra distinct nodes are allowed. `lagrange_basis(field,nodes)` exposes the
full basis polynomials; the tests reconstruct every monomial of degree below
the node count as a formal polynomial. Zero is allowed here because there is
no division by a leading parameter power. Degeneration recovery still requires
nonzero nodes. Neither operation reduces encoded nodes modulo the characteristic.

Pass this generated source to `ThreeSectorConstruction.generate_degeneration()`
or `recover_power()` using `source_scheme=...`. For a=2,h=2 over F_8, the source
C(2,7) has 8 generated terms instead of the coordinate decomposition's 14;
power-two recovery has 3*8^2=192 terms and descent to F_2 has 3^2*192=1728.
These are certified decomposition lengths, without a timing or speed claim.
The source generator needs at least a+b-1 field elements; it does not search
for an extension if the supplied field is too small.

For generic tensor products, each axis is paired independently: `(i1,i2)` is
encoded as `i1*dim2+i2`. The zeroth power is the scalar unit of shape `(1,1,1)`.
This differs from the matrix-specific product API's row/column reindexing.
The existing `BilinearScheme` and `PolynomialDegeneration` constructors and
return types remain compatible. `as_tensor_scheme()` and
`as_tensor_degeneration()` expose their generic forms. Conversion back with
`to_matrix_scheme(n,m,k)` is allowed only when the target equals the exact
requested matrix tensor; matching axis lengths alone are insufficient.

`descend()` also accepts a `TensorScheme`. It requires all target coefficients
to be base-field scalars before descending the unchanged tensor through the
fixed extension. For a target with genuine extension coefficients, explicitly
call `scheme.project_to_prime_field()` to get a decomposition of **pi(target)**.
This changes the target by applying the constant-coordinate linear functional;
pi is not a field homomorphism and does not generally commute with tensor
products. Both operations use the same dimension-squared basis formula.

## Proof-derived three-sector construction

`ThreeSectorConstruction(field, a, h)` implements the finite formulas in
`Sector/Weights.lean` and `Sector/Degeneration.lean`. Its source tensor is the
actual convolution `C(a,B)`, with `B=3*h+a-1`, in the ambient spaces of shape
`(a,B,a+B-1)`. `convolution_tensor(field,a,b)` has coefficient one exactly
when `i+j=k`. Natural-number subtraction is truncated at zero as in Lean.
Zero parameters are supported; three nonempty branches require positive a,h.

Write `s=2*h+a-1` and `m=h+a-1`. The signed weights and shifted polynomial
maps are:

| Leg | Unshifted integer weight | Diagonal polynomial map |
| --- | --- | --- |
| X, coordinate i | 0 | 1 |
| Y, coordinate j | 1 on `h <= j < s`, otherwise 0 | t on that interval, otherwise 1 |
| Z, coordinate k | -1 on `m <= k < s`, otherwise 0 | 1 on that interval, otherwise t |

The third weight is shifted by **+1** before forming the polynomial maps.
`local_maps` exposes the full matrices in `[output][input]` order, with zero
off-diagonal entries and degree bounds `(0,1,1)`. The first map is the identity
on the original a-dimensional input space.

On source support, the unshifted total weight is always zero or one.
`retained` (R) keeps weight-zero terms and `erased` (E) keeps weight-one terms.
`polynomial_coefficient(i,j,k)` applies the diagonal maps to the source, without
using R or E to construct that polynomial. `verify()` checks, at every ambient
coordinate, the formal identity

```text
P(t) = t*R + t^2*E
source = R + E
R = left + middle + right
```

The branches have disjoint support, checked as an integer count before field
addition so that overlaps cannot be hidden by characteristic two. Their
coordinate embeddings from `C(a,h)` are:

| Branch | Ambient coordinates of a local `(u,r,v)` |
| --- | --- |
| left | `(u, r, v)` |
| middle | `(a-1-u, h+v, m+r)` |
| right | `(u, s+r, s+v)` |

Nonzero local coefficients satisfy `v=u+r`. `branch_embedding()` also accepts
off-support local coordinates, permitting checks of the whole local tensor.
The middle embedding reverses X and exchanges the last two legs; it is not
an ordinary translated copy with the same leg order. All three branch tensors
are stored in the **same** ambient spaces, retaining a first input of dimension
a, not 3*a. In particular the bilinear middle branch acts as a reversed-input
correlation when contracted in the original leg order.

```python
from reference import FiniteField, ThreeSectorConstruction

construction = ThreeSectorConstruction(FiniteField(2), a=2, h=2)
assert construction.shape == (2, 7, 8)
assert construction.verify()
assert sum(construction.source.coefficients) == 14
assert sum(construction.retained.coefficients) == 12
assert sum(construction.erased.coefficients) == 2
assert [len(construction.branch_support(name))
        for name in ('left', 'middle', 'right')] == [4, 4, 4]

degeneration = construction.generate_degeneration()
assert degeneration.verify()
assert degeneration.target == construction.retained
assert (degeneration.leading, degeneration.degree_bound) == (1, 2)
assert degeneration.terms == 14
```

`generate_degeneration(source_scheme=None)` applies those same diagonal maps
to an exact source decomposition and returns a certified `TensorDegeneration`
of R with leading order one. By default it uses the trivial coordinate-pair
decomposition, with a*B terms for positive source dimensions. A different exact
source decomposition can be supplied; its full target and coefficients are
checked, and its number of terms is preserved. This stage does not implement
an optimal-rank convolution scheme or a 9/4 matrix scheme.

This construction is checked by formal polynomial coefficients even over F_2,
where parameter evaluations alone cannot distinguish t from t^2.

### Powers, recovery, and descent of the proof-derived family

`construction.recover_power(exponent, nodes=None, source_scheme=None)` connects
the generated family to the existing generic tensor-power and interpolation
APIs. It validates the source decomposition, checks the available nodes before
allocating the power, and returns a certified `TensorScheme` for
`construction.retained.tensor_power(exponent)`. The order is:

```text
proof diagonal maps -> polynomial tensor power -> interpolation -> exact scheme
                                                                -> descend once
```

All factors use the **same** polynomial parameter. For power q,
`P(t)^tensor q = t^q * (R+t*E)^tensor q`; coefficient q is R^tensor q,
and the mixed noise coefficients remain present above it. The polynomial
degree bound is 2*q, so after dividing by t^q the normalized degree is at most
q. Thus q+1 distinct nonzero nodes suffice. With r source terms, the default
recovered family has `(q+1)*r^q` terms for positive a,h. This is a supplied
decomposition length, not a minimal rank or a speed estimate. Empty source
families can have a smaller bound. Power zero returns the scalar unit.

```python
from reference import FiniteField, ThreeSectorConstruction, descend

extension = FiniteField(2, (1, 1, 1))  # F_4, fixed before powering
construction = ThreeSectorConstruction(extension, a=2, h=2)
exact = construction.recover_power(2)  # nodes 1, alpha, alpha+1
assert exact.shape == (4, 49, 64)
assert exact.terms == 3 * 14**2 == 588
assert exact.verify()  # all 4*49*64 tensor coefficients

base = descend(exact)  # degree squared overhead: 2^2, applied once
assert base.terms == 2352
assert base.verify()
assert base.target == ThreeSectorConstruction(FiniteField(2), 2, 2).retained.tensor_power(2)
left = (1, 0, 1, 1)
right = tuple(j % 2 for j in range(49))
assert base.apply(left, right) == base.target.contract(left, right)
```

The demo checks powers 0, 1, and 2 of this particular proof-derived family
over the same F_4, then descends each recovered scheme to F_2. Its first input
dimension is a^q, because the three branches share that input. It also checks
that recovering before squaring would give 784 terms instead of 588. Smaller
integration tests compare both full orders including descent: for a=2,h=1,
power-two recovery then descent gives 768 terms, while recovering and
descending the seed before squaring gives 4096. Both compute the same tensor;
the latter repeats the interpolation and extension overheads.

The caller must supply a field with enough distinct nonzero elements. F_4 has
three, so the default degree bound supports powers through 2; power 3 raises
an error before building its tensor. F_2 alone is too small for the default
positive first-power example. Custom node order, extra nodes, and a supplied
exact source decomposition are supported. No extension is chosen automatically,
and a fixed finite field cannot provide an unbounded number of nodes for this
family. These finite examples therefore do not establish an asymptotic bound.
The result computes powers of the retained convolution tensor R; it is not yet
an extracted 9/4 matrix multiplication scheme.

## Finite Fourier separation and square weights

`FiniteSeparation(branches, period=None, root=None)` implements the finite
maps from `Separation/Basic.lean`. The branches must share their field and
ambient shape `(X,Y,Z)`, with M>0 branches. Its source is
`shared_first_tensor(branches)` of shape `(X,M*Y,M*Z)`. The maps start from
L complete independent copies of that source. The default period is 5*M,
or 5*M+1 when the former is zero in the characteristic. A supplied period
must be at least 5*M, be invertible in the field, and admit the supplied root
of exact order L. Missing roots raise an error in the supplied field.

Output coordinates are encoded as `(x,g)`, `((h,y),u)`, `((h,z),v)`, with
the last coordinate varying fastest. The labels g,h,u,v are numbered 1,...,M
when computing the integer expressions. At copy index r, the three local
Fourier factors are:

```text
first:   root^(r*2*g)
second:  root^(r*(u-h))
third:   (L in the field)^(-1) * root^(r*(-v-h))
```

`projection_maps` exposes the full independent matrices. `projection_scheme()`
forms the actual direct sum of source decompositions, applies those maps,
and checks its target against the independent phase-zero support formula.
Only the third map contains the normalization. The bounded integer phase
`u-v+2*(g-h)` cannot wrap modulo L, so averaging selects exactly phase zero.

The signed weights on the three legs are `g^2`, `h*u-h^2`, and `-h*v`.
Their sum on Fourier support is `(g-h)^2`, which vanishes exactly when
g=h and u=v. `weighted_maps(nonzero_node)` exposes the unshifted scalar
maps, including negative parameter powers, and rejects zero.

For formal polynomial coefficient families, the implementation shifts the
three weights by `(0,M*(M-1),M^2)`. All local exponents are then nonnegative.
Writing S for the sum of these shifts, `generate_degeneration()` certifies:

```text
all coefficients below S vanish
coefficient S = separationTarget
summed polynomial degree <= S+(M-1)^2
```

It also compares every full tensor polynomial to the independently computed
Fourier average with integer weights. `separationTarget` has matching sector
labels on all legs and a matching auxiliary pair u=v. `as_direct_sum(scheme)`
reorders the first axis so this is literally the full direct sum of
`branch[h] tensor dot_M`, where dot_M has shape `(1,M,M)` and coefficients
`dot[0,u,v] = 1` exactly when u=v. The source's shared X becomes independent
branch inputs only through this verified separation construction.

`recover_scheme()` interpolates at `(M-1)^2+1` nonzero nodes. Its bound is on
the **summed tensor polynomial**, which can be smaller than the sum of local
family degrees because the Fourier terms cancel. The generic
`TensorDegeneration.recover_scheme(degree_bound=...)` accepts such a bound only
after inspecting every formal tensor polynomial; an invalid bound is rejected.
`recover_power(q)` uses the same parameter in all factors and needs
`q*(M-1)^2+1` nodes, checking availability before source-copy or power allocation.
Its supplied term bound is that node count times `(L*r)^q`, with r source terms.

```python
from reference import FiniteField, ThreeSectorConstruction, convolution_scheme, FiniteSeparation, descend

field = FiniteField(2, (1, 1, 0, 0, 1))  # F_16, fixed for this finite construction
c = ThreeSectorConstruction(field, 2, 1)
source = convolution_scheme(field, 2, c.source_width)  # 5 terms
retained = c.recover_power(1, source_scheme=source)     # 10 terms
tagged = c.tag_scheme(retained)                       # 10 terms, shape (2,12,15)
s = FiniteSeparation(c.branches)                      # M=3, L=15
assert tagged.target == s.source
exact = s.recover_scheme(source_scheme=tagged)         # 5*15*10 = 750 terms
direct = s.as_direct_sum(exact)
assert direct.target == s.direct_sum_target
assert direct.shape == (6,36,45)
base = descend(direct)                               # 4^2*750 = 12000 terms
assert base.verify()
```

Run `python3 -m reference.proof_demo` for both the 8-term convolution power
example and this complete finite separation pipeline, including direct
bilinear evaluation after descent. Full direct sums retain ambient padding;
there is no support compaction or minimum-rank search. These finite maps do
not provide the determinant, entropy, or spectral-existence steps that lead
to the 9/4 matrix multiplication bound.

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

## Coefficient extraction research

The [research checkpoint](research/README.md) records the source audit,
complete finite-search route, and a conditional powered-catalyst route.
An explicit useful catalyst would supply enough data for a convergent
coefficient procedure. Finding that catalyst remains a research problem.
Current code checks supplied catalysts, compiles their powered maps, and
generates exact decompositions. The conditional cost planner itself does not
construct or certify those arrays. A finite dual with actual row witnesses
can also be assembled into an explicit catalyst. Dense caps stop incomplete
constructions before large powers or coefficient families are allocated.

`proof_matrix_pipeline()` executes the finite separation-to-matrix connection
over F_16. It extracts each branch from the generated 750-term scheme,
restores the middle branch's order, and factors each C(2,1) boundary to two
terms. Three cyclic factors give an exact eight-term 2-by-2 matrix scheme.
The demo checks the full tensor, a matrix product, and powers followed by
one final descent. These are classical coefficients obtained through the
pipeline and carry no improved exponent claim.

The catalyst demo uses d=2,k=9, initially generating nine-term M_2 and
81-term M_4 schemes through the original absorption route. Retaining the
positive scalar gain reduces finite supplied term bounds by exact scalar
elimination. The actual powered b=2 construction retains eleven independent
units and yields a verified 70-term M_4 scheme, compared with 81 terms when
those gains are discarded. None of
these conservative control witnesses supplies the desired 9/4 gap.

The bounded automatic interface accepts exact supplied tensors and schemes:

```python
from fractions import Fraction
from reference import (TensorScheme, compile_gain_target, coordinate_catalyst)

certificate = coordinate_catalyst()  # conservative d=2,k=9 control
S = TensorScheme.from_tensor(certificate.S)
D = TensorScheme.from_tensor(certificate.D)
result = compile_gain_target(certificate, S, D, Fraction(31, 10))
assert result.status == "verified"
assert result.prime_scheme.verify()  # eight-term M_2 over F_2
```

`plan_gain_catalytic_parameters()` returns numerical costs only.
`compile_gain_target()` checks the final strict inequality against actual
arrays after one descent to the prime field. A `resource_cap` result carries
no final scheme; neither a finite parameter cap nor a dense allocation cap
means the requested witness is impossible.

`build_geometric_catalyst()` provides another explicit input interface. For
example a supplied seven-term M_2 scheme gives a d=2,k=8 certificate with
one scalar gain and no catalyst D. The geometric construction matches actual
matrix-power blocks and retains unused scalar terms. Its useful k=5 variant
still needs an unknown good matrix-power decomposition; see the
[same-field search reduction](research/geometric-catalyst.md).

The optional search experiment finds an unseeded seven-term F_2 matrix scheme;
its small JSON artifact can be verified without the solver. The witnessed demo
then builds a finite dual, assembles its maps, and generates seven- and
49-term matrix schemes after exact simplification. These are known algorithms
reconstructed through the new interfaces.

For projective convolution, C(a,b) needs a+b-2 finite distinct nodes and one
infinity term. Zero is allowed as a finite node. C(3,3) therefore has an exact
five-term scheme over F_4, while the ordinary finite-node interface requires
five field elements. Projective interpolation never enlarges a supplied field.

## Scope and remaining work

- Only explicitly represented finite fields are supported. Arithmetic is exact;
  floating-point numbers, arbitrary fields, and algebraic closures are absent.
- The user supplies an extension large enough for the roots and interpolation
  nodes. Missing roots or too few nodes raise errors rather than inventing a
  root or reusing colliding nodes. The demo uses F_16 for period 5 and F_25 for
  the adjusted period 6 in characteristic 5.
- Arbitrary finite three-leg targets and their supplied scalar/polynomial
  decompositions are supported. The proof's three-sector diagonal restriction
  is generated; general scalar local restrictions, independent direct sums,
  leg permutations, and side-label tagging are supported. Higher-order tensors
  and general polynomial local-map matrices are not implemented.
- The three-sector finite coefficient construction, powering, interpolation,
  fixed-extension descent, and finite Fourier/square-weight separation are
  connected. Determinant/character inequalities, entropy limits,
  and spectral existence argument are not translated.
- Strassen and naive decompositions are input fixtures. Neither is a 9/4
  decomposition extracted from this proof. The example scaling is deliberately
  artificial; descent is demonstrated, not proposed as an improvement to them.
- No running-time, memory-efficiency, bit-complexity, or 9/4 arithmetic-operation
  bound is claimed. Dense tensor powers and exhaustive coefficient checking
  become large quickly. Use small dimensions and powers.

The three-sector generator supplies proof-derived polynomial families to
the verified power/recovery/descent and finite separation pipelines.
Further finite constructions and turning the
current spectral existence proof into a generator for the 9/4 schemes remain
separate mathematical and implementation tasks.
