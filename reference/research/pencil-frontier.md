# The remaining `(2,3,3)` catalyst frontier

This note studies the positive-gain certificate
`D ⊕ Unit ⊕ (M₂⊗S) ≤ D ⊕ 5S` when the concise dimensions of `S` are
`(2,3,3)`. It follows the exclusions in `catalyst-obstructions.md` but does
not enlarge the Python solver or claim a catalyst witness. The central new
result is an elementary obstruction excluding **every regular matrix
pencil**, including the Jordan-three candidate. The singular pencil with arbitrary nonzero `D` remains
unresolved; the zero-catalyst case is excluded in `star-zero-catalyst.md`.

No Lean sources were edited or compiled for this note. The proof below was
independently reviewed by the main and constructive-research agents; this is
a mathematical/source audit, not a new kernel-checked theorem.

## Explicit normal forms

Write `S(x,y)=xP+yQ` for its two `3×3` coefficient slices.

| Type | Slices / support | Rank upper certificate | Status at `d=2,k=5` |
| --- | --- | --- | --- |
| Semisimple regular | `P=I₃`, `Q` diagonal and nonscalar | 3 | Excluded |
| Jordan-two regular, same eigenvalue | `P=I₃`, `Q=E₀₁` | 4 | Excluded by commuting-source argument |
| Jordan-two regular, distinct eigenvalue | `P=I₃`, `Q=E₀₁+E₂₂` | 4 | Excluded by commuting-source argument |
| Jordan-three regular | `P=I₃`, `Q=E₀₁+E₁₂` | 4 | Excluded by commuting-source argument |
| Singular | `P=E₀₀+E₁₂`, `Q=E₀₁+E₂₂` | 4 | Not excluded here |

Regular means `det(xP+yQ)` is a nonzero polynomial, not that an invertible
slice must already exist at a point of a tiny finite field. Over an algebraic
closure, normalize an invertible slice to `I₃` and use the Jordan form of the
other slice. This gives the listed semisimple/Jordan possibilities. Over the
original finite field there can be irreducible factors and different rank
upper bounds; the exclusion below applies directly to every regular pencil
and does not need diagonalization, eigenvalues, or a four-term scheme.

A four-term Jordan-three scheme, valid over every field with signs interpreted
in that field, is

```
A = ((1,0), (1,0), (0,1), (1,1))
B = ((1,0,0), (0,0,1), (0,1,0), (1,1,0))
C = ((1,-1,0), (0,0,1), (0,-1,1), (0,1,0)).
```

The independently found binary certificate is exactly this formula modulo
two. For the Jordan-two cases, the `2×2` regular block uses its three support
products, and the last diagonal coordinate uses one product with first-leg
form `x` or `x+y`. The singular form has exactly four support products, so
its coordinatewise decomposition has four terms. Rank three would force
regularity: concise other-leg dimensions three imply the three left and
right factor vectors of a three-term decomposition are bases, making the
determinant pencil a product of three nonzero first-leg linear forms.
Thus the singular rank is exactly four.

## The regular-pencil obstruction

### Repeat without changing the catalyst

A purported certificate can be additively repeated `q` times with one fixed
`D`. Its source is `D⊕5qS` and its target is
`D⊕q Unit⊕q(M₂⊗S)`. Drop the target's `D` coordinates, and merge the first
axes of the product copies while keeping their other axes independent.
The resulting target has a first-leg space consisting of `q` scalar labels
and the shared eight coordinates of `M₂⊗S`. Its other two legs have equal size

```
N = 13q.
```

At the identity slice, choose the unit coordinates to be one and the shared
matrix/pencil coordinates to be `I₂⊗I₃`. Two further slices are left
multiplication by `E₁₂⊗I₃` and `E₂₁⊗I₃`, acting on every product block.
Their normalized commutator has rank `12q`; it vanishes only on the `q`
scalar-unit coordinates. In characteristic two, `diag(1,-1)` becomes `I₂`
and still has full rank. The nonzero commutator minor and invertible target
slice remain nonzero at a generic first-leg point.

### The source can be normalized to commuting blocks

Take any exact rank upper certificate `C` for the fixed `D`; replace it by a
restriction from `C` scalar tensors. Pull the target first-leg map through
to the source. Each source `S` block now has slices

```
a(ξ) I₃ + b(ξ) J,
```

where `a,b` are linear forms in the target first-leg variables `ξ`, after
normalizing the original regular pencil once. There are only three cases:

1. The image has dimension two. The forms `a,b` are independent, so its
   determinant polynomial is nonzero. A generic slice is invertible; all
   normalized slices are rational functions of the same matrix `J` and
   commute.
2. The image has dimension one. Every slice is `λ(ξ)H` for one fixed matrix
   `H`. If `H` is singular, factor it once through its matrix rank `r≤3`.
   This gives a fixed square `r×r` block with slices `λ(ξ)I_r`. Crucially the
   row/column compression is **independent of `ξ`**. These normalized slices
   also commute.
3. The image is zero; discard the block.

The scalar source terms from `D` behave the same way: discard identically
zero linear forms; each surviving generic pivot is nonzero. Thus the whole
source, after fixed support compression, is a commuting normalized slice
family of square ambient size

```
L ≤ C+15q.
```

There are finitely many source pivot determinants and target determinant/minor
polynomials. Over an infinite extension, their product is nonzero, so a
generic point avoids all their zeros simultaneously. This is sufficient for
a lower-bound proof over the original finite field: its exact certificate
extends to that field unchanged. No finite-input identity replaces the
coefficient equations.

### A compression of commuting slices has a bounded commutator

Suppose the target slices are `T_i=P H_i Q`, where the square normalized
source matrices commute and the target pivot `T₀` is invertible. Normalize
the source pivot `H₀` as well. Set

```
U = T₀⁻¹ P H₀
V = Q
F_i = H₀⁻¹ H_i.
```

Then `UV=I_N` and `T₀⁻¹T_i=U F_i V`. The matrices `F_i` commute. The
projection `VU` has rank `N`, so `Z=I_L−VU` has rank `L−N`. In the
commutator, the terms containing `I_L` cancel, leaving two products factoring
through `Z`. Therefore

```
rank([T₀⁻¹T_i,T₀⁻¹T_j]) ≤ 2(L−N).
```

This proof uses source commutation, not a rank decomposition into diagonal
terms for `S`. It works equally with Jordan blocks. Applying it to the
target commutator gives

```
12q ≤ 2(C+15q−13q) = 2C+4q,
```

which is impossible for `q>C/4`. This excludes every regular `(2,3,3)`
auxiliary, for arbitrary finite `D` and every characteristic. The inequality
does not assume tensor-rank direct-sum additivity. It only uses an upper
certificate for one fixed `D`, actual local restrictions, and matrix ranks.

The argument also works for any square regular two-parameter pencil of size
`n`: `N=(4n+1)q`, `L≤C+knq`, and target commutator rank `4nq`. Letting `q`
grow forces `kn≥6n+1`, hence integer `k≥7`. This generalization relies on the
same fixed one-dimensional compression case; it does not automatically
apply to arbitrary multi-parameter singular slice spaces.

## Why the singular normal form is the only remaining shape candidate

For a concise singular `3×3` pencil, normal rank one is impossible: two
rank-one generators whose sums remain rank one have a common image or
kernel, contradicting another leg's conciseness. Thus normal rank is two.
The nonzero adjugate is a homogeneous degree-two polynomial matrix of rank
one. Over the polynomial UFD, factor it as `g u vᵀ`, with primitive
homogeneous left/right kernel vectors. A constant kernel vector would again
contradict conciseness. Thus both vector degrees are at least one; their sum
is at most two, so both degrees are one and `g` is constant.

Write the right kernel as `v=xv₀+yv₁`. Its coefficients are independent, and

```
P v₀ = 0
Q v₁ = 0
P v₁ = −Q v₀ = w ≠ 0.
```

Complete these to a column basis with `v₂`. Row conciseness makes
`w,Pv₂,Qv₂` independent. Choosing these row and column bases produces the
listed singular form, the Kronecker block `L₁⊕L₁ᵀ`. This derivation works
over the base field; there is no additional irreducible-polynomial singular
case. After a column permutation it is the symmetric star pencil

```
[[0,x,y],
 [x,0,0],
 [y,0,0]].
```

Every generic slice has rank two, so the source-pivot step of the regular
argument fails. Its two-dimensional slice family cannot be compressed once
to an invertible square pivot: it has no common left/right kernel to remove.
It does restrict to dual-number multiplication, but its rank-four source
upper bound only gives `5*4≥16`, so the previous algebra filter does not
exclude it. The product with `M₂` has no identity slice either; treating it
as a regular pencil would be incorrect.

## Exact finite checks and remaining obligation

Using the binary Koszul implementation in
`catalyst_obstructions_experiment.py`, the first-leg dimension-eight
flattenings of `M₂⊗S` have these ranks:

| Normal forms | p=1 | p=2 | p=3 | p=4 | p=5 | p=6 |
| --- | --- | --- | --- | --- | --- | --- |
| All three Jordan forms above | 96 | 312 | 540 | 540 | 312 | 96 |
| Singular star | 96 | 296 | 492 | 492 | 296 | 96 |
| Rank-one denominator `binomial(7,p)` | 7 | 21 | 35 | 35 | 21 | 7 |

The largest regular ratio is `540/35≈15.429`; the singular ratio is
`296/21≈14.095`. These finite lower bounds alone do not approach the
twenty-term source budget. With shared first axes, the Koszul matrix for `q`
independent other-leg blocks is block diagonal and its rank is exactly `q`
times the single-block rank. Keeping the `q` gain units gives a decisive
filter only if the ratio is **strictly greater than 19**; rounding a
single-block lower bound to twenty is not itself sufficient without checking
this ratio. None of the displayed tests satisfies it.

The constructive-research agent independently verified four-term binary
decompositions and reviewed the commuting-source proof. The singular case
needs either an actual catalyst certificate or a new lower obstruction that
handles its varying kernels. Neither has been obtained here. The next
specific search target is the singular star, not the already excluded
Jordan-three/truncated-linear-multiplication family.

## The tied-product boundary for the singular star

After permuting coordinates, the product `M₂⊗S_star` is the bilinear map

```
(A,B), (X₀,X₁,X₂) ↦ (A X₀, B X₀, A X₁+B X₂),
```

where all five input labels represent `2×2` matrices. Thus the first leg has
size eight and the other legs have size twelve. It ties two rectangular
multiplication tasks: vertically stack `(A,B)` for the first two outputs,
and horizontally concatenate `(A,B)` for the last output. The two tasks
share the eight first-leg variables, while their other inputs and outputs
are disjoint. Treating this as the direct sum of two independent rectangular
multiplication tensors would incorrectly duplicate those first-leg variables.

For generic `A,B`, the map from the three `X` matrices to the three outputs
has rank eight. Indeed the first two outputs depend only on `X₀` and have
rank four when `A` is invertible; the last output depends only on `X₁,X₂`
and has rank four under the same condition. With `q` scalar gains the target
normal rank is therefore `9q`. The pulled-back source has normal rank at
most `C+10q`: each nonzero singular-star block has normal rank two, and the
fixed catalyst is replaced by `C` rank-one terms. This necessary inequality
allows `k=5`; it does not force `k≥6`.

At a generic pivot, one singular-star source block has a two-dimensional
square core whose normalized slices are scalar and hence commute. However,
its remaining row and column couple to that core as the pencil parameter
changes. For example, at `(x,y)=(1,0)` the symmetric star has kernel spanned
by the third coordinate; its `y` slice couples that coordinate to the first
one on both sides. Deleting the pivot kernel destroys these slices. A
parameter-dependent factorization does not repair the regular proof,
because the same local restriction maps must represent every slice at
once. These kernel couplings contribute additional terms to a compressed
commutator; the regular proof's cancellation is therefore unavailable.
No bound controlling those terms tightly enough to exclude the star has
been obtained in this pass.

A proposed alternative was an actual coefficient restriction

```
Unit ⊕ (M₂⊗S_star) ≥ Unit ⊕ M(3,2,4).
```

The dimensions are compatible. Such a map, if it preserved the scalar
summand separately, could be repeated `q` times and followed by merging
the rectangular tensor's first legs to obtain `q Unit⊕M(3,2,4q)`. The
external rectangular theorem audited in `catalyst-obstructions.md` gives

```
rank ≥ q + ceil((96q+2)/5).
```

Its asymptotic coefficient is `20.2q`, enough to beat the source upper bound
`C+20q`. The single-block value twenty for the rectangular part does **not**
justify replacing its repeated lower bound by `20q`: the applicable theorem
has asymptotic coefficient `19.2q` before adding the units. However, the
structural image-intersection proof in the final section rules out this map
even under general mixing. Compatible dimensions do not establish a
restriction, and this specific route is now closed.

The finite-search agent reported a bounded sixty-second exact binary map
search for the original star catalyst with no result. That timeout supplies
neither a witness nor an impossibility certificate. The useful established
frontier is therefore: all regular pencils are excluded; this one singular
orbit remains; its varying kernels and tied rectangular products are the
features a new proof or certificate must handle.

## Structural exclusion of the proposed rectangular restriction

The proposed map to `M(3,2,4)` mentioned above is ruled out
without coefficient search. This is a restriction obstruction for that
particular rectangular target, **not** an impossibility proof for the
original singular-star catalyst.

Write `T=M₂⊗S_star` in the tied-product coordinates. For a fixed first input
`(A,B)`, its slice is the map

```
H(A,B): (X₀,X₁,X₂) ↦ (A X₀, B X₀, A X₁+B X₂).
```

The first two output blocks depend only on `X₀`; the third depends only on
`X₁,X₂`. The two columns of each matrix are independent copies of the same
linear operation. Therefore, exactly over every field,

```
rank H(A,B) = 2 rank([A;B]) + 2 rank([A B]) ≤ 8.
```

If this rank is eight, both displayed rectangular matrices have rank two.
In particular `[A B]` is surjective onto `K²`, so the full slice image
contains the fixed output subspace

```
E = {(0,0,Z): Z ∈ Mat₂(K)},    dim E = 4.
```

Suppose a local restriction from `T` to `M(3,2,4)` existed. Both third-leg
spaces have dimension twelve; conciseness of the target forces its output
map `P` to be surjective, hence invertible. For each first input
`X ∈ Mat_(3,2)(K)`, there are a corresponding source pair `(A_X,B_X)` and
one fixed second-input map `Q` such that the target slice is

```
L_X = P H(A_X,B_X) Q.
```

For every rank-two `X`, the rectangular multiplication slice
`L_X: Mat_(2,4)(K) → Mat_(3,4)(K)` has rank eight. Since `H` has rank at most
eight, restricting its domain by `Q` cannot lose image dimension in this
case. Thus

```
Im(H(A_X,B_X) Q) = Im H(A_X,B_X),
P(E) ⊆ Im L_X.
```

But `Im L_X = Col(X)⊗K⁴`. Choose three rank-two `X` whose column spaces are
`span(e₀,e₁)`, `span(e₀,e₂)`, and `span(e₁,e₂)`. Their three output images
have intersection zero, over every field, including `F₂`. The fixed
four-dimensional space `P(E)` cannot lie in this intersection. This is a
contradiction. It excludes **all** mixing of the local linear maps, not
just coordinate cuts or the row-identification ansatz.

### A scalar summand cannot rescue this exact restriction

More strongly, there is no restriction

```
Unit ⊕ T ≥ Unit ⊕ M(3,2,4),
```

whether or not its local maps preserve the unit separately. The output
spaces now both have dimension thirteen, so the output map is still
invertible. Every source slice has rank at most nine. Consider target first
inputs with scalar coordinate one and rectangular component of rank two.
Their target slices have rank nine, attaining this source maximum. Thus the
source scalar coordinate is nonzero, `H` has rank eight, and restricting the
second input does not change the full source slice image. That image contains
a fixed five-dimensional subspace: the source unit-output line plus `E`.
The output isomorphism carries it into every chosen target slice image.
However, using the same three column planes, the intersection of those
target images is just the one-dimensional target unit-output line. This is
impossible.

The same argument excludes `m Unit⊕T ≥ m Unit⊕M(3,2,4)` for any finite
`m≥0`: use target scalar coordinates all equal to one. The source maximum
slice rank is `m+8`; equality forces a common source image of dimension
`m+4`, whereas the target image intersection has dimension `m`.

These proofs use exact slice ranks and the invertible third-leg map. They
do not require a tensor-rank lower bound, tensor-rank additivity, an infinite
field, or an imported theorem. They also do not apply directly to
`D⊕Unit⊕T ≤ D⊕5S_star`: its output dimensions and maximum slice ranks differ,
so neither saturation nor an invertible output map is forced. The singular
star remains an unresolved catalyst candidate; this particular rectangular
restriction route is closed.

As a bounded independent matrix check of the slice calculation, all 256
pairs `(A,B)` over `F₂` and all 6,561 pairs over `F₃` were enumerated. The
rank distributions were respectively `{0:1, 4:27, 6:36, 8:192}` and
`{0:1, 4:128, 6:384, 8:6048}`. Every slice satisfied the displayed stacked/
horizontal rank identity, and every rank-eight slice contained each basis
vector of the third output block. These checks corroborate the matrix
identity; the impossibility conclusion is the field-independent proof above.


## Follow-up: the zero-catalyst star case is excluded

The later `star-zero-catalyst.md` proves that no exact
`m Unit⊕M₂⊗S_star≤5S_star` exists for any positive natural `m`, over any
field. It uses the source hyperplane's large low-rank subspaces and the
five-dimensional bad-pair locus in the target. Thus the unresolved star
frontier now concerns **nonzero `D` only**. The regular-pencil exclusions
and the rectangular restriction obstruction above are unchanged.
