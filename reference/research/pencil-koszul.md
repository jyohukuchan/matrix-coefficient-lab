# Concise (2,3,3) pencil frontier: exact Koszul checks and source compression

This research pass examined four explicit auxiliary pencils over `F_2`.
It generated actual four-term auxiliary decompositions and exact projected
Koszul matrix ranks. It also independently audited the stronger commuting-
source obstruction discovered by research A. That argument excludes the
regular square pencils; the singular pencil remains unexcluded here.

No useful catalyst was found. The new source-compression argument is a
mathematical derivation independently checked by agents, not a Lean theorem.

## Four exact normal-form inputs

All tensors have shape `(2,3,3)`; the first leg selects the two matrices.

```
Jordan3:       I_3, J_3,
Jordan2+0:     I_3, J_2 direct_sum [0],
Jordan2+1:     I_3, J_2 direct_sum [1],

singular:     [[x,y,0],
               [0,0,x],
               [0,0,y]].
```

Here the Jordan blocks are nilpotent with ones on the superdiagonal.
Unseeded normalized-input subset searches found and fully verified four-term
schemes after respectively `1720`, `1352`, `1420`, and `1413` subsets.
Every rank-three subset was checked first. Thus the four ranks are exactly
four over `F_2`; no rank upper bound was inferred from the shape.

The nontrivial Jordan3 certificate is

```
A=((1,0),(1,0),(0,1),(1,1))
B=((1,0,0),(0,0,1),(0,1,0),(1,1,0))
C=((1,1,0),(0,0,1),(0,1,1),(0,1,0)).
```

Let `T=M_2 tensor S`, with shape `(8,12,12)`. A `d=2,k=5` positive-gain
catalyst would, after repeating its comparison `q` times, imply an upper
bound `rank(D)+20q` for the repeated unit-plus-T target.

## Exact Koszul results

For a tensor of first-axis dimension `a`, a degree-`p` Koszul matrix has shape

```
(choose(a,p+1)*dimC, choose(a,p)*dimB)
```

and each rank-one tensor contributes at most `choose(a-1,p)` to its rank.
The reported lower quantity is the exact rational ratio of matrix rank to
this rank-one factor; ceilings are taken only after repetitions are included.

For first dimension eight, the unprojected ratios at degrees one through six
are:

| Auxiliary | p=1 | p=2 | p=3 | p=4 | p=5 | p=6 |
| --- | --- | --- | --- | --- | --- | --- |
| Each regular Jordan form | 96/7 | 104/7 | 108/7 | 108/7 | 104/7 | 96/7 |
| Singular | 96/7 | 296/21 | 492/35 | 492/35 | 296/21 | 96/7 |

The experiment exhausts **all 255** one-dimensional kernels of surjective
first-axis maps `F_2^8 -> F_2^7`. Changing a quotient basis only applies an
invertible coordinate change and preserves the corresponding Koszul rank.
At degree three, the rank-one factor is twenty:

| Auxiliary | Exact matrix-rank histogram | Maximum ratio |
| --- | --- | --- |
| Jordan3 | 300:9, 312:246 | 78/5 = 15.6 |
| Jordan2+0 | 288:9, 312:246 | 78/5 |
| Jordan2+1 | 300:18, 312:237 | 78/5 |
| Singular | 276:27, 278:36, 280:192 | 14 |

A separate bounded sample of **256** deterministic seed-zero surjections
`F_2^8 -> F_2^5`, at degree two, gave:

| Auxiliary | Exact matrix-rank histogram | Maximum ratio |
| --- | --- | --- |
| Jordan3 | 84:3, 88:7, 92:47, 96:199 | 16 |
| Jordan2+0 | 76:1, 80:4, 84:5, 88:47, 96:199 | 16 |
| Jordan2+1 | 82:1, 84:6, 88:20, 92:76, 96:153 | 16 |
| Singular | 72:1, 76:1, 78:29, 80:225 | 40/3 |

The first sampled maximum used row bit masks `(217,99,195,228,108)`.
The dimension-five sample is not exhaustive, and no extension-field
projection inventory was exhausted.

Other first-axis dimensions cannot improve the required ratio by dimensions
alone: for the original `(8,12,12)` orientation, only first dimensions eight,
seven, and five admit a matrix-dimension ratio above nineteen. For the other
orientations `(12,8,12)` and `(12,12,8)`, every first-axis projection dimension
from one to twelve gives ratio at most eighteen, irrespective of coefficients.
Further projections on the other two axes cannot increase a fixed Koszul
matrix rank or change its rank-one factor.

## Shared-first amplification, including the scalar gain

Ordinary tensor-rank direct-sum additivity is not assumed. Instead, restrict
`q` independent copies of `T` to one tensor whose first axis is shared while
the other two axes retain copy labels. Its degree-`p` Koszul matrix, after
ordinary matrix row/column permutations, is block diagonal with `q` copies of
the single-block Koszul matrix. Thus its matrix rank is exactly `q*kappa`,
with the same factor `h=choose(a-1,p)`.

The positive scalar blocks can be kept independently. The special elementary
scalar-summand identity `rank(q Unit direct_sum X)=q+rank(X)` therefore gives

```
rank(repeated target) >= q + ceil(q*kappa/h).
```

Consequently **kappa/h > 19** suffices to contradict the source budget
`rank(D)+20q` for large enough `q`. A single-block ceiling must not be
multiplied by `q`: for example, `ceil(kappa/h)=20` does not imply lower bound
`20q`. The rational slope is the relevant quantity.

The standalone experiment checks actual sharing maps for `q=2`, retains both
scalar units, isolates the core, and verifies exact Koszul matrix-rank doubling
for the best seven- and five-dimensional projections of every form. No checked
ratio exceeds nineteen, so these computations give no arbitrary-D exclusion.

## Independently audited stronger obstruction for regular pencils

For the three regular forms, research A found a different obstruction that
does not need the inconclusive Koszul bounds or rank(S)=4.

Suppose a finite `D` gives the desired positive-gain comparison. Repeat it
`q` times while keeping `D` fixed, and replace the source `D` by a supplied
`C`-term decomposition into scalar units. Pull the eventual target's first
axis back through the candidate restriction. Every source auxiliary block
then has matrix family

```
H_h(z)=a_h(z)*I_3+b_h(z)*J.
```

There are exactly three cases for the image of its two first-axis linear
forms:

1. **Dimension two:** the map onto the pencil is surjective. Its generic
   determinant is a nonzero polynomial because the original pencil contains
   `I_3`. A generic pivot normalizes this block to commuting operators:
   all its matrices and the inverse pivot are rational functions of `J`.
2. **Dimension one:** all slices are `lambda(z)*H` for one fixed matrix `H`.
   If `H` is singular, factor it **once** through its rank-r space. The
   resulting square block is `lambda(z)*I_r`, `r<=3`. This factorization is
   independent of `z`; choosing a different factorization for every slice
   would be invalid.
3. **Dimension zero:** remove the zero block.

The scalar units from the source decomposition of `D` are already commuting
one-dimensional blocks. The entire normalized source therefore acts on a
commuting square ambient space of dimension

```
L <= C+15q.
```

The shared-first target, with its `q` scalar units retained, has square other
axes of size `N=13q` and a generic invertible slice. At the identity pivot,
the two matrix slices `E12 tensor I_3` and `E21 tensor I_3` give a commutator
of rank `12q`: matrix multiplication carries a spectator matrix-column factor
of dimension two. This rank remains valid in characteristic two because
`[E12,E21]=diag(1,-1)` is invertible there too.

The target determinant and a full commutator minor are nonzero at the known
identity pivot. Source determinant factors and that minor can therefore be
made nonzero **simultaneously at a generic point** over a rational-function
extension. This uses nonzero formal polynomials, not a claim that a finite
base field has enough evaluation points. Any alleged coefficient restriction
persists under that extension.

Write the compression of the normalized source as

```
AB=I_N,     target_operator_i=A*F_i*B,
```

where the `F_i` commute and act on dimension `L`. Set `Q=I_L-BA`. It has rank
`L-N`. The compressed commutator is a sum of two terms factoring through Q,
so its rank is at most `2*(L-N)`. Hence

```
12q <= 2*(L-N) <= 2C+4q,
4q <= C.
```

Choosing `q>C/4` contradicts the fixed finite catalyst. This excludes each
regular square-pencil auxiliary for **any finite D** at `d=2,k=5`.
The argument is stronger than a tensor-rank upper/lower comparison: it uses
the structure of every pulled-back source pencil.

This derivation was independently checked by main and both research agents.
It has not been added as a Lean theorem. It does not assume cancellation or
general tensor-rank additivity, and no division by the field's integer two
is used in the rank inequality.

## Why the singular pencil remains a different problem

For `L1 direct_sum L1-transpose`, every nonzero pencil slice has rank two,
but the joint row and column spans both have dimension three. A
two-dimensional pulled-back family therefore has identically zero determinant
and cannot be compressed once through one fixed rank-two row/column space.
The regular-source normalization argument does not apply to this case.

It can also be viewed as two rectangular branches sharing the first leg:
one branch has shape `(2,1,2)` and the other `(2,2,1)`. They are not a full
independent direct sum. Separating their tensor ranks additively would require
an additional valid argument or actual restriction maps.

This pass supplies the explicit singular tensor, its verified rank-four
scheme, exhaustive seven-dimensional GF2 quotient results, and bounded
five-dimensional projection results. None gives a `k=5` catalyst or excludes
that singular auxiliary for arbitrary D. It is the remaining pencil frontier
among the four forms studied here.

Run the executable checks with

```
python3 -m reference.research.pencil_koszul_experiment
```

The script uses the existing exact binary Koszul implementation, bounded
bit-matrix caps, and full coefficient verification of its maps and auxiliary
schemes. No Lean source was edited or built.

## Independent audit: the singular star cannot have an empty-D gain witness

The following new argument, proposed by main and independently audited here
and by Research A, excludes an exact restriction

```
5 S_star -> m Unit direct_sum (M_2 tensor S_star)
```

for every positive integer m, over every field. It does not exclude a
nonempty catalyst D. This is an unformalized mathematical derivation, not
a Lean theorem or an inference from the numerical Koszul bounds.

Write the target product slice at first input `(A,B)` as the tied map

```
(X_0,X_1,X_2) -> (A X_0, B X_0, A X_1+B X_2).
```

For `V=vertical_stack(A,B)` and `H=horizontal_concat(A,B)`, its exact matrix
rank is `2 rank(V)+2 rank(H)`: the two output groups use disjoint input
blocks, and each matrix has two independent columns. A nonzero star-pencil
slice has rank two over every field, including characteristic two. Thus a
source slice is a block diagonal sum of five such slices; setting any one
of the five parameter pairs to zero gives source slice rank at most eight.

Both source and target are concise on the first leg. A restriction therefore
induces an injective pullback from target first covectors into the source's
ten-dimensional first covector space: a covector whose pullback is zero
would annihilate every target coefficient. For every pulled-back covector,
target slice rank is bounded by source slice rank. This implication requires
no invertibility of the other two restriction maps.

It suffices to argue over an algebraic closure. Any hypothesized base-field
restriction remains an exact restriction after scalar extension; all linear
ranks and conciseness are preserved. Dimension arguments below refer to
algebraic dimension, so they do not confuse a finite-field point set with
the underlying determinantal variety.

### One gain

For m=1 the pullback image is a hyperplane W in the source space. Write its
one nonzero normal as five two-coordinate blocks, and choose a nonzero
block numbered ell. Choose distinct i,j different from ell. Then

```
W_i = W intersect {pair_i=0}, dimension 7,
W_j = W intersect {pair_j=0}, dimension 7,
dimension(W_i intersect W_j) = 5,
W_i + W_j = W.
```

The dimension statements hold because the normal still has a nonzero
ell-block after either or both selected pairs are set to zero. Their
preimages L_i,L_j in target first coordinates are seven-dimensional spaces
of slices of rank at most eight, and together span the target's full
nine-dimensional first covector space.

Let epsilon be the target scalar coordinate. A slice with epsilon nonzero
has rank `1+2 rank(V)+2 rank(H)`. Rank at most eight therefore implies
`rank(V)+rank(H)<=3`, which requires at least one of V,H to have rank at
most one. Each of these two rectangular rank-at-most-one loci has dimension
`4+2-1=5` in the eight-dimensional `(A,B)` space. Their union has dimension
five.

If a seven-dimensional linear L_i is not contained in `{epsilon=0}`, its
epsilon-one fiber is an affine six-dimensional space. Forgetting epsilon
embeds that fiber into `(A,B)` coordinates. It cannot lie in a union of
five-dimensional varieties. Thus L_i and L_j both lie in the eight-dimensional
hyperplane `{epsilon=0}`, contradicting that they span all nine dimensions.

### Two or more gains

For m=2, the target first dimension is ten and its pullback is an
isomorphism. Each source coordinate-pair-zero space V_i has dimension eight
and yields a distinct eight-dimensional preimage L_i of target slices of
rank at most eight.

The target scalar coordinate vector has two entries. If L_i has a nonzero
projection to these coordinates, fix any nonzero vector in its projected
image. The corresponding affine fiber has dimension at least six. At least
one scalar slice is nonzero, so the product slice must again have rank at
most six and `(A,B)` must lie in the same five-dimensional bad locus.
Projection to `(A,B)` is injective on this fiber, giving the same dimension
contradiction. Consequently every L_i lies in the eight-dimensional space
where both scalar coordinates vanish, and hence equals that space. This
would make all distinct V_i equal under the pullback isomorphism, impossible.

For m>=3 the target concise first dimension `8+m` exceeds the source's ten,
already contradicting first-leg flattening rank monotonicity.

### Elementary finite-field alternative

For a field of order q, each rectangular bad locus has
`q^5+q^4-q` points. Their intersection consists of zero together with the
nonzero decomposable 2-by-2-by-2 coefficient cubes, whose count is

```
1 + (q^2-1)^3/(q-1)^2 = q^4+2q^3-2q.
```

The union therefore has `2q^5+q^4-2q^3` points. For q>=3 this is strictly
less than q^6. A six-dimensional affine fiber has q^6 distinct `(A,B)`
points, giving the same contradiction without algebraic geometry. For q=2
the count equals q^6, so this particular counting proof must first extend
the coefficient field to F_4. Exact maps persist under this extension.

This closes the empty-D star case over all fields. The source-rank argument
used the ten-dimensional first space of exactly five star blocks; adding
D changes that geometry. No conclusion about arbitrary nonempty D follows
from this proof.
