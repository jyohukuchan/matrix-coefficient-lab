# The singular star admits no zero-catalyst positive gain at d=2,k=5

For every field `K` and every positive natural number `m`, there is no exact
local restriction

```
m Unit ⊕ (M₂⊗S_star) ≤ 5 S_star,
```

where `S_star` is the concise `(2,3,3)` pencil

```
S_star(a,b) = [[0,a,b], [a,0,0], [b,0,0]].
```

This excludes the `D=0` singular-star search, including all possible mixing
of the local linear maps. A certificate with a nonzero fixed `D` is not
excluded by this proof. The main agent discovered the hyperplane argument;
this note independently audits its dimensions, slice ranks, field extension,
and finite counting variant. Another research agent independently audited
it. No Lean files were edited or built, and the theorem is not claimed as a
new kernel-checked Lean declaration.

## Exact slice ranks

Every nonzero star slice has rank two. For instance, if `a≠0` the first two
rows and columns have determinant `−a²`; if `b≠0` use the first and third.
Thus a slice of `5 S_star`, whose first-leg space is `K¹⁰` with five pairs
`(a_i,b_i)`, has rank twice the number of its nonzero pairs. In particular,
when any one pair is zero its rank is at most eight.

Use the tied-product coordinates from `pencil-frontier.md` for the product:

```
H(A,B): (X₀,X₁,X₂) ↦ (A X₀, B X₀, A X₁+B X₂),
A,B,X₀,X₁,X₂ ∈ Mat₂(K).
```

Its first leg is `K⁸` and other legs are `K¹²`. Exactly,

```
rank H(A,B) = 2 rank([A;B]) + 2 rank([A B]) ≤ 8.
```

Consequently the product slice rank is even. It is at most six precisely
when the stacked `4×2` matrix or horizontal `2×4` matrix has rank at most
one. Denote that union of pair loci by `B_bad⊆K⁸`.

For `Unit⊕(M₂⊗S_star)`, write the first input as `(ε,A,B)`. Its slice rank is
`1_(ε≠0)+rank H(A,B)`. Hence its rank-at-most-eight locus is exactly

```
{ε=0}  ∪  { (ε,A,B): (A,B)∈B_bad }.
```

The integer counts above are dimensions or numbers of nonzero independent
scalar blocks, not casts of those integers into the field.

## A large low-rank subspace must have zero unit coordinate

Over an algebraic closure, each rank-at-most-one rectangular matrix locus
has algebraic dimension five. An elementary chart verifies this: at a
specified nonzero entry of a `4×2` rank-one matrix, the pivot entry, the other
three entries in its column, and the one other entry in its row determine
all remaining entries. There are five parameters. Finitely many such charts
and the zero matrix cover the locus. The transposed `2×4` locus has the same
parameter count. Their union `B_bad` therefore has dimension five.

Let `L⊆K⊕K⁸` be a seven-dimensional linear subspace all of whose target
slices have rank at most eight. If the unit-coordinate functional is nonzero
on `L`, its fiber at `ε=1` is a six-dimensional affine space, and projection
to `(A,B)` embeds that fiber as a six-dimensional affine subspace of `K⁸`.
Every one of its pairs must belong to `B_bad`, which has dimension five.
This is impossible. Thus

```
L ⊆ {ε=0}.
```

For a certificate over an arbitrary field, extend all its coefficient maps
to an algebraic closure. The same exact tensor identity and the same ranks
of the fixed linear maps persist. Therefore the preceding dimension
argument applies regardless of the original field's characteristic or
cardinality. It uses standard algebraic dimension, rather than a new Lean
formalization of that dimension argument. The fully finite alternative
below avoids algebraic geometry when the original field is finite.

## The hyperplane contradiction for one gain

Suppose `Unit⊕(M₂⊗S_star)≤5S_star` were an exact restriction. In bilinear-map
orientation, its first-leg map is an injection

```
F: K⁹ → K¹⁰.
```

Injectivity follows from the target's first-leg conciseness: if `F(x)=0`,
the entire target slice at `x` is zero, forcing `x=0`. Thus its image is a
hyperplane `H=ker ℓ` for some nonzero linear functional `ℓ` on the five source
pairs.

Choose one pair index `h` on which `ℓ` has a nonzero component. Choose
**distinct** other indices `i,j≠h`. Let `V_i⊆K¹⁰` be the subspace where pair
`i` is zero, and similarly for `V_j`. The surviving `h` component ensures
that the hyperplane equation remains nonzero on `V_i`, `V_j`, and their
intersection. Therefore

```
dim(H∩V_i) = dim(H∩V_j) = 7,
dim(H∩V_i∩V_j) = 5,
(H∩V_i)+(H∩V_j) = H.
```

Every source slice on `V_i` or `V_j` has rank at most eight. Restricting the
other input and applying the output map cannot increase matrix rank. Hence
the two transported subspaces

```
L_i = F⁻¹(H∩V_i),    L_j = F⁻¹(H∩V_j)
```

are seven-dimensional target-first subspaces whose slices all have rank at
most eight. By the preceding lemma, both lie in `{ε=0}`, a space of
dimension eight. But their span is all of `K⁹`, because their source spans
are all of `H`. This is a contradiction.

Only first-leg injectivity and slice-rank monotonicity were used. The other
leg maps may be arbitrary and may mix every block. No output isomorphism or
tensor-rank additivity is assumed.

## Two gains and larger gains: supplementary checks

The one-gain result alone excludes every positive `m`: discard `m−1` gain
blocks from the target to obtain a one-gain restriction with the same
source. The following checks also validate the dimension boundary directly.

For two gains, target and source first-leg dimensions are both ten, so the
first-leg injection `F` is an isomorphism. Each source `V_i` has dimension
eight and slice rank at most eight. Its preimage `L_i` is therefore an
eight-dimensional low-rank target-first subspace.

Let the projection of any such `L` to its two unit coordinates have rank
`r`. If `r>0`, choose any nonzero unit vector in its image. The fiber over
that vector is affine of dimension `8−r≥6` and embeds into the pair space.
There is at least one nonzero unit contribution to the slice rank. Since
the product contribution is even, rank at most eight forces it to be at
most six; thus every pair in that affine fiber lies in `B_bad`, again
impossible. Therefore `r=0`, and `L` must equal the eight-dimensional
units-zero subspace. Under the isomorphism `F`, all five different `V_i`
would then have the same preimage, which is impossible.

For `m≥3`, the target's first flattening rank is `8+m≥11`, while the source's
is ten. Flattening rank cannot increase under a restriction, so such a map
is already excluded without a slice argument.

## A finite counting proof

For a finite field of size `q`, the number of `4×2` matrices of rank at most
one is

```
1 + (q⁴−1)(q²−1)/(q−1) = q⁵+q⁴−q.
```

The transposed locus has the same size. Their intersection in the pair
space consists exactly of rank-at-most-one `2×2×2` coefficient cubes:
`(A,B)` is a pair of scalar multiples of one rank-one matrix. Including zero,
its size is

```
1 + (q²−1)³/(q−1)² = q⁴+2q³−2q.
```

For nonzero pairs, the two rank conditions force both a common left factor
and a common right factor, giving the displayed three-factor cube; cases
where just one matrix is nonzero are included. Thus inclusion-exclusion gives

```
|B_bad| = 2q⁵+q⁴−2q³ < q⁶    for every q≥3.
```

Indeed the difference is `q³(q³−2q²−q+2)=q³(q−2)(q²−1)>0`.

A six-dimensional affine pair space contains exactly `q⁶` points and cannot
fit inside `B_bad`. This proves the low-rank subspace lemmas above directly
over every finite field with `q≥3`. For `F₂`, extend a purported exact
coefficient certificate to `F₄` and apply the same argument. A coefficient
identity persists under this scalar extension; an identity only on binary
input points would not justify it. This distinction is essential to the
certificate model.

The bad-pair counts for `q=2,3,4` are respectively `64,513,2176`, compared
with `q⁶=64,729,4096`. Equality at `q=2` explains why the counting argument
needs extension there. The earlier complete slice inventories over `F₂`
and `F₃` independently give bad counts `1+27+36=64` and
`1+128+384=513`. The main agent owns the persistent exact experiment
`star_slice_obstruction_experiment.py`; this report does not duplicate that
implementation.

## Scope of the result

The normal-rank inequality alone allowed this candidate: ten source versus
nine target ranks for one gain. The new obstruction uses the arrangement of
large low-rank first-leg subspaces, not just their maximum rank. It closes
all exact `D=0` positive-gain star searches at `d=2,k=5` over every field.

A nonzero `D` changes both the first-leg ambient space and the low-rank
arrangement. The hyperplane and span dimensions above then no longer apply.
This proof therefore supplies no exclusion of a star auxiliary with an
arbitrary nonzero catalyst. It also supplies no positive-gain witness for
another auxiliary. The geometric `D=0` existence result in
`geometric-catalyst.md` concerns much larger, different auxiliaries and is
consistent with this exclusion.
