# Necessary conditions for a nonzero catalyst with the singular star

Let `S=S_star` be the concise `(2,3,3)` symmetric star pencil, and let
`T=M₂⊗S`. For any field `K`, suppose an exact coefficient restriction exists:

```
D ⊕ m Unit ⊕ T ≤ D ⊕ 5S,       m≥1.
```

Then the following are necessary:

1. The first flattening rank of `D` is at least three.
2. `D` has a two-dimensional linear space of first slices defined over `K`
   on which every nonzero slice, even after algebraic-closure extension, has
   matrix rank at least four.
3. Its other two flattening ranks are at least four, and they cannot both be
   four.
4. If its first flattening rank is three, its generic first-slice matrix
   rank is at least five; thus both other flattening ranks are at least five.
5. Its exact tensor rank over `F₂` or `F₃` is at least six; over every other
   field it is at least five.

Consequently the smallest remaining componentwise boundaries for the concise
shape of `D` are `(3,5,5)`, `(4,4,5)`, and `(4,5,4)`. They are **not**
catalyst witnesses.
The existence of a star certificate with arbitrary nonzero `D` remains
unresolved.

The main agent derived the argument; the existence-audit and
constructive-research agents independently checked the rank, dimension,
annihilator, and concision steps. This is a mathematical proof audit, not a
new kernel-checked Lean theorem. No Lean sources were edited or built.
Bounded SAT timeouts for diagonal catalysts are separate search outcomes;
they are not evidence proving the conditions or impossibility for all `D`.

## Reductions and slice ranks

For `m>1`, discard `m−1` output unit blocks. Therefore it suffices to derive
necessary conditions from the one-gain restriction

```
D ⊕ Unit ⊕ T ≤ D ⊕ 5S.
```

Work over an algebraic closure of `K`. An exact coefficient identity extends
to this field unchanged, and flattening ranks and ranks of the fixed local
linear maps do not change under scalar extension. The algebraic closure is
infinite, including in positive characteristic. All geometric arguments
below take place there. An equality only on the original finite-field input
points would not justify this step.

Let `h` be the first flattening rank of `D`. Replace its first leg by its
concise support of dimension `h`, leaving the other legs unchanged. More
explicitly, the map from first covectors to coefficient slice matrices has
an `h`-dimensional image. Choose a basis of that image: source slices factor
through the corresponding first-coordinate quotient, and target slices can
be lifted along a linear section. Composing the original restriction maps
with this quotient and section gives the same restriction for the first-
concise `D`. The set of slice matrices is unchanged. Thus this reduction
does not alter either the certificate implication or generic slice rank.

Let `rho` be the generic first-slice matrix rank of `D`: the rank over the
rational function field in its first variables, equivalently the maximum
rank over the algebraic closure. Every specialization has rank at most
`rho`, and some specialization achieves `rho`, since a nonzero maximal
minor polynomial is nonzero at some point of this infinite field.

The star slices are

```
S(a,b) = [[0,a,b], [a,0,0], [b,0,0]].
```

They have rank two for `(a,b)≠0` and rank zero otherwise, in every
characteristic. Therefore the source `5S` has five first-coordinate pairs
and a slice with any one pair zero has rank at most eight.

The product admits the tied-product coordinates

```
H(A,B): (X₀,X₁,X₂) ↦ (A X₀, B X₀, A X₁+B X₂),
A,B,X₀,X₁,X₂ ∈ Mat₂(K).
```

Its first dimension is eight; the other dimensions are twelve. Exactly,

```
rank H(A,B) = 2 rank([A;B]) + 2 rank([A B]).
```

The two columns of each matrix are independent copies of the same operation;
the first two output blocks depend only on `X₀`, the last only on `X₁,X₂`.
This proves the formula. The rank is even, at most eight, and at least four
whenever `(A,B)≠0`.

Its rank-at-most-six locus is

```
B_bad = {rank([A;B])≤1} ∪ {rank([A B])≤1} ⊆ K⁸.
```

Each part has algebraic dimension five. In a rank-one `4×2` matrix with a
specified nonzero pivot, the pivot, its other three column entries, and its
one other row entry determine the whole matrix: five parameters. Finitely
many pivot charts plus the zero matrix cover the locus. The transposed
`2×4` locus has the same dimension. Thus `B_bad` also has dimension five
and cannot contain a six-dimensional affine space. These are the same
slice-rank and dimension facts used in `star-zero-catalyst.md`.

## The source hyperplane and the projection lemma

In bilinear-map orientation, the first local map pulls target first inputs
back into source first inputs. Both `D` first supports have been made
concise, so the target first space has dimension `h+9` and source first
space dimension `h+10`:

```
V = K^h ⊕ K(gain) ⊕ K⁸(A,B),
W = K^h(D-source) ⊕ E₁ ⊕ ... ⊕ E₅,       E_i=K².
```

The first pullback `F:V→W` is injective. Indeed `F(v)=0` would make the entire
target slice zero. The target first flattening is concise: the independent
`D`, unit, and product blocks have first flattening ranks `h`, one, and
eight. Hence its zero slice forces `v=0`.

Identify `V` with the hyperplane `H=im(F)=ker(eta)` in `W`, for one nonzero
normal functional `eta`. For each pair index `i`, set

```
V_i = {w∈W: pair i is zero},
H_i = H∩V_i,                 dim H_i≥h+7.
```

Let `Y=K^h⊕K(gain)` and let `pi:H→Y` be the target `D`/gain projection via
`F⁻¹`. This projection is surjective on all of `H`.

**Projection lemma.** Its restriction to any `H_i` is not surjective.

To prove this, assume it is surjective. Fix a target `D` first input where
its slice has rank `rho`, and fix the gain coordinate to one. The affine
fiber in `H_i` over these fixed values has dimension

```
dim H_i−dim Y ≥ (h+7)−(h+1) = 6.
```

Under `F⁻¹`, its remaining `(A,B)` coordinates embed it as an affine space
of the same dimension in `K⁸`. On `H_i`, the source slice has rank at most
`rho+8`: the source `D` contribution is at most `rho` and one star pair is
zero. Restricting the other input and applying the output map cannot increase
matrix rank. The target slice on the selected fiber has rank

```
rho + 1 + rank H(A,B) ≤ rho+8.
```

Therefore its even product rank is at most six, so the entire affine fiber
lies in `B_bad`. This contradicts its dimension five. The argument allows
arbitrary mixing of all local maps and needs no output isomorphism.

## Four independent annihilators

For each `i`, choose a nonzero functional `ell_i∈Y*` annihilating `pi(H_i)`;
the projection lemma guarantees one. Since `pi:H→Y` is surjective, its
pullback `pi*(ell_i)` is a nonzero element of `H*`.

The subspace `H_i` is the common kernel of pair `i`'s two source coordinate
forms restricted to `H`. Its annihilator is exactly the span of those two
restrictions. Thus `pi*(ell_i)` has a representative supported on source
pair `i` in

```
H* ≅ W*/span(eta).
```

If `eta` has a nonzero component on a star pair `r`, choose the four other
indices. Their coordinate spans are independent in this quotient: a sum
supported outside pair `r` cannot be a nonzero scalar multiple of `eta`,
whose `r` component is nonzero. The representatives also have disjoint pair
supports, so one nonzero vector from each gives four independent pullbacks.
If `eta` has no star component, its nonzero component lies in `D`; all five
pair spans are independent modulo `eta`, and any four suffice.

Consequently there are four independent `ell_i` in `Y*`. Since
`dim Y=h+1`, this gives

```
h+1≥4,        h≥3.
```

In the case where `eta` has only a `D` component, all five annihilators are
independent, giving the supplementary stronger bound `h≥4` for that case.

## A constant lower-rank two-plane in D, for every h

Set

```
U = H∩K^h(D-source),         dim U≥h−1.
```

Every source pair coordinate vanishes on `U`. The four independent
annihilators above therefore all vanish on `pi(U)`. It follows that

```
dim pi(U) ≤ (h+1)−4 = h−3,
dim ker(pi|U) ≥ (h−1)−(h−3) = 2.
```

Let `L=ker(pi|U)` and take any nonzero `v∈L`. Its target first input `F⁻¹(v)`
is nonzero, but its target `D` and gain coordinates vanish. Hence `(A,B)≠0`,
and the target slice has rank at least four by the product formula. The
source slice at `v` is supported only on `D`; no star coordinates contribute.
Matrix-rank monotonicity therefore yields

```
rank D(v) ≥ rank H(A,B) ≥ 4,        every nonzero v∈L.
```

Choose any two-dimensional subspace of `L`. This establishes the stated
constant lower-rank two-plane for **all** `h`, not just `h=3`, and in
particular `rho≥4` for all `h`.

The two-plane can be chosen over the original field `K`. Before extending,
`F,H,pi,H_i` are all `K`-linear. A fixed linear map's rank is unchanged by
scalar extension, so the projection lemma's non-surjectivity gives
annihilators over `K` by ordinary linear algebra. Choose the four independent
annihilators there, and form their `K`-linear kernel `L`. Its extension is
exactly the kernel used above. Thus `D` must in particular have a rank-at-
least-four slice at an original-field point. This last fact is a stronger
corollary than merely knowing its generic rank.

Each other flattening rank of `D` must now be at least four. Concising its
other legs preserves the ranks of its slice matrices; if either support
had dimension at most three, no slice could have rank four.

If both other concise dimensions were four, restrict `D` to two independent
vectors `u,v` in this two-plane. Every nonzero combination would be an
invertible `4×4` matrix. But

```
f(x,y) = det(xD(u)+yD(v))
```

is a homogeneous degree-four polynomial. If identically zero, it immediately
provides singular nonzero slices; otherwise it has a projective zero over
the algebraic closure. Equivalently, `det(D(u)+tD(v))` has nonzero leading
coefficient `det D(v)`, positive degree four, and hence a root. The
corresponding combination is nonzero since `u,v` are independent. Either
case contradicts its required rank four. This works in every characteristic.

Thus the global necessary concise shape bounds are first axis at least
three, other axes at least four, with the `(4,4)` other-axis boundary
excluded. The next section strengthens the other-axis bound when `h=3`;
the former `(3,4,5)` and `(3,5,4)` boundaries are consequently excluded.

## Strengthening at first dimension three: generic slice rank at least five

The independently audited rank-drop argument in
[star-convolution-catalysts.md](star-convolution-catalysts.md) closes the
former `(3,4,5)` and `(3,5,4)` boundaries. This section records its general
`h=3` consequence and the distinction from tensor-rank bounds.

Assume `h=3`. The slice-plane argument above already gives generic slice
matrix rank `rho≥4`. Suppose for contradiction that `rho=4`. Over the
algebraic closure, let

```
Delta = {delta∈Kbar³: rank D(delta)<4}.
```

This is a closed homogeneous cone, defined by all four-by-four minors.
At least one minor is a nonzero polynomial because the generic rank is
four. Every nonzero vector in the necessary two-plane `L` has rank at least
four; every slice has rank at most four. Thus

```
Delta∩L = {0}.
```

It follows that `dim Delta≤1`. Indeed a homogeneous cone component of
dimension at least two projectivizes to a positive-dimensional closed
variety in projective two-space. It meets the projective line `P(L)` by
the projective dimension inequality, contradicting the displayed
intersection. Equivalently, cutting such a cone component by the linear
hyperplane `L` lowers dimension by at most one and cannot leave only the
origin. This is an algebraic-closure argument; a point sample only over the
original finite field would not establish this intersection property.

It remains to check why a rank-drop cone of dimension at most one excludes
the certificate in this first dimension. Use the target-first coordinates
`(delta,epsilon,A,B)` and the two-row source-pair maps `Q_i` obtained from
`F`. The target first space has dimension twelve and

```
J_i = ker Q_i,       dim J_i≥10.
```

On `J_i`, the source slice has rank at most `4+8=12`; hence the target slice
also has rank at most twelve. Exactly, its low-rank locus is the finite union

```
E       = {epsilon=0},                         dim E=11,
Z_D     = Delta×K×K⁸,                         dim Z_D≤10,
Z_core  = K³×K×B_bad,                         dim Z_core=9.
```

To verify the union, outside all three sets the `D`, gain, and product
contributions have ranks four, one, and eight, totaling thirteen. Inside
`E` the gain vanishes; inside `Z_D` the `D` rank drops by at least one;
inside `Z_core` the even product rank drops by at least two. Each exception
therefore gives rank at most twelve.

Each `J_i` is an irreducible linear space contained in that union of closed
sets. It must be contained in one member and cannot lie in the smaller
nine-dimensional `Z_core`. Call it type A if contained in `E`, and otherwise
type B, which must be contained in `Z_D`.

For type A, the gain-coordinate functional belongs to the row space of
`Q_i`. Since that row space has dimension at most two, it adds at most one
other form beyond the common gain form.

For type B, the projection of `J_i` to the three `D` coordinates has
dimension at most one, since its image lies in `Delta`. Its kernel under
that projection has dimension at least nine. The entire target-`D`-zero
subspace has dimension nine, so that kernel is the full subspace; moreover
`dim J_i=10`. Hence every row of `Q_i` is supported on the three target-`D`
coordinates. All type-B row spaces together lie in one three-dimensional
space, regardless of how many components `Delta` has or whether those
components are linear.

If all five indices are type A, the star-pair rows together span at most
`1+5=6` dimensions, and the three source-`D` rows bring the first pullback
rank to at most nine. If at least one is type B, there are at most four
type-A pairs. All pair rows together span at most `3+1+4=8` dimensions;
including the source-`D` rows gives rank at most eleven. Both contradict the
injective rank twelve of `F`.

Thus `h=3,rho=4` is impossible. Combined with `rho≥4`, the strengthened
necessary condition is

```
h=3 ⇒ generic first-slice matrix rank rho≥5.
```

It forces both other flattening ranks to be at least five. For `h≥4`, the
all-h two-plane condition still requires the other ranks at least four and
excludes the `(4,4)` pair; the above row-rank contradiction does not
automatically extend to `h≥4`. The revised minimal componentwise boundaries
are consequently `(3,5,5)`, `(4,4,5)`, and `(4,5,4)`.

This is a **matrix slice rank** constraint. A six-term tensor-rank scheme
for a tensor of shape `(3,4,5)` does not evade it: every slice of that shape
has matrix rank at most four, so the whole shape is excluded. Conversely a
generic-rank-five tensor may drop to rank four on the necessary plane; this
argument does not force generic rank six or supply any witness in the new
boundaries.

## A necessary tensor-rank bound over the same field

The two-plane is defined over `K`, which also yields an elementary lower
bound on the **exact tensor rank of `D` over K**. Take any supplied rank
upper certificate with `R` simple terms. Contracting its first factors with
`v` in the two-plane gives

```
D(v) = Σ_{r=1}^R ell_r(v) b_r c_rᵀ,
```

where each restricted `ell_r` is a `K`-linear functional on that plane.
Each active term has matrix rank at most one. Thus the required inequality
`rank D(v)≥4` for every nonzero `v` implies at least four nonzero evaluated
coefficient forms at every such vector.

For `K=F_q`, a nonzero linear form on the two-dimensional plane is nonzero
at exactly `q²−q` vectors; an identically zero form contributes none. Summing
the active-term counts over all `q²−1` nonzero vectors gives

```
4(q²−1) ≤ Σ_{v≠0} #{r: ell_r(v)≠0} ≤ R(q²−q).
```

Consequently every exact rank upper certificate must satisfy

```
R ≥ ceil(4(q²−1)/(q²−q)) = 4+ceil(4/q).
```

In particular,

- over `F₂` or `F₃`, the exact tensor rank of `D` is at least six;
- over every finite field of size at least four, it is at least five.

The forms and the decomposition must be over the same field used for this
count. A decomposition available only after field extension is not an
original-field rank upper certificate.

Over every field, finite or infinite, a direct argument already gives
rank at least five. Three or fewer terms cannot produce one matrix slice of
rank four. If there were four terms, restrict any one coefficient form to
the two-plane. If it is zero, at most three terms are active everywhere; if
it is nonzero, its kernel contains a nonzero vector, at which at most three
terms are active. Either option contradicts the slice rank lower bound.
Thus the finite counting argument strengthens this universal bound for
`F₂` and `F₃`.

This applies to arbitrary `D`, not just diagonal or algebra-multiplication
catalysts. It assumes no additivity property of tensor rank and uses no
external coding-theory theorem. The `K`-defined two-plane, rather than an
unrelated generic-rank estimate or a plane existing only after extension,
is the essential input.

## Generic rank must not be replaced by finite-point maximum

The projection proof uses a target `D` slice of **generic** rank `rho`, after
scalar extension. Enumeration only over the base finite field can miss that
rank. A concrete binary coefficient tensor of concise shape `(3,4,5)` is

```
D(x,y,z) = [[x, 0, 0, 0, y+z],
            [0, y, 0, 0, x+z],
            [0, 0, x+y, 0, x+y],
            [0, 0, 0, z, x]].
```

Its first four columns have determinant `x·y·(x+y)·z`, a nonzero formal
polynomial. Its generic slice rank is therefore four. Exact finite matrix
checks give flattening ranks `(3,4,5)` and these slice-rank histograms:

| Field | rank 0 | rank 2 | rank 3 | rank 4 |
| --- | --- | --- | --- | --- |
| `F₂` | 1 | 1 | 6 | 0 |
| `F₄=F₂[a]/(a²+a+1)` | 1 | 3 | 24 | 36 |

Thus a maximum over binary inputs is three, despite generic rank four.
It would be incorrect to assign generic rank three to this tensor from
that enumeration. The example is a generic-rank computation control, not a
catalyst witness.

The strengthened two-plane condition changes the rejection question: since
its required two-plane is defined over the base field, a maximum base-field
slice rank below four **does** exclude a base-field certificate. Accordingly,
this particular tensor is excluded over `F₂` by the two-plane condition,
despite passing the earlier all-h generic-rank threshold of four. Its
shape is now independently excluded over every field by the strengthened
`h=3` threshold of five. It would be misleading to call
it a counterexample to that strengthened exclusion filter. Computing generic
rank and applying a base-point lower-rank consequence are different justified
operations; the projection lemma itself still needs generic rank.

For finite fields, the projection lemma can also be justified by extending
to a finite field large enough to contain a target `D` input achieving
`rho`, and of size `q≥3`. The exact bad-pair count from
`star-zero-catalyst.md` is

```
|B_bad| = 2q⁵+q⁴−2q³ < q⁶.
```

The affine fiber contains at least `q⁶` pair points, a contradiction. A
rank-`rho` algebraic-closure input uses only finitely many algebraic
coordinates, so it lies in one finite extension of the original finite
field. Further extension can ensure `q≥3`. This variant keeps the same exact
coefficient-identity requirement and does not license checking identities
only on base-field input points.

## Remaining research boundary

The conditions exclude all `D` with first flattening rank below three, all
`D` whose second or third flattening rank is below four, and all `D` whose
last two concise dimensions are both four. They additionally require the
constant lower-rank two-plane, which is stronger than merely finding one
rank-four slice.

The strengthened `h=3` condition additionally excludes every tensor of first
flattening rank three with either other flattening rank at most four.
These conditions neither supply a certificate in the revised remaining
shapes nor exclude every nonzero `D`. Larger searches and earlier diagonal-catalyst SAT timeouts
should be reported separately from these proved structural restrictions.
The zero-catalyst exclusion and the complete-in-principle geometric
auxiliary family concern different cases and remain consistent with this
result.
