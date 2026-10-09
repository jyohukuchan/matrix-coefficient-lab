# Variable auxiliaries: calibration and complete-convolution exclusion

This independent mathematical audit concerns actual exact restrictions

```
D ⊕ m Unit ⊕ (M₂ ⊗ S) ≤ D ⊕ 5S,     m ≥ 1.
```

The principal new result is that **every complete convolution `C(a,b)`, and
any direct sum of complete convolutions sharing their polynomial-output leg,
is impossible as `S`, for every finite `D` and every coefficient field**.
This includes `C(2,4)` and `C(3,3)`, whose asymptotic rank five had previously
made them reasonable calibration candidates. The proof strengthens the
[earlier commutator filter](catalyst-obstructions.md) by retaining the output
coordinates outside an identity slice. It uses elementary coefficient
elimination and the characteristic-free commutator bound, rather than the
external algebra or rectangular lower-bound theorems. It has not been
formalized in Lean.

The proof's retained three-sector tensors differ from complete convolutions.
The smallest retained tensor `R(2,1)` survives the rank-cost bounds: its
asymptotic rank is five, but its supplied exact rank upper bound is six.
The final slice-space bounds proved below require any `R(2,1)` catalyst
to have first flattening rank at least five, generic slice rank at least
eight, and other flattening ranks at least eight. For the fixed star the
first flattening rank must be at least six, generic slice rank at least six,
and other flattening ranks at least six. Earlier first-dimension-three/four
refinements remain in the proof record but are superseded by these general
exclusions. A
mixture of three length-two dot tensors in all three singleton-leg orientations
also survives these necessary bounds. Neither is a catalyst witness.

## What the existing existence argument actually supplies

The source theorem
[`normalizedStates_nonempty_of_no_catalyst`](../../lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Spectrum/StateObstruction.lean)
assumes the absence of positive-gain certificates for **all** `D,s,m`.
The detector constraints range over every test tensor `s`. The subsequent
multiplicative-state argument uses invariance under every normalized
multiplicative translate, which is precisely where this unrestricted family
matters. Restricting the tests to one chosen convolution, star, or finite
list does not preserve those hypotheses automatically.

Over an algebraically closed coefficient field, the `9/4` character bound
contradicts a global multiplicative state detecting `M₂` at value five. Its finite inconsistency therefore yields some finite
auxiliary, assembled as a nonnegative integral mixture of the detector tests
used by the finite certificate. It does not promise a prescribed auxiliary,
nor every auxiliary whose spectral inequalities happen to be favorable.
A finite witness over the algebraic closure uses only finitely many algebraic
coefficients, hence lies in a finite coefficient extension of the original
field. The separate geometric route below establishes a same-field existence
class from the all-fields exact exponent theorem.
A finite feasible LP state only rules out duals in its selected row family.

There is also an explicit sufficient variable class, already implemented in
[`build_geometric_catalyst`](../geometric_catalyst.py) and audited in
[geometric-catalyst.md](geometric-catalyst.md). From a supplied exact `r`-term
scheme for `M_(2^N)` with `r < 5^N`, it constructs

```
S_N = ⊕_(j=0..N−1) 5^(N−1−j) M_(2^j),
D = 0,     m = 5^N−r.
```

The common middle blocks in `5S_N` and `M₂⊗S_N` are matched by identity maps;
the `5^N` scalar source blocks supply the given top matrix scheme and the
remaining gain units. This is a finite restriction construction, with no
cancellation theorem. The same-field exact-rank exponent theorem supplies
existence of a finite top block after powering and restricting an infimum
witness, as detailed in that note. It does not identify useful coefficients
or a usable initial finite block. The geometric constructor does not replace
the missing discovery stage.

## Calibration before an exact-obstruction screen

For a normalized nonnegative, additive, multiplicative,
restriction-monotone character `χ`, an actual witness implies

```
m + χ(M₂) χ(S) ≤ 5 χ(S),
χ(M₂) ≤ 5 − m/χ(S).
```

If a proved bound on powered exact schemes gives `χ(S) ≤ A`, then
`χ(M₂) ≤ 5−m/A`. Using the detecting-character and exact-rank-to-arithmetic
bridges named in [the fixed-star calibration](star-convolution-catalysts.md),
this implies `ω ≤ log₂(5−m/A)`. Here `A` must be a proved character upper
bound; a failed search is not one.

For gain one, avoiding an accidentally stronger claim than `9/4` requires

```
A ≥ 1/(5−2^(9/4)) ≈ 4.11232334197.
```

Thus `A≤4` asks for a stronger result; `A=5` gives only the weaker numerical
upper bound `log₂(4.8)≈2.26303440583`. Passing this calibration says nothing
about finite catalyst existence. Complete convolutions provide an explicit
example: they can pass this numerical check and still be excluded below.

Exact rank and asymptotic rank have separate roles. A supplied exact scheme
with `R` terms gives the source bound `rank(D)+5qR` after **additive**
repetition. A bound `rank(S^n)≤poly(n)A^n` cannot replace `R` by `A` in that
calculation: additive repetition does not take tensor powers of `S`.

## Output elimination strengthens the commutator bound

Use matrix slices as maps from the second coordinate space `Y` to the third
coordinate space `Z`. Suppose `T` is concise on `Z`, `dim Y=N`, `dim Z=z`,
and three first-leg slices `A,B,C` satisfy:

- `A:Y→Z` is injective, with image `V` of dimension `N`;
- both `B` and `C` map all of `Y` into `V`;
- identifying `V` with `Y` using `A`, their normalized commutator has rank `c`.

Then

```
rank(T) ≥ z + ceil(c/2).
```

To prove it, take any `r`-term exact decomposition. Its third-factor vectors
span `Z`: otherwise the output flattening would not have rank `z`. Their
classes therefore span `Z/V`. Select `z−N` terms whose third-factor classes
are independent in this quotient, and let `K` be the span of their actual
third-factor vectors. Independence modulo `V` gives `K∩V=0` and
`Z=V⊕K`.

Project `Z` to `Z/K`. The selected `z−N` terms vanish, leaving a decomposition
with at most `r−(z−N)` terms. The quotient identifies `V` isomorphically with
`Z/K`, so `A` becomes an invertible `N×N` slice. Because the full images of
`B,C` were in `V`, their normalized operators and commutator rank are
preserved. The elementary commutator theorem then gives

```
r−(z−N) ≥ N + ceil(c/2).
```

Rearrangement proves the bound. The quotient depends on the proposed
rank decomposition; no fixed quotient is asserted to kill terms of every
scheme. There is no cancellation of tensor summands. The proof works over
every field. The commutator theorem can pass to a rational-function extension
to choose a generic pivot slice, as in the earlier audited proof; rank over
an extension never exceeds base-field rank, and no division by two occurs.

This bound applies directly to a whole repeated target. Output conciseness
is checked for that whole tensor, so selected decomposition terms may mix
all its blocks. The argument does not require a rank-additivity theorem or a
term decomposition respecting those blocks.

## Every complete convolution is excluded

For positive `a,b`, write `L=a+b−1`. The coefficient tensor `C(a,b)` is
polynomial multiplication of inputs of lengths `a,b`, with all `L` output
coefficients retained. All three flattenings are concise. Over an algebraic
closure, evaluation/interpolation supplies an exact `L`-term scheme, and the
output flattening gives rank at least `L`.

Consider

```
T_q = q copies of (m Unit ⊕ M₂⊗C(a,b)).
```

Its second coordinate dimension is `N=q(m+4b)`, and its concise output
dimension is `z=q(m+4L)`. Choose the first slice corresponding to the identity
matrix times the constant polynomial coefficient, and value one on each
scalar unit. This slice injects `Y` onto the subspace `V` consisting of the
scalar outputs and matrix-polynomial coefficients of degrees `0,…,b−1`.

Choose the other two slices as `E12` and `E21`, again at the constant
polynomial coefficient, and zero on the scalar blocks. Their entire images
lie in `V`. After normalization, their commutator is `diag(1,−1)` repeated
over the two matrix columns and `b` polynomial coefficients in each block.
Its rank is `4bq`, including characteristic two. The strengthened bound gives

```
rank(T_q) ≥ q(m+4L+2b).
```

Swapping the two input legs gives the additional bound with `2a`; thus the
lower bound is at least `q(m+4L+2 max(a,b))`.

Now suppose a fixed finite catalyst existed. Iterating the restriction
additively while retaining one fixed `D`, then projecting away the target
`D`, would give

```
T_q ≤ D ⊕ 5q C(a,b).
```

Extend all coefficients to an algebraic closure. The witness still holds;
`D` still has some fixed finite rank upper bound `C`; and the exact `L`-term
convolution scheme is now available. Consequently

```
rank(T_q) ≤ C+5qL.
```

But `2 max(a,b)≥a+b=L+1`, so the lower bound exceeds `5qL` by at least
`q(m+1)`. Choose `q` with `q(m+1)>C` to obtain a contradiction. This excludes
all positive complete convolutions, regardless of the original field or a
finite coefficient extension used by the hypothetical maps. Missing
interpolation nodes over a small original field cannot evade the argument.

For the two previously interesting candidates:

| Auxiliary | Exact rank over closure | Product-plus-one lower bound per repeated block | Source slope |
| --- | ---: | ---: | ---: |
| `C(2,4)` | 5 | `1+4·5+2·4=29` | 25 |
| `C(3,3)` | 5 | `1+4·5+2·3=27` | 25 |

These bounds supersede the earlier inconclusive convolution/algebra screens.

The asymptotic rank of `C(a,b)` over the original field is also `L`: the
output flattening is the lower bound, and a fixed finite extension carrying
an `L`-term scheme can be descended after all tensor powers with one fixed
coefficient-algebra overhead. That statement is useful for calibration, but
the exclusion uses the exact closure scheme in the additive repetition.

## Mixtures, diagonal tensors, and dot orientations

Let `S=⊕_(i=1..t) C(a_i,b_i)` with all polynomial outputs in the same tensor
leg, and put `R=Σ_i(a_i+b_i−1)`. The closure source scheme has `R` terms.
Apply the same whole-target identity-image argument simultaneously to its
components. It gives lower bounds

```
q(m+4R+2Σ_i b_i),      q(m+4R+2Σ_i a_i).
```

Since `Σ_i a_i+Σ_i b_i=R+t`, their maximum is at least
`q(m+5R+t)`, exceeding the source slope `5qR`. Hence this entire mixture
class is excluded for arbitrary `D`.

A diagonal auxiliary `n Unit` is a mixture of `C(1,1)` and is excluded for
every `n≥1`, even though a large `n` can pass the character-rank calibration.
Dot tensors with a singleton first or second leg are `C(1,r)` and `C(r,1)`
after coordinate identification. Thus any direct sum of dots using at most
two singleton-leg orientations has a common nonsingleton leg that can be
chosen as output, and is excluded by the mixture theorem.

Three orientations escape this particular common-output argument. For
example, the direct sum of one `dot₂` in each singleton-leg orientation has
shape `(5,5,5)` and a supplied exact six-term decomposition. Compressing one
side to retain a four-dimensional identity slice while preserving all five
output coordinates gives product rank lower bound `4·5+2·4=28`; after one
unit gain this is 29, versus source slope 30. This is survival of a necessary
screen. Its supplied asymptotic rank upper is six, and flattening lower is
five; no exact or asymptotic direct-sum rank equality is asserted.

The full outer-product family `M(a,1,c)` retains the previously documented
[rectangular necessary condition](catalyst-obstructions.md)
`ac≥m+4 max(a,c)` for `a,c≥2`, conditional on that note's external
rectangular lower theorem. In particular, `min(a,c)≤4` is excluded. Larger
outer products are not promised by the existing existence theorem. The new
complete-convolution result has no dependency on that external theorem.

## Retained proof sectors remain different candidates

Let `R(a,h)` be the retained tensor of
[`ThreeSectorConstruction`](../sectors.py), with `a,h≥1`. It shares its first
axis among three branches, each equivalent to `C(a,h)`, with the middle
branch's other legs exchanged. Its concise dimensions are

```
(a, 3h+a−1, L),      L=3h+2a−2.
```

The proof degeneration from `C(a,3h+a−1)` gives asymptotic rank at most `L`;
its output flattening has rank `L`, so its asymptotic rank is exactly `L`.
This remains valid over the original field by fixed coefficient-algebra
descent after tensor powers, if interpolation coefficients require extension.
Separately, concatenating exact closure convolution schemes for its three
branches gives an **exact** rank upper bound

```
R_exact ≤ 3(a+h−1).
```

The branches share a first axis, so neither tensor-rank additivity nor
direct-sum additivity of character values is assumed here.

Retain `h` second-leg coordinates from each branch, using middle coordinates
`h+(a−1)+r` for `0≤r<h`. This compression preserves all `L` output
coordinates and gives an injective first slice of rank `3h`, at first
coefficient `e0`. Tensoring with `M₂`, the same identity and off-diagonal
matrix slices yield output dimension `4L`, input dimension `12h`, and
commutator rank `12h`. Thus the whole repeated product-plus-gain target has
rank at least `q(m+4L+6h)`. Using the supplied exact closure upper bound,
a necessary condition is

```
5·3(a+h−1) ≥ m+4L+6h,
3h ≤ 7a−7−m.
```

For `a=2,m=1`, this excludes every `h≥3`; `h=2` is an equality boundary.
The useful small distinctions are:

| Auxiliary | Asymptotic rank | Supplied exact closure upper | Repeated product-plus-one lower | `5R_exact` |
| --- | ---: | ---: | ---: | ---: |
| `R(2,1)` | 5 | 6 | 27 | 30 |
| `R(2,2)` | 8 | 9 | 45 | 45 |
| `R(2,3)` | 11 | 12 | 63 | 60 |

Thus `R(2,1)` and `R(2,2)` survive this rank-cost screen, while
`R(2,3)` is excluded. Replacing
`R_exact=6` by the asymptotic value five would incorrectly exclude
`R(2,1)`. Conversely, its asymptotic value five is the relevant value for
checking that a gain-one witness would not automatically beat `9/4`.

Mixtures containing retained sectors or differently oriented tensors require
new checks. An excluded component combined with a surviving component need
not remain excluded: the source and lower-bound margins can offset. No
closure of the exclusion under arbitrary mixtures is asserted.

## The smallest retained sector requires a nonzero catalyst

The root coordinator subsequently derived a slice-space obstruction for
`R(2,1)`, independently checked in this audit. Its pencil consists of two
column blocks `(x,y)^T` and one row block `(y,x)`, on disjoint other-leg
coordinates. Every nonzero slice has rank three.

For its product with `M₂`, let `A,B` be the two actual `2×2` matrix variables
on the first leg. Each matrix column is an independent spectator copy.
Consequently the product slice rank is

```
4 rank([A;B]) + 2 rank([A B]).
```

It is generically twelve. Its rank-at-most-eleven locus is exactly
`{rank([A;B])≤1} ∪ {rank([A B])≤1}` in the eight-dimensional pair space.
Each component has dimension five: a rank-one `4×2` matrix has five pivot-chart
parameters, as does its transpose. Thus this locus cannot contain a
six-dimensional affine space. These coefficients and dimensions hold in
every characteristic.

Discard extra gains and let `h` be the first flattening rank of `D`, after
making its first leg concise. Let `ρ` be its generic first-slice rank over
an algebraic closure. The first pullback is injective from the concise
target first space of dimension `h+9` into the source first space of dimension
`h+10`. Its image is a hyperplane `H` in
`K^h ⊕ E₁ ⊕ ... ⊕ E₅`, with each source pair space `E_i=K²`.

Put `H_i=H∩{pair i=0}`; its dimension is at least `h+7`. Project `H` via
the inverse first pullback to the target `D` first coordinates and its gain
coordinate, a space `Y` of dimension `h+1`. If this projection were surjective
on `H_i`, fix a target `D` slice of rank `ρ` and gain value one. The remaining
core pair variables form an affine fiber of dimension at least six. The
source slice rank on `H_i` is at most `ρ+12`, since four retained-sector
pairs remain. Slice-rank monotonicity then forces every product slice on
this affine fiber to have rank at most eleven. This contradicts the
five-dimensional bad-pair locus. Thus every projection image of `H_i`
is a proper subspace of `Y`.

Choose a nonzero annihilating functional for each image. Its nonzero pullback
to `H` is represented by forms supported on source pair `i`, modulo the
hyperplane normal. If that normal has a component on pair `r`, the other
four pair supports are independent modulo the normal; if the normal has
only a `D` component, all five are independent. At least four independent
functionals therefore lie in `Y*`, giving

```
h+1 ≥ 4,      h ≥ 3.
```

This excludes `D=0` and catalysts with first flattening rank below three.
The maps on the other two legs may mix all blocks; only slice-rank
monotonicity was used. The argument parallels the projection lemma and
annihilator proof in [star-nonzero-catalyst.md](star-nonzero-catalyst.md), with
source ranks `3/12` and product generic rank twelve replacing `2/8` and eight.
The additional transferred bounds below are proved with the new minimum
core rank six. Bounds involving the star's other rank formulas are not
assumed without checking them.

## Further necessary catalyst bounds for `R(2,1)`

These refinements were derived by the root coordinator and independently
audited here. They concern an actual certificate with this fixed retained
auxiliary, rather than a generic family of proof-sector tensors.

### An initial lower-rank two-plane for every first dimension

Keep the hyperplane `H`, the target `D`/gain projection `π`, and the four
independent annihilators from the preceding argument. Put

```
U = H ∩ K^h(D-source),      dim U ≥ h−1.
```

All five source pair coordinates vanish on `U`, so the four independent
annihilators vanish on `π(U)`. Thus `dim π(U)≤h−3`, and
`L=ker(π|U)` has dimension at least two. Every nonzero `v∈L` corresponds
to a nonzero target first input with only core coordinates. For a nonzero
pair `(A,B)`, both stacked and horizontal matrices have rank at least one,
so its product slice rank is at least `4+2=6`. The source slice at `v`
contains only `D`, giving

```
rank D(v) ≥ 6          for every nonzero v∈L.
```

All maps, annihilators, and kernels can be taken over the actual certificate
coefficient field `K`. Extension to an algebraic closure preserves their
linear ranks and this nonzero-slice inequality. Consequently `D` has a
`K`-defined two-plane with every nonzero slice of rank at least six even
geometrically. Its second and third flattening ranks are at least six.
They cannot both be six: on such a plane, two independent `6×6` slices
would have every nonzero combination invertible, but their degree-six
homogeneous determinant has a projective zero over the algebraic closure.

For a finite coefficient field `K=F_q`, take any `r`-term rank decomposition
and restrict its first coefficient forms to this plane. Each nonzero form
is active at `q²−q` vectors. At each of the `q²−1` nonzero vectors at least
six terms must be active, since each active term has matrix rank at most
one. Therefore

```
r(q²−q) ≥ 6(q²−1),
rank_K(D) ≥ 6+ceil(6/q).
```

This means at least nine over `F₂`, eight over `F₃,F₄,F₅`, and seven over
larger finite fields. Over any infinite field the elementary bound is also
at least seven: a decomposition with at most six terms has some nonzero
restricted coefficient form, and its nonzero kernel vector would leave at
most five active terms. If witness maps are supplied only over an extension
`E`, the plane and finite-field count are over `E`; its cardinality must not
be replaced by the original base-field cardinality.

### The source-normal cases strengthen this to a three-plane globally

The preceding two-plane estimate leaves one dimension unused. If the source
hyperplane normal has a nonzero `D` component, all five pair-supported
annihilators are independent modulo that normal, not merely four. Their
span in the target `D`/gain dual has dimension five. Then

```
dim U ≥ h−1,      dim π(U)≤h−4,      dim ker(π|U)≥3.
```

If the normal has no `D` component, the whole source `D` space lies in the
pullback image, so `dim U=h`. Four independent annihilators now give
`dim π(U)≤h−3`, again leaving a kernel of dimension at least three.
These are all possible normal cases. Thus there is a coefficient-field-defined
three-dimensional source `D` space mapping injectively to pure core, with
rank at least six at every nonzero slice, including geometrically.

If the generic slice rank of `D` were at most seven, all nonzero slices of
that core three-space would have rank exactly six: the only positive core
ranks are six, eight, ten, and twelve. It would lie in the rank-six simple
`2×2×2` tensor cone, whose maximum linear dimension is two. Consequently

```
ρ≥8        for every catalyst first dimension h.
```

Both other flattening ranks are at least eight. The temporary earlier
first-dimension-four/five generic-rank-seven bounds remain true but are
superseded by this global conclusion.

Counting active first coefficient forms on the three-plane over `F_q`
gives

```
r(q³−q²)≥6(q³−1),
rank_K(D) ≥ max(8, ceil(6(1+1/q+1/q²))).
```

Thus the global same-field tensor-rank lower bound is eleven over `F₂`,
nine over `F₃`, and eight over all other fields. The first-dimension-three
refinement below raises the generic rank and the general-field tensor-rank
lower bound to nine there. As throughout, the finite-field cardinality is
that of the actual certificate coefficient field.

### First flattening rank three forces generic slice rank at least eight

Assume `h=3`. Write the twelve-dimensional target first space as target
`D` coordinates, gain, and the eight core coordinates. Let `Y*` be the
four-dimensional span of the target `D` and gain coordinate forms. Let `W_i`
be the row space of the two pulled-back source pair coordinates, `P=Σ_iW_i`,
and `Z` the row space of the three pulled-back source `D` coordinates.

Let `ρ` be the generic slice rank of `D`, and let `Δ⊂Kbar³` be its proper
closed generic-rank-drop cone. Each `J_i=ker W_i` has dimension at least ten.
On it the source slice rank is at most `ρ+12`, so its target slices lie in

```
{gain=0}
  ∪ (Δ × Kbar(gain) × Kbar⁸(core))
  ∪ (Kbar³(D) × Kbar(gain) × B_bad).
```

Here `B_bad` is the five-dimensional rank-at-most-eleven core locus. The
last component has dimension nine and cannot contain `J_i`. As an
irreducible linear space, `J_i` lies in one of the other two closed sets.
If it lies in gain zero, `W_i` contains the gain form. Otherwise its linear
projection to the target `D` space lies in the proper cone `Δ`, hence has
dimension at most two. Some nonzero pure target `D` form annihilates that
projection and belongs to `W_i`. Thus each `W_i` contains a nonzero form
`u_i∈Y*`.

Each pair row space contributes at most one further dimension beyond its
chosen `u_i`, so

```
dim P ≤ 4+5=9,      dim Z ≤ 3.
```

First-pullback injectivity gives `dim(P+Z)=12`; all these bounds are
therefore equalities, `P∩Z=0`, and the chosen `u_i` span all of `Y*`.
In particular `Y*⊂P`. This argument uses only properness of `Δ`; it assumes
no independent linear components of that cone.

The source first space has dimension thirteen and pullback image `H`
dimension twelve. The ten source pair coordinate forms have a
one-dimensional restriction kernel because their image is `P` of dimension
nine. Their relation is the hyperplane normal, which consequently has no
source `D` component. Hence the whole three-dimensional source `D` space,
with all pairs zero, lies in `H`. Its target preimage annihilates `P` and
therefore `Y*`: it has only core coordinates. It maps injectively to a
three-dimensional linear subspace of the core pair space. Slice-rank
monotonicity proves

```
rank D(v) ≥ 6          for every nonzero source D input v.
```

The core rank-six locus requires both stacked and horizontal matrices to
have rank one. Equivalently, `A,B` are scalar multiples of the same rank-one
matrix. In coefficient-cube coordinates this is the cone of simple tensors
in `Kbar²⊗Kbar²⊗Kbar²`. Its linear subspaces have dimension at most two.
Indeed two independent simple tensors whose sum stays simple can differ in
only one tensor factor. Every third vector in a simple-only linear space
must vary that same factor: variation in another factor makes its sum with
one of the first two vectors have matrix flattening rank two. The remaining
factor space has dimension two.

The three-dimensional pure-core image cannot lie entirely in that
rank-six cone. Some vector therefore has core rank at least eight, since
`4rank([A;B])+2rank([A B])` takes only the positive values six, eight, ten,
and twelve. The corresponding source `D` slice has rank at least eight.
Consequently

```
h=3  ⇒  ρ ≥ 8,
```

and both other flattening ranks are at least eight. No assumption about
independent rank-drop plane components was needed for this conclusion.
This is an intermediate bound; the rank-drop refinement below raises it
to nine.

Over a finite certificate field, the full three-dimensional nonzero-slice
bound also strengthens the active-term count:

```
r(q³−q²) ≥ 6(q³−1),
r ≥ ceil(6(1+1/q+1/q²)).
```

Together with the refined `ρ≥9` below, this gives tensor rank at least
eleven over `F₂` and nine over every other field, for `h=3`. Again, the finite
count concerns the actual field of the certificate maps.

### Three independent rank-drop planes and the stronger bound nine

The equality `P⊕Z=Kbar¹²*` and inclusion `Y*⊂P` also impose linear components
on the target `D` generic-rank-drop cone. The quotient `P/Y*` has dimension
five. Each of the five spaces `W_i` contributes at most one dimension to it;
all five contributions must be independent. Consequently each `W_i` has
dimension two and `W_i∩Y*` dimension exactly one.

For a gain-zero pair, this intersection is precisely the gain line.
Otherwise the target `D` projection of `J_i` lies in `Δ` and has dimension
at most two. It cannot have dimension below two: its annihilator would
provide at least two independent pure target `D` forms in `W_i∩Y*`. Thus
this projection is a plane in `Δ`, whose normal spans that intersection.
Because `Δ` is proper in three-dimensional space, the plane is an
irreducible component of `Δ`. The gain type and plane type cannot coincide,
as that would again provide two independent intersection forms.

The five intersection lines span all four dimensions of `Y*`. There must
therefore be at least one gain type and at least three plane types with
normal forms spanning all three target `D` dual coordinates. In particular,
`Δ` contains three linear plane components with independent normals. These
planes follow from the restriction; they were not assumptions in the
preceding pure-core argument.

Now suppose `ρ=8`. The three-dimensional pure-core image has every slice
of rank at most eight, by source-to-target slice monotonicity. The formula
`4r_stack+2r_row` implies that every stacked `4×2` matrix in this image has
rank at most one: stack rank two would give core rank at least ten.

A linear space of rank-at-most-one `p×q` matrices has either a common image
line or a common row factor. To check this elementary classification, take
two independent matrices `u vᵀ,p qᵀ`. Their sum has rank at most one only if
`u,p` are proportional or `v,q` are proportional. If two independent image
vectors occur, their row factors must agree, and comparison with both
forces every other matrix's row factor to agree. Otherwise all images lie
in one fixed line. A common image line in `4×2` permits dimension at most
two; our dimension-three image must therefore have a common row factor.
Write

```
A = u vᵀ,      B = w vᵀ,
```

with fixed nonzero `v∈Kbar²` and linearly varying `u,w∈Kbar²`. The core's
rank-eight drop is now `det[u,w]=0`, a homogeneous quadratic in the three
source `D` variables. This quadratic is not identically zero: otherwise
the whole dimension-three image would have core rank six, contradicting
the simple-tensor linear-space bound proved above.

Since source `D` slice rank dominates core slice rank and its generic rank
is eight, its rank-drop cone `Δ` is contained in this quadratic zero set.
The three independent linear plane components just proved would make
three distinct independent linear forms divide a nonzero quadratic.
That is impossible over any field. Thus

```
h=3  ⇒  ρ≥9.
```

Both other flattening ranks are therefore at least nine. The additional first-dimension-four/five filters below further restrict
the componentwise shape boundaries. No surviving boundary is a witness.

### The first-dimension-three diagonal class is excluded

A further refinement excludes every first-concise `h=3` catalyst with a
common diagonal slice form

```
D(δ) = diag(ℓ₁(δ),…,ℓ_N(δ)).
```

Remove identically zero entries first, so the nonzero entry forms span the
three-dimensional first dual space and the generic slice rank is `N`.
The target generic-rank-drop cone is the union of the entry hyperplanes.

Choose three independent source diagonal-entry forms. Their pullbacks span
`Z`, so their common kernel `J=ker Z` in the twelve-dimensional target
first space has dimension nine. On `J`, at least three distinct source
entries vanish. The whole source slice rank is at most
`N−3+5·3=N+12`. Its target slice therefore lies in the finite union of

```
{gain=0},
{ℓ_j(target δ)=0} for j=1,…,N,
Y × B_bad.
```

Although the last component has variety dimension nine, it cannot contain
a nine-dimensional **linear** space. If it contained `J`, the linear core
projection of `J` would have dimension at least `9−4=5`. Irreducibility
places that projection inside one of the two rank-at-most-one rectangular
matrix components. By the elementary matrix-space classification above,
each component has linear subspaces of dimension at most four. This is a
contradiction. Equivalently, the largest possible linear dimension in
`Y×B_bad` is eight.

Thus the irreducible linear space `J` lies in the gain hyperplane or a
target diagonal-entry hyperplane. Some nonzero form in `Y*` annihilates
`J`, hence belongs to `Z`. But `Y*⊂P` and `P∩Z=0`, another contradiction.
The exclusion works for every diagonal length `N` and every characteristic.
It relies on the common diagonal form and the removal of three independent
entries; commuting or upper-triangular slice families are not automatically
covered.

### A rank-drop cone inside one hyperplane excludes any first dimension

The following geometric filter holds for arbitrary `h`. Suppose the
proper generic-rank-drop cone `Δ` of `D` is contained in one proper linear
hyperplane `ker ℓ` over the algebraic closure. Each pair-zero target first
space has dimension at least `h+7`. The core-bad component times target
`D`/gain coordinates has dimension only `h+6`, so irreducibility puts each
pair-zero space in the gain-zero exception or the `D` rank-drop exception.
Its pair row space therefore contains the gain form or the pure `D` form
`ℓ`.

Choose four pair indices whose supported coordinate forms are independent
modulo the source hyperplane normal, as in the first-rank lower-bound proof.
Their selected nonzero target forms must be independent, yet all lie in
`span{gain,ℓ}`, of dimension two. This contradiction excludes such `D` for
any `h`. It applies, for example, to a complete-convolution catalyst
`D=C(a,b)`: every nonzero first slice is injective of rank `b`, so its
generic rank-drop cone is just the origin. This statement concerns `D`
with fixed auxiliary `R(2,1)`, separately from the exclusion of complete
convolutions as `S` above.

## First dimensions four and five need linear rank-drop components

For fixed auxiliary `R(2,1)`, the target first space has dimension `h+9`.
Let `P` again span the five source pair row spaces and `Z` the source `D`
rows, with `dim Z≤h`. First conciseness implies `dim(P+Z)=h+9`, so
`dim P≥9`.

Suppose the generic-rank-drop cone `Δ` has no linear hyperplane component.
Every linear target `D` projection contained in `Δ` then has dimension at
most `h−2`. As before, each pair-zero space is in gain zero or the `D`
rank-drop exception: the core-bad exception has dimension `h+6`, below its
required dimension `h+7`. For a pair of the latter type, at least two pure
`D` forms annihilate its projection. The pair row space has dimension at
most two, so it consists entirely of two pure `D` forms.

If there are `t` such pure-`D` pairs, their combined span has dimension at
most `min(2t,h)`. Each of the other `5−t` gain-zero pairs contains the common
gain form and contributes at most one additional row. For `t<5`,

```
dim P ≤ min(2t,h) + 1 + (5−t).
```

For `t=5`, `dim P≤h`. At `h=4` and `h=5`, the maximum of these bounds is
eight. This contradicts `dim P≥9`. Thus both first dimensions require a
linear hyperplane component of `Δ` over the algebraic closure.

At `h=4`, one hyperplane component is also insufficient. If its normal is
`ℓ`, non-gain pairs either contain `ℓ` or consist of two pure `D` forms
coming from a smaller linear projection. With `t` of the second kind,

```
dim P ≤ min(2t+1,4) + 1 + (5−t) ≤ 8
```

when `t<5`; the all-pure-`D` case is smaller. The extra one inside the minimum
accounts for `ℓ`; the other one accounts for the gain form. Consequently
`h=4` requires at least two distinct linear hyperplane components, whose
normals are independent. A generic irreducible determinant is not sufficient.

In particular, `h=4` or `h=5` cannot have generic slice rank six. The
necessary rank-at-least-six two-plane would avoid the generic-rank-drop
cone except at the origin. Every linear hyperplane component intersects
that two-plane in a nonzero vector, contradicting the required slice rank.
Thus

```
h∈{4,5}  ⇒  ρ≥7.
```

Both other flattening ranks must be at least seven. These first-dimension-specific rank-seven bounds are weaker than the
subsequently proved global rank-eight bound. With that global bound and the
first-dimension-three rank-nine refinement, the componentwise boundaries
of these shape inequalities are `(3,9,9)` and `(4,8,8)`. They are
necessary-screen boundaries, not witnesses.

This excludes the whole proposed binary band family

```
D(a,b,c,d) = L₆(a,b) + cC + dE,      shape (4,6,7),
```

regardless of the two additional `6×7` matrices. Its `c=d=0` plane has rank
six at every nonzero closure point: the first and last maximal minors of
the complete band pencil are `a⁶` and `b⁶`. Its generic rank is six because
there are only six rows. Random additional directions cannot repair this
obstruction. If its first flattening rank drops below four, the earlier
first-rank-three or first-rank-at-least-three conditions exclude it instead.
The transposed `(4,7,6)` boundary is excluded as well. Testing whether binary
rank-drop inputs span all first coordinates would miss this stronger
necessary condition.

## A convolution plus two units is also excluded as catalyst

Another proposed catalyst was

```
D = C(2,b) ⊕ Unit ⊕ Unit,      b≥1,
```

with first inputs `(a,beta,c,d)` for the convolution pair and the two units.
For `b=5`, it has concise shape `(4,7,8)`, generic slice rank seven, and a
closure-defined two-plane `(a,beta,c,d)=(x,y,x,y)` whose nonzero slices have
rank at least six. It passes those basic checks, but the following source
row geometry excludes it. Use `beta` here to distinguish a first input
coordinate from the integer convolution length.

Its generic-rank-drop cone is exactly

```
{c=0} ∪ {d=0} ∪ {a=beta=0}.
```

Each pair-zero target subspace lies in the gain-zero exception, one of the
first two hyperplanes, or the last codimension-two plane exception. In the
last case its row space must be exactly `span{a,beta}`. All such row spaces
are therefore the same two-dimensional space. If none have this last type,
`P` is generated by at most three common forms (gain, `c`, `d`) and five
additional rows, hence has dimension at most eight. If there are `t≥1`
last-type pairs, then

```
dim P ≤ 2+3+(5−t)=10−t.
```

The required `dim P≥9` forces exactly one such pair and `dim P=9`.
Equality makes the selected gain/`c`/`d` forms span all three of their
coordinates. Together with `a,beta`, the entire five-dimensional target
`D`/gain form space `Y*` lies in `P`. The four source `D` rows span a
complementary `Z`, since injection requires rank thirteen.

The ten source pair forms have image of dimension nine, so their one
relation supplies the source hyperplane normal with no source `D`
component. The whole source `D` first space consequently lies in the
pullback image and maps to pure target core, because `Y*⊂P`. Every nonzero
source `D` slice would then have rank at least six. But a unit-only input
has slice rank one. This contradiction excludes every length `b`, even
though the `b=5` candidate passes the basic shape and generic-rank screens.

## First-rank-four diagonal catalysts are excluded

Here diagonal means a common diagonal slice form
`diag(ℓ₁(δ),…,ℓ_N(δ))`, with nonzero entry forms spanning a four-dimensional
first space. This exclusion permits repeated proportional forms.

Let `Y*` be the five-dimensional target `D`/gain form space and retain
`P,Z` as above. In the source first space, write the pullback hyperplane
normal as `η=(η_D,η_pairs)`. If `η_D≠0`, all five nonzero pair-supported
annihilators are independent modulo this normal: a pair-supported sum
cannot equal a nonzero multiple of a normal with a source `D` component.
They span `Y*`, giving `Y*⊂P`. If `η_D=0`, at least four remain independent,
so `dim(P∩Y*)≥4`.

When `η_D=0`, `dim P=9`, `dim Z=4`, and `P∩Z=0`; the image of `Y*` modulo
`P` has dimension at most one, hence `dim(Z∩Y*)≤1`. When only `η_D` is
nonzero, `dim P=10`, `dim Z=3`, and `P∩Z=0`, so the intersection is zero.
When both normal components are nonzero, `dim P=10`, `dim Z=4`, and
`dim(P∩Z)=1`; since `Y*⊂P`, the intersection again has dimension at most one.
In all cases

```
dim(Z∩Y*) ≤ 1.
```

The necessary six-rank two-plane and exclusion of two other dimensions
both equal to six imply `N≥7`. Moreover any proportional class of diagonal
entry forms has multiplicity at most `N−6`: its hyperplane intersects the
necessary two-plane in a nonzero vector, where all those entries vanish
but the remaining slice rank must be at least six.

Kill four independent source diagonal-entry forms, whose pullbacks span
`Z`. Their common kernel `J=ker Z` has dimension at least nine. The source
`D` slice is identically zero on `J`, so the complete source slice has rank
at most fifteen, from the five retained-sector blocks. Thus target slices
on `J` also have rank at most fifteen.

Outside the core-bad exception, core rank is twelve, so target `D` rank
is at most three, regardless of gain. At least `N−3` diagonal entries
vanish. They cannot all belong to one proportional entry-form class,
because its multiplicity is at most `N−6`. Therefore two independent
entry forms vanish. This is the needed correction for repeated forms:
merely counting two vanished entries would not suffice.

It follows that `J` is contained in the finite union of pairwise kernels
of two independent target entry forms and the component `Y×B_bad`.
Irreducibility puts it in one member. The pairwise-kernel case puts two
independent forms of `Y*` in `ann J=Z`, contradicting the intersection bound.
In the core-bad case its core projection has dimension at least
`dim J−5≥4`. The maximum linear dimension in `B_bad` is four, so necessarily
`dim J=9`, the projection has dimension four, and its kernel is all of
`Y`. Consequently `J=Y⊕L` for a four-dimensional linear `L⊂B_bad`.

A maximal four-dimensional linear space in the stacked rank-one component
has a common row factor, with freely varying `u,w`, and generic core rank
eight. In the horizontal rank-one component it has a common column factor
and generic core rank ten. Thus `L` contains a core input of rank at least
eight. Independently choose target `D` coordinates of generic rank `N`
and gain one in the whole `Y` factor. Its slice has rank at least `N+9>15`,
contradicting the bound on `J`.

This excludes every first-rank-four common diagonal catalyst, for all
lengths and all characteristics. It uses the actual source bound fifteen
and the six-plane multiplicity bound; a weaker generic-loss count would
leave an unjustified gap for repeated diagonal forms. No analogous
first-rank-five diagonal exclusion is claimed.

## The three-plane argument also strengthens the fixed-star catalyst bounds

This section explicitly supersedes the earlier weaker catalyst bounds in
[stage4.md](stage4.md), [star-nonzero-catalyst.md](star-nonzero-catalyst.md),
and [star-convolution-catalysts.md](star-convolution-catalysts.md). Their
older necessary inequalities were valid; the new argument makes more
catalysts impossible. It does not produce or exclude all fixed-star
certificates, and does not change the fixed-star target's stronger-exponent
calibration.

For the singular star auxiliary, the same source-normal cases give a
coefficient-field-defined three-plane in `D` mapping injectively to pure
product core. Its slice-rank formula is

```
2 rank([A;B]) + 2 rank([A B]),
```

with positive values four, six, and eight. Every nonzero slice in the source
three-plane has rank at least four geometrically. Generic `D` rank at most
five would force the whole core three-space into the rank-four cone, again
the simple `2×2×2` tensor cone of maximum linear dimension two. Therefore

```
ρ_D≥6         for every fixed-star catalyst first dimension.
```

Both other flattening ranks are at least six. Active-term counting on the
three-plane gives the global same-field tensor-rank bound

```
rank_K(D) ≥ max(6, ceil(4(1+1/q+1/q²)))    over F_q.
```

It is seven over `F₂` and six over every other field. The general-field
rank lower bound follows directly from `ρ_D≥6`.

At first flattening rank three, the previously audited equality of pair
row-space dimensions gives the whole source `D` three-space a pure-core
image and three independent linear rank-drop plane components. Suppose
its generic rank were six. Its core image has rank at most six, hence
lies in the union of the stacked and horizontal rank-at-most-one
rectangular matrix cones. Irreducibility places the linear image in one
component. A dimension-three stacked rank-one space must have a common
row factor; the common-image alternative has maximum dimension two. The
horizontal case similarly has a common column factor.

The generic core rank is six, since a rank-four-only linear space has
dimension at most two. Its rank-four drop is a nonzero quadratic determinant
of the two varying vector factors. Since `D` rank dominates core rank and
has generic rank six, its generic-rank-drop cone is contained in that
quadratic zero set. Three independent linear plane components cannot divide
a nonzero quadratic. Thus

```
h=3  ⇒  ρ_D≥7.
```

Both other flattening ranks are at least seven in that first dimension.
The same-field tensor-rank lower bound there is seven over every field,
combining generic rank and the finite count. In particular the historical
binary `(4,4,5)` random catalyst is now mathematically excluded by the
other-axis bounds. Its older SAT timeout remains a historical solver result,
superseded for viability by this proof.

The generic-rank bounds do not assert that a finite base field contains a
point attaining the generic rank. A base-field slice enumeration with
maximum below six is not automatically a generic-rank computation. The
coefficient-field-defined three-plane guarantees pointwise rank at least
four on its nonzero vectors; exact finite coefficient counts and geometric
rank claims must keep these different assertions separate.

## General first-dimension exclusions supersede the small-catalyst analyses

This argument treats arbitrary `D`, including non-diagonal tensors. It
uses the source-`D`-zero subspace rather than killing a few source entries.
It applies to both fixed auxiliaries considered here. Put `r=2` for the
star and `r=3` for `R(2,1)`: a nonzero auxiliary slice has rank `r`, the
product core has generic rank `4r`, and the necessary source `D` three-plane
has every nonzero slice of rank at least `2r`.

### The intersection of source and target `D` coordinate forms

Retain the target `D`/gain form space `Y*` of dimension `h+1`, pair row span
`P`, and source `D` row span `Z`. The source hyperplane normal has the
following cases.

- If its `D` component is zero, `dim P=9`, `dim Z=h`, and `P∩Z=0`.
  Four independent pair-supported annihilators give `dim(Y*∩P)≥4`.
- If its `D` component is nonzero, all five pair-supported annihilators
  are independent, giving `dim(Y*∩P)≥5`. Here `dim P=10` and
  `dim(P∩Z)≤1`.

In either case, the projection of `Y*∩Z` into the quotient by `P`, allowing
for `P∩Z`, proves

```
dim(Y*∩Z) ≤ h−3.
```

Let `J=ker Z` in the target first space. Since `dim Z≤h`, its dimension is
at least nine. The source `D` slice vanishes on `J`, so every target slice
there has rank at most `5r`.

### The low-target-`D` exception is impossible at h≤5

Off the core-bad locus `B_bad`, the core rank is `4r`. Thus the target `D`
slice rank is at most `r`, ignoring the possibly nonzero gain. Hence `J`
is contained in the union of the two closed sets

```
{target D rank≤r} × gain × core,
Y × B_bad.
```

By irreducibility of the linear space `J`, it lies in one set. If it lies
in the first, its linear projection `V` onto the target `D` first coordinates
consists of slices of rank at most `r`. It meets the necessary rank-at-least-
`2r` three-plane only at zero. Linear dimension counting gives
`dim V≤h−3`. At least three independent pure target `D` forms therefore
annihilate `J`, putting three independent elements of `Y*` in
`ann J=Z`. This contradicts the preceding intersection bound when `h≤5`.

No rank-additive state inequality is used here. The closed rank condition
is the vanishing of ordinary matrix minors; its linear projection dimension
is bounded by the explicitly proved source `D` slice plane.

### The core-bad exception excludes h≤4

If `J⊂Y×B_bad`, its linear core projection has dimension at least
`dim J−(h+1)≥8−h`. The maximum linear dimension in `B_bad` is four.
For `h=3` this is already impossible. Smaller first dimensions were
excluded earlier.

For `h=4`, equality is forced: `dim J=9`, the core projection has dimension
four, and its kernel is all five-dimensional `Y`. Thus `J=Y⊕L`, for a
four-dimensional linear `L⊂B_bad`. Such a maximal rectangular rank-one
space has generic core rank at least six for the star, or at least eight
for the retained tensor. Independently choose target `D` of its generic
rank (at least six or eight) and gain one. The resulting target slice rank
exceeds the source bound ten or fifteen. This contradiction excludes
**every `h=4` catalyst**, with no diagonal hypothesis.

Combining the two exceptions proves that both fixed auxiliaries require
`h≥5`.

### The star also excludes h=5

For the star with `h=5`, take `K=J∩Y`. The core projection `L` has dimension
at least three and at most four, and `dim K≥5`. Its generic core rank is
at least six: a rank-four-only core space has maximum linear dimension two.

If the gain form vanishes on `K`, this five-dimensional space is the whole
target `D` coordinate space. Fix a core input of rank at least six and
vary its fiber by `K` to achieve target `D` generic rank at least six. Its
target slice rank exceeds ten regardless of its gain value.

Otherwise the gain is nonzero on `K`. Fix a core input of rank at least
six and gain value one; this is possible by adding an element of `K`.
The remaining target `D` coordinates form an affine fiber `δ₀+V`, where
`V` is the `D` projection of `K∩{gain=0}` and has dimension at least four.
Its intersection with the necessary `D` three-plane has dimension at least
two. Choose a nonzero direction `v` there with `rank D(v)≥4`.

Some four-by-four minor of `D(v)` is nonzero. The same minor of
`D(δ₀+t v)` is a polynomial with that nonzero leading coefficient at
degree four. Over the algebraic closure it is nonzero for some `t`.
Thus the fixed-gain affine fiber contains a target `D` slice of rank at
least four. Adding gain one and core rank at least six gives rank at least
eleven, again exceeding the source bound ten.

This excludes star `h=5`. For the retained tensor the corresponding estimate
is only `6+1+8=15`, equal to its source bound; this argument leaves `h=5`
unresolved.

The final first-axis and generic-slice-rank conclusions are therefore:

| Fixed auxiliary | Catalyst first flattening | Catalyst generic slice rank | Both other flattenings |
| --- | ---: | ---: | ---: |
| `R(2,1)` | at least 5 | at least 8 | at least 8 |
| `S_star` | at least 6 | at least 6 | at least 6 |

These are necessary conditions, not existence results. In particular all
first-dimension-three/four shape boundaries discussed earlier in this note
are excluded by the more general argument. A surviving dimension triple
must still satisfy every applicable geometric and exact-rank screen.

## Low-core normal forms impose stronger other-axis bounds

The global three-plane in source `D` maps injectively to a three-dimensional
linear core pair space `L`. Restrict the first input to this plane and
project target other coordinates onto the pure core blocks. This gives an
actual restriction from that source `D` subtensor to the corresponding core
subtensor. Its other flattening ranks therefore cannot exceed those of
`D`. This uses ordinary flattening monotonicity, not tensor-rank additivity.

For the retained auxiliary, generic `D` slice rank eight or nine forces
stacked rank at most one throughout `L`; stacked rank two would already
give core rank at least ten. Its dimension-three matrix space has a common
row factor, so `A=u vᵀ,B=w vᵀ` with fixed `v` and a dimension-three linear
space of pairs `(u,w)∈K²⊕K²`. For generic rank ten or eleven, the image is
instead contained in the union of the stacked and horizontal rank-one
components. Irreducibility selects one, giving this common-row form or the
common-column form `A=u vᵀ,B=u wᵀ` with fixed `u`.

These normal forms have the following exact other-flattening ranks.
Each operation includes the two independent matrix columns as spectators.

| Core operation | Common row factor | Common column factor |
| --- | --- | --- |
| Column block `(AX,BX)` | input 2, output 6 | input 4, output 4 |
| Row block `AX+BY` | input 4, output 4 | input 6, output 2 |

For example, with common row factor the column input only sees the two
entries of `vᵀX`; its output coefficient span is the three-dimensional
space of `(u,w)` times two columns, giving six. In the row block the two
scalar inputs per column are independent and the union of all `u,w` spans
`K²`: otherwise their dimension-three pair space would lie in a
fixed-line pair space of dimension only two. Thus its input and output
flattenings are four. The common-column calculation is the transpose
version: the column operation sees the full two-dimensional row-factor
span, while the row operation has coefficient-pair span dimension three
and fixed one-dimensional output image.

The blocks occupy independent other-coordinate subspaces, so their
ordinary matrix flattening ranks add even though the first axis is shared.
The retained core has two column blocks and one row block; the star core
has one of each. Consequently every dimension-three image in these normal
forms has other flattenings

| Fixed auxiliary core | Common row factor | Common column factor |
| --- | --- | --- |
| `R(2,1)` product | `(8,16)` | `(14,10)` |
| Star product | `(6,10)` | `(10,6)` |

All dimensions concern the original order of the second and third legs.
The retained auxiliary is asymmetric in those legs.

It follows that retained catalysts satisfy the additional conditions

```
ρ∈{8,9}   ⇒ third flattening ≥16,
ρ∈{10,11} ⇒ third flattening ≥16 or second flattening ≥14.
```

Star catalysts with `ρ∈{6,7}` require second flattening at least ten or
third flattening at least ten. These are generic-rank implications: no
finite-point maximum is substituted for `ρ`. Since generic slice rank is
at most each other flattening rank, they also provide immediate necessary
shape filters without computing it.

Combining the first-axis exclusions, generic-rank bounds, and these
normal-form bounds leaves the following componentwise dimension boundaries
for these particular necessary inequalities:

- Retained `R(2,1)`: `(5,8,16)`, `(5,14,10)`, or `(5,12,12)`.
- Fixed star: `(6,6,10)`, `(6,10,6)`, or `(6,8,8)`.

A successful catalyst must dominate at least one listed triple in the
stated leg order. These boundaries do not assert existence, and other
component or rank-drop obstructions can exclude tensors having these
shapes.

An independent temporary coefficient experiment checked 32 instances of
these flattening identities over `F₂,F₃,F₄,F₅`, for both auxiliaries and both
normal forms, including the rank-one and rank-two hyperplane normal types
for the dimension-three pair space. Every expected flattening rank matched.
These checks are performed evidence rather than a tracked control artifact;
the span calculations above prove the general normal-form dimensions.

## The outer-product plus dot-tail candidate fails a component screen

A proposed retained-sector catalyst was the direct sum of a full outer
product of shape `(4,6,24)` and a singleton-first-leg length-two dot tensor.
Its concise shape is `(5,8,26)`. Writing its first inputs as `(v,t)` with
`v∈K⁴`, its slice rank is six whenever `v≠0`, plus two whenever `t≠0`.
Its generic rank-drop cone is

```
{t=0} ∪ {v=0}.
```

It has a constant-rank-six source three-plane and passes those elementary
checks, but cannot be a catalyst. A pair-zero target first space lies in
gain zero or the target `D` rank-drop exception; the smaller core-bad
component cannot contain it. Its linear target `D` projection into the
rank-drop cone lies in one of the displayed components. If it lies in
`{v=0}`, four independent pure target `D` forms annihilate it and must
belong to its pair row space, whose dimension is at most two. This is
impossible. Thus every pair row space contains the gain form or `t`.
Their whole span has dimension at most `2+5=7`, contradicting `dim P≥9`.

The same argument excludes any constant-nonzero-slice-rank first space of
dimension at least three, direct-summed with one singleton-first dot or
scalar tail: the zero-first-space component has codimension at least three
and cannot be the projection exception for a two-row pair map. This
applies to either fixed auxiliary. No SAT run is needed to establish these
exclusions; passing dimension and finite-point rank-drop-span checks alone
would miss them.

## A concrete constructive relation: `M₂≤R(2,1)²`

The root coordinator found actual selector maps from the square of the
retained auxiliary to `M₂`, independently checked here at every coefficient.
The square has shape `(4,16,25)`. In zero-based lexicographic tensor-product
coordinates, the selected source indices for the three target legs are

```
first:  (1,0,3,2),
second: (1,13,2,14),
third:  (2,17,7,22).
```

Equivalently, for standard matrix coordinates, the pullbacks are

```
first  (i,j) ↦ (i,1−j),
second (j,k) ↦ (3k,1+j),
third  (i,k) ↦ (i+3k,2).
```

The retained tensor's left and right branches obey
`R(i,3k,i'+3k')=1` precisely when `i=i'` and `k=k'`. Its middle branch obeys
`R(1−j,1+j',2)=1` precisely when `j=j'`. Thus the pulled-back square
coefficient is exactly

```
R(i,3k,i'+3k') R(1−j,1+j',2)
  = 1_(i=i') 1_(j=j') 1_(k=k'),
```

the matrix multiplication tensor coefficient. The selected tensor words
are left-middle and right-middle. Only coordinate selections occur, with
integer coefficients zero and one; no cancellation, interpolation, or
characteristic assumption is involved.

Independent exact checks succeeded over `F₂,F₃,F₄,F₁₁`. The displayed
coefficient identity proves the restriction over every field. This is an
actual legal ordinary state row, `f(M₂)≤f(R(2,1)²)`, and can be tensored with
any supplied context maps. It connects previously separate power levels
and can invalidate a proposed finite state that assigned inconsistent
values to them.

This construction is not a positive-gain catalyst or a useful near-`9/4`
rank scheme. Selecting from the coordinate six-term auxiliary decomposition
squared leaves the naive eight matrix products. The numerical upper bound
obtained merely from the auxiliary's asymptotic rank five would charge two
auxiliary copies per matrix factor; it supplies no new competitive matrix
exponent. Its value for discovery is the new exact restriction relation
and the resulting finite-state constraints.

## Two retained copies supply matrix multiplication: `M₂≤2R(a,h)`

For every `a≥2,h≥1`, there is also an actual restriction from two
**direct-sum** retained copies to `M₂`. The portable exporter is
[`export_retained_two_copy_matrix_order`](../sector_matrix.py), with
coefficient, copy-placement, and resource-cap tests in
[`test_sector_matrix.py`](../test_sector_matrix.py). This construction
requires no tensor square or multiplicative state identity.

Set `μ=h+a−1` and `r=2h+a−1`. In each source copy select first coordinates
`0,a−1`, and send the target second coordinates `(B₀₀,B₀₁,B₁₀,B₁₁)` to
source coordinates `(μ,0,h,r)`. For selected first entries `(x₀,x₁)`, the
middle output at `μ` is `x₀B₀₀+x₁B₁₀`. The sum of the left output at zero
and right output at `r+a−1` is `x₀B₀₁+x₁B₁₁`. The branch embeddings give
these four coefficients as one; every cross-selected coefficient is zero
by its branch range. This is an integer coefficient identity, valid over
every field. The two copies use independent first entries for the two
rows of `A`, while both receive the same `B`. Their two output pairs are
the four entries of `AB`.

In the concatenated first input, the target selectors are
`(0,a−1,a,2a−1)`. The second map duplicates the displayed four selections
in both source copies. The third map retains each copy's middle output
and its sum of the indicated left and right outputs. Thus there are two
positive source block occurrences, even in characteristic two; no field
scalar representing their count occurs in the restriction.

Independent inspection of the support formula proves this relation for
all the stated parameters. The exporter was checked at every coefficient
in 16 agent cases over `F₄,F₁₁,F₁₆,F₃₁`, for `(a,h)=(2,1),(3,1),(2,2),(3,2)`,
and in 20 additional root controls. Its ordinary additive state row is

```
f(M₂) ≤ 2 f(R(a,h)).
```

Consequently a detector state with `f(M₂)=5` must have `f(R(a,h))≥5/2`.
This excludes earlier assignments with retained value two even without
cross-level multiplicative identities. It is neither a positive-gain
catalyst nor a new `9/4` rank upper bound: the small coordinate scheme
compiled through the maps gives the naive eight matrix products.

## A signed all-field restriction: `M₂≤S_star²`

The fixed star also has an actual square-to-matrix restriction. Unlike the
retained-sector selector above, this construction uses signed linear maps.
It holds over the integers and therefore over every field, including
characteristic two. The root coordinator found the maps; the coefficient
identity and interpretation below were independently checked here.

Take the star support to be

```
(0,0,1), (0,1,0), (1,0,2), (1,2,0).
```

Its square has shape `(4,9,9)`, with zero-based lexicographic coordinates.
The first map is the identity on the four entries of `A`. The second map's
rows, in target order `(B₀₀,B₀₁,B₁₀,B₁₁)`, have the following nonzero source
coefficients:

```
{1: 1}, {0: 1, 6: 1}, {2: 1}, {3: −1, 4: 1, 8: 1}.
```

The third map's rows, in target order `(C₀₀,C₀₁,C₁₀,C₁₁)`, are

```
{3: 1}, {2: −1, 4: 1, 8: 1}, {6: 1}, {0: 1, 1: 1}.
```

For a direct bilinear explanation, write the second source input as a
`3×3` array with scalar `β` in position `(0,0)`, vector `u` in positions
`(0,1),(0,2)`, vector `v` in positions `(1,0),(2,0)`, and matrix `W` in the
remaining `2×2` block. The star square outputs `Au`, `Aᵀv`, `βA`, and
`tr(AᵀW)` in its corresponding positions. For the target second input
`B=[u w]`, with `w=(w₀,w₁)ᵀ`, set

```
β=w₀,    v=(−w₁,w₀)ᵀ,    W=w₁ I₂.
```

If `A=[[a,b],[c,d]]`, its output `y=Aᵀv` has coordinates
`y₀=−a w₁+c w₀`, `y₁=−b w₁+d w₀`. The third map recovers the second
column of `AB` through

```
−y₁+tr(βA)       = a w₀+b w₁,
 y₀+tr(AᵀW)     = c w₀+d w₁.
```

The first column is the supplied `Au`. These identities prove every target
coefficient without dividing by a scalar. In particular, reducing the
signed maps modulo two preserves the construction.

Independent controls checked all 64 coefficients over the integers and
over `F₂,F₃,F₄,F₅,F₁₁`. This supplies the ordinary monotonicity row
`f(M₂)≤f(S_star²)`, and its tensor-context versions when the context maps
are supplied. It gives a concrete relation between power levels; an
assignment treating those levels independently must respect this row.
It supplies no positive-gain catalyst or new near-`9/4` rank upper bound.

## Verification and use in finite state experiments

The durable coefficient controls are
[`convolution_commutator_experiment.py`](convolution_commutator_experiment.py)
and [`test_convolution_commutator.py`](../test_convolution_commutator.py).
They cover 144 complete-convolution cases over `F₂,F₃,F₄,F₅`, lengths
`a,b=2,…,4`, repetition counts one and two, and both input orientations,
plus four explicit mixed-complement controls. The nine focused tests also
include Strassen/projective-interpolation rank upper controls, characteristic
two, full tensor restrictions, invalid inputs, and resource caps.

In addition, an independent temporary experiment completed 16 retained-sector
cases over those four fields, for `(a,h)=(2,1),(2,2),(2,3),(3,1)`. It checked
concise flattenings, the actual middle-branch coordinate compression,
preservation of full output concision, and the `3h`-dimensional identity
image. Those retained-sector controls were performed but are not themselves
a tracked reproducible artifact.

The controls test the coefficient constructions underlying the proofs;
they do not enumerate arbitrary decompositions or prove a rank theorem by
finite sampling. The proofs above supply the arbitrary-decomposition output
elimination, arbitrary-`D` repetition, and first-space geometry arguments.
No Lean edits, Lake build, or new kernel theorem are claimed.

A tensor-rank lower bound is **not** a normalized-state lower constraint:
`rank(T)≥r` does not imply `f(T)≥r` for additive monotone states. These rank
arguments screen candidate exact catalysts. Finite LP lower rows must instead
carry actual restriction maps, such as `r Unit≤T` proving subrank at least
`r`. Rank lower bounds and subrank witnesses must remain distinct.
