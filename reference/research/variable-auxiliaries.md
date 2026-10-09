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
initially survived these necessary bounds. The later
[replicated rank audit](mixed-dot-rank-catalysts.md) excludes it for every
finite catalyst, including unequal positive dot lengths. The intermediate
slice-space results below retain their narrower scope. The retained and
star families are not catalyst witnesses. The subsequent
[integer Koszul proof](star-koszul-catalysts.md) also excludes the fixed star
for every finite catalyst.

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

## Three oriented dot tensors: exclusions for three small catalysts

Let `S` be the direct sum of one length-two dot tensor in each singleton-leg
orientation. Its concise shape is `(5,5,5)` and its supplied exact
decomposition has six terms. Although the previous common-output
convolution argument does not
apply, the following new slice-space argument excludes

```
D ⊕ Unit ⊕ (M₂⊗S) ≤ D ⊕ 5S
```

for **each of `D=0`, `D=Unit`, and `D=M₂`, over every field**. This also
excludes any positive gain for these same catalysts, since projecting a
larger gain to one unit preserves the other target blocks. It does not
exclude an arbitrary larger catalyst.

Work over the algebraic closure of the coefficient field. An exact tensor
restriction remains exact there. Label the singleton-first block of `S`
by its scalar first input `s`, and label its other first inputs by two
vectors `u,v∈K²`. Its slice rank is

```
2·1_(s≠0) + 1_(u≠0) + 1_(v≠0).
```

Consequently each of the five source copies has slice rank at most four.
Let `s₁,…,s₅` denote its five scalar coordinates, pulled back as linear
forms on the target first space through the putative first map. Setting
two of these forms to zero reduces the source rank upper bound by four;
setting all five to zero bounds the source contribution of `5S` by ten.
These statements use distinct source block occurrences even if some of
the pulled-back forms happen to coincide.

The core first space of `M₂⊗S` has three disjoint parts: a `2×2` matrix
`A`, a `4×2` matrix `B`, and a `2×4` matrix `C`. Its slice rank is

```
4 rank(A) + 2 rank(B) + 2 rank(C).
```

For example, the first core block is `M(2,2,4)`; the other two are
`M(4,2,2)` and `M(2,4,2)`. Their generic slice ranks are eight, four, and
four. There is also the target gain form `g`. In the three cases the
target first dimension, generic target slice rank, and source rank after
killing two scalar forms are

| `D` | Target first dimension | Generic target slice rank | Source upper bound on a pair kernel |
| --- | --- | --- | --- |
| `0` | 21 | 17 | 16 |
| `Unit` | 22 | 18 | 17 |
| `M₂` | 25 | 21 | 20 |

For `D=M₂`, write its independent first matrix as `E`; its slice rank is
`2 rank(E)`. For `D=Unit`, write its scalar form as `d`.

Each pair kernel `Jᵢⱼ=ker(sᵢ,sⱼ)` therefore lies in the target's generic
rank-drop cone. A linear space over an infinite field is irreducible, so
it lies in one of that cone's finitely many components. Those components
are `g=0`, also `d=0` for a unit catalyst, also `rank(E)≤1` for a matrix
catalyst, and the three core matrix rank-drop components.

The maximal dimension of a linear subspace in the rank-one cone of a
`2×2` matrix space is two; for a `4×2` or `2×4` matrix space it is four.
This follows from the common-image or common-row classification used
above. A pair kernel has dimension at least 19, 20, or 23 in the three
cases. The `B` and `C` bad components have maximal linear dimensions only
17, 18, or 21, so neither can contain a pair kernel. If `Jᵢⱼ` lies in an
`A` or `E` rank-one component, equality of dimensions forces its two
normal forms to be independent, pure forms in that matrix block. Their
normal two-plane is itself a rank-one matrix space: the annihilator of
a maximal common-row or common-image two-plane in a `2×2` matrix space
has the other common factor. If it lies in a scalar component, the
corresponding scalar form belongs to `span(sᵢ,sⱼ)`.

Here is the resulting bound on the whole span of the five forms. Let
`Y` be the pure scalar row space: `span(g)` for `D=0` or `M₂`, and
`span(g,d)` for `D=Unit`. If the images of the `sᵢ` modulo `Y` span at
most one dimension, their full span has dimension at most `dim(Y)+1`.
Otherwise choose two forms with independent quotient images. Their pair
cannot contain a scalar normal, so both belong to the same pure matrix
block. Any further form with nonzero quotient image is independent of
at least one of this pair, forcing it into that same block. Their matrix
span has dimension at most two: every independent pair spans a rank-one
space, so the determinant quadratic vanishes on the whole span, including
in characteristic two. The remaining forms lie in `Y`. Thus

```
dim span(s₁,…,s₅) ≤ 3   for D=0 or M₂,
dim span(s₁,…,s₅) ≤ 4   for D=Unit.
```

For `D=Unit`, attaining four forces the span to be the direct sum of a
pure `A` rank-one normal two-plane and all of `span(g,d)`.

To finish, set `J=ker(s₁,…,s₅)`. The following elementary dimensions
bound linear spaces of uniformly low target slice rank:

| `D` | Source rank upper bound on `J` | Maximum dimension of a target linear space with this rank bound |
| --- | --- | --- |
| `0` | 10 | 17 |
| `Unit` | 11 | 18 |
| `M₂` | 14 | 21 |

For completeness, let `a,b,c∈{0,1,2}` be the generic ranks of the three
matrix projections of such a space. Their dimension bounds are
`f(a)+t(b)+t(c)`, where `f=(0,2,4)` and `t=(0,4,8)`; its generic core
slice rank is `4a+2b+2c`. Add zero or one for each scalar projection;
for the matrix catalyst add `f(e)` to dimension and `2e` to rank.
The dimension of the space is at most the sum of its projection
dimensions. Enumerating these three-valued ranks proves the displayed
maxima. Each maximum requires `a=0`, `b=c=2`, and all remaining target
blocks unrestricted.

For `D=0`, the scalar-span bound gives `dim J≥18`, contradicting the
maximum 17. For `D=M₂`, it gives `dim J≥22`, contradicting the maximum
21. For `D=Unit`, a scalar span of dimension at most three gives
`dim J≥19`, contradicting the maximum 18. The only remaining possibility
is a four-dimensional scalar span. Then

```
J = (an A rank-one two-plane) ⊕ (the whole B space) ⊕ (the whole C space),
```

with `g=d=0`. This space has generic target slice rank `4+4+4=12`,
contradicting the source upper bound eleven. This completes all three
exclusions without a tensor-rank additivity theorem or a SAT result.

An independent finite enumeration checked the last dimension maximization;
its unique maximizing rank profiles are `(a,b,c)=(0,2,2)` with unrestricted
scalar or catalyst blocks. The proof uses algebraic-closure linear spaces,
not finite-point rank sampling. Larger `D` remains unknown from this
argument; the surviving three-orientation auxiliary is still a search
candidate, with these three catalyst choices removed.

## Every diagonal catalyst is excluded for the mixed-dot auxiliary

The small-catalyst argument strengthens to **`D=n Unit` for every finite
`n≥0`, over every field**. The proof below concerns actual restrictions
and slice ranks; it does not assert an optimal tensor-rank value over a
finite field. It also excludes any positive gain by projection to one
gain unit.

Work over the algebraic closure as before. The target first space has
dimension `n+21`; its independent scalar-coordinate forms are the `n`
target diagonal forms and the gain form. Denote their full row space by
`Y`, of dimension `n+1`. The remaining first coordinates are the core
matrices `A,B,C` from the preceding argument. Let `c₁,…,c₅` be the five
source singleton-first scalar forms, pulled back to this target space.
The generic target slice rank is `n+17`. On each pair kernel the source
slice rank is at most `n+16`. The same dimension argument as above forces
each pair normal span to contain a target scalar-coordinate form or to
be a pure `A` rank-one normal two-plane. In particular, the larger scalar
row space does not permit an arbitrary scalar combination as the first
alternative: it must contain an individual coordinate form.

Let `J=ker(c₁,…,c₅)`. All five source singleton-first inputs vanish on
`J`, so every source slice there has rank at most `n+10`.

First suppose the images of the `cᵢ` modulo `Y` span at most one
dimension. The projection of `J` to each core matrix block then has
codimension at most one: its annihilator is the intersection of the
source normal span with that block's pure row space, which has dimension
at most one. The projected matrix spaces therefore attain ranks two,
two, and two, so the generic core slice rank on `J` is sixteen. At most
five independent target scalar-coordinate forms vanish identically on
`J`, since its normal span has dimension at most five. Hence generic
target slice rank is at least `16+n+1−5=n+12`, a contradiction.

Otherwise choose two forms with independent images modulo `Y`. They
are pure `A` forms spanning a rank-one normal two-plane. Every other
form with nonzero quotient image is likewise pure `A`, by pairing with
one of these two. All such forms span at most two dimensions. Every
remaining form lies in `Y`; pairing it with a pure `A` form forces it to
be a multiple of one individual target scalar coordinate. A zero form
cannot meet either pair alternative. Let `r` be the number of distinct
scalar coordinates among these remaining forms. There are at most
three remaining forms, so `r≤3`. Their common kernel is

```
J = (an A rank-one two-plane) ⊕ K⁸_B ⊕ K⁸_C ⊕ K^(n+1−r),
```

with generic target slice rank `12+n+1−r=n+13−r`. The source upper bound
`n+10` forces `r=3`. Thus `n≥2`, the five source forms are independent,
and the only possible equality case has

```
J = (an A rank-one two-plane) ⊕ K⁸_B ⊕ K⁸_C ⊕ K^(n−2),
generic target slice rank on J = n+10.
```

It remains to use the source diagonal forms. Let `t₁,…,tₙ` be their
pullbacks restricted to `J`. None can vanish identically there: losing
one source diagonal block would bound every source slice on `J` by
`n+9`, below its generic target rank. For each nonzero `tⱼ`, the
hyperplane `ker(tⱼ)` must lie in the target rank-drop cone on `J`, since
source rank on it is at most `n+9`.

The only linear hyperplanes contained in that cone are the `n−2`
remaining target scalar-coordinate hyperplanes. Indeed, the core `A`
two-plane has constant rank one away from zero, so its rank-drop locus
has codimension two in `J`. The `B` and `C` rank-drop cones have maximal
linear subspaces of codimension four in `J`. None can contain a
hyperplane. Irreducibility of a linear hyperplane then forces
`ker(tⱼ)` to equal one of the remaining scalar-coordinate hyperplanes.

Moreover, two distinct source forms cannot choose the same hyperplane.
On that hyperplane at least two source diagonal blocks vanish, giving
source rank at most `n+8`. Exactly one independent target scalar block
is lost there, so the generic target rank is `n+9`, again a contradiction.
Thus the `n` source forms would require `n` distinct target hyperplanes,
while only `n−2` remain. This proves the exclusion for all `n`; the
first case and the impossibility of `r=3` cover the smaller `n` as well.

The root derivation and a separate proof-family audit were independently
checked here. The argument uses direct-sum block counts as integers and
formal coefficient restrictions, so it remains valid in characteristic
two. The earlier `D=0` and `D=Unit` exclusions are now special cases;
the `D=M₂` argument remains a separate result.

## A broader mixed-dot filter from the catalyst rank-drop cone

Keep the three-orientation dot auxiliary of the preceding section. Let
`D` be first-concise with first dimension `h` and generic first-slice rank
`ρ`, both computed over the algebraic closure. Define

```
Δ_D = {x : rank(D_x)<ρ}.
```

If every linear subspace contained in `Δ_D` has dimension at most `h−3`,
then **no positive-gain catalyst with this `D` exists**. Equivalently, a
potential catalyst must have a linear rank-drop slice subspace of
codimension one or two. This necessary subspace is asserted over the
algebraic closure; it need not have a basis over the original coefficient
field. Padding away nonconcise first coordinates does not evade the
condition: the tensor and its concise version restrict to each other,
and the putative certificate can be composed with these maps.

Here is an independent audit of the general proof. As above, pull back the
five singleton-first scalar forms `sᵢ` to the target first space, whose
dimension is `h+21` and whose generic slice rank is `ρ+17`. On every pair
kernel `ker(sᵢ,sⱼ)`, source rank is at most `ρ+16`, so this kernel lies in
the union of the five generic rank-drop sets: `Δ_D`, gain zero, or one of
the three core matrix rank-drop sets. The kernel has dimension at least
`h+19`. Its projection to the `D` first space has dimension at least
`h−2`, so it cannot lie in `Δ_D`. The `B` and `C` bad sets have maximal
linear dimensions `h+17`, so they cannot contain the kernel either.

Thus each pair's normal span contains the gain form `g` or is an
independent pure `A` rank-one normal two-plane. Let `Y=span(g)`.
If the images of the five forms modulo `Y` span at least two dimensions,
the same pair argument as above forces all their nonzero quotient forms
to be pure `A` forms, spanning at most two dimensions. The remaining
forms are multiples of `g`. Their common kernel therefore retains the
whole `D`, `B`, and `C` first spaces and at least two `A` dimensions.
It has generic target slice rank at least `ρ+4+4+4=ρ+12`.

If their quotient images span at most one dimension, their total span has
dimension at most two. Their common kernel then has codimension at most
two. Its `D` projection has dimension at least `h−2`, and hence attains
generic slice rank `ρ` by hypothesis. Its `A` projection has dimension
at least two, so it attains rank at least one. Its `B` and `C` projections
each have dimension at least six, exceeding the maximal dimension four
of their rank-one linear spaces; each therefore attains rank two.
These rank conditions hold simultaneously on a nonempty open subset of
the common kernel, giving generic target slice rank at least `ρ+12`
again. In either case all five source scalar forms vanish there, so
source rank is at most `ρ+10`, a contradiction.

This proof uses exact slice-rank inequalities and linear spaces over an
infinite field, with no tensor-rank additivity or finite-point inference.
It extends to every original field by scalar extension. It also applies
with any chosen catalyst leg first: the mixed-dot auxiliary is invariant
up to coordinate changes under leg permutations, as is `M₂`.

Some concrete applications need only the elementary form of the filter.
If `D` has constant positive slice rank on every nonzero first input and
`h≥3`, its rank-drop cone is just zero, and it is excluded. This includes
an outer-product tensor with at least three first coordinates. For
`D=C(a,b)`, multiplication by any nonzero first-input polynomial is
injective on the second-input polynomial space, so its slice rank is
always `b`. Thus any complete convolution catalyst with `a≥3`, or with
`b≥3` after exchanging the two input legs, is excluded for this mixed-dot
auxiliary.

Matrix-multiplication catalysts have a further application invoking a
standard bounded-rank matrix-space theorem. For `D=M(a,b,c)`, a first
slice is `c` copies of an `a×b` matrix, so `Δ_D` is the rank-at-most
`min(a,b)−1` matrix cone. Over an infinite field, Flanders' theorem bounds
the dimension of a linear subspace of this cone by

```
max(a,b)·(min(a,b)−1) = ab−max(a,b).
```

The reference is H. Flanders, “On spaces of linear transformations with
bounded rank,” *Journal of the London Mathematical Society* 37 (1962),
10–16, [doi:10.1112/jlms/s1-37.1.10](https://doi.org/10.1112/jlms/s1-37.1.10).
The algebraic closure satisfies the theorem's field-size requirement.
Consequently the filter excludes this catalyst whenever `max(a,b)≥3`.
Applying it on another leg excludes every `M(a,b,c)` with any parameter
at least three. This application invokes the cited theorem; the filter
itself does not require it. In the special case `M(3,2,2)`, the relevant
rank-one `3×2` spaces have dimension at most three by the elementary
common-row or common-image classification, so this example needs no
general bounded-rank theorem.

## The mixed-dot exclusions extend to every length `t≥2`

Let `S_t` be the direct sum of a length-`t` dot tensor in each of the three
singleton-leg orientations. It has shape `(1+2t)³` and a supplied exact
`3t`-term coordinate decomposition. The previous three kinds of exclusion
extend to **every `t≥2`, over every field**:

- Every finite diagonal catalyst `D=n Unit` is excluded.
- The catalyst `D=M₂` is excluded.
- Any first-concise `D` whose generic slice-rank-drop cone contains no
  linear subspace of codimension at most two is excluded.

The last condition may be checked on any leg, using the same leg symmetry
as above. These results concern actual restrictions with any positive
gain; projection reduces that gain to one unit. They do not assert
existence for the catalyst families remaining after these filters.

Here is the complete parameter audit. A source `S_t` slice has rank

```
t·1_(s≠0) + 1_(u≠0) + 1_(v≠0),
```

where `u,v∈K^t`. If `D` has generic slice rank `ρ`, the source slice
upper bound after killing two singleton-first scalar forms is
`ρ+3t+10`; after killing all five it is `ρ+10`. The three core first
matrix spaces have shapes `2×2`, `(2t)×2`, and `2×(2t)`, dimensions
`4,4t,4t`, and slice contribution

```
2t rank(A) + 2 rank(B) + 2 rank(C).
```

This follows directly from the matrix-column spectators in the product
with `M₂`. The target generic slice rank, including one gain unit, is
`ρ+4t+9`. It exceeds the pair-kernel source upper bound by `t−1`, which
is positive for the stated range. The target first dimension is
`h+8t+5`. A pair kernel has codimension at most two. A rank-one linear
subspace of either rectangular core matrix space has dimension at most
`2t`, so every linear subspace in its rank-drop cone has codimension at
least `2t≥4`.
Thus those two bad components cannot contain a pair kernel. The `A`
rank-drop component still forces a pure rank-one normal two-plane.

For the general catalyst rank-drop hypothesis, the proof above is
unchanged up to its final rank count. The pair kernels cannot project
into `Δ_D`, since their `D` projections have codimension at most two.
If the scalar forms modulo the gain form span at least two dimensions,
their common kernel retains full `D,B,C` and at least two `A` dimensions.
Its generic target rank is at least `ρ+2t+8`. Otherwise their common
kernel has codimension at most two; its `D` projection attains rank
`ρ`, its `A` projection attains rank at least one, and its `B,C`
projections have dimension at least `4t−2>2t`, hence rank two.
The same lower bound follows. Since `ρ+2t+8>ρ+10` for `t≥2`, both cases
contradict the source upper bound.

For `D=n Unit`, if the five scalar forms modulo the target scalar row
space span at most one dimension, every core projection has codimension
at most one and attains its full generic rank. At most five target scalar
coordinates vanish identically on the common kernel. Hence its generic
target slice rank is at least

```
4t+8+n+1−5 = n+4t+4 > n+10.
```

In the other case the five forms are pure `A` forms spanning a rank-one
normal two-plane, together with at most three pure target scalar-coordinate
forms. If these eliminate `r≤3` distinct scalar axes, the common kernel
has generic target slice rank

```
2t+8+n+1−r = n+2t+9−r.
```

The source upper bound `n+10` would require `r≥2t−1`. This is impossible
for `t≥3`. For `t=2` it forces exactly the equality case `r=3`, excluded
by the preceding `n` source hyperplanes versus `n−2` target hyperplanes
argument. Thus no diagonal catalyst survives at any stated length.

For `D=M₂`, retain its independent first matrix `E`, with source and
target slice contribution `2 rank(E)≤4`. The pair-kernel alternatives
are gain zero or pure rank-one normal two-planes in `A` or `E`.
If the five scalar forms modulo the gain form span at least two
dimensions, their nonzero quotient forms all lie in just one of these
matrix blocks, spanning at most two dimensions. The common kernel
therefore has generic target slice rank at least `2t+12` if that block
is `A`, and at least `4t+10` if it is `E`.

Otherwise the common kernel has codimension at most two. Its `A,E`
projections each have dimension at least two and attain rank at least
one; its `B,C` projections attain rank two. The two former projections
cannot both be rank-one matrix spaces: that would require at least two
pure `A` normals and two independent pure `E` normals, contradicting
codimension at most two. At least one therefore attains rank two.
The generic target rank is again at least
`min(2t+12,4t+10)=2t+12`. Every case exceeds the source common-kernel
upper bound `4+10=14` for `t≥2`.

All simultaneous generic-rank statements use intersections of nonempty
open subsets of a linear space over the algebraic closure. The rank sums
are ordinary matrix slice ranks in disjoint blocks, not an additivity
assumption for tensor rank. The coefficient identities and block counts
therefore give all-field exclusions, including characteristic two.

## Exact single-copy subrank of the mixed-dot family

For every `t≥1`, the ordinary single-copy subrank is **`Q(S_t)=3` over
every field**. This is a single-copy statement, with no asymptotic-subrank
claim.

Order the three direct-sum blocks by singleton first, second, and third
leg, using the original dot support `(0,j,j)`, `0≤j<t`, in each permuted
block. Zero-based coordinate selectors for three independent units are

```
first:  (0,1,t+1),
second: (0,t,t+1),
third:  (0,t,2t).
```

They select exactly one supported triple in each block, and every mixed
selected triple has coefficient zero. Thus the selected tensor is the
three-unit diagonal tensor over the integers. For `t=2`, the portable
[`export_mixed_dot_subrank_order`](../mixed_dot_orders.py) exports these
actual maps, with controls in
[`test_mixed_dot_orders.py`](../test_mixed_dot_orders.py).

For the upper bound, suppose four independent units restrict from `S_t`,
and extend the maps to the algebraic closure. Pull back the source
singleton-first scalar form to the four-dimensional target diagonal
first space. Its kernel has dimension at least three. On this kernel
the source singleton-first dot block disappears, while each other block
has slice rank at most one. Hence every target diagonal slice on this
kernel has rank at most two. The rank-at-most-two cone of a four-entry
diagonal matrix is the finite union of its coordinate two-planes.
An irreducible linear space of dimension at least three cannot be
contained in that union, a contradiction. A restriction to more than
four units would project to four, so the upper bound follows.

Independent temporary controls checked the displayed lower restriction
for `t=1,2,3` over `F₂,F₃,F₄,F₅`, twelve cases in total. The integer
coefficient identity and linear-space proof establish the all-field
statement. The lower maps supply the legitimate ordinary state row
`3≤f(S_t)` for normalized additive monotone states. The subrank upper
bound is not a state upper row and does not imply `f(S_t)≤3`.

## Two same-orientation dot-two catalysts are excluded

This is a stepping stone to the
[all-finite-copy theorem](#every-finite-same-orientation-dot-catalyst-is-excluded)
below, which supersedes this two-copy screen.

For `S=S_2`, the catalyst consisting of **two copies of a singleton-first
length-two dot tensor** is impossible over every field. Its shape is
`(2,4,4)` and its source slice rank is at most four. By leg permutation,
the result also holds for two dots both having the singleton second leg,
or both having the singleton third leg. This proof does not cover a
mixture of different dot orientations.

Write its two independent target first scalars as `d₁,d₂`, whose slice
contributions are two each, and retain the gain scalar `g` of weight one.
As before, any pair of the five source scalar normal forms contains one
of these individual scalar-coordinate forms or spans a pure `A` rank-one
normal two-plane. If their images modulo `span(d₁,d₂,g)` span at most
one dimension, their common kernel attains the full generic core rank
sixteen. This exceeds the source bound `4+10=14` there.

Otherwise the nonzero quotient forms span a pure `A` rank-one normal
two-plane; every other form is a pure scalar-coordinate form. Let `r≤3`
be the number of distinct eliminated target scalar axes. The surviving
core has generic rank twelve. For `r≤1`, surviving scalar weight is at
least three, giving target rank at least fifteen, again above fourteen.

If `r=2`, exactly one target scalar axis remains. The common kernel has
generic target rank fourteen if that axis is a dot scalar, and thirteen
if it is the gain. Let `t₁,t₂` be the two source catalyst scalar forms
restricted to this kernel. Neither can vanish identically, since that
would give source rank at most twelve throughout. Killing either reduces
source rank by two, to at most twelve. Its hyperplane must therefore lie
in the target rank-drop cone, and the only available linear hyperplane
is the one remaining scalar-coordinate hyperplane. The core `A` drop
has codimension two and every linear subspace in a `B,C` bad cone has
codimension at least four. Thus both `t₁,t₂` have that same kernel. Killing both
gives source rank at most ten while the target core still has generic
rank twelve, a contradiction.

It remains to exclude `r=3`, when all target scalar axes vanish and the
common kernel is exactly the whole `B,C` first spaces plus an `A`
rank-one two-plane. Here ordinary tensor flattening concision supplies
the missing obstruction. On the source, all five singleton-first dot
blocks disappear. Each remaining copy is the direct sum of a shape
`(2,1,2)` outer-product block and a shape `(2,2,1)` dot block, whose
second and third concise dimensions are three each. Adding the two
source catalyst dots gives the other-leg upper bounds

```
5·(3,3) + (4,4) = (19,19).
```

For the target `A` block, a rank-one plane is equivalent to one of the
following actual coordinate restrictions of `M(2,2,4)`:

```
A=[[x,0],[y,0]]:  A-first selectors (0,2), other flattenings (4,8),
A=[[x,y],[0,0]]:  A-first selectors (0,1), other flattenings (8,4).
```

In the first case the four second-input coordinates in the first matrix
row produce eight independent output coordinates `(xw,yw)`. In the
second case all eight second-input coordinates are used but only four
output coordinates remain. Common-row or common-image coordinate
changes preserve these concision values for every rank-one two-plane.
The full `B` and `C` core blocks have other flattenings `(4,8)` and
`(8,4)`, respectively. Because all three blocks have disjoint coordinates
on each leg, their actual flattening ranks add to

```
(16,20) or (20,16).
```

One target flattening is therefore twenty, while the corresponding
source flattening is at most nineteen. Flattening ranks cannot increase
under the remaining restriction maps, so the last case is impossible.
These are matrix flattening ranks, without a tensor-rank additivity
assumption. Eight independent temporary controls checked the two `A`
coordinate concision calculations over `F₂,F₃,F₄,F₅`; the coordinate
span calculation proves the field-independent values.

## A second-leg exclusion for one larger-dot-catalyst boundary

The [all-finite-copy theorem](#every-finite-same-orientation-dot-catalyst-is-excluded)
below now subsumes this boundary calculation.

For the mixed-dot auxiliary `S_2`, consider `D=n DotX₂`, the direct sum
of `n` singleton-first length-two dots. The argument here excludes one
specific boundary case for every `n≥3`; it does **not** establish an
all-`n` exclusion for this whole catalyst family.

Suppose the five source singleton-first scalar normal forms span a pure
`A` rank-one normal two-plane together with three distinct target `D`
scalar-coordinate axes, while leaving the gain scalar uneliminated.
This is the `r=3` case in which all three eliminated scalar axes belong
to `D`. On their common first kernel the target is the direct sum of
the three core blocks, with `A` restricted to a rank-one two-plane, one
gain unit, and `n−3` remaining singleton-first dots. The source first
restriction removes all five singleton-first auxiliary blocks. The
remaining source is a restriction of

```
n DotX₂ ⊕ 5(DotY₂ ⊕ DotZ₂).
```

Analyze slices on the second leg of this new restriction. The source
second-slice rank is at most `n+15`: each source catalyst dot contributes
at most one, while each auxiliary copy contributes two from its
singleton-second block and one from its singleton-third block. Pull
back the five source singleton-second scalar forms as `α₁,…,α₅` on the
target second space. Killing any two gives source rank at most `n+11`;
killing all five gives at most `n+5`.

The target second-slice contributions of the core blocks are as follows:

| Block | Generic second-slice rank | Minimum codimension of a linear rank-drop subspace |
| --- | --- | --- |
| `A` common-row plane | 2 | 4: a four-vector must vanish |
| `A` common-image plane | 2 | 4: a `2×4` matrix drops to rank one |
| `M(4,2,2)` block | 8 | 2: a `2×2` matrix drops to rank one |
| `M(2,4,2)` block | 4 | 4: a `4×2` matrix drops to rank one |

In the second table column for `M(4,2,2)`, the slice rank is four times
the rank of its second matrix. For `M(2,4,2)` it is twice that rank.
Each remaining target dot has generic second-slice rank one and drops
only when its second two-vector vanishes. Including the gain, the
generic target second-slice rank is therefore

```
2+8+4+1+(n−3) = n+12.
```

The pair kernels of the `αᵢ` have codimension at most two and must lie
in this generic rank-drop cone, because source rank there is at most
`n+11`. The `A` and `M(2,4,2)` bad components cannot contain such a
linear space. Irreducibility leaves only three pair-normal alternatives:
the span contains the gain form; it is a pure rank-one normal two-plane
in the `M(4,2,2)` second matrix; or it is the whole pure two-coordinate
normal plane of one remaining dot's second vector.

If the five `αᵢ` modulo the gain form span at least two dimensions,
choose a pair with independent quotient images. All further nonzero
quotient forms must belong to the same pure block as this pair; forms
with zero quotient image are multiples of the gain. In the matrix case
their span is a rank-one normal two-plane; in the dot case it is that
dot's full two-coordinate plane. Killing all five leaves generic target
second-slice rank at least

```
2+4+4+(n−3) = n+7        in the matrix case,
2+8+4+(n−4) = n+10       in the remaining-dot case.
```

Both exceed the source bound `n+5`. The remaining-dot case is possible
only when there is at least one remaining dot.

Otherwise their quotient span has dimension at most one, so the common
second kernel has codimension at most two. The `A` block still attains
second-slice rank two: its common-row four-vector projection cannot
vanish, and its common-image matrix projection has dimension at least
six, greater than the maximum four of a rank-one `2×4` matrix space.
The `M(4,2,2)` second projection has dimension at least two, so it
attains rank at least one and contributes at least four. The
`M(2,4,2)` second projection has dimension at least six and contributes
four. At most one remaining dot's second vector can vanish identically,
because each needs two independent pure normals and the whole kernel
has codimension at most two. The generic target rank is consequently
at least `2+4+4+(n−4)=n+6`, again above `n+5`.

This exhausts the second-leg possibilities and excludes the stated
first-kernel boundary. Both orientations of the `A` rank-one plane are
covered explicitly. The proof extends to every coefficient field via
its algebraic closure. The first-leg case in which the five scalar
forms have quotient span at most one remains outside this argument;
no conclusion for all finite `n DotX₂` is made here.

## The three same-orientation dot-two catalyst is also excluded

The [all-finite-copy theorem](#every-finite-same-orientation-dot-catalyst-is-excluded)
below now supersedes this three-copy screen.

For `S=S_2`, the previous partial arguments and one final flattening
step prove that **`D=3 DotX₂` is impossible over every field**. Here
`D` is exactly three copies of the singleton-first length-two dot, with
concise shape `(3,6,6)` and generic first-slice rank six. Its source and
target shapes before restricting the first space are `(28,31,31)` and
`(24,27,27)`, respectively. Leg symmetry gives the same exclusion for
three dots all singleton on the second leg or all singleton on the
third leg. This result does not assert an exclusion for arbitrary
larger numbers of dots or for mixtures of their orientations.

Let the three target catalyst scalars be `d₁,d₂,d₃`, each of slice weight
two, and let the gain form `g` have weight one. Put
`Y=span(d₁,d₂,d₃,g)`. The five pulled-back source scalar forms are
`c₁,…,c₅`. Killing any pair gives source first-slice rank at most
`6+16=22`, below the generic target rank `6+17=23`. The pair-normal
alternatives are therefore the familiar individual target scalar axis
or pure `A` rank-one normal two-plane. Their common kernel `J` has
source slice-rank upper bound `6+10=16`.

First consider the case where the `cᵢ` modulo `Y` span at least two
dimensions. Their matrix forms span a pure `A` rank-one normal
two-plane; their remaining forms eliminate `r≤3` individual scalar
axes. The restricted target core has generic rank twelve.

- If `r≤1`, its surviving scalar weight is at least five, giving target
  rank at least seventeen, above sixteen.
- If `r=2` and two catalyst axes disappear, the gain and one catalyst
  axis remain, giving generic target rank fifteen. If instead one
  catalyst axis and the gain disappear, two catalyst axes remain,
  giving generic rank sixteen. In either situation all three source
  catalyst forms restricted to `J` are nonzero. Killing any one loses
  two source rank units, so its hyperplane must be a remaining target
  scalar-coordinate hyperplane; the core bad sets cannot contain a
  hyperplane. Two source forms choosing the same axis would lose four
  source rank units but at most two target rank units, a contradiction.
  Three source forms would therefore need three distinct axes, while
  only two remain.
- If `r=3` and the gain disappears, exactly two catalyst axes also
  disappear. The remaining catalyst dot adds `(2,2)` to the restricted
  core other-leg flattenings `(16,20)` or `(20,16)`. Thus a target
  flattening is twenty-two, whereas the source after all five `cᵢ`
  vanish has other flattenings at most `5·3+3·2=21`.
- If `r=3` and the gain remains, all three catalyst axes disappear.
  This is exactly the preceding second-leg exclusion, with `n=3`.

It remains to analyze quotient span zero or one. Every core projection
of `J` has codimension at most one and retains its full generic slice
rank, so generic core rank on `J` is sixteen. Since source rank there
is at most sixteen, every target catalyst and gain scalar must vanish
identically on `J`. Hence the five-form normal span `W` contains all
four dimensions of `Y`. There are two possibilities:

```
W=Y,                 J=the whole 20-dimensional core first space;
W=Y⊕span(λ),         J=ker(λ) in that core, dimension 19,
```

where in the second line `λ` is a nonzero pure core form. Let `t` be
any of the three source catalyst scalar forms restricted to `J`. It is
nonzero: otherwise source rank is at most fourteen throughout `J`,
contradicting its generic target rank sixteen. On `H=ker(t)` the source
rank is at most fourteen, so `H` lies in the core generic rank-drop cone.
The `B,C` bad cones have maximal linear dimension sixteen in the full
twenty-dimensional core first space, and the `A` bad cone has maximal
linear dimension eighteen.

If `W=Y`, then `dim H=19`, so irreducibility and those dimension bounds
give an immediate contradiction. In the other case `dim H=18`. It
must lie in the `A` bad cone, and equality of dimensions forces

```
H = (an A rank-one two-plane) ⊕ K⁸_B ⊕ K⁸_C.
```

In particular, `λ` and the new hyperplane normal are pure `A` forms
spanning a rank-one normal two-plane; the possible triangular-matrix
normal form is justified by this argument, rather than assumed.

The target on `H` has other flattenings `(16,20)` or `(20,16)`, by the
actual coordinate concision calculation above. The source has all five
auxiliary singleton-first blocks zero and at least one of its three
catalyst dot blocks zero. Its corresponding other flattenings are at
most

```
5·(3,3) + 2·(2,2) = (19,19).
```

A flattening twenty cannot restrict from a flattening of rank at most
nineteen. This excludes the last possibility and proves the full
three-copy result. The source block counts, first-slice inequalities,
and normal-form concision values are all coefficient constructions over
the algebraic closure, so the proof applies to every original field and
every positive gain, after projection to one gain unit.

## Every finite same-orientation dot catalyst is excluded

For the mixed-dot auxiliary `S=S_2`, **no finite catalyst `D=n DotX₂`
exists, for any `n≥0`, over any field**. Here `DotX₂` is the length-two
dot tensor with its first leg singleton, so `D` has concise shape
`(n,2n,2n)` and generic first-slice rank `2n`. Leg symmetry gives the
same theorem for any finite number of dots all singleton on the second
leg or all singleton on the third leg. Different orientations mixed
inside `D` are outside this theorem. Any positive gain projects to one
gain unit. This is independently audited slice-space mathematics over
the algebraic closure, with no Lean formalization claimed.

The earlier two-copy, three-copy, and second-leg boundary arguments are
retained above as stepping stones. The following branch analysis closes
the entire finite family. The zero-copy case is already excluded by the
small-catalyst theorem, so assume `n≥1` below.

Use the previous core matrices `A,B,C`, with generic slice rank sixteen.
The target first dimension is `n+21`; the source first dimension is
`n+25`. Its target scalar row space `Y` has `n` independent catalyst
coordinate forms of slice weight two and one gain form of weight one.
Let `c₁,…,c₅` be the pulled-back scalar forms of the five source
singleton-first auxiliary blocks, and set

```
W=span(c₁,…,c₅),       J=ker(W).
```

Source rank on `J` is at most `2n+10`. Killing any pair of the `cᵢ`
gives source rank at most `2n+16`, below the generic target rank
`2n+17`. Thus every pair normal span contains an individual target
scalar-coordinate form or is a pure `A` rank-one normal two-plane.
The `B,C` bad components have linear codimension at least four and
cannot contain a pair kernel.

If the images of the `cᵢ` modulo `Y` span at least two dimensions,
the same pair argument forces all nonzero quotient forms into one pure
`A` rank-one normal two-plane. Forms with zero quotient image are pure
individual scalar-coordinate forms. Let `r_D` count distinct eliminated
catalyst axes and `r_g∈{0,1}` record whether the gain disappears. Then
`r_D+r_g≤3`, and the restricted target core has generic rank twelve.
The complete possibilities are:

| Eliminated axes | Generic target rank or decisive bound | Exclusion |
| --- | --- | --- |
| At most one | At least `2n+11` | Source upper bound is `2n+10` |
| Two catalyst axes | `2n+9` | Source catalyst forms cannot fit the remaining `n−1` scalar hyperplanes |
| One catalyst axis and gain | `2n+10` | Source catalyst forms cannot fit the remaining `n−1` scalar hyperplanes |
| Two catalyst axes and gain | A target flattening is `2n+16` | Source other flattenings are at most `2n+15` |
| Three catalyst axes, gain retained | Second-leg recursion | The preceding boundary proof applies |

For clarity, in either two-axis case every source catalyst scalar form
is nonzero on `J`. Killing it reduces source rank by two and forces its
hyperplane into the generic target rank-drop cone. The restricted `A`
plane drops only at zero, a codimension-two locus, and `B,C` bad linear
spaces have codimension at least four. Hence it must be a remaining
target scalar-coordinate hyperplane. Two source forms cannot select the
same axis: source would lose four rank units, while target loses at
most two and its initial source surplus is at most one. The `n` source
forms would need `n` distinct axes, but only `n−1` remain. The other
rows are exactly the concision and second-leg proofs already given.

It remains to handle quotient span zero or one. The core projection of
`J` has codimension at most one and retains generic core rank sixteen.
Define `p` as the number of **actual target catalyst coordinate forms**
vanishing identically on `J`, and `q∈{0,1}` according to whether the
gain vanishes identically. These count individual axes, not dimensions
of a possibly noncoordinate scalar normal space. Its generic target
rank is therefore

```
G = 16+2(n−p)+1−q = 2n+17−2p−q.
```

Since source rank is at most `2n+10`, one needs `2p+q≥7`.
If the quotient span is one, `dim(W∩Y)≤4`; if it is zero,
`dim W≤5`. This gives only the following branches:

| Quotient dimension | `(p,q)` | Generic target rank `G` | Source surplus `2n+10−G` |
| --- | --- | --- | --- |
| One | `(3,1)` | `2n+10` | 0 |
| One | `(4,0)` | `2n+9` | 1 |
| Zero | `(3,1)` | `2n+10` | 0 |
| Zero | `(4,0)` | `2n+9` | 1 |
| Zero | `(4,1)` | `2n+8` | 2 |
| Zero | `(5,0)` | `2n+7` | 3 |

Let `t₁,…,tₙ` be the source catalyst scalar forms restricted to `J`.
Whenever they all vanish, the whole source has slice rank at most ten.

First take quotient dimension zero and source surplus zero or one.
Then `J` is the direct product of the **whole core** with a scalar
subspace `L`. This `L` may contain an additional noncoordinate relation;
different surviving scalar-coordinate restrictions may even be
proportional. Each `tⱼ` is nonzero on `J`, since losing two source rank
units would put source rank below `G` throughout `J`. Its hyperplane
must lie in the target rank-drop cone. The full core `A` bad set has
linear codimension at least two and `B,C` have at least four, so none
can contain a hyperplane in this product. Thus `tⱼ` is a pure scalar
form on `L`, proportional to a surviving coordinate restriction.
Killing all the `tⱼ` retains the whole core, of generic rank sixteen,
while source rank is at most ten. No independence or distinctness of
the surviving scalar axes is used.

Next take quotient dimension zero and source surplus two or three.
Killing any two source catalyst forms reduces source rank to at most
`2n+6`, below `G`. A pair kernel must therefore lie in the target
rank-drop cone. The pair normals contain a nonzero actual scalar
coordinate restriction or span a pure `A` rank-one normal two-plane;
the `B,C` components again cannot contain a codimension-two linear
space. Work modulo the **full scalar row space `L*`**, allowing all
relations among the coordinate restrictions. If the `tⱼ` images modulo
that space span at most one dimension, their common kernel retains a
core projection of codimension at most one and generic rank sixteen.
Otherwise a pair with independent quotient images forces all nonzero
quotient forms into the same pure `A` rank-one normal two-plane; the
remaining forms are pure scalar. Their common kernel retains at least
two `A` dimensions and the whole `B,C` spaces, giving core rank at
least twelve. Both bounds exceed the source rank ten after all source
catalyst forms vanish.

Finally take quotient dimension one. There are exactly four eliminated
coordinate axes, and their span is precisely `W∩Y`. After removing
them, the ambient first space `U` is the direct product of the core
and all surviving scalar coordinates. Its remaining normal is of the
form

```
w=λ(core)+y(scalars),       λ≠0,       J=ker(w) inside U.
```

**Counting alone does not force `y=0`.** The fifth normal may mix core
and surviving scalar coordinates. Nevertheless, the surviving scalar
coordinate restrictions on `J` are independent, because any relation
among them would enlarge `W∩Y` beyond the four eliminated axes.

For either surplus zero or one, every `tⱼ` is nonzero: an identically
zero form would lose two source rank units against a surplus smaller
than two. Its hyperplane lies in a generic rank-drop component. A scalar component
gives a pure scalar-coordinate restriction. If an `A` component
contains it, equality of maximal linear dimensions forces that
hyperplane to be the product of an `A` rank-one two-plane with the
whole `B,C` spaces and **all** surviving scalar coordinates. Its
annihilator in `U` is consequently pure `A`. In particular, this
alternative forces `y=0` and `λ` to be a rank-one pure `A` normal.
The `B,C` alternatives are too small.

In branch `(p,q)=(3,1)`, an `A`-type hyperplane has a target other-leg
flattening

```
20+2(n−3)=2n+14,
```

whereas the source after killing its source catalyst form has other
flattenings at most `15+2(n−1)=2n+13`. Thus all `tⱼ` must be scalar
types. Killing all retains at least `ker(λ)` in the core, which has
codimension one and generic core rank sixteen. Source rank is ten,
giving a contradiction even when the fifth normal was mixed.

In branch `(p,q)=(4,0)`, scalar-type source forms cannot repeat the
same surviving axis: two copies would lose four source rank units
against at most two target rank units, with initial surplus only one.
There are `n−4` surviving catalyst axes and the gain, so at most `n−3`
source forms are scalar types. At least three must therefore be
`A` types. Their existence forces the remaining normal pure `A` and
rank one, as just proved. By coordinate changes the core `A`
hyperplane is the upper-triangular three-space

```
A=[[a,b],[0,d]],       det(A)=ad.
```

Its only linear rank-drop hyperplanes are `a=0` and `d=0`. Thus the
three source `A`-type forms span at most two pure `A` normal directions.
Killing any three of them loses six source rank units, giving source
rank at most `2n+4`. It leaves a nonzero `A` direction, the whole
`B,C` spaces, and all surviving scalar coordinates. Its generic
target rank is at least

```
4+4+4+2(n−4)+1 = 2n+5,
```

a final contradiction.

This completes every branch. The scalar-coupling distinction is essential:
quotient-zero branches retain a full core and permit arbitrary scalar
relations, while quotient-one branches permit a mixed fifth normal
until an `A`-type maximal linear space forces that normal pure `A`.
The proof uses algebraic-closure linear spaces, formal restrictions, and
ordinary matrix slice or flattening ranks. It gives no assertion about
an optimal finite-field tensor-rank decomposition.

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
