# Shared-dot catalyst obstructions and coefficient controls

Checkpoint: 2026-10-09. For every field, every q>=1, every positive gain m,
and every finite catalyst D, the shared tensor

```
C_q = sum_i(e0 tensor ei tensor ei + ei tensor e0 tensor ei + ei tensor ei tensor e0)
```

**cannot** satisfy `D + m Unit + M2 tensor C_q <= D + 5 C_q`.
These independently audited arguments are not Lean theorems. Integral
certificates and polynomial identities below are replayed by Python; the
unbounded catalytic iteration is a mathematical argument.

| q | Supplied polynomial source slope | Whole-target lower control | Reason |
| --- | ---: | ---: | --- |
| 1 | 10 | m+12 | Two-term W degeneration and normalized commutator |
| 2 | 20 | 124/6 | Unimodular integer Koszul transformation |
| 3 | 25 | 153/6 | Determinant-one integer Koszul minor |
| >=4 | 5(q+2) | m+6(q+1) | q+2-term degeneration and normalized commutator |

The rational-function transfer and shared-first repetition proved below
compare n times the lower control against one fixed finite D cost plus
n times the source slope. Every row has a strict gap for positive gain.
The integer coefficients and unimodular operations work in every
characteristic; no inference from testing several primes is needed.

## Closing q=1 and q=2

For q=1, `(e1+t e0)^tensor3-e1^tensor3` has leading coefficient `t C1`
and two terms. The invertible slice `(1,1)` has determinant minus one.
The two normalized matrix-unit slices have commutator rank eight, so
n whole-target copies have rank lower n(m+12). The polynomial source
upper is C_D+10n. Thus C_D>=n(m+2), impossible for unbounded n.
`shared_dot_w_degeneration` checks this actual two-term family formally.

For q=2, the binary five-row map in `shared_dot_two_koszul.json` projects
`Unit + M2 tensor C2` to a five-dimensional first leg. Rebuilding its
signed exterior-degree-two matrix gives a 130 by 130 integer matrix.
The 903 stored elementary integer row/column operations are all swaps
or additions of an integer multiple of a distinct row/column. They are
unimodular, hence invertible over every field. Replaying them gives a
triangular leading 124 by 124 block with diagonal +/-1 and determinant
minus one. Independent fraction-free and symbolic determinants agreed.
This is a unit minor of the transformed matrix; no claim that a single
original minor is a unit is made.

A rank-one tensor contributes exterior-matrix rank at most six. For
n target copies merge only their projected first coordinates, preserving
both other copy indices: the exterior matrix becomes block diagonal,
with rank at least 124n. Its polynomial approximation has at most
C_D+20n terms. Consequently

```
124n <= 6(C_D+20n),   4n <= 6C_D,
```

which again contradicts a fixed finite C_D. Larger gains restrict to one.
The portable `verify_shared_dot_two_integer_transform` reconstructs all
49 target support terms and every signed matrix entry, replays the stored
unimodular operations, and checks the unit triangular block. It requires
only the standard library. The original binary-map field ranks are 126
in F2,F4,F5,F11 and 124 in F3; those samples are supplementary controls.

## Independent slice-space and transfer arguments

The following derivations also establish all-field D=0 geometry and the
commutator obstruction. The C3 unit-minor argument closes the equality
case that the commutator alone leaves unresolved.

## Source and target slices

A first slice of `C_q`, with scalar coordinate `s` and vector `v∈K^q`, is

```
M(s,v) = [[0, vᵀ], [v, s I_q]].
```

It has rank at most `q+1`, and rank at most two when `s=0`. Its
determinant is the integer polynomial

```
det M(s,v) = −s^(q−1) Σ_i v_i².
```

In characteristic two the vector quadratic becomes `(Σ_i v_i)²`.
This characteristic dependence matters in the geometry below.

After separating the two matrix-column spectators of `M₂`, a first
slice of `M₂⊗C_q` is equivalent to `I₂⊗B`, where

```
B = [[0,   V₁, V₂, …, V_q],
     [V₁, A,  0,  …, 0  ],
     [V₂, 0,  A,  …, 0  ],
     [⋮,   ⋮,  ⋮,  ⋱, ⋮  ],
     [V_q,0,  0,  …, A  ]].
```

Here `A,V₁,…,V_q` are arbitrary `2×2` matrices. Thus the core slice rank
is `2 rank(B)`. The core first space has dimension `4(q+1)` and its
generic slice rank is `4(q+1)`: for example, `A=I₂,V₁=I₂,V_i=0` for
`i>1` gives invertible `B` in every characteristic.

For `q≥2`, let

```
W = Σ_i V_i adj(A) V_i,       F_q = det W.
```

Schur complementation when `det A≠0` gives

```
det B = (det A)^(q−2) F_q.
```

Both sides are polynomials over the integers. Their equality on the
generic invertible locus proves the formal identity in every
characteristic. The polynomial `F_q` has degree two in `A` and degree
four in the `V` coordinates.

## All-characteristic small-codimension lemma

For every field and `q≥2`, the cone `F_q=0` contains **no linear
subspace of codimension at most two** in the full `(A,V)` space.
Work over its algebraic closure.

For fixed invertible `A`, normalize to `A=I₂` using the covariance
proved below. Suppose a matrix-tuple subspace `N` of codimension at most
two had `det(Σ_i V_i²)=0` identically. Intersect `N` with the
`3q`-dimensional space of traceless tuples

```
V_i=[[a_i,b_i],[c_i,−a_i]].
```

This intersection has dimension at least `3q−2`. On that space

```
Σ_i V_i² = Q I₂,       Q=Σ_i(a_i²+b_i c_i),
```

so determinant zero forces `Q=0` identically on the intersection.

In odd characteristic, `Q` has nondegenerate polar form on its
`3q` coordinates. A linear space on which `Q` vanishes is totally
isotropic for that polar form, and has dimension at most
`floor(3q/2)`. For `q≥2`, this is smaller than `3q−2`.

In characteristic two, the polar form is the symplectic form on the
`2q` coordinates `(b,c)` and has radical the `q`-dimensional `a` space.
On that radical, `Q=(Σ_i a_i)²`; its zero subspace has dimension `q−1`.
For a linear zero space `L`, its radical intersection therefore has
dimension at most `q−1`, and its projection to the symplectic quotient
has dimension at most `q`. Thus `dim L≤2q−1`, again smaller than
`3q−2` for `q≥2`. This proves the fixed-invertible-`A` lemma in every
characteristic.

The full-space adaptation is also characteristic-free. Given
`H⊂{F_q=0}` of codimension at most two, its kernel in the `V` variables
has codimension at most two. If its `A` projection contains an
invertible `A₀`, the highest quartic coefficient on the affine fiber
at `A₀` gives precisely the fixed-`A₀` contradiction. Otherwise its
`A` projection is a singular matrix space of dimension at most two;
dimension equality forces the entire `V` space into that kernel. At
any nonzero rank-one `A`, the two explicit matrices below make
`Σ_i V_i adj(A)V_i` invertible, a contradiction. Normalization and
the explicit affine-fiber calculation are given in the historical
proof below and require no invertibility of two.

## Historical scalar-shift proof in odd characteristic

Assume `char K≠2` and `q≥3`, with `K` algebraically closed. Then
**the cone `F_q=0` contains no linear subspace of codimension at most
two** in the full `(A,V₁,…,V_q)` space.

First fix an invertible `A`. It is enough to prove the assertion for

```
f(V₁,…,V_q) = det(Σ_i V_i²).
```

The normalization is legitimate. For arbitrary invertible `P,Q`, the
changes `A'=P A Q`, `V_i'=P V_i Q` give

```
W' = det(P)det(Q) P W Q,
F_q' = (det(P)det(Q))³ F_q.
```

This follows from `adj(PAQ)=adj(Q)adj(A)adj(P)`. Taking `P=A⁻¹,Q=I`
normalizes the fixed `A` to the identity without changing its zero set.

Suppose a linear subspace `N` of the `4q` matrix coordinates has
codimension at most two and `f` vanishes identically on it. The scalar
tuple space `{(α₁I,…,α_qI)}` has dimension `q`, so its intersection with
`N` has dimension at least `q−2≥1`. Choose a nonzero scalar tuple
`αI` in that intersection. Evaluating `f` there gives
`(Σ_i α_i²)²=0`; hence `κ=Σ_i α_i²=0`.

For arbitrary `V∈N`, the whole line `V+tαI` lies in `N`. Set

```
Q_V=Σ_i V_i²,       L_V=Σ_i α_i V_i.
```

Because scalar matrices commute,

```
Σ_i (V_i+tα_iI)² = Q_V+2t L_V+t²κI = Q_V+2t L_V.
```

The coefficient of `t²` in its determinant is `4 det(L_V)`. Odd
characteristic therefore forces `det(L_V)=0` for every `V∈N`.
The linear map `V↦L_V` is surjective onto the four-dimensional matrix
space, since `α≠0`. Its image on `N` has dimension at least two.
A linear space of singular `2×2` matrices has dimension at most two,
by the common-row or common-image classification. Dimension equality
forces `N` to contain the whole kernel of `V↦L_V`.

The standard dot quadratic `Q(β)=Σ_i β_i²` is nondegenerate in odd
characteristic. Since `α` is a nonzero isotropic vector, its restriction
to `α⊥` has radical exactly `Kα` and rank `q−2≥1`. Choose
`β∈α⊥` with `Q(β)≠0`, and choose an invertible matrix `U`. The tuple
`V_i=β_i U` belongs to the kernel just found, but

```
f(β_i U) = Q(β)² det(U)² ≠ 0,
```

a contradiction. This proves the fixed-invertible-`A` lemma.

Now suppose a full linear space `H⊂{F_q=0}` has codimension at most
two. Let `L` be its `A` projection and `N` its kernel in the `V`
coordinates. Then `dim N≥4q−2`. If `L` contains an invertible `A₀`,
its fiber is an affine translate of `N`. For `V∈N`, take the highest
quartic coefficient in `F_q(A₀,V₀+tV)=0`; it gives
`F_q(A₀,V)=0`. The fixed-invertible-`A` lemma contradicts this.

Otherwise `L` is a singular matrix space and has dimension at most two.
Since `dim H≥4q+2`, equality forces `dim L=2` and `N` to be the whole
`V` space. Pick a nonzero rank-one `A∈L` and normalize it to
`diag(1,0)` using the covariance already proved. With

```
V₁=[[0,1],[1,0]],       V₂=diag(0,1),       V_i=0 for i>2,
```

one gets `Σ_i V_i adj(A)V_i=I₂`, so `F_q=1`. This contradicts the
whole affine fiber: its highest quartic coefficient would force
`F_q(A,V)=0` for all `V`. The full lemma follows.

This historical scalar-shift step uses the invertibility of four and
does not apply in characteristic two. The preceding traceless-quadratic
argument supplies the characteristic-two proof.

## All-characteristic `D=0` exclusion

Assume a putative gain-one restriction for `D=0,q≥3`. Pull back the
five singleton-first source scalars as forms `c₁,…,c₅` on the target
first space, and let `g` be its gain form. The target first dimension
is `4q+5`. A pair kernel has codimension at most two. Its source slice
rank is at most

```
3(q+1)+2·2 = 3q+7 < 4q+5,
```

below the generic target rank. Hence it lies in the generic rank-drop
cone. By the determinant formula this cone is the union of gain zero,
`det A=0`, and `F_q=0`. Irreducibility of a linear space chooses one
component. The small-codimension lemma excludes the `F_q` component.
If it lies in `det A=0`, equality of maximal linear dimensions forces
its normals to span a pure `A` rank-one two-plane. Otherwise its normal
span contains `g`.

If the images of the five forms modulo `Kg` span at least two
dimensions, choose a pair with independent quotient images. Both are
pure `A` forms. Every other nonzero quotient form is forced pure `A`
by pairing with one of them. These forms span at most a rank-one
normal two-plane; the remaining forms are multiples of `g`. Their
common kernel therefore retains all `V` variables and an `A` rank-one
two-plane.

For any fixed nonzero rank-one `A`, the generic rank of the reduced
matrix `B` with unrestricted `V` is exactly `q+4`. To verify this,
normalize `A=diag(1,0)`. Its kernel equations are

```
Σ_i V_i y_i=0,       V_i x+A y_i=0 for every i.
```

Generically the rows `e₂ᵀV_i` span the two-dimensional row space, forcing
`x=0`. Then `y_i=z_i e₂` and `Σ_i z_i V_i e₂=0`, with a generic rank-two
map from `K^q` to `K²`. The kernel dimension is `q−2`, giving rank
`2q+2−(q−2)=q+4`. The two conditions occur simultaneously: `V₁=I₂`
and `V₂=[[0,1],[1,0]]`, with the rest zero, already achieve them.
Thus the core rank on the common kernel is generically `2q+8>10`,
while its source rank is at most ten because all five `cᵢ` vanish.

If their quotient images span at most one dimension, their full span
has dimension at most two. The common kernel projects to a core
linear space `H` of codimension at most two. Source rank at most ten
would force `rank(B)≤5` on `H`, and hence `det B=0` there. The sextic
component is impossible by the lemma, so `H` lies in `det A=0`.
Dimension equality forces `H` to be an `A` rank-one two-plane together
with all `V` variables. Its generic reduced rank `q+4>5` is a final
contradiction.

This excludes `D=0` for every coefficient field and
every `q≥3`. Larger gains project to one gain unit.

Ten independent coefficient-matrix controls checked the rank `q+4`
for the displayed rank-one-`A` witness, for `q=3,4` over
`F₂,F₃,F₄,F₅,F₁₁`. These finite checks support the normalization
calculation; the kernel equations prove its field-independent rank.

## All-characteristic rational-function commutator obstruction

The portable `reference/shared_dot_orders.py` implements an actual
`q+2`-term polynomial degeneration with leading coefficient `t³ C_q`:

```
t Σ_i (e₀+t e_i)⊗³ − (e₀+t²Σ_i e_i)⊗³ +(1−qt)e₀⊗³.
```

Its full formal coefficient identity is checked by
`require_shared_dot_identity`; no exact `q+2` rank assertion follows
merely by setting a nonzero interpolation parameter.

Nevertheless this explicit family transfers the exact commutator bound
over a rational-function field. Suppose any finite `D` supplied a
certificate. Its exact catalytic iteration gives actual restrictions

```
T_n = n(m Unit ⊕ M₂⊗C_q) ≤ D ⊕ 5n C_q
```

for every positive integer `n`. Fix any finite coordinate rank scheme
for `D`, of size `C_D`, independent of `n`. Multiply this scheme by
`t³`, use the displayed degeneration for each of the `5n` source
copies, and apply the constant iterated restriction maps. Divide one
factor of the resulting tensor by `t³`. The resulting rational tensor
`Q_n(t)` has coefficients regular at zero, satisfies `Q_n(0)=T_n`, and
has exact rank over `K(t)` at most

```
C_D+5n(q+2).
```

This is an exact rank upper bound for the parameter family over `K(t)`,
not an unsupported substitution of border rank for exact rank.

The tensor `T_n` has square second and third dimensions

```
N_n=n(m+4(q+1)).
```

Choose the first slice using `M₂` input `I₂` and `C_q` input
`(1,1,0,…,0)` on each core block, and coefficient one on all gain
units. The shared-dot slice has determinant minus one, so this
identity slice is invertible over every field. Use first slices
`E₁₂` and `E₂₁` with that same shared-dot input, and gain coordinates
zero, for the other two matrices.

After normalizing by the identity slice, their commutator on each
core block is

```
[L_E12,L_E21] ⊗ I_(q+1) = L_diag(1,−1) ⊗ I_(q+1),
```

of rank `4(q+1)`. In characteristic two `diag(1,−1)` is the identity,
so this rank is unchanged. The gain blocks contribute zero. Thus
the limiting normalized commutator has rank `4n(q+1)`.

The corresponding slice of `Q_n(t)` is invertible over `K(t)`, because
its determinant is regular with nonzero constant term. Its inverse
and normalized commutator are likewise regular at zero. A maximal
nonzero minor of the limiting commutator has nonzero constant term,
so remains nonzero over `K(t)`. The exact characteristic-free Strassen
commutator inequality therefore gives

```
rank_(K(t))(Q_n(t)) ≥ N_n+½·4n(q+1) = n(m+6(q+1)).
```

Combining the upper and lower bounds yields

```
C_D ≥ n(m+q−4)       for every n≥1.
```

The left side is fixed and finite. This is impossible when `m+q>4`.
All steps are formal coefficient constructions or matrix-minor
identities over the original field and its rational-function field;
no assumption on its characteristic or cardinality is needed.

For `q=3,m=1`, source rate `5(q+2)=25` equals commutator lower rate
`m+6(q+1)=25`. This argument leaves that case open. For `q=4,m≥1`,
the lower rate is strictly greater than the source rate and excludes
every finite catalyst, including in characteristic two.

## Exact Koszul-minor transfer and shared-first repetition

For the remaining `q=3,m=1` case, let `P` be any **supplied actual
first-leg projection** from `Unit⊕M₂⊗C₃` to a first space of dimension
five. Keep the second and third spaces of dimension seventeen. Let
`κ_K` be the exact matrix rank over `K` of its degree-two Koszul matrix

```
Λ²(K⁵)⊗(K¹⁷)* → Λ³(K⁵)⊗K¹⁷.
```

If `κ_K>150`, then this projector excludes **every finite catalyst `D`
over that field**. The same conclusion holds for any extension over
which the given original coefficient restriction would remain exact.

Indeed, after catalytic iteration use `P` on each of the `n` target
blocks, merging all projected first coordinates into one common
five-dimensional space. The first map is the concatenation
`(P,…,P)` on the direct-sum first input, while the second and third
copy labels remain distinct. This is an actual tensor restriction.
Its Koszul matrix is block diagonal, up to coordinate permutations,
with `n` copies of the single-block matrix. Its rank is exactly `nκ_K`;
this is ordinary matrix-rank additivity.

Apply the same fixed projection to the regular rational family `Q_n(t)`
constructed above. The Koszul matrix entries are linear in tensor
coefficients, so are regular at zero. Choose a nonzero rank-`κ_K`
minor of the single-block limiting matrix, and its block-diagonal
repeat. That limiting minor has nonzero determinant, so the rational
Koszul matrix has rank at least `nκ_K` over `K(t)`.

For one rank-one tensor, its Koszul matrix is the tensor product of
`x↦a∧x` with a rank-at-most-one matrix in the other legs. For nonzero
`a∈K⁵`, the wedge map `Λ²→Λ³` has rank `binomial(4,2)=6`, in every
characteristic. Matrix-rank subadditivity therefore bounds the Koszul
rank of an exact `r`-term tensor by `6r`. The rational family's exact
rank upper bound gives

```
nκ_K ≤ 6(C_D+25n),
C_D ≥ n(κ_K−150)/6.
```

This contradicts fixed finite `C_D` when `κ_K>150`.

The criterion needs exact coefficient matrices with the proper exterior
signs. In characteristic two signs collapse; in odd characteristic they
must be retained. An exact rank computation over one finite field is
a theorem for that field and its extensions, not automatically for all
fields. A common all-field certificate can instead provide an explicit
integer minor of determinant one or minus one, or provide a nonzero
integer minor and separately handle all prime divisors of its determinant.
Merely testing several prime fields is insufficient to assert that.

The following actual integer certificate now fulfills the criterion
over every field.

## Universal unit-minor certificate for `C₃`

The five-coordinate first map, in target order
`Unit⊕(M₂⊗C₃)` with row-major matrix coordinates and lexicographic
tensor-product coordinates, has rows

```
0 1 0 0 1 0 1 0 1 0 0 0 0 0 1 1 1
0 0 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1
1 1 0 1 0 0 1 0 0 0 0 1 0 0 1 0 1
0 1 0 0 1 1 1 1 0 1 0 1 0 1 1 1 1
0 1 1 0 1 0 0 0 0 0 1 0 1 1 0 1 1
```

The exact signed degree-two Koszul matrix has shape `170×170`. The
compact certificate in
[shared_dot_three_koszul.json](../shared_dot_three_koszul.json)
supplies the first map and 153 distinct row and column indices.
[`shared_dot_koszul.py`](../shared_dot_koszul.py) rebuilds the integer matrix.
The corresponding minor has determinant one in the listed order.
It therefore has rank at least 153 over every field.

The finite-auxiliary-row agent found the sparse unit-pivot certificate
and independently checked the original minor by Bareiss elimination.
This audit independently reconstructed the integer target from all
eight matrix-multiplication support triples and nine shared-dot support
triples, plus the gain unit: 73 supported triples. It reconstructed the
projected coefficients and the signed exterior matrix and checked every
entry against the supplied artifact. It then extracted the specified
original minor with both index sets sorted and ran a separate
fraction-free Bareiss computation, including row-swap signs and exact
divisibility checks. Its determinant was **minus one**; sorting changes
the certificate's sign and preserves its unit value. This verifies an
integer determinant, not a collection of finite-field rank samples.

The original binary numerical rank 166 is stronger for that field, but
the unit minor alone supplies the universal rank 153 needed here.
Under the shared-first repetition, its block-diagonal `n`-fold minor
has determinant `±1`, so the rational family has Koszul rank at least
`153n` in every characteristic. The exact rank-one factor is six and
the family rank upper bound is `C_D+25n`; hence

```
153n ≤ 6(C_D+25n),       3n ≤ 6C_D.
```

The fixed finite number `C_D` cannot satisfy this for all positive
integers `n`. Thus `C₃` admits no finite catalyst for gain one over any
field. Larger positive gains project to one gain unit. The universal
minor, formal degeneration, and actual shared-first restriction together
give the arbitrary-catalyst theorem; the numerical rank alone is not
used as an all-field claim.

## Executable constructions and their limits

`shared_dot_orders.py` supplies literal shared tensors, the quotient
`C_q <= S_q` from the three-oriented length-q dot direct sum, and actual
three-unit maps for q>=3. Its exact `shared_dot_scheme` uses 2q+1 terms
in characteristic two; in odd characteristic it averages the plus/minus
cubes and subtracts q times the scalar cube. If q casts to zero that last
term vanishes, giving 2q terms. These are certified upper lengths, not
optimal-rank assertions. All coefficients use the represented field.

The q+2-term polynomial family is checked at every formal degree, including
degree-four and degree-six noise. `recover_shared_dot_power` connects it to
the existing power/interpolation engine: exponent e pays 3e+1 distinct
nonzero nodes and (3e+1)(q+2)^e terms. It rejects insufficient nodes before
construction and inserts no field extension. The C3 square over F11 was
verified with 175 terms at all 16³ coefficient positions.

`shared_dot_matrix.py` adds actual star, square-to-M2 and two-copy-to-M2
orders. For the two-copy construction each source first input is
(a0,a1,0,...), the second input is (B11,B10,B00,B01-B10,0,...), and
source outputs Z0+Z2 and Z1+Z3 give that matrix row. The two copies have
independent first rows and share B. Square maps compose coordinate star
selectors with the independently checked signed star-square construction.
These ordinary coefficient maps do not provide a positive gain.

`shared_dot_catalyst_screen.py` recognizes literal C_q and records
the all-field exclusion for every q>=1 and any finite D. It checks the full polynomial identity and actual
invertible slices and commutators; finite controls are not proofs of the
unbounded iteration or generic geometry. Unsupported coefficients and
resource caps remain distinct from exclusion. None of these new
mathematical arguments has been formalized in Lean.

The optional `--shared-dot-connections` pass adds only ordinary exact
quotient, subrank, matrix, rank-upper and selected context rows to the
finite-state experiment. Analytic screens do not become state rows.
A feasible finite state says nothing about an arbitrary unsearched
catalyst; an actual negative-unit dual would still require exact replay.
