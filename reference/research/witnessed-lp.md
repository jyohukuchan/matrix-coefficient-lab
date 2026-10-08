# Research B: a witnessed finite LP can assemble a catalyst directly

This report derives an executable interface from the current Lean sources.
It is a proposed coefficient construction, checked at the level of finite
tensor identities, not a Lean-verified new theorem. It modifies no Lean source.

The key observation is that a dual certificate can retain **positive and
negative direct-sum blocks**, rather than discard them into an additive group
completion. Those blocks give an explicit catalyst and three local matrices.
No search for the localization witness used by the abstract proof is needed.

## Source facts

Paths below are inside
`lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/`.

- `Spectrum/StateObstruction.lean:34–52` defines the six state-constraint kinds
  and their integer vectors.
- `StateObstruction.lean:109–151` separates ordinary restriction differences
  from the distinguished detector difference.
- `StateObstruction.lean:162–180` interprets a finite natural combination
  equal to `-m*e_unit` as a positive-gain catalyst.
- `Convex/Farkas.lean:65–74` supplies the finite normalized alternative;
  `Convex/RationalCone.lean:162–205` rationalizes and clears denominators.
- `Spectrum/Catalyst.lean:28–56` selects an additive-localization witness at
  line 49. The construction below replaces that choice by explicit retained
  finite blocks.
- `Tensor/Semiring.lean:46–48` defines an order witness as three local matrices
  whose restriction gives exactly the target tensor.

## Concrete registry and relation records

Use one explicitly represented finite field throughout one certificate. A
registry maps opaque IDs to exact `FiniteTensor` arrays over that field. An ID
always denotes one fixed tensor, with fixed axes and coordinate order.

Do not identify IDs merely because dimensions agree, evaluations agree, or
an approximate invariant agrees. If a representative changes, either register
the new exact tensor under a new ID, or supply checked equivalence maps and
translate the relation to the fixed representative.

For a block list `L`, write `DS(L)` for its full three-axis direct sum;
`DS([])` is the empty tensor with shape `(0,0,0)`. `unit` is the exact scalar
tensor with shape `(1,1,1)` and coefficient one. Set `T=M_d`.

Each non-detector relation record contains:

```
kind, integer vector v,
positive opaque-ID list P,
negative opaque-ID list N,
three matrices W with DS(P).restrict(W) == DS(N).
```

Thus `v` must equal `count(P)-count(N)`, where counts are ordinary integers.
The required data for each relation are:

| Kind | Vector | P | N | Concrete source-to-target maps |
| --- | --- | --- | --- | --- |
| lower(x) | `e_x` | `[x]` | `[]` | delete all coordinates |
| upper(x,r) | `r*e_unit-e_x` | `r` units | `[x]` | transpose an exact r-term scheme's A/B/C families |
| order(x,y) | `e_y-e_x` | `[y]` | `[x]` | supplied checked restriction `x <= y` |
| addPos(x,y,w) | `e_w-e_x-e_y` | `[w]` | `[x,y]` | checked map `w -> x ⊕ y` |
| addNeg(x,y,w) | `e_x+e_y-e_w` | `[x,y]` | `[w]` | checked map `x ⊕ y -> w` |

For the additive relations, `w` must be the tensor class of `x ⊕ y`, with
actual maps establishing the direction used. If `w` is literally the exact
direct-sum array, associativity/reindexing maps are sufficient. A general
mutual-restriction representative requires its explicit equivalence witnesses.
Each direction can be a separate record; no unverified equality is assumed.

A detector record instead contains

```
x_id, tx_id, v = e_tx-k*e_x,
checked identity registry[tx_id] == T.tensor_product(registry[x_id]).
```

It has positive list `[tx]` and negative list consisting of `k` copies of `x`,
but it deliberately supplies **no restriction from the negative to the
positive side**. The detector is the hypothetical state inequality whose
failure is being certified. Assembling it as an ordinary witnessed order
relation would assume the desired conclusion.

`upper` may use any supplied exact rank upper bound, including naive schemes.
No minimum-rank oracle is required. One need not impose a single computable
class-wide rank function on this executable finite registry.

## Exact rational LP and integer certificate

Select finitely many records, including every ID mentioned by their vectors.
For vector columns `v_i`, the normalized primal is

```
f(unit)=1,
dot(v_i,f)>=0 for every selected i.
```

An exact dual target is

```
q_i >= 0 rational,
sum_i q_i*v_i = -e_unit.
```

Clear denominators with `m=lcm_i denominator(q_i)` and `n_i=m*q_i`. Then

```
m > 0, n_i >= 0 integers,
sum_i n_i*v_i = -m*e_unit.                  (1)
```

The certificate verifier requires ordinary integer equality at **every ID**;
it does not reduce counts modulo the field characteristic. An arbitrary
solver's approximate infeasibility status is insufficient. The verifier must
check the rational solution or the cleared integer certificate independently.
Dividing all `n_i` and `m` by their common gcd is a sound optional reduction.

## Explicit assembly without group-completion cancellation

Partition records into ordinary records `O` and detectors `H`. Repeat each
record's lists/maps `n_i` times. Define the following literal tensor block lists:

```
D_list = concatenation over i in O of n_i copies of P_i,
N_list = concatenation over i in O of n_i copies of N_i,
S_list = concatenation over i in H of n_i copies of [x_i].

D = DS(D_list), N = DS(N_list), S = DS(S_list).
```

The block-diagonal sum of the ordinary maps gives `N <= D` explicitly.
Let `P_H` denote the opaque list of `T*x_i` blocks repeated `n_i` times.
The dual identity (1) says exactly

```
count(D_list) + count(P_H) + m*e_unit
  = count(N_list) + k*count(S_list).        (2)
```

The two sides are the same multiset of exact registered tensor blocks.
Match equal-ID occurrences deterministically; copy their within-block axis
coordinates using three independent permutation matrices. This gives a
literal reindexing between their full direct sums, without cancelling any
tensor factor or direct summand.

Distributivity reindexings identify

```
DS(P_H) ==_reindex T.tensor_product(S),
DS(k repetitions of S_list) ==_reindex k copies of S.
```

Compose, in source-to-target order:

1. Expand/group the source `D ⊕ k*S` into ordinary source blocks followed by
   `k` copies of `S_list`.
2. Apply the block-diagonal ordinary restrictions, and identity on `k*S`,
   yielding `N ⊕ k*S`.
3. Reorder opaque blocks using (2), yielding
   `D ⊕ m units ⊕ DS(P_H)`.
4. Regroup `DS(P_H)` by inverse product-distribution reindexing, yielding
   `D ⊕ m units ⊕ (T tensor S)`.

The composite consists of three finite local matrices and proves

```
D ⊕ m*unit ⊕ (M_d tensor S) <= D ⊕ k*S.
```

This is precisely the input of `reference.catalysts.CatalyticCertificate`.
The actual assembled source-to-target maps must undergo its full coefficient
check. Retaining the positive additive blocks is what supplies the catalyst:
the abstract completion proof omits those blocks because their group images
cancel, then must recover an unspecified localization witness.

### Coordinate recipe for product distribution

For one tensor axis of `T` of size `a`, and block `x_j` of axis size `b_j`
at offset `o_j` inside `S`, the product-coordinate ordering is

```
T tensor S coordinate:      i*sum_j b_j + o_j + z,
sum_j (T tensor x_j):       sum_(h<j) (a*b_h) + i*b_j + z.
```

Use this bijection independently on all three axes. Reverse it when packing
the final target. Occurrence matching for (2) uses each block's independent
axis offsets; it must not assume all blocks have equal dimensions.
Zero-length axes require empty index lists instead of division by zero.

## Small complete test: d=2, k=9

Take registry IDs `u=unit` and `t=M_2`; `M_2 tensor unit` has exactly the
same array as `M_2`. Supply its naive 8-term decomposition. Use just

```
upper(t,8):       8*e_u-e_t,
detector(u):     e_t-9*e_u.
```

Weights `(1,1)` give `-e_u`. The assembler obtains

```
m=1, S=unit, D=8 units,
source = 17 units,
target = 8 units ⊕ 1 unit ⊕ M_2.
```

The ordinary stage maps the source's first eight units to `M_2` using the
naive scheme. The source's remaining nine units are reallocated to the target's
eight catalyst units and one scalar-gain unit. Thus all maps are explicit,
with target shape `(13,13,13)` and source shape `(17,17,17)`.

This is a genuine matrix-tensor catalyst certificate, but **not a useful 9/4
catalyst**: `log_2(9)>3`. It tests complete LP-to-matrix assembly without
claiming progress on the efficient witness search. A variant using Strassen's
7-term scheme and `k=8` gives another positive certificate.

For `d=2,k=5`, the same naive or Strassen upper-bound columns provide no
positive dual certificate. Additional explicitly witnessed tensors/relations
are needed. Inventing a smaller upper bound or treating numerical character
samples as restriction witnesses would be unsound.

## Suggested implementation tasks

1. **Registry and witnessed records:** strict field/ID/shape validation;
   construction of integer vectors from the relation kind; exact checking of
   each ordinary map and each detector's product identity.
2. **Certificate verifier:** accept integer `(n_i,m)` first, reject negative
   counts, zero gain, unknown IDs, and any balance mismatch. Rational solver
   integration can follow after this verifier works.
3. **Assembler:** implement repeated block lists, block-diagonal maps, opaque
   occurrence matching, and product distribution; return a checked
   `CatalyticCertificate`. Preserve a construction record for inspection.
4. **Complete toy and negative controls:** `d=2,k=9` over both odd and even
   characteristics; tampered tensor coefficients, reversed order maps,
   missing additive direction, swapped block offsets, and certificates whose
   balance holds only modulo the characteristic must be rejected.
5. **Useful witness acquisition:** populate the registry with actual finite
   constructions and certified restrictions, then solve exact LPs. This is
   the separate mathematical search; the assembler cannot create relations
   absent from its input.

Dense arrays can become large even for small dual weights, so explicit resource
caps should return incomplete construction status. An integer dual certificate
alone must never be labelled a checked tensor catalyst. Valid positive-gain
data cannot have zero `S`; the assembled certificate must check nonzeroness,
and a claimed certificate with no detector contribution must fail a concrete
validity check rather than be used as a useful witness.

One benefit of this interface is that every selected relation is checked at
the coefficient level. Completeness of the registry for arbitrary tensor
classes is unnecessary for soundness. It affects only whether the LP search
can find the useful catalyst guaranteed by the existence argument.

## Sound next constraint generators for d=2, k=5

These generators are proposals for populating the registry; this report has
not found a useful dual certificate for that target.

1. **Translate existing restrictions by products.** If checked maps establish
   `x <= y`, their Kronecker products with identity maps establish
   `x tensor z <= y tensor z` for any registered `z`. Register both actual
   product arrays and that ordinary order row. Close a bounded inventory under
   several products rather than silently imposing scalar multiplicativity on
   the additive-state variables. Detector rows for these translated tensors
   are separate rows with exact `M_2 tensor z` identities.
2. **Use finite interpolation as a witnessed order.** For diagonal sector
   maps `A_t,B_t,C_t` and recovery weights `w_t/t^leading`, concatenate their
   matrices horizontally over the evaluation nodes and scale the C block by
   the recovery weight. The resulting maps restrict the direct sum of source
   copies to the retained target. The same construction for a power gives
   `retained^n <= (n+1) copies of source^n` for the sector's normalized degree
   bound one, assuming enough nonzero nodes in the supplied field. This order
   row carries the interpolation overhead; it must not use the zero-overhead
   character inequality as if it were an exact finite restriction.
3. **Use finite Fourier separation as a witnessed order.** Concatenate the
   actual Fourier-source and weighted local maps over interpolation nodes.
   The output is the separated direct sum, whose permutation to branch/dot
   blocks must also be retained. Charge every Fourier copy and interpolation
   node. Register branch extraction as additional witnessed restrictions.
4. **Import determinant filtration maps when implemented.** The explicit
   adapted-basis tensor and associated graded degeneration occur in
   `Determinant/Filtration.lean:46–139`; character-only use occurs in
   `Determinant/Character.lean:152–180`. They can supply genuinely different
   finite relation families after their basis and polynomial maps have been
   implemented and the finite interpolation overhead retained. The
   concavity inequality itself is not an ordinary linear order row.

An efficient first search should keep tensor factors and block provenance
symbolically while solving the **integer/rational scalar LP**, then materialize
and verify tensor maps only for selected rows or a proposed dual support.
Map verification remains authoritative; symbolic provenance is an optimization,
not a substitute for the exact tensor certificate.

As a focused executable check during this pass, horizontal concatenation of
the generated sector diagonals over nodes `(1,2)` in `F_4`, with C scaled by
`weight/node`, was evaluated for `a=2,h=1`. The exact full coefficient check
confirmed a restriction from two source copies, shape `(4,8,10)`, to the
retained tensor, shape `(2,4,5)`. This validates the proposed finite-order
generator on that example; it is not a d=2,k=5 dual certificate.

The current small matrix pipeline gives cubic matrix schemes. Strassen supplies
`M_2 <= 7 units`; neither this row nor its tensor powers are a direct
`M_2 <= 5 units` witness. New finite relations can provide catalysts without
providing a direct rank-5 decomposition, but finding them
requires a richer inventory. No sampled real-valued character, numerical
rank estimate, or rank bound inferred from a dimension can replace its maps.
