# Independent review: witnessed LP assembly and catalyst compiler

Review scope: `reference/witnessed_constraints.py`,
`reference/catalyst_compiler.py`, and their two test modules as supplied during
this research pass. This reviewer did not edit either implementation. The
mathematical derivations here remain separate from Lean verification.

## Correctness findings

No coefficient-correctness counterexample was found in the reviewed paths.
Both focused test modules passed: **18 tests**, approximately 4.2 seconds.

The following central details are implemented correctly:

- `WitnessedOrder` represents `negative <= positive`; its local maps go from
  the positive source direct sum to the negative target direct sum.
- Detector rows are hypothetical scalar constraints. They are accepted only
  with an exact registered `M_d tensor test` identity, and are not treated as
  ordinary restriction witnesses.
- Integer dual balance is ordinary integer balance. For example, two scalar
  gain copies in characteristic two remain two independent tensor blocks.
- `assemble_catalyst` uses the positive ordinary blocks as `D`, the detector
  test blocks as `S`, restricts `D` to the negative ordinary blocks, and then
  matches the literal tensor-ID occurrence multiset. No cancellation of `D`
  or `S` is used.
- The final product packing is the inverse distributivity bijection. For
  one axis, block `j` has width `b_j`; output order is matrix index first,
  then auxiliary block and its coordinate. This agrees with the axiswise
  tensor-product encoding even when the three axes have different sizes.
- `_lift_core` preserves full cross-slot entries of the old catalyst maps;
  it does not assume that `D` is merely passed through. The matrix pair index
  splits combined row/column digits into the spectator and old matrix block.
- `compile_absorption` reuses the same catalyst between successive seed
  terms, retaining its potentially mixed coordinates. It does not cancel the
  catalyst or allocate an independent catalyst for every term.
- Final matrix coefficients are verified independently, and final descent
  applies the one fixed supplied field's degree-squared factor. There is no
  automatic repeated descent inside the absorption recurrence.
- `find_nonnegative_dual` checks independent supports up to the coordinate
  dimension. Its stated finite-family completeness is valid: a minimal
  nonnegative cone representation has linearly independent support. An
  exhausted finite family is not reported as global impossibility.

### Additional independent unequal-axis check

An in-memory example used `F_4`, two detector tests, and a tensor `x` of shape
`(1,2,3)` with coefficients `(2,0,3,0,1,0)`. The ordinary rows were the naive
matrix upper bound, its tensor translation by `x`, and extraction of a scalar
unit from a nonzero coefficient of `x`. With integer weights all one and gain
two, assembly gave

```
D.shape      = (17,26,35),
S.shape      = (2,3,4),
source.shape = (35,53,71),
target.shape = (27,40,53).
```

The full coefficient check passed. This covers unequal axis widths and
extension coefficients beyond the symmetric ragged-block test. It is a
conservative toy, not a useful `d=2,k=5` catalyst.

## Resource-preflight issues identified

These findings were sent to main and corrected during this review. The
reviewer inspected the resulting early exponent guard, sparse integer balance,
used-coordinate reduction, and dense constraint-entry cap. The two focused
resource regression tests passed after those corrections. The original issues
below concerned resource behavior, not acceptance of invalid coefficients.

### 1. Huge powered-catalyst exponent before any cap

The reviewed `compile_powered_catalyst` computes `k**exponent` and loops over
all `exponent` terms of its cost estimate before applying the declared dense
limits. A value such as `10**9` can therefore require substantial time and
memory merely to discover that the requested construction cannot fit.

There is a cheap guaranteed lower bound. Since `d>=2` and `S` is nonzero, all
three auxiliary axes are positive. The final target has at least

```
d^(6*exponent) >= 2^(6*exponent)
```

tensor coefficients. If the cap is the integer `cap`, reject before powers or
loops when

```
6*exponent >= cap.bit_length().
```

Indeed `2^cap.bit_length() > cap`. This safe early rejection does not replace
the existing detailed preflight; it handles obviously impossible exponents
without constructing huge integers. A similar bounded-power rule can protect
other APIs that allow arbitrarily large exponents.

### 2. Uncapped dense state-vector matrix

The reviewed `WitnessedStateSystem.vectors` allocates a dense integer matrix
of `number_of_rows * number_of_registered_IDs` entries. The tensor, map, and
expanded-block limits do not bound this scalar LP allocation.
`find_nonnegative_dual` accesses the matrix before its subset cap check, even
when the caller requests zero subsets. `verify_balance` and rational-dual
clearing also access it.

Many unused registered IDs combined with repeated tiny valid detector rows
can make this matrix enormous while every tensor is small. Add a scalar
matrix-entry limit before dense-vector access, or perform certificate balance
with sparse counters and cap dense construction only in the finite LP solver.
Keep the distinction between incomplete resource-limited search and finite
exhaustion.

## Proposed sector-to-ordinary-row export API

The current proof-derived sector degeneration can supply ordinary witnessed
constraints without first generating a rank decomposition. Suggested API:

```
export_sector_order(construction, source_key, retained_key,
                    *, exponent=1, nodes=None, limits=...)
```

Return a `WitnessedOrder` plus the exact inventory tensors and node metadata.
Here `source_key` identifies **one powered source atom**, not a whole repeated
direct sum. Define

```
S_e = construction.source.tensor_power(e),
R_e = construction.retained.tensor_power(e).
```

The generated row has

```
positive = (source_key,)*len(nodes),
negative = (retained_key,),
maps: len(nodes) copies of S_e -> R_e.
```

The sector identity is `P(t)=t*R+t^2*E`; consequently
`P(t)^e/t^e` has degree at most `e`, and its constant coefficient is `R^e`.
Require at least `e+1` distinct nonzero field nodes before allocating powers.
Use the existing exact constant-coefficient interpolation weights.

For each powered axis coordinate `I`, decode its `e` seed-axis digits and
define the evaluated diagonal at node `t` by

```
diag_e_axis(I,t) = product over digits i of diag_axis(i).evaluate(t).
```

Each output row is a horizontal concatenation of evaluated diagonal rows over
the nodes. In the C map only, multiply the node block by `weight(t)/t^e`.
The resulting direct-sum restriction equals the interpolation sum because
all three local maps meet inside the same source-copy label. All powers use
the same parameter. The source shape is `len(nodes)*S_e.shape` axiswise;
the output shape is `R_e.shape`.

Validate the construction, preflight node availability and dense dimensions,
materialize the three matrices, and check the entire tensor coefficient
identity. Do not use finite-field evaluations as a substitute for the formal
sector identity; this export is derived from the already checked polynomial
construction.

The zero power is a scalar unit and needs one nonzero node. Zero or empty
source axes should either receive explicit empty-axis handling or be excluded
with a documented finite implementation boundary. Arbitrarily large exponent
arguments need the same bounded-power preflight discussed above.

This adds a genuine finite restriction row to the LP inventory. It retains
the `e+1` source-copy overhead; the zero-overhead character inequality is
asymptotic and is not this row. A focused `e=1`, `a=2,h=1` check over `F_4`
already succeeded during the earlier research pass for source shape `(4,8,10)`
and retained target shape `(2,4,5)`.

## Second pass: retained gains, automatic planning, and finite-type rows

The independent second pass read `compile_gain_powered_catalyst`,
`compile_gain_absorption`, `reference/gain_planning.py`, and
`reference/finite_type_constraints.py`. The compiler, gain-planning, and
finite-type modules passed 34 focused tests together. A new metadata-cap
regression passed separately after main's correction below. No implementation
files were edited by this reviewer, and no Lean build was run.

### Positive gains and matrix coefficients

The powered maps implement `g_next=k*g+m*u`, where `u=d^b`, by retaining k
copies of the previous gain and restricting each new gain's `M_u` factor to
u diagonal scalar units. All three diagonal selectors use coordinate
`i*u+i`. This preserves ordinary integer copy counts in characteristic two;
the gain is never converted into one field coefficient. The subsequent
product packing uses the same row/column digit split as the gain-free path.

The absorption loop keeps D as one mixed block and composes every full
catalyst map, including entries crossing D and S. Its current-block and
spectator offsets account for all previously retained gain units. It then
puts matrix coordinates first and gain units last, matching the exact
scalar-output pivot elimination interface. The cost is

```
C_b + (k^b*R_s-g_b)*r_seed.
```

The 70-term empty-D and 87-term cross-mixed-D size-four schemes are now
actual generated arrays with exact coefficient checks, rather than only
predicted costs. An additional independent F_4 test used D shape `(1,2,1)`
and S shape `(2,1,3)` with non-prime coefficients: the powered map verified
from `(175,107,256)` to `(56,53,72)`, and its size-four scheme verified with
185 terms, equal to `34+(81*2-11)`. This independent check used
`CompilerLimits(max_tensor_coefficients=5_000_000)` explicitly; its source
exceeds the default two-million-coefficient cap. No coefficient counterexample
was found.

### Automatic target planning

The planner uses strict integer comparisons `Q_b^B < d^(b*A)`, not
floating-point logarithms. Accepting this fixed-block condition when
`k^B>=d^A` is valid: retained positive gains change the actual finite
coefficient Q_b. The scalar control case `d=2,k=9,m=1` has Q_1=8 and passes
the `31/10` target despite failure of the k-based sufficient condition.
This does not claim 9/4 for that conservative witness.

The integer recurrence includes the direct catalyst cost once per step.
The extension degree squared is charged once in the final target comparison
and once in actual final descent. Dense final-array preflight precedes
powered compilation; a numeric plan is not presented as generated arrays.
The reviewed exponent and iteration defaults are bounded search limits,
not a complete search guarantee.

### Exact-type and finite-separation row directions

For exact type, the source is `shared(B)^n`; the target shares all products
whose literal word positions have the requested counts. The first map is
identity. The other maps decode each word and map its coordinate to base
`number_of_branches*axis_size` digits `label*axis_size+digit`. This gives a
restriction from source to target without replacing differently ordered
tensor products by equal atom IDs. The zero power is handled as one empty
word and a scalar unit.

For finite separation, every node and Fourier label is a separate source
copy. Horizontal concatenation of the actual weighted maps realizes the
sum over those copies. Lagrange weights occur in the third leg only; negative
parameter powers are evaluated only at nonzero nodes. The final first-leg
permutation converts coordinate-major Fourier order to the branch-major
direct sum `sum_h B_h tensor dot_M`. Full tensor verification checks both
exports before they can become ordinary witnessed rows. Neither export
silently removes interpolation-copy overhead or substitutes entropy/state
inequalities for finite restrictions.

### Found and corrected metadata-cap gap

Originally, zero branch axes could bypass every dense tensor/map cap while
word metadata remained large. Two branches of shape `(1,0,0)` and counts
`(4095,1)` have n=4096 and M=4096; they satisfy a 4096-block cap, but retain
over 16 million word positions and repeatedly form long empty-axis products.
The reviewer supplied this reachable example. Main added the preflight
`n*M <= limits.max_constraint_entries` before powers or word enumeration,
and `test_empty_axes_do_not_bypass_the_word_position_cap` passes.

An initially suspected zero-size matrix-seed path was not reachable through
the current `BilinearScheme` constructor, which requires positive matrix
dimensions. It is not a confirmed defect and requires no implementation
change. No further correctness or preflight defect was identified in these
reviewed paths.
