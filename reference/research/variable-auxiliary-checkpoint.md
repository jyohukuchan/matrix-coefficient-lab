# Variable auxiliaries and actual matrix connections

Checkpoint: 2026-10-09. **No useful `d=2,k=5` positive-gain certificate or
near-`9/4` matrix decomposition has been found.** The coefficient compiler
still requires a discovered finite witness. The work below supplies actual
new restriction maps and stronger candidate filters; it does not finish
stage 4. No Lean sources were changed or rebuilt.

## Actual coefficient maps crossing matrix power levels

[`sector_matrix.py`](../sector_matrix.py) exports four independently
coefficient-checked constructions:

- `export_retained_matrix_order`: `R(a,h)^2 -> M(a,a,2)` for positive `a,h`.
  It selects the left-middle and right-middle tensor words, preserving their
  shared first matrix coordinates. At `(a,h)=(2,1)` the source is `(4,16,25)`;
  its three selectors are `(1,0,3,2)`, `(1,13,2,14)`, `(2,17,7,22)`.
- `export_star_matrix_order`: `S_star^2 -> M2`. Signed maps use the identity
  `adj(A)=tr(A)I-A` to correct a transposed matrix-vector block. The full
  coefficient identity holds over the integers, hence over every field.
- `export_retained_two_copy_matrix_order`: `2 R(a,h) -> M2` for `a>=2,h>=1`.
  Each copy supplies one matrix row. Its middle output is one dot product;
  adding selected left/right outputs gives the other. The second map shares
  the same matrix input between the two independent source copies.
- `export_matrix_column_sum_order`: `c` independent copies of `M(a,b,1)`
  restrict to `M(a,b,c)`. First-leg maps merge the copies' matrix coordinates;
  the other two maps preserve the separate output columns. Copy counts stay
  integers even when they vanish as field scalars.

Every API checks tensor/map allocation bounds before materializing arrays,
then verifies every coefficient including zeros. Tests cover extension-field
encodings, characteristic two, rectangular products, corruption, and caps.
Restricting the supplied six-term retained scheme squared gives the classical
eight-term `M2` scheme. The new maps are useful ordinary state rows, not a
new low-rank matrix algorithm.

The existing shared left-middle seed is the star after actual coordinate
changes. Its padded square therefore also restricts to `M2`; the experiment
runner composes and verifies those maps. The distinct shared LM/ML exact-type
target has two rectangular word components. Column stacking connects two
copies of each rectangle, or two copies of the shared target, to `M2`.
At detector value `f(M2)=5`, these rows give lower values `5/2` without
asserting a three-unit subrank restriction.

This matters because a proposed assignment that simply scaled first
flattening ranks to make `f(M2)=5` assigned value four to the retained and
star squares. The actual square-to-matrix maps refute that assignment.
An exact finite feasible state may still assign different values to a tensor
and its powers: it is an additive monotone state for selected rows, and is
not automatically multiplicative.

## Expanded exact finite experiments

[`variable-proof-experiments.py`](variable-proof-experiments.py) preserves
literal finite tensors and local maps, reconstructs them on replay, and
checks normalized rational states or integer negative-unit balances exactly.
No supplied rank-seven matrix scheme is included. Its selected fields are
`F4,F11,F16,F31`; no extension is inserted implicitly.

The completed historical constructor families were:

| Family | F4 atoms/rows | F11 atoms/rows | F16 atoms/rows | F31 atoms/rows |
| --- | --- | --- | --- | --- |
| Strong lower maps and Fourier splitting | 115/401 | 157/556 | 139/491 | 157/556 |
| Added mixed retained subrank maps | 115/409 | 157/564 | — | — |
| Adaptive coordinate subrank cuts | 115/501 | 157/683 | — | — |
| Filled existing product detectors | 115/519 | 157/701 | — | — |
| Added both square-to-matrix rows | 115/521 | 157/703 | — | — |
| Completed word columns, retained rows, M2 detector and fixed context | 115/529 | 157/711 | — | — |
| Mixed-dot ordinary maps, 15-unit square and actual context swaps | 125/555 | 167/737 | — | — |
| Shared C3 quotient, exact upper schemes and matrix connections | 132/574 | 174/756 | — | — |

Every listed completed family has an exactly verified feasible state and
replayed ordinary maps. No negative-unit dual was found. The original F16
construction-window-limited run is superseded by the completed 139/491 run.
Adaptive families at a prescribed round cap are not closure results.
Even a zero-new-coordinate-cut round leaves node-capped searches and general
linear restrictions unresolved.

Before the mixed/shared pass, the F4/F11 states have `R(2,1)=R(2,2)=5/2`, `R(3,1)=3`, the retained
square `25/4`, the shared LM seed square five, and `M2=5`. Individual LM/ML
words and their shared target have value `5/2`. Retained `M2` product detectors
are tight, including value `125/4` on the product with the retained square.
The actual M2 detector pins its square to 25. The fixed retained context lift
is the ordinary restriction `2 R(2,1)^2 -> M2 tensor R(2,1)`, with existing
atoms and checked local maps. These are selected finite-state values, not
rank values or a character of the full tensor semiring.

Both latest files passed fresh standalone replay, with all 574 and 756 ordinary
or detector rows checked. [The compact result record](variable-auxiliary-results.json)
contains saved artifact hashes and budgets; it does not replace the literal
coefficient certificates. The complete reference suite passed 405 tests in 353.641 seconds;
focused controls and the preceding 401-test run are also recorded there.
The eleven protected files and frozen Challenge model still match the exact
specification baseline. Portable F4 integration tests reject weakened states,
check actual product-coordinate swaps, and verify repeated passes add no rows.

## Mixed/shared tensor constructions and node-free extraction

[`mixed_dot_orders.py`](../mixed_dot_orders.py) gives actual three-unit,
star-plus-dot, 15-unit square, square-to-M2 and two-copy-to-M2 maps.
[`mixed_dot_lengths.py`](../mixed_dot_lengths.py) extends coefficient/profile
controls to arbitrary positive dot length; square selectors give 3+6t units,
without an optimality assertion. The connected finite states assign three to
the mixed auxiliary, fifteen to its square, three to C3, and 49 to the C3
square. Such states need not be multiplicative. No analytic rank lower bound
is inserted as an ordinary state row.

[`shared_dot_orders.py`](../shared_dot_orders.py) implements shared Cq
quotients, actual 2q+1-term upper schemes (2q in odd characteristic when
q casts to zero), and a q+2-term formal polynomial degeneration. All formal
noise coefficients are checked. Its powered node recovery verifies a 175-term
C3 square over F11; insufficient nodes are rejected. The separate
[`shared_dot_matrix.py`](../shared_dot_matrix.py) connects shared tensors to
the star and to actual matrix multiplication through squares or two copies.

[`coefficient_extraction.py`](../coefficient_extraction.py) now collects
degree triples from a supplied polynomial decomposition directly, with no
evaluations, interpolation nodes or field extension. C3 and its square yield
9- and 135-term exact schemes over F2,F3,F4 even when interpolation cannot
run. These raw lengths are controls, not optimal-rank claims. This resolves
node shortage for a supplied finite family; discovering the useful family
or positive-gain restriction from the existence proof remains open.

The construction parameters include sectors `(a,h,e)=(2,1,1),(2,1,2),
(3,1,1),(2,2,1)`, determinant filtrations `(d,e)=(1,1),(2,1),(1,2),(2,2)`,
and two exact-type pairs with counts `(1,1)`. Actual Fourier periods are ten
over F11/F31 and fifteen over F16; this family has no available period over
F4. Larger copied-source constructions exceeding their caps are omitted
explicitly.

## Mathematical screening changed the candidate frontier

The independently audited [derivations](variable-auxiliaries.md) are not
Lean theorems. They apply to actual exact restrictions

```
D + m Unit + M2 tensor S <= D + 5 S,  m >= 1.
```

Complete polynomial convolutions, including direct sums sharing their
polynomial-output leg, are excluded as fixed auxiliaries for every finite
`D`. Output elimination strengthens the commutator lower bound and additive
iteration absorbs the fixed catalyst cost. Passing spectral calibration
alone therefore does not make a prescribed convolution viable.

For the retained `R(2,1)` auxiliary, necessary catalyst flattening boundaries
are `(5,8,16)`, `(5,14,10)`, or `(5,12,12)` in the stated leg order. For the
fixed star they are `(6,6,10)`, `(6,10,6)`, or `(6,8,8)`. These are necessary
shape inequalities, not witnesses. Further rank-drop component conditions
exclude some tensors having those dimensions. The former random `(4,4,5)`
star catalyst, smaller band tensors, and the proposed outer-product plus
dot-tail retained catalyst are now excluded mathematically.

[`star_catalyst_screen.py`](../star_catalyst_screen.py) implements selected
exact consequences of those unformalized lemmas. It uses actual concise
flattening ranks, a basis of the first slice image, and coefficient-field
rank counting. Finite slice maxima are never substituted for generic ranks;
an interrupted enumeration does not exclude a candidate.

## Bounded direct experiments and their limits

The mixture of one length-two dot in each singleton orientation is now
excluded for every finite catalyst, every positive gain and every field.
The [rank audit](mixed-dot-rank-catalysts.md) proves a replicated rectangular
rank lower bound, applies first-leg substitution, and iterates a hypothetical
restriction until the fixed catalyst cost is absorbed. It also covers unequal
positive dot lengths mathematically. The earlier slice-space proofs remain
in [the mixed-dot audit](mixed-dot-catalysts.md) and [the general
note](variable-auxiliaries.md), but their unresolved boundary is superseded.
[`mixed_dot_rank.py`](../mixed_dot_rank.py) checks actual projections from
repeated cores to independent M(2,2,2t) blocks and first concision.
[`mixed_dot_catalyst_screen.py`](../mixed_dot_catalyst_screen.py) recognizes
literal S2 coefficients and reports the all-catalyst analytic obstruction;
it computes a sufficient repetition count without constructing those copies.
These are not ordinary normalized-state rows or new Lean theorems.

Three historical fixed `D=0,F2` searches with 512 MB/120 seconds,
2048 MB/120 seconds, and 2048 MB/300 seconds reached memory or native time
caps. The last used independently verified source automorphisms to remove
some GL2 gauge freedom. None found maps or established impossibility.
Later unit, matrix and dot-sum searches were stopped or omitted after the
independent mathematical exclusions; they were not solver UNSAT results.

The shared Cq family is now excluded for every q>=1, every finite D and
every positive gain over every field. [The proof](shared-dot-catalysts.md)
uses rational-function commutator transfer for q=1 and q>=4, and integral
Koszul certificates for q=2,3. The C3 certificate is an original 153-unit
minor with determinant one. C2 uses 903 unimodular integer row/column
operations giving a triangular 124-unit minor; its original minor is not
claimed to be a unit. Both certificates are rebuilt and replayed by the
standard-library [`shared_dot_koszul.py`](../shared_dot_koszul.py).

The earlier gauged C3/F2 search timed out after 913.358 seconds, with
900 native seconds and 2048 MB. The C4/F2 search was stopped after 765.380
seconds when its all-field obstruction was audited. Neither returned maps
or UNSAT. Independent coefficient automorphism controls remain saved.

The mixed-dot tensor itself has border rank exactly six in every field.
[Its adjugate proof](mixed-dot-border-rank.md) has a single obstruction
coefficient minus one; [`mixed_dot_border.py`](../mixed_dot_border.py)
checks the actual cofactor and six-term upper. No five-term degeneration
can improve that source budget. The adjugate identity alone does not exclude
arbitrary D; the separate replicated-rank argument above does.

The compact LM/ML shared-word tensor is `(4,4,4)` with eight unit coefficients.
Coordinate search for three units exhausted 30 nodes. A separate complete
coefficient SAT encoding returned UNSAT over F2 in about 0.88 seconds;
the F4 encoding reached a 30-second cap. The F2 result concerns this fixed
tensor and field only. Neither outcome proves an arbitrary-field subrank
upper bound. A separate `2 S_star -> M2` F2 encoding was UNSAT in about
1.17 seconds; this does not replace the all-field square construction.

Large literal replay artifacts and generated logs remain outside the
repository. Status summaries are experiment records, not certificates.
Reproduction creates fresh literal artifacts, and replay checks their actual
coefficients and exact rational/integer identities.

## Reproduce

From the repository root, choose a local output folder. SciPy is optional
for LP proposals; `z3-solver` is optional for direct SAT experiments. Exact
coefficient verification uses the standard library.

```sh
python reference/research/variable-proof-experiments.py --fields 4 --strong --output /tmp/variable-proof/base
python reference/research/variable-proof-experiments.py --refine /tmp/variable-proof/base/f4/finite-replay.json --output /tmp/variable-proof/refined
python reference/research/variable-proof-experiments.py --cut /tmp/variable-proof/refined/finite-replay.json --fill-existing-detectors --wall-seconds 600 --output /tmp/variable-proof/cuts
python reference/research/variable-proof-experiments.py --cross-level /tmp/variable-proof/cuts/finite-replay.json --output /tmp/variable-proof/connected
python reference/research/variable-proof-experiments.py --mixed-dots /tmp/variable-proof/connected/finite-replay.json --shared-dot-connections --output /tmp/variable-proof/shared
python reference/research/variable-proof-experiments.py --replay /tmp/variable-proof/shared/finite-replay.json
python3 -m unittest discover -s reference -t . -v
```

Hard caps are two million tensor entries and one million map and constraint
entries. The original family's separate 100,000-entry operational cap applies
to detectors, context lifts and direct-sum identities; swaps use the hard caps.
Coordinate proposals use 100,000 nodes, support 2048, diagonal 64 and one
second per search. Native LP proposals use ten seconds. Construction windows
are checked between operations; full validation and replay are timed
separately and may extend beyond a search window.

The next useful discovery target is an actual restriction connecting the
surviving proof tensors and matrix contexts that eliminates the remaining
finite states, together with a verified negative-unit dual. More scalar
lower rows or unsupported multiplicativity equations do not supply that
witness. The unrestricted existence theorem still guarantees some finite
auxiliary, without guaranteeing any of these prescribed small candidates.
