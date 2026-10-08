# Bounded coefficient search and the smallest catalyst candidates

This investigation concerns actual finite coefficients and restriction maps,
not the conditional numerical cost planner. It produced an unseeded exact
seven-term 2 by 2 matrix scheme, an unseeded four-term F4 C(2,3) scheme,
and a complete F2 tensor-cube inventory. It
has **not** found a useful `d=2,k=5` positive catalytic restriction.
The tested singular-star case with `D=0` is now independently excluded over
every field and for every positive scalar gain. The same auxiliary with
nonzero `D` remains unresolved and was not searched here.

The standalone executable is `reference/catalyst_search_experiment.py`.
Its default mode and inventory mode use only the Python standard library.
The optional SAT mode uses a solver installed outside the repository; it
does not add a package dependency to the ordinary reference implementation.

## An actual unseeded matrix coefficient result

The square 2 by 2 coefficient equations were solved without importing a
known scheme, supplying any seed coefficients, or constraining terms to a
known formula. The solver used 84 unknown coefficient bits and all 64 formal
tensor coefficient equations. Its symmetry conditions are complete for
rank at most seven:

* delete zero terms;
* merge equal nonzero input-factor pairs by adding their output vectors;
* sort distinct active input-factor pairs;
* pad with zero terms at the end.

Over F2 there is no nontrivial scalar normalization. Any at-most-seven-term
scheme can be represented in this form, so an UNSAT result would concern
the full rank-at-most-seven problem, not a planted support template.

With `z3-solver 5.1.0.0`, `random_seed=0`, and a 30,000 ms timeout, the run
found an actual seven-term family in approximately 0.63 seconds, including
Python startup and independent checking. This is an observation for this
run, not a runtime guarantee. After SAT supplied the input factors, output
coefficients were independently recovered by GF2 Gaussian elimination.
The complete tensor checker verified all 64 coefficients. A separate run
verified all 256 pairs of F2 input matrices and all 64 coefficients after
embedding the scheme into F4.

The resulting coefficient arrays and origin metadata are preserved in
[`unseeded-f2-rank7.json`](unseeded-f2-rank7.json). The booleans in that JSON
record performed checks; consumers should reconstruct the scheme and run
the coefficient checker themselves rather than trust those flags.

This is a meaningful coefficient-finding result, but seven-term 2 by 2
matrix multiplication is already known. Its finite-block exponent is
`log_2(7)`, which is above 9/4. It is not a decomposition extracted from
the repository's 9/4 proof.

The command used for the search was:

```sh
python3 -m reference.catalyst_search_experiment --mode sat \
  --size 2 --rank 7 --solver-path /path/to/optional-z3 --timeout-ms 30000
```

The solver path denotes a user-specified optional installation outside the
repository; replace it with the directory containing an installed z3
Python package. The solver is not a repository dependency. Without z3,
this mode reports `status="unavailable"`.

## Solve the third local map exactly

For fixed source `U`, target `V`, and first two local maps `A,B`, define

```
H[(x,y),w] = sum_(u,v) A[x,u] * B[y,v] * U[u,v,w].
```

The unknown third local map `C` must satisfy

```
H * C^T = V_flat,
V_flat[(x,y),z] = V[x,y,z].
```

`solve_third_leg` performs one row reduction with all target output slices
as right-hand sides. It sets free variables to zero and then checks the
entire restricted tensor using the existing independent implementation.
It supports explicitly represented F2 and F4, using bit-packed elimination
for F2 and exact field arithmetic for F4. This works for general local
restrictions, not only a rank-one source.

`bounded_restriction_search` enumerates the first two map matrices over F2.
It filters impossible flattening ranks and impossible map ranks before
solving the third leg. It makes the result distinction explicit:

| Status | Meaning |
| --- | --- |
| `found` | Actual maps were returned and the complete coefficient identity passed |
| `exhausted` | Every map pair for this exact fixed source and target was examined |
| `filtered` | A necessary flattening rank inequality rules out this fixed restriction |
| `resource_cap` | The declared search budget ran out; no nonexistence claim |
| `unavailable` | An optional solver dependency is absent |

The default executable demonstration produced the following results:

| Fixed restriction problem | Result |
| --- | --- |
| Four independent scalar tensors to `2 by 1` times `1 by 2` | Found actual maps at pair 15,451 of 65,536 possible first-two-map pairs |
| Two independent scalar tensors to F2 cube code 22 | Exhausted all 256 map pairs; no restriction exists |
| The four-scalar problem with a one-pair cap | Incomplete resource-cap result after one pair |
| One scalar tensor to `2 by 1` times `1 by 2` | Rejected by exact flattening ranks without enumeration |

Cube code 22 has support at `(0,0,1),(0,1,0),(1,0,0)` and flattening ranks
`(2,2,2)`. Thus the exhaustion example demonstrates an obstruction beyond
the elementary flattening filter.

Run this complete small demonstration with:

```sh
python3 -m reference.catalyst_search_experiment
```

Raw two-leg enumeration remains impractical for ordinary catalytic targets.
Even `D=0,S=Unit,d=2,k=5,m=1` has 5 by 5 matrices on both first legs, or
`2**50` first-two-map pairs. Sound algebraic obstructions should eliminate
such families before enumeration.

## Complete inventory of the 256 F2 tensor cubes

All coefficients of shape `(2,2,2)` were enumerated. For each tensor, the
search examined every subset of normalized nonzero input-factor pairs at
increasing ranks zero through three and solved the output coefficients.
There are nine input-factor pairs: three nonzero two-dimensional vectors
on each leg. Exhausting all smaller ranks establishes minimum rank;
the first found certificate establishes the matching upper bound.

Every found scheme was independently coefficient-checked and tested on all
16 pairs of two-dimensional F2 inputs. The result was:

| Exact rank | Number of tensors |
| --- | --- |
| 0 | 1 |
| 1 | 27 |
| 2 | 162 |
| 3 | 66 |

Exactly 174 cubes are concise: all three flattening ranks equal two. Of
these, 108 have exact rank two and 66 have exact rank three. Every concise
cube has an invertible output slice already over F2.

The experiment also enumerated all six elements of `GL_2(F2)` on each leg
and partitioned all 256 cubes into eight actual invertible-local-map orbits.
Code `c` means coefficient number `i` is bit `i` of `c`, with z fastest.

| Minimum orbit code | Orbit size | Exact rank | Flattening ranks | Invertible nonzero output slices |
| --- | --- | --- | --- | --- |
| 0 | 1 | 0 | 0,0,0 | 0 |
| 1 | 27 | 1 | 1,1,1 | 0 |
| 6 | 18 | 2 | 1,2,2 | 0 |
| 18 | 18 | 2 | 2,1,2 | 0 |
| 20 | 18 | 2 | 2,2,1 | 2 |
| 22 | 54 | 3 | 2,2,2 | 2 |
| 24 | 108 | 2 | 2,2,2 | 1 |
| 107 | 12 | 3 | 2,2,2 | 3 |

For the concise orbits, these correspond to the split diagonal/GHZ tensor,
the dual-number/W tensor, and the nonsplit quadratic-field tensor. The
finite orbit partition and coefficient ranks above were actually computed;
using these identifications for arbitrary catalysts requires the additional
mathematical obstructions investigated in `catalyst-obstructions.md`.

```sh
python3 -m reference.catalyst_search_experiment --mode inventory
```

That complete run took approximately 0.62 seconds in this environment.
Again, the observation is not a speed guarantee. It narrows the useful
auxiliary-tensor inventory without running large blind catalytic SAT
instances. Zero padding does not create new concise tensor types.

## C(2,2): weak filters pass, a stronger obstruction rules it out

For `d=2,k=5`, simple scalar or diagonal auxiliary families have strong
commutator obstructions. The independent audit develops those arguments
and their arbitrary-catalyst extension in `catalyst-obstructions.md`.
This search investigation did not assume rank additivity for arbitrary
direct sums.

The first candidate beyond the cube inventory was the convolution tensor
`S=C(2,2)`, shape `(2,2,3)`. Its polynomial product is

```
(a0+a1*t)*(b0+b1*t)
  = a0*b0 + (a0*b1+a1*b0)*t + a1*b1*t**2.
```

The standard Karatsuba coefficient identity supplies a rank-three scheme
over F2, and the output flattening rank three proves that this is exact.
For `D=0,m=1`, a positive catalyst would be the fixed restriction

```
source = five copies of C(2,2),              shape (10,10,15),
target = Unit direct_sum (M_2 tensor C(2,2)), shape (9,9,13).
```

All three elementary flattening inequalities pass. Two-leg map enumeration
has `2**180` pairs before solving the third leg. A direct unrestricted
coefficient SAT formulation has 375 local-map bits and 1,053 cubic XOR
equations; each equation receives up to 20 source support products. These
counts are substantial despite the tiny auxiliary tensor.

Small GF2 Koszul-flattening calculations initially failed to exclude this
candidate. With target axes ordered `(9,13,9)`, the computed ranks were:

| Koszul degree p | Matrix dimensions | Computed matrix rank | Rank-one factor | Tensor-rank lower bound |
| --- | --- | --- | --- | --- |
| 1 | 324 by 117 | 116 | 8 | 15 |
| 2 | 756 by 468 | 420 | 28 | 15 |
| 3 | 1,134 by 1,092 | 824 | 56 | 15 |
| 4 | 1,134 by 1,638 | 970 | 70 | 14 |

The lower bound is the ceiling of matrix rank divided by
`choose(dimA-1,p)`, the rank of that flattening on one nonzero simple tensor.
This does not use tensor-rank additivity. Source rank is at most 15, so
these computed bounds reach but do not exceed the source budget. Equality
of a lower bound and an upper budget does not prove a restriction exists.

The independent algebra audit then found a stronger obstruction. Selecting
output coefficients zero and one from `C(2,2)` gives the multiplication
tensor of `K[t]/(t**2)`. The executable verifies this actual selection map
and the restricted Karatsuba scheme coefficient by coefficient.

For q independent copies, the target restricts to multiplication in the
algebra

```
(K direct_product Mat_2(K[t]/(t**2)))**q,
```

which has dimension `9*q` and `2*q` maximal two-sided ideals. The classical
Alder–Strassen lower bound gives rank at least `16*q`. The source has an
explicit upper bound `rank(D)+15*q`. Repeating a putative positive
comparison with the **same fixed D** would imply `q <= rank(D)` for every
q, a contradiction. This excludes `C(2,2)` for any catalyst `D`, not just
`D=0`, assuming that classical algebra-rank theorem. Its dependency and
verification status are documented in `catalyst-obstructions.md`; this is
not a new Lean-checked theorem.

This change demonstrates why passing Koszul filters is not evidence that a
catalyst exists. We avoided running a large blind SAT problem for this
already-obstructed candidate.

## Field-sensitive convolution ranks and surviving boundaries

The standalone families mode now performs actual rank and restriction
checks and reproduces the small exterior-flattening calculations:

```sh
python3 -m reference.catalyst_search_experiment --mode families
```

For `C(2,3)` over F2, all `choose(21,4)=5,985` normalized four-pair input
subsets were exhausted without a solution. A five-term coefficient scheme
was then found without a seed: pair subset number 637 at rank five, or
6,622 total subsets including rank four. All coefficients passed the
independent checker.

A faster exact lower-bound filter inspects the output slice space `V`.
If the tensor had rank exactly `dim(V)`, every input outer product in a
decomposition would lie in `V` and those products would span `V`.
Enumerating all input-factor pairs shows:

| Convolution over F2 | dim(V) | Rank-one input pairs in V, encoded as bit masks | Their span dimension | Exact rank |
| --- | --- | --- | --- | --- |
| C(2,2) | 3 | (1,1), (2,2), (3,3) | 3 | 3 |
| C(2,3) | 4 | (1,1), (2,4), (3,7) | 3 | 5 |
| C(3,3) | 5 | (1,1), (4,4), (7,7) | 3 | 6 |

For `C(3,3)`, the deficient rank-one span proves rank at least six; a
six-product diagonal/cross coefficient identity supplies a matching exact
upper bound, independently checked. Thus no million-subset search was
necessary to settle that finite rank.

These ranks depend on the field. Actual F4 coefficient schemes were also
constructed and checked:

* C(2,3) has four terms by interpolation at all four finite field elements.
* C(3,3) has five terms by using those four finite nodes and the leading
  coefficient as the projective infinity node. Subtract its `t**4`
  contribution before interpolating the remaining degree-three polynomial.

The output flattening gives the matching lower bounds, four and five. This
does not conflict with the F2 lower bounds because extending the field
provides additional rank-one directions in the same output slice space.

The presently computed target filters remain inconclusive for the larger
convolutions:

| Auxiliary S | Best computed unprojected Koszul lower bound for Unit plus M_2 tensor S | Source bound 5*rank(S) over F2 | Source bound over F4 |
| --- | --- | --- | --- |
| C(2,3) | 20 | 25 | 20 |
| C(3,3) | 24 | 30 | 25 |

These are fixed-`D=0` bounds, not nonexistence statements for arbitrary D.
The algebra audit may yield stronger restrictions; a passed table entry
still does not assert a catalyst.

Another candidate surviving these elementary computations was the full outer-product tensor
`S=M(2,1,2)`, shape `(2,2,4)`, exact rank four. Then
`M_2 tensor S` is reindexed as `M(4,2,4)`. The target `Unit plus M(4,2,4)`
has shape `(9,9,17)` and source rank budget 20. Unprojected Koszul degrees
one through four gave ranks `152,524,960,1054` with factors `8,28,56,70`,
respectively: the largest resulting tensor-rank lower bound was 19.

Twelve deterministic full-rank first-leg projections from dimension nine
to eight, using random seed zero, improved the lower bound to 20 in ten
cases. For example row bit masks

```
323, 209, 488, 453, 266, 63, 14, 95
```

give Koszul degree-one matrix rank 135 and rank-one factor 7, hence lower
bound 20. This matches the source budget; it does not exceed it. These are
exact GF2 matrix ranks of actual projected tensors. No guessed general
rectangular-matrix rank formula was used.

The independent audit subsequently ruled out this outer-product family
using an inspected rectangular-multiplication lower-bound theorem. Repeating
the catalyst and merging first inputs restricts the targets to `M(4,2,4q)`,
whose cited lower bound is `24q`, while the source budget is
`rank(D)+20q`. This excludes every fixed D; it is stronger than the computed
Koszul bound. The external theorem, coefficient-level merge, and verification
scope are recorded in `catalyst-obstructions.md`. We retain the numerical
checks above as reproducible weaker filters, not as evidence of a surviving
outer-product catalyst.

A slot-only Kronecker ansatz for the local maps with `D=0` cannot provide
the desired `k=5` comparison: maps that act independently on the copy/slot coordinates
and leave S unchanged would restrict five scalar tensors to `M_2`, after
extracting a unit from nonzero S. That would be an ordinary rank-five
2 by 2 matrix scheme, inconsistent with its rank-seven lower bound. A
useful catalyst must mix the auxiliary coordinates nontrivially. This
ansatz exclusion is an algebraic statement about that restricted map
family, not nonexistence of unrestricted catalysts.

The current work has therefore found an actual unseeded coefficient
generator for small known schemes, useful exact rejection filters, and
field-sensitive candidate inventories. It has not supplied the large
spectral-obstruction catalyst needed for the near-9/4 extraction route.

## Reusable projective convolution construction

The previously experimental finite-plus-infinity interpolation is now
available independently as

```python
from reference.projective_convolution import projective_convolution_scheme
```

For positive input lengths a,b it requires exactly `a+b-2` distinct finite
nodes and adds one infinity term, for a total of `a+b-1` terms. Default
nodes use canonical field encodings and may include zero. It does not
silently enlarge the field, accept duplicate nodes, or pretend that F2
supplies enough points for C(2,3). Empty axes return a zero-term scheme
with no nodes; scalar multiplication uses infinity alone.

Nine focused tests passed, including formal reconstruction of every
monomial, the nontrivial lower coefficients in the infinity output vector,
exhaustive small-input checks over F2/F3/F4, non-base F16 nodes, rejection
controls, and solver-free regression of both saved unseeded coefficient
schemes.
No Lean theorem was added or checked for this API.

The projective source also connects to the existing proof pipeline:
`a=2,h=1` over F4 has a five-term C(2,4) source. Its second sector power
recovers in 75 terms and descends to F2 in 300 terms, with complete target
coefficient checks. Ordinary finite-node interpolation would require five
finite points for that source; adding infinity permits the smaller F4
coefficient field. Recovery itself still uses the existing distinct nonzero
finite interpolation nodes, and the field is fixed before the power.

## Independent checks of the experiment machinery

In addition to the result-specific full coefficient checks, the GF2
all-right-hand-side solver was compared with a complete brute-force oracle
on all 256 pairs of two by two coefficient/RHS matrices. Returned solutions
and inconsistency results matched in every case. The Koszul rank-one factor
was verified on all 63 nonzero simple tensors of shape `(3,2,2)` at exterior
degrees zero, one, and two. These finite controls support the computations;
they do not replace the classical algebra-rank theorem used in the stronger
candidate exclusions.

## Bounded full-map searches at the (2,3,3) boundary

Two actual F2 auxiliary tensors were tried with `D=0,d=2,k=5,m=1`:

1. The regular pencil with slices `(I_3,J_3)`, where `J_3` has ones directly
   above the diagonal.
2. The singular pencil

   ```
   [[x,y,0],
    [0,0,x],
    [0,0,y]].
   ```

Both auxiliaries have exact rank four. Unseeded normalized-input searches
found four-term schemes after 1,720 and 1,413 subsets, respectively,
including all exhausted rank-three subsets. Complete output-space rank-one
inventories separately proved that rank three is impossible, and the
returned schemes passed all coefficient checks.

For both fixed catalytic problems, source shape is `(10,15,15)` and target
shape is `(9,13,13)`. The unrestricted local maps have 480 unknown bits and
must satisfy 1,521 complete coefficient equations. The regular pencil's
source support gives 38,025 cubic product occurrences; the singular pencil
gives 30,420.

The optional `bounded_restriction_sat` API uses all these coefficient
equations. Its only symmetry restriction sorts first-leg map blocks after
checking that the source consists of five identical full direct-sum
blocks. Any witness can be put in that order by jointly permuting source
blocks on all three legs, so this does not select a speculative support
pattern. Before encoding, a cubic-product cap limits formula construction;
the solver also receives a 60-second timeout and a 512 MB memory limit.

Both runs returned `status="resource_cap"`, reason `timeout`, after their
60-second solver budgets. No actual catalytic certificate was returned.
The timeout outcomes themselves do not establish nonexistence, even for
the fixed F2 cases. Separate algebraic arguments have since excluded both
tested restrictions with the scopes stated below. The search outcomes do
not supply conclusions about other auxiliary tensors or a nonzero catalyst
for the singular star. If SAT returns a model, the implementation
checks the complete restricted tensor, re-solves the third leg by an
independent linear system, and invokes `CatalyticCertificate.require_valid`
before reporting an actual positive catalyst.

The commands were:

```sh
python3 -m reference.catalyst_search_experiment --mode catalyst \
  --pencil-kind regular --solver-path /path/to/optional-z3 --timeout-ms 60000
python3 -m reference.catalyst_search_experiment --mode catalyst \
  --pencil-kind singular --solver-path /path/to/optional-z3 --timeout-ms 60000
```

The generic SAT encoder was checked against complete two-leg map enumeration
for all 256 F2 tensor cubes as targets of two independent scalar tensors.
Every status matched: 190 feasible restrictions and 66 exhaustive negative
results. Every returned SAT restriction was coefficient-checked and its
third map independently recovered. This gives a finite soundness control
for the encoder and block-ordering symmetry. The solver decisions remain
incomplete; the independent algebraic arguments resolve these particular
fixed problems without using the timeout results as premises.

The regular pencil is now excluded by a separate commuting-pencil
compression obstruction, independently checked by the constructive-route
and existence-audit investigations. Its argument and scope are recorded in
`catalyst-obstructions.md`; the SAT timeout was not a premise. The singular
candidate lacks a regular slice and was analyzed separately.

The subsequent first-map obstruction in
[`star-zero-catalyst.md`](star-zero-catalyst.md), independently checked by
the main, constructive-route, and existence-audit investigations, excludes

```
5*S_star -> m*Unit direct-sum (M_2 tensor S_star)
```

over **every field and every positive integer m**. Equivalently, the
`D=0,S=S_star,d=2,k=5` positive-gain restriction is impossible even with
arbitrary mixing of the three local maps. This supersedes the original
fixed F2 search question; it is a separate proof, not a reinterpretation
of SAT's `resource_cap` status. A nonempty catalyst `D` changes the
dimensions and rank-locus argument and is **not** excluded by this proof.
No nonempty-D singular-star search or impossibility proof was obtained
in this investigation. The generic restriction SAT tool remains usable
for other supplied finite tensors.

One additional bounded search targeted an actual lower-bound embedding
rather than a positive catalyst:

```
M_2 tensor S_singular -> M(3,2,4).
```

Source shape is `(8,12,12)` and target shape `(6,8,12)`. This has 288 map
bits, 576 complete coefficient equations, and 18,432 cubic product
occurrences. It likewise reached its 60-second solver timeout, reporting
`resource_cap`; no embedding was found or disproved. No further dense SAT
rounds were started.

The proposed motivation was to append a separate scalar identity to a
verified embedding and merge q rectangular blocks by sharing their first
input. The inspected n=3 rectangular lower theorem would then give

```
q + ceil((96*q+2)/5) <= rank(D)+20*q.
```

This is impossible once `q+2 > 5*rank(D)`. However, a separate image-space
obstruction now proves that the required core embedding does **not** exist
over any field; this particular route cannot establish the desired
catalyst exclusion. The source slice associated to two 2 by 2 matrices
A,B has rank `2*(rank([A;B])+rank([A B]))`. A rank-eight target slice
forces both source ranks to be two, and its source image therefore
contains a fixed four-dimensional third-output block. The output local
map must be invertible because both output dimensions and the target
output flattening rank are twelve. Thus images of every rank-eight target
slice would share that transported four-dimensional subspace. The three
target slices associated to the column planes 01, 02, and 12 have zero
common image, a contradiction. The existence-audit report records the
full argument in the final section of [`pencil-frontier.md`](pencil-frontier.md).
That report also rules out adding the same finite number of scalar
summands on both sides, including arbitrary mixing of those summands.
Saturation forces a common source image of dimension `m+4`, while the
target plane-image intersection has dimension only `m`.

An independent check here exhausted all 256 F2 first-leg slices. The slice
rank counts were `{0:1,4:27,6:36,8:192}`; every rank-eight image contained
the four source output basis vectors at indices 2, 5, 8, and 11. Complete
enumeration of the three target plane images verified their dimensions
eight and their common image dimension zero. These finite checks support
the coordinate argument; the every-field conclusion comes from the proof.

This argument excludes the proposed rectangular restriction, rather than
all positive catalysts using the singular auxiliary. It also explains why
the earlier solver timeout must remain separate from an impossibility
proof. The experiment wrapper remains an explicitly requested historical
search, with full coefficient checking if a solver ever claims SAT.

```sh
python3 -m reference.catalyst_search_experiment --mode embedding \
  --solver-path /path/to/optional-z3 --timeout-ms 60000
```

## Extension-field coefficient search with an actual result

The generic optional SAT API also supports the explicit field
`F4 = F2[alpha]/(alpha^2+alpha+1)`. Each map coefficient is two canonical
power-basis bits. Every cubic product is expanded using exact F4 basis
multiplication, including arbitrary nonbase source coefficients. Its two
coordinate equations are formal tensor identities; there is no reduction
modulo identities of field-valued functions. Construction caps count both
field products and expanded Boolean product occurrences before building
the formula. The API never selects an extension field automatically.

For the rejected rectangular embedding, a hypothetical F4 encoding would
have 288 field variables, 576 map bits, 1,152 bit-coordinate equations,
and 202,752 Boolean cubic-product occurrences. **That dense F4 solver run
was not started**: the field-independent image obstruction arrived before
launch. Its absence is not another timeout or a SAT result.

Instead, small controls demonstrated a real field-sensitive extraction.
Without a planted scheme or projective-interpolation factors, the generic
SAT equations found a restriction of four scalar tensors to C(2,3) over
F4. The independently recovered output map passed all 24 formal tensor
coefficients and all 1,024 input pairs. The actual four-term arrays, field
modulus, solver version, and deterministic seed are saved in
[`unseeded-f4-c23.json`](unseeded-f4-c23.json). A solver-free regression
reconstructs and verifies that artifact. Earlier complete F2 searches
proved that this same convolution tensor has rank five over F2, so the
extension-field mode changes an actual finite coefficient problem.

The same fixed unseeded coefficient problem can be reproduced with:

```sh
python3 -m reference.catalyst_search_experiment --mode convolution-sat \
  --field F4 --a 2 --b 3 --rank 4 --solver-path /path/to/optional-z3 \
  --timeout-ms 10000
```

The command reports its actual arrays and complete coefficient verification;
`--output` can save that report. The archived JSON is a separate portable
certificate consumed by the solver-free regression.

The exact F4 Gaussian solver matched exhaustive four-element solution
enumeration on 128 small two-variable systems. Twelve scalar controls
covered every nonzero F4 source coefficient and all four target
coefficients. A rank-two diagonal source was also correctly rejected for
the W tensor over F4. Finally, all 256 F2 cube controls were rerun with the
unified encoder: the original 190 feasible and 66 exhaustive negative
results were preserved, and every returned map passed independent checks.

This establishes an executable extension-field coefficient finder and
actual small witnesses. It does not extract the 9/4 decomposition or
provide a useful positive catalyst.
