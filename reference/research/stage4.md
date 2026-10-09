# Stage 4 checkpoint: bounded discovery and target calibration

**No useful `d=2,k=5` certificate or new exponent-2.4 matrix decomposition
was found.** This checkpoint adds executable bounded searches and independently
audited mathematical restrictions. It completes a reporting checkpoint,
not the stage's discovery target. No Lean sources were edited or rebuilt;
the new mathematical arguments and Python checker are not Lean theorems.

This report preserves the initial experiment results. The later
[variable-auxiliary checkpoint](variable-auxiliary-checkpoint.md) supersedes
its small-catalyst frontier and the viability of its random `(4,4,5)` candidate.
The expanded families and new square-to-matrix restrictions are documented
there; historical row counts below describe their original source versions.

## The principal finding changes the search priority

The singular `(2,3,3)` star has an actual three-term polynomial degeneration
`P(t)=t*S+t²*E`, verified at every formal coefficient. Since each first factor
is constant, extracting `[t^n]P(t)^n` gives an exact decomposition over the
original field with at most `(n+1)*3^n` terms. No interpolation nodes or field
extension are needed. `reference.star_degeneration` constructs this identity
and bounded small powered coefficient schemes.

Consequently every normalized nonnegative, additive, multiplicative,
restriction-monotone character has `phi(S)<=3`. A real fixed-star certificate

```text
D ⊕ m Unit ⊕ (M₂ tensor S) <= D ⊕ 5S,       m>=1
```

would imply `phi(M₂)<=5-m/phi(S)<=14/3`. Using the existing detecting-character
and exact-rank-to-arithmetic bridges, this would give
`omega<=log₂(14/3)≈2.222392421`, strictly better than `9/4`.
The strict comparison is already integer-exact: `14**4 < 3**4 * 2**9`.

Thus the existing `9/4` theorem does **not** guarantee a certificate for this
fixed small star. It guarantees some auxiliary, not every prescribed auxiliary.
This is not a proof that the star target is impossible. A prescribed direct
sum of star copies is also not automatically guaranteed by the theorem.
The [full argument](star-convolution-catalysts.md) names the existing Lean
bridge declarations and the field-descent qualification. It is separate from
the current dense compiler's supplied four-term star cost recurrence.

The main discovery priority should therefore vary auxiliaries generated from
the proof and their mixtures. Fixed-star searches remain a limited research
branch; increasing only their catalyst size would chase a stronger result.

## Independently audited star obstructions

The [nonzero-catalyst argument](star-nonzero-catalyst.md) derives necessary
conditions from actual exact restrictions, allowing arbitrary mixing of the
three maps. Concise dimensions of `D` must satisfy:

- First flattening rank at least three; other ranks at least four.
- The other two ranks cannot both be four.
- A coefficient-field-defined two-dimensional slice space has every nonzero
  matrix slice of rank at least four, including after algebraic closure.
- Same-field tensor rank is at least six over `F₂/F₃`, and at least five
  over larger finite fields.

The [rank-drop refinement](star-convolution-catalysts.md) adds:

- At first rank three, generic slice rank and both other ranks are at least
  five. Every nonzero slice in the whole concise first space has rank at
  least four.
- That rank-drop cone must have at least three linear plane components with
  independent normals. For a concise `(3,5,5)` target its determinant must
  have corresponding independent linear factors over the witness field.
- A generic rank-drop cone contained in one proper linear hyperplane excludes
  a catalyst for any first rank.
- A first-concise three-dimensional diagonal catalyst is excluded regardless
  of its diagonal length.

These exclude the complete `(3,4,6)` and truncated `(3,4,5)` convolution
catalysts, the `(4,4,5)` band pencil `L₄(a,b)+cE₀₀+dE₁₁`, and the tested
`(3,6,6)` diagonal catalyst. These originally left small dimension boundaries
unresolved. The later [slice-space audit](variable-auxiliaries.md) excludes
every first-concise star catalyst dimension below six. Its necessary shape
boundaries are now `(6,6,10)`, `(6,10,6)` and `(6,8,8)`; they are not witnesses.
The whole fixed-star problem remains unresolved.

`reference.star_catalyst_screen` implements selected exact consequences.
Its scope fixes the restriction coefficient field. A base-field enumeration
does not compute generic rank or automatically exclude extension-field maps.
Its first-leg basis removes the flattening kernel before slice enumeration;
sampling caps never imply impossibility.

## Connected finite row experiments

`reference.discovery_inventory` adds actual scalar lower maps, checked rank
upper maps, both direct-sum identity directions, tensor-context lifts, and
detectors on multiple proof-derived and convolution auxiliaries. Equal tensors
share state coordinates. Context lifts use `T tensor M₂`, while detectors use
`M₂ tensor T`; explicit coordinate-swap maps connect those arrays. No rank
additivity or tensor cancellation is assumed.

`reference.finite_state_solver` optionally uses SciPy to propose either a dual
or a normalized state. Rational reconstruction must satisfy every equality or
inequality exactly. A dual is cleared over the integers; all ordinary maps and
detector products are checked. Numeric solver statuses are never certificates.
An exact feasible state proves absence of a negative-unit dual **only for this
selected finite row family**, rather than absence of arbitrary catalysts.

| Field/family | Atoms | Rows | Exact result |
| --- | ---: | ---: | --- |
| `F₁₁`, proof rows and three contexts, k=5 | 46 | 135 | Normalized feasible rational state |
| `F₄`, available proof rows and contexts, k=5 | 32 | 91 | Normalized feasible rational state |
| `F₂`, supplied rank-seven control, k=8 | 11 | 28 | Gain-one dual and actual catalyst |

The last dual uses the supplied known-rank row and unit detector only. It is a
control, not a proof-derived discovery. Missing Fourier/node availability is
recorded; no field extension is inserted silently. Construction, LP proposals,
map verification and dense allocation have distinct reported budgets. The
native LP/SAT time limits are advisory; total wall times include encoding and
verification. These families are finite constructor experiments, not a
complete implementation of the spectral or entropy existence step.

## Direct search and remaining SAT outcomes

`reference.flip_search` performs seeded flips and two-to-three-term moves,
each an exact tensor identity. It independently checks improved and final
schemes. Over `F₂`, seed zero turns the naive eight-term matrix decomposition
into seven terms in 45 steps. A separate unseeded SAT rank-seven matrix control
also produces actual maps and independently recovered coefficients. The `F₄`
flip control reaches seven; the two tested `F₁₁` controls hit their five-second
caps at eight. No known seven-term scheme is supplied to these direct searches.

The star-product searches do not improve their supplied decompositions: four
60-second searches started from the known 28-term product; subsequent searches
started from the naive 32-term product. These caps are not rank lower bounds.
The direct rank-19 product goal and the catalytic restriction are different
problems: finding 19 product terms alone would not give a restriction from
five star copies. No performance superiority of either route is established.

Each of the four subsequently excluded catalysts was also tested with a
120-second, 512-MB SAT limit. Those runs timed out; their exclusion instead
comes from the independent mathematics above. One additional binary `(4,4,5)`
candidate avoids the recorded hyperplane obstruction:

```text
D(a,b,c,d)=L₄(a,b)+cC+dE,
C mask = 306381, E mask = 47558, row-major 4×5 bit order,
seed=20261008, first candidate chosen from a cap of 32.
```

Its `c=d=0` plane has rank four at every nonzero closure point. Its binary
rank-drop inputs `[4,7,8,12,13]` span all four first coordinates, so its generic
rank-drop locus cannot lie in one proper linear hyperplane. The single SAT
run timed out after 120 seconds (133.16 seconds including encoding), with
865 map bits, 3,978 equations and 186,966 cubic-product occurrences. This
candidate is now excluded by the stronger other-flattening bounds in
[the later audit](variable-auxiliaries.md). Its SAT timeout is a historical
solver outcome; the mathematical exclusion does not come from that timeout.

[stage4-results.json](stage4-results.json) records finite states, budgets,
selected coefficient tensors, exact-control attribution and search outcomes.
It is an experiment record; a metadata flag alone is not a verified witness.
Large generated logs and construction DAGs are not committed.

## Reproduce and independently replay

The ordinary reference API and mathematical controls use the standard library.
SciPy/NumPy and `z3-solver` are optional search dependencies. If absent, the
respective proposal search reports `unavailable`; empty finite families can
still supply exact states without them.

```sh
python3 -m reference.discovery_experiment --field 11 --output /tmp/discovery-f11
python3 -m reference.discovery_experiment --field 4 --core-seconds 20 --output /tmp/discovery-f4
python3 -m reference.discovery_experiment --field 2 --no-contexts --known-control --k 8 --core-seconds 10 --sat-seconds 30 --output /tmp/discovery-control
python3 -m reference.discovery_experiment --replay /tmp/discovery-f11/finite-replay.json
python3 -m unittest discover -s reference -t . -v
```

Replay reconstructs the graph and all tensor/map representatives, then checks
the normalized rational state or integer dual. It does not trust saved status,
row-vector, rank or shape flags. A recorded SAT case can be reconstructed from
its `D_coefficients`/`S_coefficients` and run through
`bounded_restriction_sat(D⊕5S,D⊕Unit⊕M₂⊗S,...)`; retain its saved construction
caps as well as its native time/memory limits. Solver versions may change the
search path, but proposed witnesses still need exact coefficient checks.

## Next bounded target

Increase and vary connected finite sector/type/determinant constructions and
their available tensor contexts, allowing the dual to choose auxiliary
mixtures. Use the exact feasible states to identify missing finite relations
before increasing SAT budgets. Keep direct decomposition controls as limited
baselines and filter the now-excluded random catalyst before new searches.
Continue recording gain,
auxiliary decomposition costs and eventual coefficient allocation together.
Scalable symbolic elimination and Lean formalization of the catalytic and
obstruction arguments remain separate work; neither supplies the missing
useful witness by itself.

## Checkpoint verification

On 2026-10-08, all **275 reference tests passed in 175.650 seconds**, including
56 new tests. The connected `F₁₁` feasible-state artifact and the `F₂` control
dual were independently replayed from their serialized actual tensors/maps.
The eleven protected Lean files and frozen Challenge model match the exact
baseline. No Lean build or new kernel theorem is claimed for this checkpoint.
