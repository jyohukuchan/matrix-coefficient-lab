# Coefficient extraction: research checkpoint

This directory records an investigation of the gap between the proved exact
rank exponent bound and an executable generator of useful coefficients.
The derivations below have been checked independently by research agents, but
are **not new Lean theorems**. No decomposition achieving the 9/4 bound has
been generated, and no speed claim is made.

## Concrete extraction target

Fix a represented finite field, initially F_2, and a positive rational slack
epsilon. Write tau = 9/4 + epsilon = A/B in reduced form. A successful output
consists of n >= 2, r, and three finite coefficient families giving an exact
r-term decomposition of the n-by-n matrix multiplication tensor, with

```text
r**B < n**A.
```

Every tensor coefficient must be verified. Agreement on sampled inputs is
insufficient. Positive slack is essential: an infimum bound need not be
attained at exactly 9/4.

The existing all-fields **exact coefficient rank** theorem implies that such
an output exists. Exhaustively enumerating finite coefficient families over
the fixed field is consequently a terminating algorithm in principle. This
does not provide a useful advance bound on the matrix size or search cost.
The finite-search report gives a tiny complete rectangular experiment and
an input-family/output-linear-system reduction. Even the F_2, 2-by-2,
seven-term input-pair subset search has 5,271,289,966,800 candidates.

## Where the source proof stops supplying coefficients

The sector maps, interpolation, Fourier separation, tensor powers, and descent
are finite transformations that have executable counterparts. The spectral
argument then reasons about real-valued characters rather than carrying
rank-one coefficient arrays through to the conclusion.

There are two infinite compactness steps in the state/character construction;
even the result named a finite multiplicative slice quantifies over the whole
tensor semiring. The final exponent-infimum argument selects an unspecified
finite block. Translating character inequalities alone will therefore not
automatically yield coefficients. This is a specific missing construction,
not an impossibility claim based on Lean's `noncomputable` annotation.

The [source audit](existence-audit.md) identifies the relevant declarations,
including `Character/Existence.lean`, `Spectrum/StateObstruction.lean`,
`Spectrum/MultiplicativeStates.lean`, and
`Arithmetic/RankExponent.lean:219–234` under `AuxiliarySeparation`.

## Conditional route through an explicit catalyst

An explicit finite certificate would provide tensors D,S, with S nonzero,
and three actual local linear maps proving

```text
D ⊕ m*unit ⊕ (M_d tensor S) <= D ⊕ k*S,
```

where d >= 2 and m >= 1. The symbol <= means an exact tensor restriction
from the right side to the left side. Integers m and k count independent
direct-sum copies; they are not field scalars reduced modulo the characteristic.

The existing character bound implies existence of such a positive-gain
certificate over the algebraic closure when the integer k > d**(9/4).
This follows by contraposing the normalized-state obstruction implication:
a normalized state would yield a character with value at least k on M_d,
contradicting the universal character bound. Over the algebraic closure of
a finite field, all coefficients of a finite certificate lie in one finite
extension. This supplies a search target, without a useful dimension or
extension-degree bound. Searching that particular existence route would
vary both extensions and tensor/map dimensions. A stronger same-field
reduction from the exact exponent theorem is recorded below.

The main agent derived, and the research agents independently checked, the
following way to turn such a certificate into a convergent coefficient
procedure. First select away the m scalar gain blocks; do not cancel D.
Let T=M_d and define, for b >= 1,

```text
D_b = ⊕_(i=0)^(b-1) k**(b-1-i) copies of (D tensor T**i).
D_b ⊕ (T**b tensor S) <= D_b ⊕ k**b*S.
```

The second restriction follows by sums, tensor products, and composition
using D_(b+1) = (D tensor T**b) ⊕ k*D_b. It keeps the same S.
Given exact supplied decompositions of S and D with R_S and C_D terms, a
conservative supplied decomposition bound for D_b is

```text
C_b = C_D * sum_(i=0)^(b-1) (d**3)**i * k**(b-1-i).
```

Select a nonzero coefficient of S to obtain an explicit unit <= S map.
Repeatedly absorb the same D_b using r_t copies of the powered restriction
and the supplied r_t-term decomposition of M_(u_t). This gives the recurrence

```text
u_0 = r_0 = 1,
u_(t+1) = u_t * d**b,
r_(t+1) = R_S * (k**b*r_t + C_b).
```

Tensoring the catalyst with M_(u_t) instead would pay its cost repeatedly;
the absorption construction above is needed. For tau=A/B, choose b with

```text
(R_S*k**b)**B < d**(b*A).
```

With Q=R_S*k**b>1, r_t is bounded by
Q**t * (1 + R_S*C_b/(Q-1)), whereas u_t=d**(b*t). Thus the strict target
inequality eventually holds. All steps stay in one fixed extension of
degree e; descend only the final scheme and test

```text
(e**2*r_t)**B < u_t**A.
```

For any tau>9/4, sufficiently large d allows an integer k strictly between
d**(9/4) and d**tau. The route therefore avoids assuming a near-optimal seed
decomposition once an explicit catalyst is supplied. The map compositions
now have an executable coefficient compiler. The open practical problem is
finding a useful certificate within feasible dimensions and resource costs.
See the [constructive-route report](constructive-route.md) for the source
connections and the [finite-search report](finite-search.md) for integer
planning and search contracts.

## What is implemented at this checkpoint

- `reference.extraction` extracts all three branches from a generated
  separated scheme, restores the middle branch's coordinate order, factors
  singleton-leg boundary cases by exact Gaussian elimination, and connects
  three cyclic boundary factors to a square matrix scheme. The default
  F_16 example starts with 750 separated terms and yields the classical
  eight-term 2-by-2 scheme. It verifies the connection, not the 9/4 claim.
- `reference.catalysts.CatalyticCertificate` verifies every coefficient of a
  supplied positive-gain certificate and can select away the gain blocks.
  Tests include a real F_2 certificate with d=2,k=9. Its k is too large for
  the desired 9/4 sufficient gap.
- `plan_conditional_costs` evaluates the integer recurrence from numerical
  assumptions. `plan_catalytic_parameters` first verifies an actual catalyst
  and supplied decompositions of D,S. Both return a **conditional numerical
  plan**, not an iterated matrix decomposition. Resource limits report an
  incomplete search, and failure of the sufficient gap is not a no-solution
  result.
- `reference.catalyst_compiler` materializes the powered restrictions and
  absorption maps, then produces exact matrix coefficient families. It
  verifies the maps and arrays independently, reuses one catalyst even when
  its coordinates genuinely contribute to the matrix output, and retains
  one field until final descent. Dense caps report incomplete construction.
- `reference.witnessed_constraints` searches a small finite rational cone,
  checks the dual exactly, clears denominators over the integers, and
  assembles a catalyst from ordinary restriction witnesses. Positive blocks
  supply D and detector inputs supply S; balanced tensor-ID occurrences
  determine explicit permutations. No unspecified group-completion witness
  needs to be searched for. See [witnessed-lp.md](witnessed-lp.md).
- `reference.constraint_generators` exports actual interpolation restrictions
  from sector tensor powers and lifts ordinary rows by a tensor context.
  This keeps every source copy needed by the finite construction; it does
  not turn an asymptotic character inequality into a finite order for free.
- `reference.finite_type_constraints` exports exact-type coordinate pulls and
  finite Fourier/square separation as ordinary witnessed rows. Tensor words
  retain their literal coordinate order. Fourier periods and interpolation
  copies remain in the positive block list; no entropy-limit overhead is
  silently removed, and an additive state is not treated as multiplicative.
- `reference.determinant_filtration` implements the determinant construction's
  actual quotient/kernel bases, output-dual change, and formal identity
  `P(t)=t*graded+t**2*cross`. It generates coefficient families from a supplied
  source decomposition and exports a finite interpolation row. Characteristic
  two retains the correct signs and independent branch blocks. This finite
  implementation does not supply the determinant entropy limit.
- `reference.geometric_catalyst` turns a supplied exact matrix-power rank gap
  into an actual same-field positive-gain certificate with D=0 and a specified
  direct sum of matrix powers. It checks all local coefficients and returns
  an exact conservative auxiliary decomposition.
- `projective_convolution_scheme` adds a leading-coefficient/infinity term,
  allowing C(3,3) to have an actual five-term scheme over F_4. A finite search
  finds C(2,3) and C(3,3) ranks five and six over F_2. These field differences
  matter when selecting candidate auxiliaries.

A synthetic arithmetic experiment in the finite-search report uses
d=4,k=23,R_S=2,C_D=1,e=2, inner exponent 23/10 and outer exponent 231/100.
It finds b=14 and the first successful final inequality at t=50, with matrix
size bit length 1401 and descended term-bound bit length 3234. No tensors or
maps with those costs are supplied. These numbers illustrate how large the
conditional construction can become; they do not establish a new algorithm.

Run the implemented connection with:

```sh
python3 -m reference.extraction_demo
python3 -m reference.catalyst_demo
python3 -m reference.witnessed_demo
python3 -m unittest reference.test_extraction reference.test_catalysts -v
```

The witnessed demo loads a supplied seven-term artifact originally found by
an unseeded SAT search, checks it, derives a finite dual and catalyst, and
generates matrix arrays. The 2-by-2 result compacts from 15 terms to seven;
the next step compacts from 63 to 49. This is a reconstruction of known small
matrix algorithms, not an extraction of the theorem's 9/4 decompositions.
The optional solver is needed only to repeat the search, not to check the
artifact or run the ordinary reference package.

## Retaining the positive gain

The original compiler first selected away all scalar gains. A stronger
finite construction retains them. Its powered comparison has gain

```text
g_b = m * sum_(i=0)^(b-1) k**(b-1-i)*d**i,
D_b ⊕ g_b*unit ⊕ (M_(d**b) tensor S) <= D_b ⊕ k**b*S.
```

Each matrix gain restricts to its independent diagonal scalar coordinates.
Repeated absorption keeps r*g_b gain blocks, decomposes D_b directly, and
eliminates these scalar summands by exact linear substitution. The improved
cost recurrence is

```text
r_next = C_b + (k**b*R_S - g_b)*r.
```

Scalar elimination removes only independent units, not a general catalyst.
For any exact scheme of T ⊕ unit, choose a nonzero final C coordinate and
solve for its term; projecting away the last X/Y coordinates gives an exact
scheme of T with one fewer term. This supplies actual arrays and proves
scalar-summand additivity without assuming general tensor-rank additivity.
See [gain-preserving-route.md](gain-preserving-route.md).

The bounded `compile_gain_target()` interface now selects a sufficient block
and iteration count using this exact recurrence, constructs the arrays, and
descends only the final scheme. It accepts a fixed block when
`(k**b*R_S-g_b)**B < d**(b*A)`, even if the original `k**B < d**A` sufficient
condition fails. Such a plan remains numerical until all generated arrays
and the final strict inequality have been checked. Dense caps return an
incomplete result rather than an uncertified rank claim.

## Narrowing the useful-catalyst search

The concrete target d=2,k=5 lies above 2**(9/4), while log_2(5)<12/5.
Its existence argument does not imply that a tiny auxiliary will work.
We derived necessary flattening, commutator, algebra-restriction, and
rectangular-sharing conditions. The [obstruction report](catalyst-obstructions.md)
specifies which results use an inspected external rank theorem, including
[Yang's arbitrary-field rectangular lower bound](https://arxiv.org/abs/2609.14393v2).
The reports and Python checks do not add those theorems to Lean.

In particular, candidates with a flattening rank one are excluded, and
concise shapes (2,2,2), (2,2,3), and (2,2,4) cannot work with any finite D.
The last case uses an actual restriction from repeated rectangular products
to a shared-input rectangular product. Outer-product auxiliaries M(a,1,c)
must satisfy ac >= 1+4*max(a,c), so min(a,c)<=4 is excluded.

At the next shape (2,3,3), a new commuting-source compression argument
excludes regular matrix pencils. It replaces D by a supplied diagonal
decomposition, compresses each pulled-back source pencil into commuting
square slices, and contradicts the target's nonzero matrix commutator as
the number of copies grows. It uses one catalyst throughout, rational-function
generic pivots, and no general rank-additivity assumption. The remaining
singular pencil has varying kernels and falls outside this proof. See
[pencil-frontier.md](pencil-frontier.md) and the independent
[Koszul/audit report](pencil-koszul.md).

Bounded full-map searches on specified F_2 auxiliaries have hit time caps.
Those solver results alone are incomplete searches. A separate new argument
now excludes the singular star **with D=0**, over every field and for every
positive gain: two seven-dimensional first-leg subspaces would have to lie
in the scalar-zero hyperplane, although they span all nine dimensions.
The nonzero-scalar fibers instead have dimension six and cannot fit inside
the union of the two dimension-five bad-matrix loci. See
[star-zero-catalyst.md](star-zero-catalyst.md) for the independently reviewed
proof and exact finite-count alternative. A nonzero arbitrary D remains
unresolved. No useful d=2,k=5 witness has been obtained.

An attempted shortcut from the singular product to `M(3,2,4)` is now excluded
over every field. Saturating its maximum slice rank would force all target
slice images to contain a fixed four-dimensional subspace, whereas three
target coordinate planes have zero image intersection. The proof also handles
equal scalar augmentations and arbitrary coordinate mixing. It closes that
particular rectangular route, while leaving the original catalyst unresolved.

Checkpoint validation on 2026-10-08: all **196 reference tests** passed in
151.656 seconds. Both catalyst demos and the independent obstruction/Koszul
experiments passed. The star-slice experiment checks every matrix pair over
F_2, F_3, and F_4 and all 1,023 binary first-leg hyperplanes. The eleven
protected source files and frozen Challenge model match the exact baseline.
No Lean sources were changed or rebuilt for these Python/research checkpoints.

## Revised stages 1–3 checkpoint

The [cost/DAG/integration checkpoint](stages-1-3.md) implements the revised
first three stages: integer allocation screening, checked structured
construction graphs, and automatic finite proof-row integration with exact
search and family-removal analysis. The standard-library-only control pipeline
reconstructs the known seven- and 49-term schemes. Its four proof-derived row
families are validated but unused by the found control dual; removing the
known-rank control exhausts this selected finite problem. The same selected
d=2,k=5 family also exhausts. Neither is a global impossibility result.

The graph replays a known scheme's 25th tensor power in 18 nodes, queries
specified coefficients, and counts support without dense expansion. General
restriction identities still require small fully checked seeds; interpolation/
Fourier formulas and scalar-gain elimination are not yet fully symbolic.
Numerical forecasts remain separate from checked construction certificates.
No new useful catalyst or exponent-2.4 decomposition is claimed.
Checkpoint validation: all **219 reference tests passed** in 158.442 seconds;
the eleven protected files and frozen Challenge model still match the baseline.
The 23 new tests cover forecasts, construction replay/expansion, corruption
rejection, resource caps and finite pipeline attribution. No Lean rebuild was
performed for these Python changes.

## Stage 4 discovery checkpoint

The [bounded discovery report](stage4.md) records connected finite LP experiments,
exact feasible-state certificates, seeded direct coefficient search, and further
independently audited [nonzero-star](star-nonzero-catalyst.md) and
[rank-drop](star-convolution-catalysts.md) obstructions. No useful d=2,k=5 witness
was found. The fixed star target would imply an exponent about 2.222, stronger
than 9/4; the original theorem therefore does not guarantee this auxiliary.
Varying proof-derived auxiliaries and their mixtures takes priority over
increasing only the fixed-star catalyst dimensions. The new arguments and
Python checker have not been formalized in Lean.
Checkpoint validation: **275 tests passed**; protected sources remain unchanged.

## Next research and implementation tasks

The [geometric reduction](geometric-catalyst.md) shows that the all-fields
exact-rank theorem already implies a complete-in-principle search class over
F_2 with d=2,k=5,D=0: search matrix powers for r<5**b and convert the first
successful scheme. No extension or arbitrary auxiliary support is necessary
for this existential class. This does not provide a practical discovery
method. The direct conversion requires b>=4 by a cited algebra-rank lower
bound, and its first possible dense source already has 6,280,426,125 entries.
The smaller geometric auxiliaries' arbitrary mixed maps are a separate
question; the report distinguishes them from this supplied-scheme bound.

The main research priority is to find a low-cost explicit catalytic
restriction, or a finite rational LP whose order constraints carry actual
restriction witnesses and whose output produces one. A finite input-pair
search with exact linear output recovery is a separate route for testing
small direct decompositions. Neither route currently has a feasible cost
bound for a 9/4-plus-slack witness.

The supplied-certificate-to-coefficients path and a bounded integrated control
pipeline are implemented. The next discovery priority is richer connected
witnessed row families, assessed against direct decomposition search with
explicit budgets and switching criteria. Further construction work includes
symbolic interpolation/Fourier rules and scalar elimination. Formalizing the
central catalytic argument and obstruction arguments in Lean remains a
separate verification task.

The research was divided into main-agent derivation, independent existence
audit, constructive-route review, finite-search experiments, and peripheral
implementation. The detailed reports are durable research notes,
including limitations and failed shortcuts, rather than published proofs.
