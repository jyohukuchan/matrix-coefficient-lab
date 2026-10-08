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
extension-degree bound. A fair search must vary both the extension and the
tensor/map dimensions.

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
decomposition once an explicit catalyst is supplied. The open practical
problem is finding a useful certificate and compiling these map compositions.
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

A synthetic arithmetic experiment in the finite-search report uses
d=4,k=23,R_S=2,C_D=1,e=2, inner exponent 23/10 and outer exponent 231/100.
It finds b=14 and the first successful final inequality at t=50, with matrix
size bit length 1401 and descended term-bound bit length 3234. No tensors or
maps with those costs are supplied. These numbers illustrate how large the
conditional construction can become; they do not establish a new algorithm.

Run the implemented connection with:

```sh
python3 -m reference.extraction_demo
python3 -m unittest reference.test_extraction reference.test_catalysts -v
```

Checkpoint validation (2026-10-08): the full reference suite passed all 112
tests, the extraction demo passed, and the protected-specification comparison
passed. Python compilation and whitespace checks passed. No Lean sources were
changed or rebuilt for this Python/research checkpoint.

## Next research and implementation tasks

The main research priority is to find a low-cost explicit catalytic
restriction, or a finite rational LP whose order constraints carry actual
restriction witnesses and whose output produces one. A finite input-pair
search with exact linear output recovery is a separate route for testing
small direct decompositions. Neither route currently has a feasible cost
bound for a 9/4-plus-slack witness.

The implementation task is a coefficient compiler for the powered-catalyst
and absorption recurrence, initially exercised on the existing d=2,k=9
certificate at a looser exponent. It must produce the arrays, verify all
coefficients independently, retain one coefficient field throughout, and
distinguish tested finite examples from the convergence derivation. Dense
allocation limits will matter well before the synthetic example's size.
Formalizing the recurrence in Lean is a separate verification task.

The research was divided into main-agent derivation, independent existence
audit, constructive-route review, finite-search experiments, and peripheral
implementation. The three detailed reports are durable research notes,
including limitations and failed shortcuts, rather than published proofs.
