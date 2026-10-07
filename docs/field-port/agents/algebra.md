# Algebra agent: all-fields matrix multiplication port

## Current status

The finite algebraic constructions have been generalized from `ℂ` to a field
`K`. The characteristic-sensitive Fourier normalization uses a period whose
image in `K` is nonzero. Interpolation uses distinct nonzero elements of an
infinite field rather than natural-number casts.

All original owned algebra and character-dependent targets listed below have passed Lean 4.34.1 against
the repository's pinned Mathlib commit
`d13f23b723b8a846827a245b89c10fc7d3f11612`. The final batch completed successfully
after the generic character/rank interfaces compiled. This document does not
claim the complete all-fields theorem is verified.

Additional integration ownership: `Determinant/Character.lean` and
`Sector/Character.lean` were assigned after that successful batch. The root's
central build subsequently compiled both, including their dependencies.

## Ownership

All paths in this list are relative to
`lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/`.

- `Separation/Fourier.lean`
- `Separation/FiniteProjection.lean`
- `Separation/Basic.lean`
- `Polynomial/Interpolation.lean`
- `Character/Degeneration.lean`
- `Character/FiniteSeparation.lean`
- `Convolution/Basic.lean`
- `Convolution/Rank.lean`
- `Convolution/RankLowerBound.lean`
- `Convolution/Symmetry.lean`
- `Determinant/Filtration.lean`
- `Determinant/Character.lean`
- `Sector/Degeneration.lean`
- `Sector/Branches.lean`
- `Sector/Character.lean`
- `Tensor/TagInequality.lean`

`Determinant/Basis.lean` and `Determinant/Kernel.lean` were inspected and already
have the required generic ring/field parameters; they have not been changed.

## Mathematical changes

For `M` sectors, `separationPeriod K M` is `5*M` when its field cast is nonzero,
and `5*M+1` otherwise. Its cast is always nonzero; when `M>0`, its integer value
lies between `5*M` and `6*M`. Algebraic closedness supplies a primitive root of
this period. The Fourier/square-weight construction otherwise retains exactly
the same terms and has the same degree bound `(M-1)^2`.

Projection APIs now take an explicit period `L` and proofs
`hL : 5*M ≤ L`, `hLcast : (L : K) ≠ 0`, and
`hζ : IsPrimitiveRoot ζ L`. In particular, the second argument to
`projectedBranchPolynomial` now denotes the period itself, rather than a sector
count that was internally multiplied by five. `separationPolynomial ζ L B`
and `projectedTensor ζ L B` also expose this period.

The finite character separation bound changes its constant from `5` to `6`.
The type-counting limit removes that constant, preserving the final entropy
inequality. `Character/FiniteSeparation.lean` and `Tensor/TagInequality.lean`
contain that change. The root agent has also changed the two Character wrappers
in `Entropy/Tag.lean` to constant six.

`interpolationNode` embeds natural indices into the complement of zero in `K`.
The rank-after-powering proof performs coefficient extraction with these nodes
directly, avoiding the upstream `PolynomialApproximation.rank_power` theorem's
`CharZero` premise. Convolution rank upper bounds now require `[Infinite K]`.

The determinant filtration and sector constructions require only `[Field K]`.
The determinant basis changes have unit pivots, so no modular representation
theory or division by dimensions is needed. The scalar field must be supplied
explicitly in closed tensor statements when inference has no other source.

The original topological separation lemmas and the convolution border-rank
corollary remain specialized to `ℂ`, matching their upstream Euclidean-topology
definitions. The generic main argument uses the polynomial/character lemmas.

## Verified targets

These targets have passed individually or in the indicated build batches:

- `Separation.Fourier`
- `Separation.NoWrap` (unchanged dependency)
- `Separation.SquareWeights` (unchanged dependency)
- `Separation.FiniteProjection`
- `Separation.Basic`
- `Polynomial.Interpolation`
- `Convolution.Basic`
- `Convolution.Rank`
- `Convolution.RankLowerBound`
- `Convolution.Symmetry`
- `Character.Degeneration`
- `Character.FiniteSeparation`
- `Determinant.Kernel` (unchanged)
- `Determinant.Basis` (unchanged)
- `Determinant.Filtration`
- `Determinant.Character`
- `Sector.Weights` (unchanged dependency)
- `Sector.Degeneration`
- `Sector.Branches`
- `Sector.Character`
- `Tensor.TagInequality`

The full prefix for each target is
`OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.`.

## Build environment and exact commands

The canonical private source checkout is
`/Users/selanavot/Documents/Codex/2026-10-06/ope/work/matrix-multiplication-all-fields`.
**Run all Lake commands in its `lean/` directory.** The setup agent copied the
focused Lake configuration and installed dependencies there. The old
`openai-math` checkout is a read-only reference; do not edit it.

Earlier component checks used a focused sibling harness at `../lean-focus`
whose `OAI` symlink points into this checkout. Do not run additional builds
there: mixing both package paths can invalidate traces and produce transient
missing-olean errors during concurrent builds. Use canonical `lean/` now and
serialize the final aggregate build.

```sh
lake build OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Separation.Fourier OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Polynomial.Interpolation OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Convolution.Basic OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Determinant.Filtration OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Sector.Degeneration
lake build OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Separation.Basic OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Convolution.Rank
```

The final successful character-dependent batch was:

```sh
lake build OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Character.FiniteSeparation OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Convolution.RankLowerBound OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Convolution.Symmetry OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Sector.Branches OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Tensor.TagInequality
```

## Remaining integration work

No owned target has a remaining Lean error. The newly assigned
`Determinant/Character.lean` and `Sector/Character.lean` now supply explicit
`K` on closed branch/tensor expressions. Determinant's second `Character`
namespace redeclares `{K} [Field K]` because the preceding namespace closed
their scope. Their polynomial-degeneration comparisons require only
`[Infinite K]`; actual binary/three-sector entropy-tag inequalities require
`[IsAlgClosed K]`. Pure support, nonzero, sum, and reindexing lemmas retain
`[Field K]` alone. The root's central build verified both modules; its only
remaining error was an unrelated missing explicit field argument in
`Polynomial/Inequalities.lean`, which the root fixed before building Main.

The root coordinator is integrating
the generic character, profile, spectral, exact-rank, and arithmetic chain.
The complete generic theorem still needs its aggregate build and axiom check.
No `sorry`, new axioms, or admitted proof gaps have been added to these files.

Mechanical caution when repairing later integration calls: use bare definition
names in simp unfold lists. A partially applied definition such as
`(scalarTensor (K := K))` is not a simp theorem. Supply `(K := K)` to closed
tensor theorem applications where no scalar-valued argument determines the
field.

Three pure helper lemmas in `Tensor/TagInequality.lean` unnecessarily inherited
`[IsAlgClosed K]`. Their assumption scope has been reduced, and the canonical
`lean/` command `lake build
OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Tensor.TagInequality`
completed successfully (3663 jobs), with no remaining owned-module warnings.
This canonical build also rebuilt Fourier, FiniteProjection, Basic,
Interpolation, Degeneration, and FiniteSeparation successfully.

Build coordination: do not start another Lake process until the root agent
assigns one. The root has the single aggregate build running. Several parallel
Lake processes caused resource oversubscription and shared-artifact races;
component proof editing can continue separately from compilation.

## Final algebra source review

Reviewed the final owned diff against the preserved baseline, including the
non-mechanical Fourier and interpolation changes. No mathematical or trusted
specification issue was found.

- The normalized root sum uses the explicit premise `(L : K) ≠ 0` everywhere
  it is needed. Increasing the period to `5*M+1` preserves the integer no-wrap
  inequality and the exact square-weight constant coefficient. Polynomial
  exponents and sector equations remain integer/natural data, even when the
  scalar field has characteristic two, three, or five.
- The at-most-`6*M` copy count and later factors two and three are real
  character/dimension counts. They are not field scalars being inverted.
  Exact-type cancellation uses positivity of `(M : ℝ)`, not `(M : K)`.
- Interpolation nodes are genuinely distinct nonzero field elements; their
  existence is guarded by `[Infinite K]`. The powered interpolation proof
  retains overhead `((Lx+Ly+Lz)*n+1)*r^n`; it does not exponentiate a
  single-copy interpolation loss.
- After normalizing scalar-field substitutions, explicit field arguments,
  and whitespace, the determinant filtration, sector branches, and
  convolution definitions have no substantive algebraic changes. The sector
  degeneration gains only the required infinite-field premise on its rank
  extraction corollary. Determinant pivots remain units, and the degeneration
  is a filtration rather than a representation-theoretic direct sum.
- The shared-input tag inequality keeps its original final statement. Only
  its finite integral comparison changes from constant five to six; the
  existing entropy limit removes that fixed factor.
- These files do not redefine `Tensor.RankAtMost`, arithmetic correctness, or
  `Arithmetic.omega`, and contain no new axiom, `sorry`, or `admit` declaration.
  Complex topological comparator lemmas remain isolated specializations.
