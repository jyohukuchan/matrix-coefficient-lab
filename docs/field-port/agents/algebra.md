# Algebra agent: all-fields matrix multiplication port

## Current status

The finite algebraic constructions have been generalized from `ℂ` to a field
`K`. The characteristic-sensitive Fourier normalization uses a period whose
image in `K` is nonzero. Interpolation uses distinct nonzero elements of an
infinite field rather than natural-number casts.

The independent algebra targets listed below have passed Lean 4.34.1 against
the repository's pinned Mathlib commit
`d13f23b723b8a846827a245b89c10fc7d3f11612`. The character-dependent targets are
implemented and await a successful build of the generic character/rank
interfaces. This document does not claim the complete all-fields theorem is
verified.

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
- `Sector/Degeneration.lean`
- `Sector/Branches.lean`
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
contain that change. The two Character wrappers in `Entropy/Tag.lean` also need
constant six; they are owned by the root agent.

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
- `Determinant.Kernel` (unchanged)
- `Determinant.Basis` (unchanged)
- `Determinant.Filtration`
- `Sector.Weights` (unchanged dependency)
- `Sector.Degeneration`

The full prefix for each target is
`OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.`.

## Build environment and exact commands

The original checkout's Lake configuration pulls many unrelated projects.
The setup agent created a focused harness with an `OAI` symlink and the exact
Mathlib/fixed-point dependency pins. The current harness is
`/Users/selanavot/Documents/Codex/2026-10-06/ope/work/lean-focus`.
The canonical private source checkout is
`/Users/selanavot/Documents/Codex/2026-10-06/ope/work/matrix-multiplication-all-fields`.
The harness symlink now points into this checkout. The old `openai-math`
checkout is a read-only reference; do not edit it.
Run builds from that directory, using `lake build`, rather than invoking the
original checkout's full Lake project.

```sh
lake build OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Separation.Fourier OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Polynomial.Interpolation OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Convolution.Basic OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Determinant.Filtration OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Sector.Degeneration
lake build OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Separation.Basic OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Convolution.Rank
```

After changes to upstream character/rank interfaces, continue with:

```sh
lake build OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Character.FiniteSeparation OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Convolution.RankLowerBound OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Convolution.Symmetry OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Sector.Branches OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Tensor.TagInequality
```

## Remaining integration work

The last character-dependent build was blocked before these owned targets ran
by mechanical generic-port errors in root-owned `Character/Basic.lean` and
`Arithmetic/RankExponent.lean`. Their simp lists contained partially applied
definitions such as `(scalarTensor (K := K))`, which Lean rejects as simp
theorems; bare definition names should be used when unfolding. RankExponent
also had one closed generic theorem call needing its explicit field argument.
The root agent is repairing these interfaces.

After those dependencies build, fix any newly exposed errors in the owned
targets, rerun their batch, and finally rebuild the complete generic theorem.
No `sorry`, new axioms, or admitted proof gaps have been added to these files.

The migration to the canonical private checkout is complete. Source ownership
continues unchanged, and the focused harness/build commands remain valid.
