# Spectral and field-descent agent

Last updated: 2026-10-06. The task is active; the final all-fields theorem is
not yet verified.

## Ownership

- `lean/OAI/LinearAlgebra/MatrixMultiplication/Arithmetic/FieldDescent.lean`
- `lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Arithmetic/FieldExtension.lean`
- This note. Do not edit other agents' files without coordination.

## Environment and commands

The initial checkout is `work/openai-math`, branch
`codex/matrix-multiplication-all-fields`. Temporary focused compilation runs
from the sibling `work/lean-focus`, with `OAI` symlinked into that checkout.
The setup agent is preparing a canonical private checkout; await its migration
notice and update these paths before continuing there.

Lean is `leanprover/lean4:v4.34.1`, Mathlib is
`d13f23b723b8a846827a245b89c10fc7d3f11612`. The fixed-point dependency has the
upstream compatibility patch; its dirty-repository warning is expected.

```sh
lake build OAI.LinearAlgebra.MatrixMultiplication.Arithmetic.FieldDescent
lake build OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Arithmetic.FieldExtension
```

## Verified results

`Arithmetic.FieldDescent` builds successfully, with no source warnings and no
`sorry` or new axioms. Its main results are:

1. `Tensor.RankAtMost.project_basis`: applying an `F`-linear coefficient
   functional to a tensor over a finite-dimensional commutative `F`-algebra
   increases rank by at most the square of the basis size.
2. `Tensor.exists_unit_linear_functional`: a finite nontrivial commutative
   `F`-algebra has an `F`-linear functional taking its unit to one.
3. `Tensor.RankAtMost.descend_finite`: descent through a fixed finite algebra
   costs at most its dimension squared.
4. `Tensor.RankAtMost.descend_algebraic_powers`: a fixed finite rank scheme
   over an algebraic extension generates one finite coefficient algebra;
   its dimension is chosen before the tensor power varies. For a rank `r`
   scheme of a scalar-extended `T`, it concludes
   `∃ d > 0, ∀ n, RankAtMost (Tensor.power T n) (r^n*d*d)`.
5. `Tensor.RankAtMost.descend_algebraicClosure_powers`: specialization to
   every arbitrary field's own `AlgebraicClosure F`.

The implementation uses `Algebra.adjoin F (Set.range coeff)` and
`Algebra.finite_adjoin_of_finite_of_isIntegral`, and thus does not require a
separability assumption or construction of an intermediate field.

`#print axioms` for `descend_algebraicClosure_powers` reports exactly:

```text
[propext, Classical.choice, Quot.sound]
```

`#check @...descend_algebraicClosure_powers` confirms arbitrary universe
polymorphic `F` with only `[Field F]` and finite coordinate-type assumptions;
no hidden characteristic, closure, or finite-dimension assumptions remain.
The scratch check is `work/lean-focus/CheckDescent.lean`.

## Drafted next result and latest blocker

`AuxiliarySeparation.Arithmetic.FieldExtension` is drafted but has not yet
compiled. It contains:

- `rankAtMost_matrixMultiplication_of_power`: relabel tensor-power coordinates
  to one matrix dimension.
- `exactRankExponent_le_algebraic_extension_logb`: for every fixed block size
  `u ≥ 2`, remove the fixed dimension-squared overhead from all powers to
  prove `ν(F) ≤ log_u(rank_E(T_u))`.
- `exactRankExponent_le_algebraic_extension`: take the infimum of these block
  bounds to prove `ν(F) ≤ ν(E)` for algebraic `E/F`.
- `exactRankExponent_le_algebraicClosure`: specialize to arbitrary `F`.

The latest build stopped in the root-owned dependency
`AuxiliarySeparation/Arithmetic/RankExponent.lean`, before compiling the new
file. The coordinator has accepted responsibility to fix it immediately:

- Lines 106 and 124 use instantiated `(matrixMultiplicationTensor (K := K))`
  as simp lemmas; the bare definition name is required for unfolding.
- Line 221 needs explicit `(K := K)` on `exists_exactMatrixRank_lt_rpow`.

After those fixes, rerun the FieldExtension build and repair its own
elaboration errors. Inspect `#print axioms` for the closure-exponent theorem,
then report the exact successful target to the coordinator.

## Mathematical audit already completed

The published detecting-character appendix is characteristic independent.
Its real compactness, fixed-point, and rational Farkas arguments operate on
numerical states over the tensor semiring. Integers there count direct tensor
summands, so `p` scalar-unit summands have rank `p` even in characteristic `p`.
The catalyst argument retains one auxiliary tensor factor and does not assume
rank cancellation or additive rank. No field limitation was found in that
appendix.

The catalogue and comparator explicitly state the published `9/4` theorem
only over `ℂ`, but no original evaluation prompt or author explanation for
the limited scope was found. Do not present a motive as established fact.
