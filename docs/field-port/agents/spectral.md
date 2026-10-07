# Spectral and field-descent agent

Last updated: 2026-10-06. The task is active; the final all-fields theorem is
not yet verified.

## Ownership

- `lean/OAI/LinearAlgebra/MatrixMultiplication/Arithmetic/FieldDescent.lean`
- `lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Arithmetic/FieldExtension.lean`
- `lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Spectrum/Obstruction.lean`
- `lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Character/Existence.lean`
- `lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Arithmetic/CharacterRounding.lean`
- This note. Do not edit other agents' files without coordination.

## Environment and commands

Canonical private checkout is
`/Users/selanavot/Documents/Codex/2026-10-06/ope/work/matrix-multiplication-all-fields`,
branch `codex/matrix-multiplication-all-fields`. Run all lake commands from
its `lean/` directory. The old public clone and `work/lean-focus` are now
read-only references; the old shared-cache symlink invalidated build traces
when different package paths were used concurrently.

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

## Verified field-extension exponent comparison

`AuxiliarySeparation.Arithmetic.FieldExtension` compiled successfully after
migration (3103 jobs, no source warnings). It contains:

- `rankAtMost_matrixMultiplication_of_power`: relabel tensor-power coordinates
  to one matrix dimension.
- `exactRankExponent_le_algebraic_extension_logb`: for every fixed block size
  `u ≥ 2`, remove the fixed dimension-squared overhead from all powers to
  prove `ν(F) ≤ log_u(rank_E(T_u))`.
- `exactRankExponent_le_algebraic_extension`: take the infimum of these block
  bounds to prove `ν(F) ≤ ν(E)` for algebraic `E/F`.
- `exactRankExponent_le_algebraicClosure`: specialize to arbitrary `F`, giving
  `exactRankExponent F ≤ exactRankExponent (AlgebraicClosure F)`.

The root-owned RankExponent elaboration errors were repaired and compiled.
Both owned descent files were rebuilt successfully from the canonical
`lean/` directory after migration (3103 jobs, exit 0). The scratch axioms
check for FieldExtension has not yet completed: attempts in the retired
`lean-focus` harness encountered missing dependency object files due to trace
invalidation. This was not a source/proof failure.

## Verified generic spectral existence

Root's coordinated canonical build confirmed successful compilation of
`Spectrum.Obstruction` and `Character.Existence`, after repairing their
root-owned Scalar dependency. `Arithmetic.CharacterRounding` also compiled
successfully. `exists_detecting_character` now concludes existence of a
`Character K` detecting each integer below `d ^ exactRankExponent K`, with
only `[Field K]`; no infinitude, closedness, or characteristic assumption.

Spectrum/Obstruction has explicit `(K := K)` at rank-lemma applications to
resolve otherwise-stuck scalar metavariables. The local natural-number rank
budget is named `budget`, avoiding collision with scalar field `K`.

Build coordination is now centralized at root: **do not start a new lake build
or lake env lean without coordination**. Multiple early Lake processes raced
over shared artifacts and oversubscribed resources. All of this agent's
sessions have finished. Final signature/axioms checks will be run centrally.

## Independent implementation audit

Read the entire two new descent modules and the spectral/Character/rank
foundations diff against private baseline `d2336fc` (which preserves upstream
`adc7f1241b42e322a6451854ab7e4b4c146bf78a`). No specification weakness or hidden
field hypothesis was found:

- `RankAtMost`, `Tensor.power`, and the existing arithmetic specification are
  unchanged. Generic matrix coefficients retain the identical 0/1 equality
  predicate; the complex tensor is the original specialization.
- Projecting two factors in a rank-one decomposition gives exactly dimension
  squared overhead; `π(e_i*e_j*c(z))` is one base-field third-factor function.
  The proof does not expand all three factors or require trace/separability.
- Every coefficient of the one chosen extension-field scheme lies in one
  finite integral coefficient algebra. This algebra and its basis are fixed
  before the tensor-power parameter, giving `r^j*d*d`, not `(r*d*d)^j`.
- The exponent comparison first bounds the actual base-field rank at `u^j`,
  removes the fixed constant as `j` varies, and only then takes the infimum
  over `u≥2`. Both the logarithm positivity and inequality direction are valid.
- `Character K` retains every original normalization, additivity,
  multiplicativity and actual K-linear restriction-monotonicity axiom.
  `TensorClass K` still quotients mutual actual restriction; no relation has
  been weakened or replaced.
- `rank_natCast` identifies m copies with an m-dimensional diagonal tensor,
  so integer copies do not collapse in characteristic p. Catalyst bounds
  retain the same finite tensor throughout and never assume rank additivity
  or cancellation.
- Real compactness/fixed-point machinery is unchanged and already generic in
  its ordered semiring. It uses numerical real states, not scalar-field
  topology. Finite coordinate presentations use `Fin` while scalar universes
  remain arbitrary.

No `sorry`, `admit`, or axiom declaration occurs in the new descent files.
The final closure comparison and detecting-character axioms check is pending
root's single-process validation; the earlier descent-power check reported
only `propext`, `Classical.choice`, and `Quot.sound`.

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
