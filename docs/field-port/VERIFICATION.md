# Verification record

Date: 2026-10-06. Lean: `leanprover/lean4:v4.34.1`.

## Checked conclusion

The public declaration elaborates to:

```lean
OAI.MatrixMultiplication.omega_le_nine_quarters :
  ∀ (F : Type u_1) [inst : Field F],
    OAI.MatrixMultiplication.Arithmetic.omega F ≤ 9 / 4
```

`u_1` is arbitrary, and `[Field F]` is its only hypothesis.
The original `Model.lean` is byte-for-byte unchanged from baseline tag
`openai-baseline-adc7f12`, as enforced by `scripts/check-proof.sh`.

The explicit cost theorem, `AuxiliarySeparation.matrix_multiplication_cost_le`,
also proves that for every positive real epsilon there is one positive real
constant C such that, at every positive matrix size n, a correct original
`Arithmetic.MatrixAlgorithm F n n n` has original cost at most
`C * n ^ (9/4 + epsilon)`.

## Commands and outcomes

All Lake commands ran from the canonical repository's `lean/` directory.
The public checks did not alter upstream rectangular or dual-exponent results.

| Command | Result |
| --- | --- |
| `lake build OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Main` | Passed, 9053 jobs. |
| `lake build OAI.LinearAlgebra.MatrixMultiplication.Main` | Passed, 9436 jobs. |
| `bash scripts/check-proof.sh` (from repository root) | Passed, 9437 jobs, including the permanent audit. |
| `git diff --check` | Passed. |

`AllFieldsAudit.lean` checks a general `F : Type u`, fields `ZMod 2`,
`ZMod 3`, `ZMod 5`, and the rationals, reals, and complexes, together with
the explicit cost theorem. The characteristic-five example supplies a local
`Fact (Nat.Prime 5)` instance proved by `norm_num`. Its initial omission was
an audit-example setup error, repaired before the final successful run.

## Axiom audit

The following six declarations have guarded `#print axioms` checks. Every
one reports exactly `[propext, Classical.choice, Quot.sound]`:

- `OAI.MatrixMultiplication.omega_le_nine_quarters`
- `OAI.MatrixMultiplication.complex_omega_le_nine_quarters`
- `OAI.MatrixMultiplication.AuxiliarySeparation.matrix_multiplication_cost_le`
- `OAI.MatrixMultiplication.AuxiliarySeparation.exactRankExponent_le_algebraicClosure`
- `OAI.MatrixMultiplication.AuxiliarySeparation.exists_detecting_character`
- `OAI.MatrixMultiplication.Foundation.Tensor.RankAtMost.descend_algebraicClosure_powers`

The guards make additional axioms, including `sorryAx`, a build failure.
An earlier independent read-only axiom audit of the auxiliary conclusion and
bridges returned the same standard axioms.

## Source review

Independent agent reviews covered the field-sensitive algebra, spectral
argument and algebraic descent, and the arithmetic specification/bridge.
No specification weakening or hidden characteristic/separability premise
was found. The existing primary arithmetic definitions and builders are
unchanged; only new `Arithmetic.FieldDescent` and `Arithmetic.Growth` files
were added in that subtree.

The dependency manifest contains ten exact Git revisions and no local-path
dependencies. Local builds used the matching downloaded Mathlib cache and
the recorded upstream fixed-point compatibility patch. A separate fresh
checkout build was not run. See the bootstrap script for reproduction.

## Interpretation

This verifies the all-fields extension in the repository's formal model.
It supplies an asymptotic upper bound, without determining the exact exponent
or useful numerical constants/crossover sizes. Historical novelty and the
authors' reasons for their original complex-only statement were not established.
