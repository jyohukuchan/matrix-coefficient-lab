# Field-generalization status

Last updated: 2026-10-06. This is an active implementation, not a completed
all-fields Lean theorem.

## Current specification

Prove `OAI.MatrixMultiplication.Arithmetic.omega F ≤ (9 : ℝ) / 4` for
`(F : Type u) [Field F]`, retaining the existing arithmetic definition and
complex specialization. All changes must live in a private GitHub repository
owned by `selanavot`; the user explicitly authorized creation and pushes.

## Baseline and environment

- Upstream: `https://github.com/openai/math`
- Baseline SHA: `adc7f1241b42e322a6451854ab7e4b4c146bf78a`
- Initial branch: `codex/matrix-multiplication-all-fields`
- Lean: `leanprover/lean4:v4.34.1`
- Mathlib: `d13f23b723b8a846827a245b89c10fc7d3f11612`
- fixed-point-theorems: `770940ddf9878cf61952ed53d910b92bca841838`
- The upstream fixed-point compatibility patch is required.
- Setup agent is preparing a private standalone snapshot of the complete
  MatrixMultiplication subtree with these exact dependencies. This avoids
  unrelated manuscript and proof-package downloads. Baseline and changes
  will remain separately reviewable.

Temporary local compilation currently uses sibling `../lean-focus`, whose
`OAI` symlink points into the working checkout's `lean/OAI`. All 8908 matching
Mathlib cache files have downloaded. Example command from that harness:

```sh
lake build OAI.LinearAlgebra.MatrixMultiplication.Arithmetic.FieldDescent
```

The setup agent will record the canonical private checkout and portable build
commands when the migration is complete.

## Verified so far

The following generic targets have compiled with the pinned toolchain:

- `Arithmetic.FieldDescent`: finite-algebra rank projection and fixed-overhead
  descent of all tensor powers from every algebraic extension, including
  `AlgebraicClosure F`. No separability or characteristic restriction.
- `AuxiliarySeparation.Separation.Fourier`: invertible period selection and
  generic roots-of-unity averaging.
- `AuxiliarySeparation.Polynomial.Interpolation`: arbitrary distinct nonzero
  nodes over an infinite field.
- `AuxiliarySeparation.Convolution.Basic`.
- `AuxiliarySeparation.Determinant.Filtration`, including its basis/kernel.
- `AuxiliarySeparation.Sector.Degeneration`.

These are component checks. The final all-fields theorem has not yet compiled.

## Source interfaces

- New `AuxiliarySeparation.matrixMultiplicationTensor (K := K) a b c`
  has the original coefficient definition over generic `K`; at `ℂ` it is
  definitionally the upstream tensor. The upstream complex-only constant
  remains unchanged to avoid invalidating unrelated modules.
- `exactRank` infers the field from its tensor.
- `exactMatrixRank K n`, `exactRankExponentSet K`, `exactRankExponent K`.
- `TensorSemiring.FiniteTensor K`, `TensorSemiring.TensorClass K`.
- `Character K` with the original character axioms, now over `K`.
- `omega K` and `ArithmeticBound K` abbreviate the unchanged primary
  `Arithmetic.omega K` and `Arithmetic.AdmissibleExponent K`.

## Ownership

- Root coordinator: tensor semiring, Character.Basic/Dot/Permutation/
  Symmetrization/Existence, character-to-semiring infrastructure, spectral
  obstruction specialization, entropy wrappers, determinant/sector character
  inequalities, profile and rank-bound conclusion, final theorem and review.
- `audit_algebra`: Fourier, FiniteProjection, Separation.Basic,
  Polynomial.Interpolation, Character.Degeneration/FiniteSeparation,
  Tensor.TagInequality, Convolution subtree, Determinant algebra files,
  Sector.Branches/Degeneration.
- `audit_spectral`: new `Arithmetic.FieldDescent` and new
  `AuxiliarySeparation.Arithmetic.FieldExtension`.
- `audit_extension`: environment/private-repository setup, new primary
  `Arithmetic.Growth`, `AuxiliarySeparation.Arithmetic.Exponent`.

## Immediate next steps

1. Finish private repository setup; coordinate a pause before changing the
   canonical checkout and preserve every agent's current source.
2. Compile the generic rank and character foundations; repair elaboration
   errors before compiling their dependents.
3. Change the two character entropy wrappers in `Entropy/Tag.lean` from
   fixed constant 5 to 6 (their underlying real limit lemmas already allow
   arbitrary positive constants).
4. Propagate `[Infinite K]` and `[IsAlgClosed K]` only to results that need
   interpolation and roots of unity.
5. Compile spectral existence and the profile chain; derive the closed-field
   rank theorem.
6. Compile algebraic-extension rank descent and combine with arithmetic
   conversion for arbitrary fields.
7. Build the original complex theorem, inspect final theorem axioms, document
   exact validation commands, commit/push, and create a private PR.

## Known implementation cautions

- The initial mechanical generalization adds explicit `K` applications to
  closed numerical/tensor statements. Some elaboration errors are expected.
- `Spectrum/Obstruction.lean` used a natural-number local `K`; it was renamed
  `budget` to avoid shadowing the scalar field.
- Complex topological degeneration lemmas remain specialized to `ℂ`. The
  generic proof uses polynomial interpolation and does not depend on them.
- Existing upstream comparator files contain intentional `sorry` placeholders;
  they are not proof implementations. The new final theorem must not depend
  on them or introduce new `sorry`/axiom declarations.
