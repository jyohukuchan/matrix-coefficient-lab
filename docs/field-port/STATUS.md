# Field-generalization status

Last updated: 2026-10-06. The complete all-fields theorem in
`AuxiliarySeparation.Main` has compiled. Public-entry-point regression and
final axiom checks are still pending; do not mark the project complete yet.

## Current specification

Prove `OAI.MatrixMultiplication.Arithmetic.omega F ≤ (9 : ℝ) / 4` for
`(F : Type u) [Field F]`, retaining the existing arithmetic definition and
complex specialization. All changes must live in a private GitHub repository
owned by `selanavot`; the user explicitly authorized creation and pushes.

## Baseline and environment

- Upstream: `https://github.com/openai/math`
- Baseline SHA: `adc7f1241b42e322a6451854ab7e4b4c146bf78a`
- Canonical checkout: `/Users/selanavot/Documents/Codex/2026-10-06/ope/work/matrix-multiplication-all-fields`
- Private repository: https://github.com/selanavot/matrix-multiplication-all-fields
- Private draft PR: https://github.com/selanavot/matrix-multiplication-all-fields/pull/1
- Branch: `codex/matrix-multiplication-all-fields`
- Focused upstream baseline on `main`: `d2336fc`.
- Initial partial source checkpoint: `826d7bbe8b757f98618349b4c3657220b36f57fe`.
- Lean: `leanprover/lean4:v4.34.1`
- Mathlib: `d13f23b723b8a846827a245b89c10fc7d3f11612`
- fixed-point-theorems: `770940ddf9878cf61952ed53d910b92bca841838`
- The upstream fixed-point compatibility patch is required.
- The private standalone package contains the complete MatrixMultiplication
  subtree with these exact dependencies. The baseline and work branch are
  separately reviewable in the private draft PR.
- Visibility was verified as `PRIVATE` / `isPrivate: true` before the first
  push and again before PR creation. No public fork or upstream PR was made.

Run all new Lake commands from the canonical checkout's `lean/` directory.
All 8908 matching Mathlib cache files have downloaded. Local ignored cache
links reuse that download. The old `../openai-math` clone is a read-only
reference; `../lean-focus` is a retired build harness. Do not mix build roots,
because their differing paths can invalidate shared build traces.

Example command from canonical `lean/`:

```sh
lake build OAI.LinearAlgebra.MatrixMultiplication.Arithmetic.FieldDescent
```

For a fresh checkout, run `lake update` and `lake exe cache get` from `lean/`.
The focused Lake post-update hook applies the included upstream fixed-point
compatibility patch. Detailed environment and privacy evidence is in
`docs/field-port/agents/setup.md`.

## Verified so far

The following generic targets have compiled with the pinned toolchain:

- `AuxiliarySeparation.Main`: the full chain, including
  `omega_le_nine_quarters (F : Type*) [Field F]` and the explicit uniform
  arithmetic-program cost bound. Canonical build succeeded with 9053 jobs.
- Generic detecting-character existence, determinant and sector character
  inequalities, and the polynomial-profile/exact-rank conclusion.
- `Arithmetic.FieldDescent`: finite-algebra rank projection and fixed-overhead
  descent of all tensor powers from every algebraic extension, including
  `AlgebraicClosure F`. No separability or characteristic restriction.
- `Arithmetic.Growth`: generic exact-rank recursion, padding, and scalar
  arithmetic complexity bounds using the existing program specification.
- `AuxiliarySeparation.Arithmetic.Exponent`: for every field `K`,
  `omega K ≤ exactRankExponent K` and the explicit epsilon-complexity bridge.
- `AuxiliarySeparation.Arithmetic.FieldExtension`: fixed coefficient-algebra
  overhead disappears from tensor-power rank exponents; in particular,
  `exactRankExponent F ≤ exactRankExponent (AlgebraicClosure F)`.
- `AuxiliarySeparation.Separation.Fourier`: invertible period selection and
  generic roots-of-unity averaging.
- `AuxiliarySeparation.Polynomial.Interpolation`: arbitrary distinct nonzero
  nodes over an infinite field.
- `AuxiliarySeparation.Convolution.Basic`.
- `AuxiliarySeparation.Determinant.Filtration`, including its basis/kernel.
- `AuxiliarySeparation.Sector.Degeneration`.
- `AuxiliarySeparation.Arithmetic.Exponent`: the generic rank-to-arithmetic
  exponent bridge, using the unchanged arithmetic model.
- The generic tensor semiring and `Tensor.Characters`, including
  `Character.Basic`, `Character.Dot`, and character permutation/symmetrization.
- `Character.FiniteSeparation`, `Entropy.Tag`, and `Tensor.TagInequality`
  with the constant-six finite bound.
- All algebra-agent owned convolution, determinant-filtration, and sector
  algebra targets; see its note for the complete target list.

The full auxiliary theorem is checked. The public `Main` wrapper and retained
complex/rectangular results are currently building. `AllFieldsAudit.lean`
has been drafted for generic-universe, small-characteristic and axiom checks,
but has not run yet.

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
  Symmetrization, character-to-semiring infrastructure, entropy wrappers,
  profile and rank-bound conclusion, final theorem and review.
- `audit_algebra`: Fourier, FiniteProjection, Separation.Basic,
  Polynomial.Interpolation, Character.Degeneration/FiniteSeparation,
  Tensor.TagInequality, Convolution subtree, Determinant algebra files,
  Sector.Branches/Degeneration, and the determinant/sector character
  inequalities.
- `audit_spectral`: new `Arithmetic.FieldDescent` and new
  `AuxiliarySeparation.Arithmetic.FieldExtension`, plus spectral Obstruction,
  Character.Existence, and Arithmetic.CharacterRounding integration.
- `audit_extension`: environment/private-repository setup, new primary
  `Arithmetic.Growth`, `AuxiliarySeparation.Arithmetic.Exponent`.

## Immediate next steps

1. Root alone coordinates Lake. Public `Main` is building in session 10954.
   Auxiliary Main has succeeded (9053 jobs).
2. Build `AllFieldsAudit`, inspect final theorem axioms and exact generic type,
   and make axiom checks into guarded regression checks.
3. Independent algebra and descent/spectral reviews found no specification
   weakening or hidden field assumptions; arithmetic model review also passed.
4. Agent setup is preserving the baseline with tag `openai-baseline-adc7f12`
   before any private PR merge. The user explicitly authorizes private merges.
5. Save checked proof, final validation and review guide, update/merge private
   PR1 after verification. Preserve both baseline and proof checkpoints.

## Known implementation cautions

- The initial mechanical generalization adds explicit `K` applications to
  closed numerical/tensor statements. Some elaboration errors are expected.
- `Spectrum/Obstruction.lean` used a natural-number local `K`; it was renamed
  `budget` to avoid shadowing the scalar field.
- Complex topological degeneration lemmas remain specialized to `ℂ`. The
  generic proof uses polynomial interpolation and does not depend on them.
- Do not run concurrent Lake build processes against shared artifacts. This
  caused resource contention and transient missing-olean errors even after
  the canonical build root was standardized. Root coordinates aggregate builds.
- Existing upstream comparator files contain intentional `sorry` placeholders;
  they are not proof implementations. The new final theorem must not depend
  on them or introduce new `sorry`/axiom declarations.
