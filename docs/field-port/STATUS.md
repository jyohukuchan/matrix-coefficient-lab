# Field-generalization status

Last updated: 2026-10-06. **The all-fields proof is merged; an additional
adversarial audit is in progress.** Three fresh source reviews found no proof
counterexample. Isolated rebuilding and stronger mechanical checks are still
running; see [ADVERSARIAL.md](ADVERSARIAL.md) and the root handoff note.

## Result

```lean
theorem OAI.MatrixMultiplication.omega_le_nine_quarters
    (F : Type u) [Field F] :
    OAI.MatrixMultiplication.Arithmetic.omega F ≤ (9 : ℝ) / 4
```

The field universe is arbitrary. The original `Arithmetic.omega`, program
correctness and operation-count definitions in `Model.lean` are unchanged.
The explicit epsilon-cost theorem and original complex specialization also
compile. No characteristic, finiteness, separability or algebraic-closedness
assumption is required on F.

## Verification

- Auxiliary Main passed (9053 jobs).
- Public Main, including the retained upstream results, passed (9436 jobs).
- `bash scripts/check-proof.sh` passed (9437 jobs).
- Arbitrary-universe and characteristic-2/3/5 examples passed, as did rational,
  real and complex examples and the original explicit correctness/cost claim.
- Six guarded axiom checks passed with only `propext`, `Classical.choice`,
  and `Quot.sound`. No `sorryAx` or additional mathematical axioms occur.
- Independent algebra, spectral/descent, and arithmetic-specification reviews
  found no weakened definitions or hidden field premises.

See [VERIFICATION.md](VERIFICATION.md) for exact commands and qualifications,
and [REVIEW.md](REVIEW.md) for the proof route. No build process is active.

## Repository and environment

- Canonical checkout: `/Users/selanavot/Documents/Codex/2026-10-06/ope/work/matrix-multiplication-all-fields`
- Private repository: https://github.com/selanavot/matrix-multiplication-all-fields
- Development PR: https://github.com/selanavot/matrix-multiplication-all-fields/pull/1
- Implementation branch: `codex/matrix-multiplication-all-fields`
- Preserved baseline tag: `openai-baseline-adc7f12`, at `d2336fc`.
- Exact public source: `openai/math` commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a`.
- Initial partial checkpoint: `826d7bb`; infrastructure checkpoint: `2e23d75`.
- Checked auxiliary proof checkpoint: `d56ef0a`; guarded audit/review checkpoint: `ae8fdd5`.
- Lean: `leanprover/lean4:v4.34.1`.
- Mathlib: `d13f23b723b8a846827a245b89c10fc7d3f11612`.
- fixed-point-theorems: `770940ddf9878cf61952ed53d910b92bca841838`, with the included upstream compatibility patch.

The user explicitly authorizes private repository creation, pushes, PRs and
merges. Visibility was repeatedly verified as `PRIVATE` / `isPrivate: true`.
No public fork, upstream PR or publication is authorized. PR1 is the
canonical record of the final merge state. The baseline tag remains the
review reference after merging.

For reproduction, run `bash scripts/bootstrap.sh` followed by
`bash scripts/check-proof.sh` from the repository root. The manifest has ten
exact Git revisions and no local-path dependency entries. Local ignored cache
links reuse the downloaded Mathlib cache. The old `../openai-math` clone is a
read-only reference; `../lean-focus` is a retired harness. Use only the
canonical checkout's `lean/` directory for direct Lake commands.

## Proof architecture

1. Parameterize actual finite tensors, characters, the rank exponent, and the
   spectral argument by the scalar field K.
2. In an algebraically closed field, choose Fourier period 5M or 5M+1 so its
   scalar cast is nonzero. The resulting finite overhead is at most 6M.
3. Interpolate at arbitrary distinct nonzero points; no integer-node or
   characteristic-zero assumption is needed.
4. Apply the original determinant/sector and real growth argument to obtain
   the 9/4 exact-rank bound over algebraically closed fields.
5. Descend all powers of one finite decomposition through one fixed finite
   coefficient algebra: rank grows by at most its dimension squared, a
   constant independent of the tensor power.
6. Remove that constant in the exponent and use the original generic program
   model to conclude the bound for each field F via its own algebraic closure.

## Continuation cautions

- Preserve the proven theorem and original model. Future source changes must
  pass `scripts/check-proof.sh`; `sorry` or extra mathematical axioms are
  unacceptable in the final proof dependencies.
- Only one coordinator starts Lake builds. Concurrent builds against shared
  artifacts caused contention and transient missing-olean errors.
- Do not rerun the external initial port scripts; they are not idempotent.
- Some unchanged upstream comparator files intentionally contain placeholders;
  the final theorem's guarded axiom audit confirms it does not depend on them.
- Complex topological degeneration lemmas remain specialized. The all-fields
  argument uses algebraic interpolation instead.
- Human review of the formalization is still useful. Historical novelty,
  author intent, exact omega, and practical constants were not established.

Agent ownership/history is retained in `agents/`. All assigned components
are finished; the root note and this file are authoritative for final status.
