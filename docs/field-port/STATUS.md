# Field-generalization status

Last updated: 2026-10-06. **The all-fields proof and adversarial audit passed.**
The README now explains the OpenAI source fork, premise, theorem-statement diff,
proof changes, exact conclusions, and verification limits. The repository
remains private; a future-publication goal is not authorization to publish.

## Verified result

```lean
theorem OAI.MatrixMultiplication.omega_le_nine_quarters
    (F : Type u) [Field F] :
    OAI.MatrixMultiplication.Arithmetic.omega F ≤ (9 : ℝ) / 4
```

The original arithmetic model is unchanged. The field universe is arbitrary;
no characteristic, finiteness, separability, perfectness, or algebraic-closedness
hypothesis is imposed on F. The explicit positive-epsilon cost theorem and
original public complex specialization also passed.

Over finite fields, the original program predicate is correctness on
field-valued inputs, not formal polynomial equality. The separate exact
coefficient-rank theorem supplies the algebraic certificate used by this proof.
The expanded audit checks that statement and records the distinction on F₂.

## Verification completed

- Core proof `9bd2a64fa47efe678664b5401b6a8ea6d93238ad`: auxiliary all-fields
  Main and its full OAI import chain rebuilt from committed source in an
  isolated checkout with initially empty OAI build outputs.
- Fresh `leanchecker --fresh` replay of the auxiliary Main and every imported
  declaration passed with exit 0. The isolated source tree remained unchanged.
- Expanded audit `704431ea67e20489dba45b47c70813a6075c554c`:
  `scripts/check-proof.sh` passed with exit 0. It includes the public entry
  point, arbitrary-universe and representative-field examples, explicit
  program cost, exact coefficient-rank witness, and ten axiom guards.
- Every guarded declaration uses exactly `propext`, `Classical.choice`, and
  `Quot.sound`. Scratch controls with an added axiom and `sorry` were both
  rejected by their guards; the compiler exited 1 with exactly those errors.
- Eleven original specification/builder files match the immutable baseline.
  All ten dependency source pins and the exact compatibility patch passed.
- Three fresh adversarial source reviews found no fatal defect. A stale
  five/six comment was corrected. Two new test-proof elaboration errors were
  repaired without changing their statements or the main theorem.
- Subsequent notices on 41 modified upstream files change only comments:
  removing those exact headers restores every pre-notice source byte.

The isolated fresh build covers the core dependency chain, not the entire
public entry point: the broader build was stopped before unrelated retained
results finished. The canonical public build was checked separately using
its cache. Third-party pinned caches were reused. The replay uses Lean's own
kernel, not an independent kernel implementation. No verification process
remains active. See [VERIFICATION.md](VERIFICATION.md) and
[ADVERSARIAL.md](ADVERSARIAL.md) for evidence and coverage limits.

## Stable references and continuation

- Repository: https://github.com/selanavot/matrix-multiplication-all-fields
- Original proof integration: [PR 1](https://github.com/selanavot/matrix-multiplication-all-fields/pull/1),
  tag `all-fields-proof-v1` at `9bd2a64`.
- Audit and README integration: [PR 2](https://github.com/selanavot/matrix-multiplication-all-fields/pull/2).
- Immutable baseline: `d2336fc571f1f8cdabf0c6d3a2d3ef1ec3327653`, tag
  `openai-baseline-adc7f12`, preserving OpenAI source commit
  `adc7f1241b42e322a6451854ab7e4b4c146bf78a`.
- Lean: `leanprover/lean4:v4.34.1`; dependencies are pinned in the manifest.

For new proof changes, run `scripts/check-proof.sh`; `scripts/check-kernel.sh`
adds fresh kernel replay. Run those scripts from the repository root and
direct Lake commands from `lean/`. Use one worker by default and only one
coordinator for compiler processes. Do not rerun the
non-idempotent initial port scripts or use the retired focus harness.

Preserve the baseline and original model. Do not change visibility, contact
upstream, or publish without explicit authorization. Historical novelty,
author intent, exact omega, practical constants, and an effective algorithm
generator remain outside the established result.
