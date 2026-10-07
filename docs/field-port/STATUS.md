# Field-generalization status

Last updated: 2026-10-06. **The all-fields proof is merged; the additional
adversarial audit is still running.** Three fresh source reviews found no
proof counterexample. The core all-fields theorem and its complete OAI import
chain have now rebuilt in an isolated checkout. Fresh kernel replay is active;
the expanded public audit and scratch negative controls remain pending.
See [ADVERSARIAL.md](ADVERSARIAL.md) and [agents/root.md](agents/root.md).

## Result already proved in the initial audit

```lean
theorem OAI.MatrixMultiplication.omega_le_nine_quarters
    (F : Type u) [Field F] :
    OAI.MatrixMultiplication.Arithmetic.omega F ≤ (9 : ℝ) / 4
```

The original arithmetic model is unchanged. The field universe is arbitrary;
no characteristic, finiteness, separability, perfectness, or algebraic-closedness
hypothesis is imposed on F. The explicit epsilon-cost theorem and original
public complex specialization passed the initial public build and six guarded
axiom checks. That build was incremental and reused shared OAI artifacts;
its Lake job totals must not be described as counts of fresh compilations.

The new review identified a semantic distinction worth making explicit:
`MatrixAlgorithm.Correct` means equality on field-valued inputs, which is
weaker than polynomial identity over finite fields. The separate exact
coefficient-rank theorem supplies the algebraic certificate used in this proof.
The extended audit adds that certificate, an infinite characteristic-two
field example, and four further axiom guards. These new checks have not yet run.

## Current verification scope

- Isolated checkout: detached at proof commit
  `9bd2a64fa47efe678664b5401b6a8ea6d93238ad`, initially no OAI build outputs.
- Its auxiliary all-fields Main and full OAI dependency chain rebuilt
  successfully. Only exact-pinned third-party checkouts/caches are shared.
- The broader fresh public-entry-point build was stopped after the relevant
  core compiled; remaining modules concern other retained upstream bounds.
  Do not claim a completed fresh build of the whole public entry point.
- The bundled `leanchecker --fresh` is replaying the auxiliary Main and
  all imported declarations using Lean's own kernel, not a second kernel
  implementation. Only root starts compiler/replay processes.
- All ten installed dependency source revisions and the exact allowed
  compatibility patch pass the new source validator.
- After replay: run canonical `scripts/check-proof.sh`, then the combined
  scratch type checks and deliberately failing axiom/sorry guard controls.

## Repository and next deliverable

- Repository: https://github.com/selanavot/matrix-multiplication-all-fields
- Visibility remains **private**; no publication is authorized yet.
- Merged proof: [PR 1](https://github.com/selanavot/matrix-multiplication-all-fields/pull/1),
  commit `9bd2a64fa47efe678664b5401b6a8ea6d93238ad`, tag `all-fields-proof-v1`.
- Active audit: [draft PR 2](https://github.com/selanavot/matrix-multiplication-all-fields/pull/2),
  branch `codex/adversarial-audit`.
- Immutable baseline: `d2336fc571f1f8cdabf0c6d3a2d3ef1ec3327653`, tag
  `openai-baseline-adc7f12`, preserving OpenAI source commit
  `adc7f1241b42e322a6451854ab7e4b4c146bf78a`.
- Lean: `leanprover/lean4:v4.34.1`; dependencies are pinned in the manifest.

The user requests a public-reader README after verification: provenance,
premise, small theorem-statement diff, technique, and exact statements.
A reviewed draft is prepared outside Git, with an explicit pending-verification
placeholder. Apply it only after the remaining checks finish. Do not convert
future publication intent into permission to change repository visibility.

Read [REVIEW.md](REVIEW.md) for the proof route and
[VERIFICATION.md](VERIFICATION.md) for the initial verification record.
Preserve the model, proof and baseline. Do not rerun the non-idempotent initial
port scripts. Use one verification worker on this host; see `AGENTS.md`.
Historical novelty, author intent, exact omega, practical constants, and an
effective algorithm generator have not been established.
