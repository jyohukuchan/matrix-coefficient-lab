# Root coordinator handoff

The requested proof is complete. Read `../STATUS.md`, `../VERIFICATION.md`,
`../REVIEW.md`, and repository `AGENTS.md` for the final result and workflow.
No Lean processes or mathematical proof obligations remain active.

The final public theorem is
`OAI.MatrixMultiplication.omega_le_nine_quarters (F : Type u) [Field F]`,
using unchanged `Arithmetic.omega F`. Public Main passed with 9436 jobs, and
`bash scripts/check-proof.sh` passed with 9437 jobs. Six guarded axiom checks
accept only `propext`, `Classical.choice`, and `Quot.sound`.

The permanent audit checks arbitrary universes, fields of characteristics
2/3/5, rationals/reals/complexes, and the explicit original program-cost
statement. The initial characteristic-five example lacked a local
`Fact (Nat.Prime 5)` instance; adding its `norm_num` proof completed the audit.
No change to the universal theorem was needed.

Private PR1 is the authoritative workflow/merge record. The user explicitly
authorizes merging private development PRs. Preserve the upstream baseline
at tag `openai-baseline-adc7f12` when merging. Do not publish anything.
The repository has exact dependency pins, bootstrap and check scripts,
review documentation, and individual agent notes for continuation.

Only one coordinator may start aggregate Lake builds. Shared-source editing
can be parallel, but concurrent builds caused memory contention and transient
missing artifacts. Always build in the canonical private checkout's `lean/`.
The old focus harness is retired. External initial port scripts are not
idempotent and must not be rerun.

The explanatory LaTeX note at the outer workspace's
`outputs/matrix-multiplication-all-fields.tex` predates the completed Lean
verification. The checked Lean files and final verification record are the
source of truth for the finished result.
