# Root coordinator handoff

## Active adversarial audit (2026-10-06)

The user requested three fresh adversarial reviewers. Agents `adversarial_spec`,
`adversarial_proof`, and `adversarial_repro` are auditing merged proof commit
`9bd2a64fa47efe678664b5401b6a8ea6d93238ad` read-only. Their owned reports are
outside Git at `../adversarial-{spec,proof,repro}-report.md`. The active branch
is `codex/adversarial-audit`; preserve the completed proof while investigating.
Private draft PR2: https://github.com/selanavot/matrix-multiplication-all-fields/pull/2.
All three source reviews are complete and copied to `../adversarial/`.
See `../ADVERSARIAL.md` for findings and pending mechanical checks.

Root started a clean OAI build in `../adversarial-clean/lean`, a local clone
detached at that exact commit, with an initially empty project build directory.
Only the exact-pinned dependency checkouts/caches are shared. Its log is
`../adversarial-clean-build.log`. This is a deliberate isolated verification
exception to the normal canonical-directory rule. Only root starts builds.
After that build, run bundled `leanchecker --fresh` on AuxiliarySeparation.Main
to replay all imported declarations in a fresh Lean kernel environment.
It is not an independent kernel implementation.

Preliminary findings: no proof counterexample; Model.Correct is extensional
field-valued correctness and differs from formal polynomial equality over
finite fields. The actual proof uses exact coefficient tensor rank, which
appears to supply the stronger guarantee. Investigate and document precisely.
The original build used shared OAI artifacts, so the isolated build is needed.
No final mechanical-audit verdict has been issued yet. Expanded audit examples
and ten axiom guards are drafted but not yet compiled; do not call them passed.
All ten dependency source revisions and the exact compatibility delta pass
the new `scripts/verify-dependencies.py` in both canonical and isolated trees.

## Completed proof checkpoint (before this audit)

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
