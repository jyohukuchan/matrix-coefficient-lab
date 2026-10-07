# Root coordinator handoff

## Active adversarial audit (2026-10-06)

Latest process state: the fresh kernel replay in session 75218 **passed with
exit 0**. Its log contains the correct auxiliary Main module name, and the
isolated checkout still matches `9bd2a64` exactly. The sole active build is
canonical `bash scripts/check-proof.sh`, execution session **47121**, log
`../adversarial-expanded-check.log`, source commit `00595b5`. Auxiliary Main
and public Main have rebuilt; expanded AllFieldsAudit is still running.
After success run `../AdversarialControls.lean` directly (expected exit 1 from
exactly two axiom-list guard mismatches), then finish the README and docs.
Draft final verification record: `../VERIFICATION-public-draft.md`, with
explicit placeholders for still-pending results. A concise attribution script
`../prepare-upstream-notices.py` is prepared but must be applied only after
Lean checks, with an external snapshot and `--verify` byte-preservation check.
The earlier process narrative below is historical; this paragraph is current.

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
`../adversarial-clean-build.log`. At 9377/9437 jobs, memory pressure from nine
concurrent Lean workers required stopping that process (exit 143) and resuming
with `LEAN_NUM_THREADS=2`; continuation log `../adversarial-clean-resume.log`.
After 22 more modules, root stopped the two-worker build (exit 143) to avoid
overlapping large imports and resumed with `LEAN_NUM_THREADS=1` in execution
session 96424, logging to `../adversarial-clean-single.log`. That compatibility
build was then stopped (exit 143, completed modules retained) to prioritize
the already rebuilt core theorem's kernel replay. The **active process** is
execution session 75218:

```sh
LEAN_NUM_THREADS=1 lake env leanchecker --fresh --verbose \
  OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Main
```

It runs in `../adversarial-clean/lean`, logging to `../adversarial-kernel.log`.
Do not start another Lean compiler or replay until it finishes. The focused
auxiliary all-fields Main and its entire OAI import dependency chain have
already compiled successfully there. Root confirmed that the remaining fresh
public-build modules concern other retained bounds, outside that dependency
chain. The audit will finish with the kernel replay plus a normal canonical
public build/expanded AllFieldsAudit and direct scratch controls. Do not claim
a complete fresh public-entry-point build: that broader rebuild was stopped.
This is a deliberate isolated verification
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

Latest user instruction: after verification, rewrite the README for public
readers with OpenAI provenance, premise, the small theorem-statement diff,
technique, and precise statements. The reviewed draft is outside the repository
at `../README-public-draft.md`; do not apply its verification placeholder.
The draft's diff uses the same-name AuxiliarySeparation theorem, and explicitly
retains the public complex specialization. Publication is future intent, not
permission to change repository visibility. The repository stays private.
The reproducibility reviewer is preparing (outside Git, not yet applied)
comment-only modification notices for changed upstream files; preserve their
proof bodies and original license headers when adding attribution.

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
