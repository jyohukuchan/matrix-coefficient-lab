# Root coordinator handoff

## Current result — 2026-10-06

The all-fields proof, three adversarial source reviews, mechanical audit, and
requested public-reader README are complete. There are no active Lean/Lake
processes or unfinished mathematical obligations from this audit. The
repository remains **private**. User intent to make it public soon does not
authorize changing visibility. Private pushes, PRs, and merges are authorized.

Read `../STATUS.md`, `../VERIFICATION.md`, `../ADVERSARIAL.md`, and repository
`AGENTS.md` before continuing. PR2 is the audit/README integration record:
https://github.com/selanavot/matrix-multiplication-all-fields/pull/2.
Inspect its actual GitHub state and local Git status before any new action;
do not restart completed verification merely because a context reset occurred.

## Exact verification evidence

- Merged proof source: `9bd2a64fa47efe678664b5401b6a8ea6d93238ad`, tag
  `all-fields-proof-v1`; original private PR1.
- Isolated checkout at `../adversarial-clean`, with no initial OAI outputs
  and only exact-pinned dependency checkouts/caches shared. Source unchanged
  before/after. Auxiliary all-fields Main and its entire OAI dependency chain
  rebuilt. The broader public fresh build was intentionally stopped; do not
  report it as completed.
- Fresh kernel replay on that core passed, session 75218 exit 0. It used
  Lean's own `leanchecker --fresh`, not an independent kernel implementation.
- Expanded canonical audit passed at
  `704431ea67e20489dba45b47c70813a6075c554c`, session 69489 exit 0. Public Main
  and AllFieldsAudit compiled; ten standard-axiom-only guards passed.
- Combined scratch controls in `../AdversarialControls.lean` printed the
  expected fully elaborated types and original definitions. Session 48326
  exited 1 with exactly the two intended guard errors: `adversarialInjected`
  and `sorryAx`. Neither scratch declaration is in the repository proof.
- Two test-proof errors in the first expanded run were repaired (addition
  lemma resolution and explicit `Polynomial.coeff_X` simplification). The
  main theorem and test statements did not change.
- After Lean checks, 41 changed upstream files received short attribution
  comments. `../prepare-upstream-notices.py --verify` with external
  `../upstream-notices-snapshot.json` confirmed that removing exactly those
  headers restores every pre-insertion source byte. Do not conflate this
  check with a second full compilation of the comment-only revisions.

Compact evidence is committed in `../adversarial/mechanical-checks.txt`.
Large local logs are outside Git: `../adversarial-kernel.log`,
`../adversarial-expanded-check.log`, and `../adversarial-controls.log`.
The source-review reports are committed under `../adversarial/`. Their
machine-specific path prefixes were shortened; findings were preserved.

## Interpretation and continuation cautions

The original model is unchanged. Its finite-field program correctness is
functional, while the separately verified exact coefficient-rank theorem
supplies the stronger algebraic certificate. See README and verification
record; do not attribute a symbolic-evaluation predicate to the original
program type or claim an equivalence of the two exponent definitions.

One fixed coefficient algebra is chosen before all tensor powers, so its
square-dimension overhead disappears in the exponent. Do not recursively
reuse a single descended block and charge that overhead at every level.

Only root starts compiler/replay processes. Scripts default to one worker
because several whole-Mathlib imports exhausted available memory during the
fresh audit. Use the canonical repository's `lean/`; the old focus harness
is retired. The initial external port scripts are not idempotent.

No fatal source-review finding required a proof change. The concrete source
correction was a stale five/six comment. Historical priority and author intent
are unestablished. Earlier explanatory drafts outside this repository predate
some verification; the checked Lean sources and current records take priority.
