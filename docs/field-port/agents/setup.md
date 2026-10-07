# Setup and arithmetic bridge

> **Archived development note.** Written by an AI coding agent on 2026-10-06, during development
> and while the repository was private. Status statements below are historical; current status
> and verification scope are in [STATUS.md](../STATUS.md) and [VERIFICATION.md](../VERIFICATION.md).

Owner: `audit_extension`. Updated 2026-10-06.

## Canonical checkout and privacy

- Editable checkout: `<repository checkout>`.
- Active branch: `codex/matrix-multiplication-all-fields`.
- Private repository: https://github.com/selanavot/matrix-multiplication-all-fields.
- Before the first push, `gh repo view selanavot/matrix-multiplication-all-fields --json nameWithOwner,url,visibility,isPrivate` returned `visibility: PRIVATE` and `isPrivate: true`.
- The preserved tag `openai-baseline-adc7f12` is the focused upstream baseline, commit `d2336fc`. Compare against this tag after merging development PRs into `main`.
- Private draft PR: https://github.com/selanavot/matrix-multiplication-all-fields/pull/1, attached to the app task. Its description explicitly marks the all-fields theorem unverified and the source checkpoint partial.
- Initial partial-source checkpoint pushed: `826d7bbe8b757f98618349b4c3657220b36f57fe`.
- The `upstream` remote points to the public OpenAI repository and its push URL is `DISABLED`. `origin` is the private repository.
- Old `../openai-math` is now a read-only upstream reference. Do not resume edits there.

## Exact environment and commands

- Lean installed via `elan toolchain install leanprover/lean4:v4.34.1`.
- Mathlib commit: `d13f23b723b8a846827a245b89c10fc7d3f11612`.
- fixed-point-theorems commit: `770940ddf9878cf61952ed53d910b92bca841838`.
- Apply the included upstream `lean/patches/fixed-point-theorems-lean4341.patch`; the standalone Lake post-update hook does this for a fresh setup.
- Initial setup used sibling `../lean-focus` with the same exact dependency pins. All 8908 matching Mathlib cache files were downloaded with `lake exe cache get`.
- Its `OAI` symlink now points to this canonical checkout. The canonical checkout's ignored `lean/.lake/packages` and `lean/.lake/build` links share those local caches.
- Portable fresh setup: run `lake update`, then `lake exe cache get` from `lean/`. The full upstream project's unrelated packages are intentionally absent.
- Focused check: from `lean/`, run `lake build OAI.LinearAlgebra.MatrixMultiplication.Arithmetic.Growth`.
- Generic bridge check: `lake build OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Arithmetic.Exponent`.
- Run all new builds from canonical `lean/`. The retired focus harness shares caches, and mixing build roots can invalidate build traces and cause transient missing-olean errors.

## Ownership and verified result

Owned proof files:

- New `lean/OAI/LinearAlgebra/MatrixMultiplication/Arithmetic/Growth.lean`.
- `lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Arithmetic/Exponent.lean`.

`Arithmetic.Growth` compiled successfully with Lean 4.34.1, reporting
`Build completed successfully (8935 jobs)`. It generalizes the upstream
complex recursion, padding, and exact-rank arithmetic construction to every
field using the existing generic arithmetic-program specification.

The canonical `AuxiliarySeparation.Arithmetic.Exponent` target also compiled
successfully, reporting `Build completed successfully (8942 jobs)`. This
checks `omega K ≤ exactRankExponent K` and the epsilon-complexity bridge for
an arbitrary field `K`, using the original arithmetic definition. The generic
`Arithmetic.Growth` target rechecked successfully in that canonical build.
The public entry-point build and final axiom audit were pending at this point; both later passed (see [VERIFICATION.md](../VERIFICATION.md)).

## Setup friction and durable resolution

- The full upstream Lake configuration eagerly downloads many unrelated proof packages. The standalone exact-pins package avoids those dependencies.
- Full-history Mathlib clone/update was unnecessarily expensive; initial setup used a shallow filtered clone and the pinned revision. The committed lockfile avoids silently changing dependency versions.
- `git archive` on the upstream partial clone attempted to hydrate unrelated blobs. The standalone baseline was instead copied from the materialized selected subtree, with modified sources restored from their original `HEAD` blobs. The new private repository contains only that focused history.

Infrastructure documentation checkpoint `2e23d75` is pushed to the private
draft PR. The portable bootstrap script passed `bash -n`; the canonical Lake
configuration compiled and shared caches were audited as ignored. A brand-new
dependency download through the wrapper has not been rerun, because the exact
dependencies and compiler are already installed and checked.

Next: the root coordinator runs one serialized aggregate proof build. Do not
launch concurrent Lake builds that write the same output paths. Final theorem
axioms and the complete complex specialization remain aggregate validation.

## Specification review

A separate source-review pass compared the arithmetic conversion to the
focused upstream baseline. `Model.lean` and the existing primary arithmetic
specification/builder files are unchanged. `SquareAlgorithm` is definitionally
the original square `MatrixAlgorithm`; the auxiliary exponent and admissible
bound are direct aliases of the original definitions. The generic block
recurrence retains all scalar-operation overhead, and padding uses the
unchanged correctness/cost theorem. No specification weakening was found.

The route and review checkpoints are recorded in `docs/field-port/REVIEW.md`.
The root coordinator confirmed `AuxiliarySeparation.Main` compiled (9053 jobs),
including the all-fields exponent and epsilon-cost theorems. At that point the
public entry-point build and final axiom audit were pending; both later passed. This review launched
no Lake builds, preserving the single-build coordination rule.
