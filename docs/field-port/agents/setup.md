# Setup and arithmetic bridge

Owner: `audit_extension`. Updated 2026-10-06.

## Canonical checkout and privacy

- Editable checkout: `/Users/selanavot/Documents/Codex/2026-10-06/ope/work/matrix-multiplication-all-fields`.
- Active branch: `codex/matrix-multiplication-all-fields`.
- Private repository: https://github.com/selanavot/matrix-multiplication-all-fields.
- Before the first push, `gh repo view selanavot/matrix-multiplication-all-fields --json nameWithOwner,url,visibility,isPrivate` returned `visibility: PRIVATE` and `isPrivate: true`.
- `main` is the focused upstream baseline, commit `d2336fc`; the work branch contains the extension.
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

The auxiliary exact-rank-to-arithmetic bridge is drafted but its target has
not yet been checked. It depends on the root coordinator's generalized
rank definitions. The final all-fields theorem remains unfinished.

## Setup friction and durable resolution

- The full upstream Lake configuration eagerly downloads many unrelated proof packages. The standalone exact-pins package avoids those dependencies.
- Full-history Mathlib clone/update was unnecessarily expensive; initial setup used a shallow filtered clone and the pinned revision. The committed lockfile avoids silently changing dependency versions.
- `git archive` on the upstream partial clone attempted to hydrate unrelated blobs. The standalone baseline was instead copied from the materialized selected subtree, with modified sources restored from their original `HEAD` blobs. The new private repository contains only that focused history.

Next: compile and repair the auxiliary exponent bridge, update this note with
the exact checked target, and preserve progress in the private work branch.
