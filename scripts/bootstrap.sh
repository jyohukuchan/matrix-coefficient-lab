#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd -- "$project_root/lean"

if ! command -v lake >/dev/null 2>&1; then
  printf '%s\n' 'Install elan first: https://github.com/leanprover/elan' >&2
  exit 1
fi

# The committed configuration pins every direct dependency. Its post-update
# hook applies the included upstream fixed-point compatibility patch.
lake update
lake exe cache get

printf '%s\n' 'Dependencies and Mathlib cache prepared. Run focused proof checks from lean/.'
