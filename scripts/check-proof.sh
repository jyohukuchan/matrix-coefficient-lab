#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd -- "$project_root"

# The theorem must use the upstream correctness, operation-count, and omega
# definitions verbatim. Fail if that trusted specification has changed.
git diff --exit-code openai-baseline-adc7f12 -- \
  lean/OAI/LinearAlgebra/MatrixMultiplication/Model.lean

cd -- lean
# This target imports the public entry point, preserves the complex theorem,
# checks arbitrary universes and representative fields, and audits axioms.
lake build OAI.LinearAlgebra.MatrixMultiplication.AllFieldsAudit
