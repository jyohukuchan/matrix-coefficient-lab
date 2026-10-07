#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd -- "$project_root"

# Keep the original specification and supporting arithmetic/tensor definitions.
# Use an immutable commit rather than a movable tag for the trusted baseline.
baseline=d2336fc571f1f8cdabf0c6d3a2d3ef1ec3327653
spec_root=lean/OAI/LinearAlgebra/MatrixMultiplication
git diff --exit-code "$baseline" -- \
  "$spec_root/Model.lean" \
  "$spec_root/Arithmetic/Complexity.lean" \
  "$spec_root/Arithmetic/Programs.lean" \
  "$spec_root/Arithmetic/RecursiveBlockPrograms.lean" \
  "$spec_root/Arithmetic/Padding.lean" \
  "$spec_root/Arithmetic/LowerBound.lean" \
  "$spec_root/Arithmetic/Exponent.lean" \
  "$spec_root/Tensor/ComplexTensor.lean" \
  "$spec_root/Tensor/ComplexTensorFlattening.lean" \
  "$spec_root/Tensor/ComplexMatrixTensor.lean" \
  "$spec_root/Polynomial/ExpressionFamily.lean"

python3 scripts/verify-dependencies.py "$project_root"

cd -- lean
# This target imports the public entry point, preserves the complex theorem,
# checks arbitrary universes and representative fields, and audits axioms.
lake build OAI.LinearAlgebra.MatrixMultiplication.AllFieldsAudit
