#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd -- "$project_root"
bash scripts/check-proof.sh
cd -- lean

# Replay the auxiliary conclusion and all its imported declarations in a fresh
# environment. This uses Lean's own kernel, not an independent implementation.
lake env leanchecker --fresh --verbose \
  OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Main
