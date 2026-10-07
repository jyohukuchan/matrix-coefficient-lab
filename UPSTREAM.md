# Upstream baseline

The `lean/OAI/LinearAlgebra/MatrixMultiplication` sources are extracted from
OpenAI's public mathematics repository at commit
`adc7f1241b42e322a6451854ab7e4b4c146bf78a`:
https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a

The baseline commit preserves the upstream matrix multiplication subtree.
The standalone Lake package keeps the same Lean 4.34.1 toolchain, Mathlib
commit, fixed-point dependency, and upstream compatibility patch. Other
mathematical projects and dependencies are omitted from this focused package.

The working branch contains an unfinished extension to arbitrary fields.
See the working progress documentation for validation status.
