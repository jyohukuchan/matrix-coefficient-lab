# Upstream baseline

The `lean/OAI/LinearAlgebra/MatrixMultiplication` sources are extracted from
OpenAI's public mathematics repository at commit
`adc7f1241b42e322a6451854ab7e4b4c146bf78a`:
https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a

The preserved tag `openai-baseline-adc7f12` points to baseline commit
`d2336fc`. It preserves the upstream matrix multiplication subtree and must
not be moved. Use this tag for comparisons after development PRs are merged
into `main`.
The standalone Lake package keeps the same Lean 4.34.1 toolchain, Mathlib
commit, fixed-point dependency, and upstream compatibility patch. Other
mathematical projects and dependencies are omitted from this focused package.

The working branch contains the extension to arbitrary fields. Its core
aggregate theorem has compiled; see the progress documentation for the public
entry-point build and final axiom-audit status.
