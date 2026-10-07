# Matrix multiplication over arbitrary fields

This private Lean project extends OpenAI's matrix multiplication proof from
complex scalars to arbitrary fields. The intended conclusion uses the existing
`OAI.MatrixMultiplication.Arithmetic.omega F` definition and is
`omega F ≤ (9 : ℝ) / 4`.

**The complete all-fields theorem is under construction.** Several components
have compiled; this repository does not yet claim a checked final theorem.
See [the current status](docs/field-port/STATUS.md) and
[the private draft PR](https://github.com/selanavot/matrix-multiplication-all-fields/pull/1).

## Source and review baseline

The `main` branch preserves the complete upstream MatrixMultiplication subtree
from `openai/math` commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a`.
The work branch contains the extension. [UPSTREAM.md](UPSTREAM.md) records the
source and the focused package layout. The original Apache 2.0 license is
retained in [LICENSE](LICENSE).

## Build setup

Install [elan](https://github.com/leanprover/elan) if it is not already available,
then run from the repository root:

```sh
bash scripts/bootstrap.sh
```

The script resolves the exact declared dependencies and downloads the matching
Mathlib cache. The toolchain file selects Lean 4.34.1. The package pins Mathlib
to `d13f23b723b8a846827a245b89c10fc7d3f11612` and fixed-point-theorems to
`770940ddf9878cf61952ed53d910b92bca841838`. Its Lake post-update hook applies
the included upstream Lean 4.34.1 compatibility patch. A subsequent
dirty-repository warning for that patched dependency is expected.

Run proof checks from `lean/`. These component targets have compiled:

```sh
lake build OAI.LinearAlgebra.MatrixMultiplication.Arithmetic.FieldDescent
lake build OAI.LinearAlgebra.MatrixMultiplication.Arithmetic.Growth
lake build OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Arithmetic.FieldExtension
lake build OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Separation.Fourier
lake build OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Polynomial.Interpolation
```

Additional verified algebra components and the remaining integration work are
listed in the status document. Prefer one aggregate Lake invocation rather than
several simultaneous builds writing the same outputs. Build caches and local
dependency checkouts are ignored and are not part of the source checkpoint.

The final validation must build the all-fields theorem, retain the complex
specialization, and inspect the theorem's `#print axioms` output. A component
build alone does not establish the intended final conclusion.

## Privacy

This work is authorized for the private `selanavot/matrix-multiplication-all-fields`
repository only. Do not publish it or create an upstream PR without separate
authorization. [AGENTS.md](AGENTS.md) records the working agreements.
