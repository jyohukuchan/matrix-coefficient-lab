# Matrix multiplication over arbitrary fields

This private Lean project proves **`Arithmetic.omega F ≤ 9/4` for every field
`F`**, extending OpenAI's complex-scalar result with the original arithmetic
complexity definition unchanged.

The public theorem in namespace `OAI.MatrixMultiplication` is:

```lean
theorem omega_le_nine_quarters (F : Type*) [Field F] :
    Arithmetic.omega F ≤ (9 : ℝ) / 4
```

The field universe is arbitrary. There is no characteristic, finiteness,
algebraic-closedness, or separability premise on `F`. The original
`complex_omega_le_nine_quarters` theorem remains as a specialization.

The full public entry point and permanent audit passed on 2026-10-06 with
Lean 4.34.1. The final audit reports 9437 successful build jobs and checks
six declarations against exactly `propext`, `Classical.choice`, and
`Quot.sound`; no `sorryAx` or extra mathematical axiom occurs. See the
[verification record](docs/field-port/VERIFICATION.md).

## Review

Start with [the private PR](https://github.com/selanavot/matrix-multiplication-all-fields/pull/1)
and [the reviewer guide](docs/field-port/REVIEW.md). The proof changes the
Fourier period and interpolation nodes, then descends tensor powers through
one fixed finite coefficient algebra. Its constant overhead disappears from
the exponent.

The exact upstream MatrixMultiplication subtree is preserved by tag
`openai-baseline-adc7f12`, at private commit `d2336fc`, corresponding to
`openai/math` commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a`.
This tag remains the comparison baseline after private development merges.
[UPSTREAM.md](UPSTREAM.md) records provenance and package layout. The upstream
Apache 2.0 [LICENSE](LICENSE) is retained.

## Reproduce the checks

Install [elan](https://github.com/leanprover/elan), then run from this repository:

```sh
bash scripts/bootstrap.sh
bash scripts/check-proof.sh
```

The bootstrap script resolves the pinned dependencies and downloads the
matching Mathlib cache. Mathlib is pinned to
`d13f23b723b8a846827a245b89c10fc7d3f11612`; fixed-point-theorems is pinned to
`770940ddf9878cf61952ed53d910b92bca841838`. The Lake post-update hook applies
the included upstream Lean 4.34.1 compatibility patch. Its subsequent
dirty-dependency warning is expected.

The audit first checks that `Model.lean` is unchanged from the baseline,
then builds the public theorem, examples over arbitrary universes and
representative fields, the explicit program-cost statement, and guarded
axiom checks. Run only one Lake build process at a time. Build caches and
local dependencies are ignored by Git.

## Scope and privacy

This is an asymptotic upper bound in the upstream arithmetic-operation model.
The project does not establish historical novelty or practical algorithm
constants, and it does not determine the exact value of the exponent.

Work is authorized only in the private
`selanavot/matrix-multiplication-all-fields` repository. Publication or an
upstream PR requires separate authorization. [AGENTS.md](AGENTS.md),
[STATUS.md](docs/field-port/STATUS.md), and the agent notes support continuation
after a reset.
