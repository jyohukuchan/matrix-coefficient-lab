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

The full public entry point and initial permanent audit passed on 2026-10-06
with Lean 4.34.1. The initial incremental run reported a successful dependency
graph of 9437 jobs, including six axiom guards; this did not mean 9437 fresh
source compilations. See the [verification record](docs/field-port/VERIFICATION.md)
and the subsequent [adversarial audit](docs/field-port/ADVERSARIAL.md).

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
# Additional replay of all declarations imported by the auxiliary conclusion:
bash scripts/check-kernel.sh
```

The bootstrap script resolves the pinned dependencies and downloads the
matching Mathlib cache. Mathlib is pinned to
`d13f23b723b8a846827a245b89c10fc7d3f11612`; fixed-point-theorems is pinned to
`770940ddf9878cf61952ed53d910b92bca841838`. The Lake post-update hook applies
the included upstream Lean 4.34.1 compatibility patch. Its subsequent
dirty-dependency warning is expected.

The audit protects eleven original specification and builder files against
the immutable baseline commit, verifies dependency revisions and the exact
allowed patch, then builds the public theorem, arbitrary-universe examples,
explicit program-cost and coefficient-rank statements, and guarded axiom
checks. Python 3's standard library and Git are used for dependency validation.
The kernel script uses Lean's own kernel in a fresh environment; it is not an
independent kernel implementation. Run only one build or replay process at a
time. Build caches and local dependencies are ignored by Git.

## Scope and privacy

This is an asymptotic upper bound in the upstream arithmetic-operation model.
The project does not establish historical novelty or practical algorithm
constants, and it does not determine the exact value of the exponent.

Over finite fields, the upstream `Correct` predicate means equality of
functions on field-valued inputs, which is weaker than formal polynomial
equality. The proof separately establishes the stronger exact coefficient-
tensor rank exponent bound. That identity supports the algebraic upper bound
without exploiting finite-field identities such as `x²=x` over F₂. We have
not added a symbolic-evaluation correctness predicate to the exported program
witness. Constants may depend on the field and positive exponent slack.

Work is authorized only in the private
`selanavot/matrix-multiplication-all-fields` repository. Publication or an
upstream PR requires separate authorization. [AGENTS.md](AGENTS.md),
[STATUS.md](docs/field-port/STATUS.md), and the agent notes support continuation
after a reset.
