# Root coordinator handoff

Read `../STATUS.md` and repository `AGENTS.md` first.

## Completed this session

- Established ownership with three continuing agents.
- Added generic matrix tensor in `AuxiliarySeparation/Tensor/MatrixMultiplication.lean`.
- Generalized `Arithmetic/RankExponent.lean`, `Tensor/Semiring.lean`, and
  `Character/Basic.lean` to scalar `K`.
- Applied an initial field-parameter pass to the character/semiring, entropy,
  profile, and spectral specialization files listed in STATUS ownership.
- Wrote durable project/restart documentation.

## Validation and remaining work

The complete `AuxiliarySeparation.Main` build succeeded (9053 jobs), including
all-fields exact-rank, arithmetic cost and `Arithmetic.omega` conclusions.
All generic algebra, spectral and arithmetic bridge components are checked.
Independent agent source reviews found no hidden assumptions or weakened
specification. The public `Main` regression is running (session 10954).

`AllFieldsAudit.lean` is drafted with an arbitrary-universe `[Field F]`
example, fields of characteristics 2/3/5, rational/real/complex examples,
the original correctness/cost statement, and six `#print axioms` checks.
It must be built after public Main. A separate read-only auxiliary audit
has passed: arbitrary universe, only `[Field F]`, standard axioms
`propext`, `Classical.choice`, `Quot.sound` for omega, explicit cost,
closure descent and detecting-character existence. The permanent audit now
uses guarded messages for six declarations; not yet run. Three comment-only
source fixes changed stale "complex" descriptions to arbitrary-field ones;
the next build will refresh their traces.

The user now explicitly authorizes private PR merges. Agent setup is saving
the pristine baseline as tag `openai-baseline-adc7f12`. Once all verification
passes, update documentation and PR1, commit/push, and merge privately.

Only root may start aggregate Lake builds. Multiple concurrent Lake processes
caused excessive memory use and transient missing-olean errors. All agents
now use the private canonical checkout's `lean/` directory. Source repair
remains parallel across owned files.

External one-off scripts `work/port_core.py` and `work/port_remaining.py`
performed the initial replacements. **Do not rerun them** on modified files;
they are not idempotent. Future repairs should edit source directly.

The original workspace also contains a compiled explanatory note at
`outputs/matrix-multiplication-all-fields.tex`. Its mathematical audit is
background; compilation of that LaTeX document is not Lean verification.
