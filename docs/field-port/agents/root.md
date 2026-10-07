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

Root-owned foundations have not yet compiled. Agent spectral is trying
`Arithmetic/FieldExtension.lean`, which first builds generic RankExponent.
Expect missing explicit `(K := K)` arguments in closed theorem applications.
Repair the first failing foundation before dependent files.

`Entropy/Tag.lean` still has two wrappers specialized to constant5; update
to6. Then add appropriate infinite/algebraically-closed assumptions to the
character/profile chain. Final `AuxiliarySeparation/Main.lean` is still the
old complex theorem and must be updated after components compile.

External one-off scripts `work/port_core.py` and `work/port_remaining.py`
performed the initial replacements. **Do not rerun them** on modified files;
they are not idempotent. Future repairs should edit source directly.

The original workspace also contains a compiled explanatory note at
`outputs/matrix-multiplication-all-fields.tex`. Its mathematical audit is
background; compilation of that LaTeX document is not Lean verification.
