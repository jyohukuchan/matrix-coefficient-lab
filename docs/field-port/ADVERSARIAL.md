# Adversarial audit of the all-fields theorem

Date: 2026-10-06. Audited proof: `9bd2a64fa47efe678664b5401b6a8ea6d93238ad`
(`all-fields-proof-v1`). Three fresh reviewers received separate assignments
without the authoring conversation. Their conclusions are source-review
findings, not independent implementations of Lean's kernel.

**Audit in progress.** No proof counterexample has been established. The root
coordinator is performing isolated source re-elaboration, expanded statement
and axiom checks, and fresh-environment kernel replay before a final verdict.

## Confirmed findings

1. **Finite-field semantic distinction.** The unchanged upstream
   `MatrixAlgorithm.Correct` asks for equality on field-valued inputs.
   Over F₂, `x²=x` is true as a function but false as a polynomial identity.
   The predicate must not be described as formal polynomial correctness.
   This does not invalidate the present upper bound: the proof separately
   establishes `exactRankExponent_le_nine_quarters_allFields`, using exact
   equality of every tensor coefficient. `Tensor.RankAtMost.map` preserves
   that identity under arbitrary commutative-semiring homomorphisms. We are
   adding an expanded finite-rank witness and an axiom guard for this stronger
   certificate. We have not added a symbolic-evaluation predicate to the
   exported program witness.
2. **Shared-cache verification limitation.** The initial project build used
   a build-directory symlink to the retired focus harness. The new isolated
   checkout has an initially empty OAI build directory. It shares only exact
   pinned dependency checkouts and their caches; their revisions and allowed
   compatibility patch have been checked. This is not a from-source rebuild
   of Lean or all third-party dependencies.
3. **Documentation defect.** The entropy tag lemma's comment retained the
   paper's loss five although the generalized theorem uses six. Its theorem
   and proof already use six; only the comment is corrected.
4. **Guard coverage.** The original check script protected only Model.lean.
   The widened immutable-baseline guard also protects primary arithmetic
   builders, rank/coefficient definitions, and polynomial expression data.
   Dependency and negative-control checks are being added.

## Reports

- [Specification and semantics](adversarial/specification.md).
- Proof/algebra and reproducibility reports pending integration.

## Pending verification

- Isolated original-commit build of public Main and AllFieldsAudit.
- Expanded audit: exact coefficient-rank witnesses, rational-function field
  over F₂, finite-field semantic counterexamples, and ten axiom guards.
- Fresh-environment replay with Lean's bundled `leanchecker --fresh`.
- Negative control: inject an axiom in a scratch file and confirm the guard
  fails; no injected axiom belongs to repository proof sources.
- Private audit PR with final outcomes and remaining limitations.
