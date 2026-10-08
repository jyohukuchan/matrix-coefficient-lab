# Where the proof stops being a coefficient generator

This is a bounded source audit of the proof at commit `8d6d358`. It does not
change Lean sources, rerun Lean compilation, or establish a practical `9/4`
scheme. Line references below refer to that source tree. Paths beginning with
`A/` abbreviate `lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/`.

## Main finding

The proof bounds the exact-rank exponent through **dual real-valued tensor
characters**. It does not give a forward chain of restrictions taking the
convolution construction to a particular small-rank matrix multiplication
decomposition. The concrete convolution, determinant, degeneration, tagging,
and separation constructions supply inequalities for every character. A
separate compactness argument ensures that characters detect an excessively
large rank exponent. The contradiction bounds the exponent; an infimum
argument then supplies an unspecified finite rank decomposition.

Consequently, implementing more of the concrete constructions does not by
itself select the matrix size, rank, or coefficients in the final existence
statement. It remains useful for certificates and structured searches.

There is an important qualification: over a **specified computable finite
field**, the existing theorem already implies that exhaustive search for a
positive-slack exact scheme terminates. There is no mathematical impossibility
of coefficient extraction in that setting. The missing deliverable is an
effective, useful search bound or structured generator. Extracting a computable
character is neither necessary for brute-force search nor sufficient by itself
to yield rank-one coefficients.

## Exact proof route and witnesses

| Stage | Exact source | Witness or conclusion | What a generator would still need |
| --- | --- | --- | --- |
| Tensor restriction order | `A/Tensor/Semiring.lean:46–48, 163–167, 191–210` | Existence of three local matrices; quotient by mutual restriction | Explicit matrices and finite representatives; the quotient is not a coefficient data structure |
| Least exact rank | `A/Arithmetic/RankExponent.lean:36–60` | `exactRank := Nat.find (exists_rankAtMost T)`; `exactRank_spec` supplies `RankAtMost` | A decision/search procedure for rank, not just the theorem that some rank exists |
| Catalytic obstruction | `A/Spectrum/Obstruction.lean:102–173, 267–304, 356–377` | No positive scalar-gain catalyst below `d^ν`, where `ν=exactRankExponent K` | This rules out a comparison; it does not construct a low-rank decomposition |
| Finite state constraints | `A/Spectrum/StateObstruction.lean:185–229, 233–279` | Real solution of a finite integer linear system | Present finite tensors and restriction witnesses; choose rank bounds; an exact rational LP implementation is possible for a supplied finite system |
| All additive state constraints | `A/Spectrum/StateObstruction.lean:283–308` | A function on the entire tensor restriction semiring | Compact infinite intersection; no finite cutoff or approximation certificate for all constraints is specified |
| Multiplicativity | `A/Spectrum/MultiplicativeStates.lean:117–190` | A multiplicative state on that semiring | Infinite-dimensional fixed points and another infinite intersection; no finite constructive extraction is supplied |
| Detecting character | `A/Character/Existence.lean:39–57` | A character `χ` with `k ≤ χ(T_d)` for `k < d^ν` | The character is a real-valued functional, not `a,b,c` rank-one families |
| Character growth bound | `A/Polynomial/Inequalities.lean:111–134` | Every character has `pX+pY+pZ ≤ 9/4` | Evaluating this argument does not select a decomposition |
| Rank exponent bound | `A/Arithmetic/RankBound.lean:32–39`; `A/Arithmetic/CharacterRounding.lean:93–105, 125–132` | `ν ≤ 9/4` | A concrete near-optimal block and its rank-one families |
| Positive-slack block | `A/Arithmetic/RankExponent.lean:219–234` | Some `n≥2` and `RankAtMost T_n R` with `R≤n^(ν+ε)` | `n`, `R`, and coefficients are obtained by infimum existence and least-rank witnesses, without an effective upper bound on `n` |
| Descent | `Arithmetic/FieldDescent.lean:21–33, 63–88, 93–128`; `A/Arithmetic/FieldExtension.lean:42–79` | One finite coefficient algebra and rank factor `d²` for all powers | Given represented algebraic coefficients, basis/linear functional arithmetic is implementable; the theorem alone does not produce those coefficients or bound `d` beforehand |
| Arithmetic program | `Arithmetic/RecursiveBlockPrograms.lean:269–366`; `A/Arithmetic/Exponent.lean:68–78` | A correct block program from rank-one families | This is downstream of the missing scheme; most operations are explicit once families are supplied |

The final theorem combines the algebraic-closure bound with descent:
`A/Main.lean:36–40`. For a base field, the proof fixes one extension-field block
before taking powers. It must not charge the descent factor separately at
every recursion level.

## Which existential steps are genuinely infinite?

### Finite state constraints are a plausible executable research component

`StateConstraint` has lower, upper, both signs of additivity, order, and
detector cases (`StateObstruction.lean:34–52`). For a supplied finite collection,
the integer constraint matrix and normalization vector are explicit at
lines 196–197. A Farkas alternative yields a normalized real linear functional
at lines 219–220. A feasible rational-coefficient polyhedron has rational
solutions, so exact LP or rational elimination can replace this isolated step.

However, `.order x y` includes a proof of `x≤y`, meaning an actual restriction
existence statement. It is not permission to guess a monotonicity inequality.
The character-existence application uses `TensorSemiring.rank`, namely exact
rank, for upper bounds. Supplying certified coarse upper bounds instead is
possible in experimental relaxations but must not silently be called the same
finite approximation of the original exact-rank state space.

### The passage to all additive constraints is not a finite LP

`normalizedStates_nonempty_of_no_catalyst` forms the product of intervals
`[0,R(x)]` over **every** semiring element and intersects one closed set for
every state constraint. It uses `isCompact_pi_infinite` at line 294 and compact
intersection at line 306. Finite satisfiability is enough for the proof, but
the source gives no number of constraints sufficient for an executable
character or a small-rank scheme.

### “Finite multiplicative slice” does not mean finite-dimensional

`exists_finite_multiplicative_slice` handles finitely many multipliers `z`, but
its conclusion still requires `∀ x, f(z*x)=f(z)*f(x)`
(`MultiplicativeStates.lean:122–125`). It retains functions on the entire
semiring and calls `compact_convex_fixedPoint_of_continuous` at lines 138–139.
The final multiplicativity theorem intersects slices for all nonzero `z` at
line 186. Merely solving a finite set of bilinear character equations is not a
faithful replacement.

The fixed-point theorem itself proves finite-coordinate approximate fixed
points using a compactness-selected finite cover and partition of unity
(`A/Convex/FixedPoint.lean:108–138`), and obtains the exact point through compact
intersection (`165–203`). Its finite-dimensional Brouwer reduction does not
provide a computable representation of the original infinite-dimensional
state space, a finite-cover construction, or a convergence rate.

### Growth estimates contain a further non-effective limit layer

The scalar growth proof chooses row slopes from limits:
`A/Growth/Profile.lean:45–80`. The limit is the infimum of a decreasing increment
sequence (`A/Growth/ConcaveSlopes.lean:85–98`). Existence and monotonicity supply
the proof but not an error bound after a specified number of indices. The
diagonal inequality and final exponent comparison (`Profile.lean:84–122`) are
again statements about character values. Quantifying this limit layer might
help finite certificates, but it still would not automatically output a rank
decomposition.

## What “noncomputable” does and does not establish

Many files use `noncomputable section` broadly. That annotation is **not** a
proof of inherent noncomputability. Explicit operations inside such sections
can be implemented for concrete fields and finite index sets.

* `exists_rankAtMost` gives an explicit coordinatewise decomposition at
  `RankExponent.lean:40–42`; this is computable for a finite tensor.
* `exactRank` uses `Nat.find` with an existential proposition. With finite-field
  arithmetic and finite enumeration, `RankAtMost T r` is decidable and the
  least rank can be computed, however expensively.
* Tensor reindexing by `Fintype.equivOfCardEq` needs a chosen enumeration; mixed
  radix indexing replaces it computationally.
* The projection formulas in `FieldDescent.lean:28–33` are explicit once a basis
  and linear functional are provided. Their current proof uses an abstract
  chosen basis, but a concrete extension presentation permits linear algebra.
* A proposition `∃ a b c, ...` checked by Lean is a valid existence result. It
  does not follow that the checked proof object is an executable coefficient
  table; propositions and classical choices are not an algorithm interface.
* Algebraic closure and arbitrary field types do not come with a computational
  presentation. Algorithms over finite fields, effective number fields, and
  general abstract fields must be discussed separately.

The true structural obstacle in this proof route is the dual compactness and
infimum reasoning, plus the absence of quantitative bounds. It is not simply
the presence of the word `noncomputable` or `Classical.choice`.

## A complete extraction procedure in principle over finite fields

For a supplied finite field of order `Q` and positive rational slack `ε=s/t`:

1. Enumerate `n=2,3,...`.
2. For this `n`, enumerate rank budgets `0≤r≤n³` satisfying the exact integer
   test `r^(4t) ≤ n^(9t+4s)`.
3. Enumerate the three coefficient families of shape `r×n²` over the field.
4. Check the identity at **every** one of the `n⁶` tensor coordinates. Return
   the first passing scheme, together with its field presentation and exact
   identity certificate.

Every fixed-size stage is finite. The all-fields exponent theorem and
`exists_rankAtMost_of_exponent_slack` imply some stage succeeds, so sequential
enumeration terminates. Searching a stricter positive slack can avoid
numerical equality concerns; the integer threshold above already needs no
floating-point logs. A budget-limited experiment must return “budget exhausted”
when stopped, not imply a proof of nonexistence.

Without normalization, a rank budget contains `Q^(3rn²)` possible tables. This
is a correctness/termination baseline, not a useful implementation plan for
near-`9/4` instances. The proof gives no estimate of the first successful `n`
and no practical bound on the required work. Fixing a finite field and
enumerating it directly avoids having to enumerate characters or extension
coefficients at all.

This establishes only a theoretical procedure. It does not assert that a
fast extractor has been derived, or that existing experiments reach the
exponent promised asymptotically.

## Actionable obligations and bounded next experiments

1. **Define the target certificate precisely:** a field presentation, `n,r`,
   three immutable coefficient arrays, full coefficient verification, and an
   exact rational exponent check. This keeps useful small schemes distinct
   from a claimed `9/4+ε` scheme.
2. **Implement a complete but bounded finite solver as a baseline:** verify
   tiny known positive and negative cases, include search exhaustion status,
   and measure the actual combinatorial growth. Solver research can then
   compare normalization, rank-one dictionaries, SAT, or polynomial methods.
3. **Record a witness ledger for forward constructions:** every local map,
   degeneration, coefficient extraction, extension presentation, and basis
   should be data. This supports later constructive variants and separates
   supplied scheme coefficients from generated ones.
4. **If following the spectral proof, formulate a new finite extraction
   theorem:** identify a finite set of tensors and comparisons, prove that its
   infeasibility implies a specified finite rank restriction, and bound the
   dimensions and degrees. The current finite LP theorem alone does not prove
   this implication or provide a complete coefficient extractor.
5. **Quantify only what is needed:** obtain a finite cutoff depending on
   rational `ε`, an effective state/character approximation with enough margin,
   and a route from its obstruction certificate to explicit rank-one arrays.
   Approximating a character without this final route does not solve the task.

The audit does not establish whether a useful constructive theorem can be
obtained by reworking the existing argument or requires a different idea.
It does establish where the current source provides witnesses, where it uses
infinite compactness, and why additional peripheral transformations alone do
not answer the coefficient-extraction question.

## Follow-up: a finite catalyst is an alternative search target

An independent follow-up audit confirms a stronger conditional route developed
by the main and constructive-research agents. This is a mathematical derivation
from inspected sources; no new Lean statement or build verifies it yet. It
narrows the earlier report's boundary: search can target one finite catalyst,
then explicitly generate schemes without a near-optimal-rank oracle.

### Existence by the correct contrapositive

Fix an algebraically closed `K`, integer `d≥2`, and integer `k>d^(9/4)`.
If the normalized state space for `TensorSemiring.rank`, `matrixClass d`, and
`k` were nonempty, `exists_monotone_semiringHom_of_normalizedStates`
(`MultiplicativeStates.lean:213–232`) would produce a character with
`χ(M_d)≥k`. All its semiring hypotheses are supplied by
`Character/Existence.lean:47–56` independently of that file's `k<d^ν`
assumption. The matrix-value formula and universal growth bound
(`Character/Dot.lean:179`, `Polynomial/Inequalities.lean:131–134`) instead give
`χ(M_d)≤d^(9/4)`, a contradiction.

Contraposing `normalizedStates_nonempty_of_no_catalyst`
(`StateObstruction.lean:283–288`) therefore gives finite tensor classes `D,s`
and natural `m` satisfying

```
m > 0
D + m + M_d*s ≤ D + k*s.
```

Strict `k>d^(9/4)` is essential. Equality need not exclude a detecting
character. The existing `exists_detecting_character` theorem has the opposite
comparison `k<d^ν`; it must not itself be cited for catalyst existence.

Natural numbers here mean **independent scalar tensor summands**, not field
elements reduced modulo characteristic. Thus `m>0` remains meaningful even
if `p|m`. `not_tensor_add_positive_nat_le_self`
(`Obstruction.lean:325–351`) forces `s≠0`, since otherwise the certificate
would assert `D+m≤D`. Any nonzero support entry of a finite tensor gives an
explicit scalar restriction `1≤s`. Selecting away the `m` block gives
`D+M_d*s≤D+k*s`, without cancellation.

For coefficient search, one may fix `m=1`: because `m≥1`, selecting one
scalar summand from the target yields `D+1+M_d*s≤D+k*s` with unchanged `D,s`.
There is no need to enumerate positive scalar gains separately.

### Searchable finite coefficient certificates

Choose finite tensor representatives of `D,s`. The quotient
order-to-restriction equivalence (`Tensor/Semiring.lean:197–210`) turns the
comparison into three finite local matrices from
`D ⊕ (k copies of s)` to `D ⊕ scalar_m ⊕ (M_d ⊗ s)`.
A certificate must include field presentation, `d,k,m`, both tensors,
source/target indexing, all three matrices, and full coefficient equality,
including mixed direct-sum coordinates. Nonzero support in `s` is useful
explicit data. Choosing quotient representatives does not compute them.

Over `AlgebraicClosure F_p`, all finitely many tensor and map coefficients in
such a witness lie in a common finite extension `F_(p^e)`. Enumerating finite
extensions, finite shapes, tensor entries, and local-map entries is complete
in principle. The search must **dovetail** degree and size bounds: enumerating
all shapes over the first field before trying another degree never finishes.
Candidate checking requires neither a character nor exact-rank oracle.
Certified coarse decompositions of `D,s` suffice. This completeness statement
does not automatically extend to an arbitrary abstract field without an
effective presentation, and it gives no practical degree/size bound.

### Powering leaves the auxiliary tensor `s` fixed

Write `T=M_d`. For `b≥1`, define

```
D_b = Σ_{i=0}^{b-1} D*T^i*k^(b-1-i).
```

Then `D_b+T^b*s≤D_b+k^b*s`. Use
`D_b=D*T^(b-1)+k*D_(b-1)`: multiply the original comparison by `T^(b-1)`,
multiply the previous comparison by `k`, add the spectator catalyst pieces,
and compose. This uses monotonicity and distributivity, **not cancellation or
a power of `s`**. The local maps are tensor products, direct sums, and
compositions of the supplied maps. `T^b` reindexes to `M_(d^b)` by
`TensorSemiring.matrixClass_pow` (`Semiring.lean:502–505`).

### Explicit coefficient recurrence

Take supplied rank upper certificates `R≥rank(s)` and `C_b≥rank(D_b)`.
Because `1≤s`, they give `D_b≤C_b*s`. For any supplied seed scheme `M_u≤r`
with `r≥1`, `catalytic_power_comparison` (`Obstruction.lean:186–213`) with
`n=r`, `j=1`, `t=M_(d^b)`, `k=k^b`, `D=D_b`, and `C=C_b` gives

```
M_(u*d^b) ≤ (r*k^b+C_b)*s ≤ R*(r*k^b+C_b).
```

Hence actual schemes can be constructed recursively with fixed constants:

```
u_(t+1) = u_t*d^b
r_(t+1) = A*r_t+B
A = R*k^b
B = R*C_b.
```

For rational target `τ` with `k<d^τ`, choose `b` such that
`A<d^(bτ)`; this exists since `R` is fixed. In the present regime `A>1`, and

```
r_t = A^t*r_0+B*(A^t-1)/(A-1)
    ≤ (r_0+B/(A-1))*A^t.
```

Thus `r_t/u_t^τ→0`, and an exact integer search over `b,t` terminates after
the catalyst is known. The scalar scheme `u_0=r_0=1` is sufficient: an already
improved matrix algorithm is not required. If the witness lies in a finite
extension, fix that extension for the whole construction, descend only the
final scheme, and enlarge `t` to absorb its one fixed dimension-squared
factor. No extension factor is charged at each recurrence step.

For positive rational slack above `9/4`, choose `d,k` with
`d^(9/4)<k<d^τ` for a strictly smaller target `τ` than the requested exponent.
Such pairs exist for sufficiently large `d`. Tests can be exact:
`k^4>d^9`, and a denominator-cleared power comparison for rational `τ`.

### Revised research boundary

A **single finite positive-gain catalyst certificate** is a sufficient
searchable witness. Once found, its powered-catalyst recurrence is an explicit
route to coefficients below any exponent strictly above `log_d(k)`. We need
not compute a global character or use the infimum's unknown near-optimal
scheme. The remaining obstacle is finding that catalyst in manageable tensor
dimensions and extension degree; the contrapositive supplies no useful
search bound. Every local-map identity, nonzero auxiliary tensor, strict
parameter interval, and fixed overhead must be verified. This is a concrete
new extraction route, not a claim that a near-`9/4` instance has been found.
