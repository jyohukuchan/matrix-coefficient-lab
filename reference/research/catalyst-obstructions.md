# Finite obstructions for positive-gain catalyst search

We seek an exact coefficient restriction over one represented finite field:

```
D ⊕ Unit ⊕ (M_d ⊗ S) ≤ D ⊕ (k copies of S).
```

`≤` means three local linear maps, not entrywise inequality. Scalar copies
are independent tensor blocks, not natural numbers reduced modulo the field
characteristic. The concrete target here is `d=2,k=5`:
`2^(9/4)<5<2^(12/5)`, verified by `5^4>2^9` and `5^5<2^12`.
The existing character argument guarantees some finite-extension witness but
does not bound its tensor shapes or coefficient field.

This report provides necessary filters, not a coefficient witness. The
elementary flattening and commutator checks are executable. The stronger
algebra-rank filters additionally invoke the classical Alder–Strassen theorem,
which is not newly formalized or reproved here.

## Three exact flattening filters

Write `f_i(T)` for the rank of the matrix flattening at leg `i`. Ordinary matrix
rank is monotone under local restrictions, multiplicative under tensor
products, and additive under **full independent direct sums**. Since
`f_i(M_d)=d²` and `f_i(Unit)=1`, a candidate must satisfy

```
f_i(D)+1+d²*f_i(S) ≤ f_i(D)+k*f_i(S)
(k-d²)*f_i(S) ≥ 1,  i=0,1,2.
```

Consequently `k>d²`, `S≠0`, and every flattening of `S` is nonzero. For
`d=2,k=5` this filter alone only requires `f_i(S)≥1`. It imposes no lower bound
on `D`: its flattening ranks cancel exactly. Unused coordinates should be
removed first; the concise leg dimensions of `S` are exactly its three
flattening ranks. Extra zero padding cannot improve a candidate.

Tensor rank is **not** used additively here. General tensor-rank direct-sum
additivity cannot be assumed. The rank upper bound from joining supplied
decompositions is always valid, and suffices in the arguments below.

## Elementary commutator obstruction

For a tensor with square other legs of size `N`, an invertible first-leg slice
`A`, and two more slices `B,C`, set `B'=BA⁻¹,C'=CA⁻¹` (or consistently use left
normalization). The tensor rank satisfies

```
rank(T) ≥ N + ceil(rank([B',C'])/2).
```

Here is a characteristic-free proof, including the zero-weight pivot issue.
For an `R`-term decomposition choose a generic invertible slice whose
coefficient on every nonzero first-leg term is nonzero. If necessary work
over a rational-function extension: the forbidden linear forms and required
nonzero determinant/minor polynomials cannot cover the generic point. A
specific nonzero commutator minor persists generically. Rank over an extension
is no larger than rank over the base field, so the lower bound applies over
finite fields too. Normalize the slice and write
`UV=I_N`, `B'=U diag(b)V`, `C'=U diag(c)V`. The projection `VU` has rank `N`,
so `Q=I_R-VU` has rank `R-N`. Since the diagonal matrices commute,
`[B',C']` is the difference of two products each factoring through `Q`.
Its rank is at most `2(R-N)`. No division by the field element `2` is used.

### Unit-sided auxiliaries are impossible for `d=2,k=5`, for any `D`

Any nonzero `S` with a flattening rank `1` is equivalent, after leg permutation
and removal of unused coordinates, to `dot_r` of shape `(1,r,r)`. Its tensor
rank is exactly `r`. Full direct sums of dot tensors in the same singleton-leg
orientation can be handled together in the same way.

Additively repeat a catalyst `q` times while keeping **one fixed `D`**:

```
D ⊕ q*(Unit ⊕ M_2⊗S) ≤ D ⊕ 5q*S.
```

This follows by monotonicity, spectator direct sums, and composition; it does
not cancel `D`. Project the target onto the `q` positive-gain/product blocks.
For `S=dot_r`, that target has an identity slice of size `N=q*(1+4r)`.
Two slices corresponding to left multiplication by `E12,E21` in each `M_2`
block have commutator rank `4qr`: their commutator is `diag(1,-1)` repeated
over the two columns and `r` dot coordinates. This remains invertible on the
matrix blocks in characteristic `2`. Scalar-unit blocks have zero commutator.

The resulting tensor-rank lower bound is `q*(1+6r)`. The source has an explicit
upper bound `C+5qr`, where `C` is any fixed supplied upper bound for `rank(D)`.
For sufficiently large `q` these contradict each other. Therefore **every
flattening rank of a successful `S` must be at least `2`**, irrespective of
`D` or extension degree. A nonzero catalyst cannot rescue `S=Unit`, any single
unit-sided tensor, or a same-orientation direct sum of these.

For square-slice `S` of other-leg size `n`, an invertible slice, and certified
rank upper bound `R`, the same argument gives the necessary inequality
`5R≥6n+1`. If the upper bound is too coarse this filter is inconclusive, not
proof of existence. Mixtures of dot tensors with different singleton-leg
orientations are not covered by this combined commutator proof.

## Classical algebra-rank filter

The external theorem used here is the Alder–Strassen lower bound:

```
rank(multiplication in A) ≥ 2*dim(A)-t(A),
```

for a finite-dimensional associative unital algebra over a field, where
`t(A)` counts maximal **two-sided** ideals. It is valid in positive
characteristic. A checked secondary statement is Alexey Pospelov,
[*Bounds for Bilinear Complexity of Noncommutative Group Algebras*](https://arxiv.org/abs/1003.4679),
Section 5, Theorem 2: `C(A×B)≥2 dim(A)-t(A)+C(B)`. Section 3 defines the
multiplicative complexity `C` and states `C≤bilinear rank`. There is no
characteristic or field-cardinality restriction in this theorem's statement.
Its reference [1] is A. Alder and V. Strassen, *On the algorithmic complexity
of associative algebras*, Theoretical Computer Science 15 (1981), 201–211.
The downloaded article was inspected; the original theorem's proof was not
rederived in this pass, and it is not added to the repository's Lean theorem.

Suppose `S` has a certified rank upper bound `R` and restricts to the
multiplication tensor of a unital algebra `A` of dimension `n` and `t`
maximal two-sided ideals. `M_d⊗A` is multiplication in `Mat_d(A)`, of
dimension `d²n`. Matrix-algebra two-sided ideals correspond to those of `A`,
so their number remains `t`. The direct sum of `q` copies of
`Unit ⊕ Mat_d(A)` is multiplication in the direct-product algebra with
dimension `q*(1+d²n)` and `q*(1+t)` maximal ideals. The lower bound is
`q*(1+2d²n-t)`.

Repeating the catalyst and using source rank upper bound `C+kqR` proves the
necessary prefilter

```
k*R ≥ 1+2d²*n-t.
```

This does not assume tensor-rank direct-sum additivity. The theorem is applied
to the whole direct-product algebra; the target is projected from an arbitrary
`D`, and the only source operation is an explicit decomposition upper bound.

### Every concise binary `2×2×2` auxiliary is excluded

Exact enumeration over `F_2` partitions all `256` cubes. There are `174`
concise tensors, falling into three disjoint `GL₂(F₂)^3` orbits:

| Representative | Nonzero entries / algebra | Orbit size | Certified `R` | `n,t` | Required `5R` |
| --- | --- | ---: | ---: | --- | --- |
| GHZ | `000,111`; `F₂×F₂` | 108 | 2 | `2,2` | `10≥15`, false |
| W | `001,010,100`; dual numbers after output swap | 54 | 3 | `2,1` | `15≥16`, false |
| Irreducible quadratic | multiplication in `F₄` over `F₂` | 12 | 3 | `2,1` | `15≥16`, false |

For dual numbers `A=F₂[ε]/ε²`, the unique maximal ideal is `(ε)`.
For `Mat₂(A)`, the unique maximal two-sided ideal is `Mat₂((ε))`, with simple
quotient `Mat₂(F₂)`. For `F₄`, `Mat₂(F₄)` is simple as an `F₂`-algebra;
its only maximal two-sided ideal is zero. In both cases it has `F₂`-dimension
`8`. Adding the unit gives dimension `9` and two maximal ideals, hence
rank at least `16` per repeated block. A three-term scheme for `S` would
provide at most `15` terms per source block.

The W tensor is exactly dual-number multiplication with the output basis
swapped: multiplication has support `000,011,101`; swapping output coordinates
gives `001,010,100`. Three-term Karatsuba coefficient arrays were independently
verified for both dual numbers and `F₄`. The split tensor has a two-term
certificate. The experiment enumerates every invertible change of basis and
checks that the three orbits cover the full concise inventory.

Thus the binary `2×2×2` search universe can be removed **with the cited
algebra-rank theorem as an explicit dependency**. The finite experiment proves
the orbit coverage and rank upper identities, not the external lower theorem.
Allowing `D` to grow cannot rescue these choices. Over an algebraic closure,
the quadratic field case splits, which also remains excluded.

The cube exclusion extends beyond binary coefficients. Over an algebraic
closure, a concise two-dimensional space of `2×2` slices contains an
invertible slice. If the determinant pencil were identically zero, two
rank-one generators would have a common image or common kernel (put one
generator into `E11` and inspect the other), contradicting another leg's
conciseness. Normalize the invertible slice to `I`; the second slice is
nonscalar. Its Jordan form is either diagonal with distinct eigenvalues or
one size-two Jordan block. Changing the first-leg basis yields respectively
the split-algebra or dual-number normal form above. A hypothetical certificate
over the original field remains a certificate on extending scalars, where
both normal forms are excluded. Thus **all concise `(2,2,2)` auxiliaries over
any field are excluded**, with the same classical lower theorem as dependency.
The binary enumeration is an independent finite check, not the sole reason
for the field-general claim.

### Convolution `C(2,2)` is excluded for arbitrary `D`

`C(2,2)` has shape `(2,2,3)` and the exact three-term Karatsuba scheme over
every field. Selecting its first two output coefficients gives multiplication
in `K[x]/x²`, a unital algebra with dimension `2` and one maximal ideal.
The algebra filter requires `5*3≥16`, which fails. This is a proof of
nonexistence of this **specified auxiliary family**, not a failed bounded
search or a denial of the catalyst existence theorem.

More generally, `C(n,n)` restricts by retaining its first `n` outputs to
`K[x]/x^n`, which has dimension `n` and one maximal ideal. With a certified
evaluation/interpolation decomposition of `R=2n-1` terms over a sufficiently
large represented field, the necessary inequality is
`5(2n-1)≥8n`, so `n≥3`. This excludes `n=1,2`; it does **not** exclude all
convolutions. `C(3,3)` with a certified five-term scheme satisfies `25≥24`,
so this filter alone is inconclusive. A small field lacking interpolation
nodes may require a larger certified upper bound; never pretend it has the
`2n-1` scheme over that field without constructing one.

### Every concise `(2,2,3)` auxiliary is excluded too

A tensor with these concise dimensions gives a surjective bilinear coefficient
map `K²⊗K²→K³`, so its kernel is one-dimensional in the four-dimensional
outer-product space. A nonzero generator of that kernel is a `2×2` matrix.
Independent input basis changes classify it by matrix rank, either `1` or `2`.
Output basis changes identify quotient maps with the same kernel. There are
exactly two normal forms over any field:

* Rank-one kernel: the three independent products `a₀b₀,a₀b₁,a₁b₀`, omitting
  `a₁b₁`. This has a three-term decomposition. Keep the first output and add
  the other two to obtain dual-number multiplication.
* Rank-two kernel: send its invertible generator by left/right basis changes
  to the antisymmetric relation `a₀b₁−a₁b₀`. The quotient is `C(2,2)`.
  Its Karatsuba scheme has three terms; selecting the first two outputs gives
  dual-number multiplication. The relation matrix remains invertible in
  characteristic `2`.

Both normal forms have `R=3` and a witnessed restriction to the dimension-two
local algebra. Thus both violate `5R≥16`, for every `D`. The classification
does not depend on rank additivity or determinant polynomials having rational
roots.

Over `F₂`, enumeration of all `4096` tables finds `2520` concise tensors.
There are `15` nonzero kernel lines and `168` output bases each. Nine kernel
matrices have rank one and six rank two, giving orbit sizes `1512` and `1008`.
The experiment enumerates `GL₂(F₂)×GL₂(F₂)×GL₃(F₂)`, checks disjointness and
coverage, and verifies each normal form's rank-three scheme and exact output
map to dual numbers.

After compression and leg permutation, a successful auxiliary must therefore
have at least two concise dimensions above `2`, or third concise dimension
at least `4`. The smallest remaining shape with two dimension-two legs is
`(2,2,4)`, equivalent to the full outer-product tensor `M(2,1,2)` of rank four.
Its product with `M₂` is `M(4,2,4)` of shape `(8,8,16)`. Adding the unit gives
flattening lower bound `17`, below the twenty-term source upper bound.
The initial elementary checks below were inconclusive; a subsequently
inspected rectangular lower theorem resolves this boundary in the follow-up.

The standalone experiment independently computes a Koszul flattening for
`Unit⊕M(4,2,4)` over `F₂`. With first-leg dimension `9`, exterior degree `5`,
the matrix has shape `1428×1134` and exact rank `960`. A rank-one tensor has
Koszul rank at most `binomial(8,5)=56`, so the tensor-rank lower bound is
`ceil(960/56)=18`, which does not exclude source budget `20`. This is an
inconclusive exact check, not an existence result. The collaborating solver
agent obtained lower bound `20` on some full-rank first-leg projections; those
additional projections are not independently rerun by this experiment. No
checked bound above `20` is available from these calculations, and a
single-block bound would need an argument for repeated blocks to rule out
arbitrary `D` by the fixed-catalyst method.

## Follow-up: rectangular lower bounds exclude the full outer-product family

The main agent identified a newer primary source:
Jason Yang, [*New lower bounds on tensor rank of `(n,2,m)` matrix multiplication
with GPT-6*](https://arxiv.org/abs/2609.14393), version 2 dated September 28,
2026. Its Theorem 1 states

```
rank(M(n,2,m)) ≥ (n+2)*m,  n≥4,
```

over arbitrary fields. This is **exact tensor rank**, not border rank. The
definition is an actual bilinear decomposition; finite-field coefficient
identities follow by evaluating on basis inputs. Section 2 passes to an
algebraic closure to make generic normalization available. This is legitimate
for a lower bound because extension cannot increase rank. There is no
characteristic restriction, division by an integer prime, or limitation on
the number of columns `m` in Theorem 1. Section 3.1 proves it using injectivity
of a finite linear map (Lemma 6) and dimension counting. The n=3 theorem is
not needed here, and no finite-set approximation to border rank is involved.

The primary proof was read, including the generic rank-one normalization and
the cancellation/injectivity argument in Lemma 6. No Lean proof is bundled
with this source or checked by this pass. There are presentational issues in
Section 2: its displayed basis-change formula needs the output transformation
consistent with the input transformations; standard basis-change covariance
supplies this normalization. Also matrix products span the output space,
rather than literally exhaust every output matrix when inner dimension is
two. The surjectivity assertion used in the proof follows from that span,
since every matrix unit is a product. These do not require a different
hypothesis for Theorem 1. The new filter explicitly depends on this published
preprint theorem, not the repository's existing Lean theorem.

### Shared-first merging avoids rank-additivity assumptions

The `q` independent copies of `M(n,2,m)` restrict to `M(n,2,qm)`:
on the first leg, sum the `q` coordinate copies of each `n×2` input entry;
on the second and third legs, concatenate the column blocks and reindex from
block-major to matrix row-major order. Mixed blocks remain zero. This is a
three-local-map restriction, not cancellation or a tensor-rank additivity
assertion. Applying Theorem 1 to the **single merged matrix tensor** gives
lower bound `(n+2)*qm` uniformly in `q`.

For `S=M(2,1,2)`, the tensor rank is exactly `4`, and `M₂⊗S` reindexes to
`M(4,2,4)`. Repeated catalyst comparisons therefore imply

```
rank(M(4,2,4q)) ≤ rank(D)+20q
rank(M(4,2,4q)) ≥ 24q.
```

Choosing `q>rank(D)/4` contradicts them, even if the positive scalar summands
are discarded. Hence the previously inconclusive `(2,2,4)` auxiliary is
excluded for every field and every finite `D`, conditional on the inspected
rectangular lower theorem. All auxiliaries with two concise dimensions equal
to two are now excluded: the third is at most four, and the `(2,2,2)`,
`(2,2,3)`, `(2,2,4)` cases have all been covered.

### Scalar-summand additivity has an elementary coefficient proof

For any tensor `T`, `rank(T⊕m Unit)=rank(T)+m`. This special additivity does
not extend by assumption to arbitrary tensor summands. To remove one scalar
block from an `r`-term decomposition, select term `j` with nonzero coefficient
`c_j[s]` at its scalar output. On old input coordinates, the scalar output
identity is `0=Σ_i c_i[s]*a_i_old⊗b_i_old`. Eliminate its j-th product from
each old output. The remaining `r-1` terms have

```
a_i' = a_i_old
b_i' = b_i_old
c_i' = c_i_old - (c_i[s]/c_j[s])*c_j_old,  i≠j.
```

These are exact coefficient arrays for `T`. A term with `c_j[s]≠0` exists
because the scalar block's coefficient is one. Repeat for `m` scalar blocks;
the converse rank upper bound concatenates decompositions. This proof works
over every field and does not divide by the integer `m`.

### General outer-product parameter obstruction

Let `S=M(a,1,c)`, the full outer-product tensor. It has shape `(a,c,ac)` and
rank exactly `ac`. The tensor product with `M₂` reindexes to `M(2a,2,2c)`.
After repeating a positive-gain catalyst, keep the `q` unit summands while
merging the rectangular blocks. For `a≥2`, Theorem 1 plus the scalar
elimination argument give

```
rank(q Unit ⊕ M(2a,2,2cq)) ≥ q + (2a+2)*2cq
                                = q*(1+4ac+4c).
```

The source rank is at most `rank(D)+5acq`. Swapping row and column roles
gives also `5ac≥1+4ac+4a`. A necessary condition is therefore

```
ac ≥ 1+4*max(a,c).
```

For `a,c≥2`, this implies `min(a,c)≥5`. Cases `a=1` or `c=1` are already
excluded by the unit-sided argument. At `min(a,c)=4`, retaining the positive
units is essential: the leading rectangular bound without them can be equal
to the source budget, whereas the extra `q` unit rank makes the inequality
strict. The condition is only necessary; `a,c≥5` does not provide a witness.

The resulting smallest concise shape not removed by these filters, allowing
leg permutations, is `(2,3,3)`; it remains a candidate, not an existence claim.

## What remains viable, and how to use these filters

* Require all three concise dimensions of `S` to be at least `2`.
* Remove the entire concise binary cube inventory and all concise `(2,2,3)`
  auxiliaries using the documented classical theorem and checked normal forms.
* With the explicitly cited rectangular theorem, remove every auxiliary with
  two concise dimension-two legs and all full outer-product auxiliaries with
  `min(a,c)≤4`.
* If a candidate `S` restricts to an algebra multiplication tensor, attach the
  actual restriction maps and rank upper certificate, then apply the integer
  inequality `kR≥1+2d²n-t` before solving catalyst maps.
* Direct sums of scalar tensors fail, even with arbitrary `D`. Matrix
  auxiliaries are **not all excluded**: for `S=M_u`, this algebra filter gives
  `5R(M_u)≥8u²`. The known seven-term `M₂` scheme passes `35≥32`; this is only
  survival of a necessary filter, not a catalyst witness.
* `D=0` requires a direct restriction from `kS` to `Unit⊕M_dS`. Nonzero `D`
  can matter for surviving families, but it cannot repair any of the proven
  auxiliary exclusions because repetition leaves its cost fixed.

## Executed checks and scope

Run `python3 -m reference.research.catalyst_obstructions_experiment`.
It completed `16` exact flattening/commutator cases over `F₂,F₃,F₄,F₅`, with
`r=1,2`, `q=1,2`. It also exhausts the binary cube inventory and verifies
the three orbit sizes and algebra representative rank upper certificates.
It further exhausts all `4096` binary `(2,2,3)` tables and certifies both normal
forms, their rank-three decompositions, and exact dual-number maps.
The rectangular boundary check constructs and ranks its full `1428×1134`
binary Koszul matrix, obtaining the inconclusive lower bound `18`.
It additionally verifies the shared-first rectangular merge at every tensor
coefficient over `F₂,F₃` for `q=1,2`. The inequality `24q` itself is the cited
primary theorem, not inferred from the numerical Koszul check.
These are exact field computations, not randomized input evaluations.

No useful near-`9/4` catalyst coefficients were produced, no Lean sources were
edited or built, and no rank-additivity conjecture was used.
