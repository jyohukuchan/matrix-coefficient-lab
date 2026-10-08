# Exact coefficient extraction over a fixed finite field

This note separates a terminating procedure in principle from a practical
coefficient generator. It reports source inspection and small Python
experiments; no new Lean theorem, 9/4 witness, minimum-rank result for square
matrix multiplication, or speed claim is established here.

## What the existing theorem guarantees

The relevant premise is the **exact coefficient-rank** theorem, not merely the
arithmetic-program exponent theorem.

* `AuxiliarySeparation/Main.lean:36` proves
  `exactRankExponent_le_nine_quarters_allFields`.
* `Arithmetic/RankExponent.lean:157–161` defines that exponent as the infimum
  of `log_n(exactMatrixRank K n)` over integers `n >= 2`.
* `Arithmetic/RankExponent.lean:219` supplies an actual finite block whenever
  a proposed exponent is strictly above this infimum.

Consequently, for every fixed finite field `F_q` and every rational
`epsilon > 0`, some `n >= 2` has an exact rank decomposition of length `r`
satisfying

```
r < n^(9/4 + epsilon).
```

Write `9/4 + epsilon = A/B` in lowest terms, with positive integers `A,B`.
The acceptance criterion is exactly

```
r**B < n**A.
```

This is integer arithmetic, with no floating logarithms or approximate
rounding. A binary search can calculate the largest allowed integer `r`
between the flattening lower bound `n**2` and the naive upper bound `n**3`.
It is also valid to accept the weak inequality; the theorem above gives a
strict witness, so the stricter search already has a termination guarantee.

The endpoint `epsilon = 0` is different. An infimum bounded by 9/4 need not
be attained by a finite block. Searching for `r**4 <= n**9` may therefore
never stop, even if the existing theorem is true. A finite generator with
positive slack is a better initial research target than an unspecified
"exact 9/4 algorithm".

This reasoning assumes the existing exact-rank theorem and its unchanged
specification are valid. This investigation read the sources and the existing
verification record; it did not rebuild or re-audit the Lean proof.

## The complete coefficient equations

Use row-major indices

```
x = i*n + j       for the left input,
y = j2*n + k      for the right input,
z = i2*n + k2     for the output.
```

The matrix multiplication tensor has coefficient

```
T[x,y,z] = 1 iff j == j2 and i == i2 and k == k2, otherwise 0.
```

An `r`-term decomposition consists of vectors `a[t], b[t], c[t]`, each of
length `n**2`, and the complete system over `F_q` is

```
sum_t a[t,x] * b[t,y] * c[t,z] = T[x,y,z]
```

for all `n**6` triples `(x,y,z)`. There are `3*r*n**2` scalar unknowns.
Every satisfying assignment can be handed to the existing `TensorScheme`
and independently checked at all tensor coordinates.

For `F_2`, each unknown is a bit; the equations become

```
XOR_t (a[t,x] AND b[t,y] AND c[t,z]) = T[x,y,z].
```

SAT encodings can introduce auxiliary product bits and parity constraints;
systems with native XOR support are natural. For odd prime fields use exact
modular arithmetic. For extension fields, encode the actual basis-coordinate
multiplication table, not integer multiplication modulo the field order.

These are formal tensor-coefficient equations. The general finite-field
identity `u**2 = u` on `F_2` must not be used to simplify unrelated polynomial
degeneration certificates. On the other hand, for a genuinely bilinear map,
testing every pair of input basis vectors *does* recover all tensor
coefficients, even over a finite field. Random input tests are unnecessary
and weaker than the available complete coefficient check.

## Eliminate the output coefficients by linear algebra

For fixed `A,B`, construct the matrix `W` with `n**4` rows and `r` columns:

```
W[(x,y),t] = a[t,x] * b[t,y].
```

Let `T_flat` be the `n**4` by `n**2` matrix whose column `z` is
`T[:, :, z]`. The remaining condition is simply

```
W * C = T_flat.
```

One row reduction with all output columns as right-hand sides decides
feasibility and produces `C`. Free variables may be set to zero: the search
needs one certificate, not every possible output family. Thus the nonlinear
search can concern only input factors `A,B`. Equivalently, it seeks `r`
rank-one matrices whose span contains the prescribed `n**2`-dimensional
matrix-multiplication output subspace in `F_q^(n**4)`.

The left and right input factors cannot themselves be eliminated this way
without handling their rank-one constraints. Solving the output system is a
useful reduction, not a solution to the central low-rank search.

## Sound normalization and a complete finite search

Remove terms with a zero input factor. Scale each remaining `a[t]` and
`b[t]` so its first nonzero coordinate is 1, and absorb the two scale factors
into `c[t]`. If two normalized input pairs coincide, add their output
vectors and replace the pair by one term; delete it if the output sum is
zero. Term permutations then allow input pairs to be sorted.

For vectors of length `N = n**2`, the number of normalized nonzero choices
is

```
P = (q**N - 1) // (q - 1),
L = P**2  # possible normalized input-factor pairs.
```

Enumerating all subsets of at most `r` pairs from these `L` candidates, and
solving `W*C=T_flat`, is complete for rank at most `r`. One can enumerate
exactly `r` distinct candidates when `r <= L`: any smaller feasible subset
can be extended with unused candidates, assigning zero output vectors to
the added columns. A more general implementation should preserve the
at-most-r interpretation explicitly.

Further matrix-multiplication stabilizer symmetries could reduce work, but
arbitrary basis changes on all three legs do not preserve the target. Any
additional canonicalization must prove it retains at least one member of
each target-preserving orbit. In particular, forcing every factor to be a
unit vector or presupposing a sparse support pattern is an incomplete
search unless justified separately.

A terminating theoretical generator is now explicit:

1. Fix `F_q` and rational positive slack; compute `A/B` exactly.
2. Visit `n = 2,3,...`.
3. Compute the largest allowed rank, capped by `n**3`; skip if below `n**2`.
4. Exhaust the finite normalized input-factor search at that block size.
5. Solve for `C`, and verify a found scheme with the independent coefficient
   checker and the integer exponent inequality before returning it.

Every fixed-size search is finite, and the exact-rank theorem guarantees a
successful size for positive slack. Thus this procedure terminates in
principle. It is not computationally practical as written.

There is no contradiction between noncomputable Lean choices and this
finite-field search. Once coefficients live in a fixed finite field,
certificate validity is decidable. The proof guarantees a witness without
supplying its efficient construction.

Be precise about bounds: a guaranteed terminating search itself defines a
computable witness-size function of `q` and rational slack. It would be
incorrect to conclude that *no computable bound can exist*. What is missing
here is an explicit useful bound available before running the search, or a
quantitative modulus extracted from the spectral/entropy argument.

## Small experiments actually performed

Only standard-library Python and the existing coefficient checker were used.
No SAT solver was installed, and no long search was run.

| Experiment over F_2 | Observed result | Scope |
| --- | --- | --- |
| Fix the existing naive 2 by 2 seed's `A,B`; solve for `C` | The original eight-term output family is recovered exactly | Tests linear elimination |
| Fix the existing Strassen test seed's `A,B`; solve for `C` | The original seven-term output family is recovered exactly | Tests elimination; this is a supplied known seed, not a discovered proof witness |
| Replace the first Strassen left vector `(1,0,0,1)` by `(0,0,0,1)` | Only two of the four output right-hand sides are solvable | Detects a broken candidate |
| Exhaust all four-pair selections for `2 by 1` times `1 by 2` | 81 feasible selections out of all 126 | A complete tiny rectangular search |
| Independently certify the 81 rectangular schemes | All full tensor coefficient checks pass; all 16 input pairs pass for each scheme | Confirms generated output coefficients, not just linear-solver residuals |

The rectangular experiment has three normalized nonzero vectors on each
input leg, nine input pairs, and `choose(9,4)=126` subsets. Its target has an
identity output flattening of rank four, so four terms are necessary and
the feasible schemes achieve that elementary bound. It does not approach
the square matrix multiplication exponent.

By comparison, square `n=2`, `r=7` already has 225 normalized input-pair
candidates and

```
choose(225,7) = 5,271,289,966,800
```

seven-pair subsets. The unreduced formulation has 84 coefficient bits.
These counts explain why merely replacing coefficient enumeration by
normalization and linear output solves does not make the square problem
easy.

Exact threshold examples at `n=2`:

| Slack | Acceptance inequality | Largest permitted integer r |
| --- | --- | --- |
| 1/4 | `r**2 < 2**5` | 5 |
| 1/100 | `r**50 < 2**113` | 4 |
| 1/1000 | `r**1000 < 2**2251` | 4 |

The supplied seven-term seed fails all these thresholds. It certifies the
solver mechanism, not the desired positive-slack 9/4 conclusion.

## Extension-field search and fixed-overhead descent

A different complete search can enumerate a finite extension field,
matrix size, and rank. Any finite family of coefficients in the algebraic
closure of a finite field lies in some finite extension, so this eventually
covers an algebraic-closure witness. Enumeration must dovetail parameters
if they are all unbounded, rather than searching all sizes over one
extension before visiting the next.

Given an actual `r`-term block of size `u` over a fixed degree-`d` extension,
the implemented descent certificate after powering has `d**2 * r**m` terms
at size `u**m`. If `r**B < u**A`, choose `m` by increasing integers until

```
(d**2)**B * r**(m*B) < u**(m*A).
```

This test uses only integer arithmetic and eventually succeeds. It keeps
the coefficient field fixed and pays descent overhead once. Descending
before powering would instead produce `(d**2*r)**m` terms and may destroy
the exponent margin. Existing `FieldExtension.lean:42` follows the correct
fixed-overhead route.

This explicit final step does not locate the original extension-field
block. Direct search over F_2 already has a theorem-backed termination
argument and avoids extending the search parameters, while an extension
search may be attractive for algebraically structured candidates.

## Concrete next experiments

1. Implement a bounded normalized input-factor search and output solver with
   independent certificates. Start with the complete nine-candidate
   rectangular experiment, plus rejection controls. Report resource limits
   as "search incomplete", never as nonexistence.
2. Implement an exact rational-slack acceptance function, explicitly reject
   zero slack for any API claiming theorem-backed termination, and check
   all inequalities using integers. Distinguish found witnesses, fixed-size
   infeasibility proofs, and open searches.
3. Encode the square `n=2,r=7` equations for a solver and require recovery of
   a known feasible family without hard-wiring its output coefficients.
   This benchmarks the engine; it still cannot establish a near-9/4 block.
4. Transfer actual structure from the separation construction into candidate
   input-factor supports. Keep the unconstrained formulation as the trusted
   specification; a failed restricted search excludes only that support
   family.
5. Investigate a quantitative finite witness from the spectral proof. This
   is the important mathematical task: find an effective size/margin
   estimate or a constructive substitute, rather than presenting the huge
   complete search as a practical answer.

## Conditional extraction from an explicit positive catalyst

The coordinator and the constructive-route investigation identified a
different proposed witness contract: an actual finite restriction

```
D + m*Unit + M_d*S <= D + k*S,
```

over one fixed finite extension `E`, with `m > 0`, `d >= 2`, and
`k < d**tau_inner`. Here `+` means tensor direct sum, multiplication means
tensor product, scalar `k` means `k` independent scalar tensors, and `<=`
means a triple of ordinary local linear maps. `M_d` is the coefficient
tensor for square matrix multiplication. This subsection independently
checks what **would** follow if such a coefficient certificate were supplied.
Existence at the desired threshold must be justified by a separate audit of
the character/state argument; it is not proved by the numerical planner.

### Required actual witness

An executable certificate must include the fixed field representation,
finite tensor coefficients for `D,S`, integers `d,k,m`, and three matrices
mapping

```
source = D direct_sum (k copies of S)
```

onto

```
target = D direct_sum (m copies of Unit) direct_sum (M_d tensor S).
```

The checker must establish equality at **every** target tensor coordinate.
It must also establish the interpretation of `M_d`, the exact exponent
inequality, and either `m > 0` or an independent nonzero-coordinate witness
for `S`. The positive gain rules out `S=0` by the finite-rank obstruction
proved in `Spectrum/Obstruction.lean:325–367`; an implementation can inspect
the supplied `S` directly and reject zero.

Any nonzero coefficient `S[x,y,z]=v` supplies a concrete restriction of `S`
to `Unit`: take the corresponding coordinate on each leg and scale one leg
by `v**(-1)`. This is the important `Unit <= S` certificate. Discard the
positive scalar summand by a coordinate projection to obtain

```
D + M_d*S <= D + k*S.
```

Give exact schemes or naive certified upper bounds for `D,S`. One does not
need their minimum ranks. For example `R_S=dimX(S)*dimY(S)` and
`C_D=dimX(D)*dimY(D)` are valid general upper bounds, and the standard
matrix scheme gives `rank(M_d) <= d**3`.

### Powering and absorption without cancellation

Let `T=M_d`. A valid block-power catalyst is

```
D_b = direct_sum_{i=0}^{b-1} (k**(b-1-i) copies of (D tensor T**i)),
D_b + T**b*S <= D_b + k**b*S.
```

This comparison follows by sequentially applying the original restriction,
not by adding equations and canceling a common tensor. For instance, at
`b=2`, use `D_2 = D*T + k*D` and the chain

```
k*D + T*(D+T*S)
    <= k*D + T*D + k*T*S
     = T*D + k*(D+T*S)
    <= T*D + k*D + k**2*S.
```

In general `D_(b+1)=T*D_b+k**b*D` gives the inductive chain. A certified
decomposition length for `D_b` is

```
C_b = C_D * sum_{i=0}^{b-1} d**(3*i) * k**(b-1-i).
```

Since `Unit <= S`, `D_b <= C_b*S`. Repeating the same catalyst comparison
for `n` copies gives

```
D_b + n*T**b*S <= D_b + n*k**b*S,
n*T**b*S <= (C_b+n*k**b)*S.
```

The first line keeps **one** catalyst `D_b`; the induction swaps one tensor
copy at a time. This is why its rank cost is additive rather than multiplied
by `n`.

If an actual rank-`r_t` scheme is known for `M_(u_t)`, then

```
M_(u_t) <= r_t*Unit,
M_(u_t*d**b) <= (C_b+r_t*k**b)*S.
```

All these comparisons can be implemented by composing the provided maps,
direct-sum maps, tensor products, and explicit matrix coordinate reindexings.
The resulting certificate has

```
u_(t+1) = u_t*d**b,
r_(t+1) = R_S*(k**b*r_t+C_b).
```

It is important not to obtain this formula by simply tensoring `D_b` with
`M_(u_t)`: that route multiplies the catalyst cost by `r_t` and yields a
worse bound. The copy-absorption argument is essential.

Start at `u_0=r_0=1`, using ordinary scalar multiplication. Therefore this
route is noncircular: it does not require finding a near-9/4 matrix scheme
before constructing the first one. It requires the actual catalyst
restriction witness instead. No additive rank equality or cancellation of
the nonzero tensor `S` is assumed.

### Integer-only parameter planning

Write positive rational `tau_inner=A_i/B_i` and
`tau_outer=A_o/B_o`, with `tau_inner < tau_outer`. Check

```
k**B_i < d**A_i.
```

Increase positive `b` until

```
(R_S*k**b)**B_i < d**(b*A_i).
```

This search ends because the strict margin absorbs the fixed factor `R_S`.
Define `Q=R_S*k**b` and `H=R_S*C_b`. For `Q>1`, the recurrence is exactly

```
r_t = Q**t*r_0 + H*sum_{j=0}^{t-1} Q**j
    = Q**t*r_0 + H*(Q**t-1)//(Q-1).
```

The last quotient is exact. The empty sum is zero at `t=0`. The simple
bound `r_t <= Q**t*(r_0+H/(Q-1))` proves eventual success. If a purely
numerical planner accepts `Q=1`, use `r_t=r_0+t*H` rather than dividing by
zero; genuine witnesses for `d>=2` are constrained by the existing
quadratic rank lower bound and cannot have `k<d**2`.

Let `e` be the degree of the fixed coefficient field over the prime field.
Increase nonnegative `t`, computing the recurrence, until

```
(e**2*r_t)**B_o < (u_0*d**(b*t))**A_o.
```

The inner margin and `tau_outer>tau_inner` guarantee that this search ends.
This final inequality pays for descent once; field `E` stays fixed across
all iterations. Even if `tau_inner=tau_outer`, a strict inequality for `Q`
still suffices to absorb the fixed constants, but allocating separate inner
and outer slack makes the argument easier to organize.

A planner returning these integers is **not** a rank certificate. Only
the actual map composition, resulting coefficient decomposition, and its
independent full coefficient check establish the rank claim.

### Synthetic planner experiment, not a tensor witness

An integer-only script was run with artificial numerical parameters

```
d=4, k=23, R_S=2, C_D=1, e=2,
tau_inner=23/10, tau_outer=231/100.
```

There are no supplied `D,S` coefficients or catalyst maps for these numbers,
and no claim that these costs are realizable. Here `23**4 > 4**9`, so the
putative catalyst's exponential gain lies above the proved 9/4 threshold,
while `23**10 < 4**23` leaves a strict inner margin. The script found `b=14`,
then the first successful recurrence step at `t=50`. It checked the
recurrence against the geometric sum at every step and checked that `t=49`
fails the exact final inequality. The resulting **numerical bound** has a
matrix-size bit length of 1401 and a descended term-bound bit length of 3234;
`C_b` has bit length 79 and `Q` has bit length 65. Reporting bit lengths
avoids printing huge coefficients or relying on Python's decimal conversion
limits.

This confirms the arithmetic planning logic and demonstrates potentially
enormous overhead. It discovers no matrix multiplication algorithm and
provides no empirical evidence for a feasible catalyst.

### Enumerating catalyst certificates versus raw coefficients

Over a fixed finite field, once `d,k,m` and the three dimensions of `D,S`
are fixed, enumeration of their coefficients and the three restriction
matrices is finite and exact. If target axis lengths are `a_i` and source
axis lengths are `b_i`, the restriction matrix portion alone has

```
|E|**sum_i(a_i*b_i)
```

candidates before constraints and symmetries. The large `M_d*S` target
makes even small-looking catalysts expensive. Enumerate field extension
degree, tensor dimensions, and the scalar parameters by a fair dovetail if
all are unbounded. Enumerating `m=1` is enough to cover positive-gain
witnesses: any `m>0` scalar direct sum restricts to one scalar summand.

A finite coefficient and local-map witness over the algebraic closure of
F_p lies in some finite extension of F_p. This means extension enumeration
can cover a witness **if** the separate character-contrapositive analysis
establishes one for the chosen `d,k` and threshold. It does not itself prove
that such a witness exists. Matching the quantifiers in that analysis is a
necessary remaining mathematical check.

The catalyst search has more unknown objects than direct `A,B` search and
is not obviously faster. Its potential advantage is alignment with the
proof's finite spectral obstruction: it may expose structured local maps
whose coefficients are accessible without solving a large matrix rank
problem directly. Without such structure or quantitative bounds, it is
another potentially enormous complete search. The best immediate work is
to define and verify the actual witness contract, test composition on
genuine small examples without an exponent claim, and research how to
produce the required spectral-obstruction certificate.
