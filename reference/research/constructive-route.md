# Research B: constructive routes through the spectral proof

This is a source-grounded research report, not a new proved theorem or a
generated `9/4` decomposition. The Lean files were read; no Lean sources were
modified or rebuilt. Paths below are relative to the repository root. Line
references refer to the current source tree.

## Finding

For a specified finite field and a positive rational exponent slack, there is
already a mathematically terminating coefficient-extraction algorithm:
enumerate candidate finite decompositions and check their exact coefficients.
The all-fields exact-rank theorem guarantees that some finite candidate meets
the bound. This gives neither a usable running time nor a size bound.

The more interesting proof-guided route is to search for a **finite rational
state-inconsistency certificate with explicit tensor-restriction witnesses**.
Finite Farkas problems themselves are computable. The current proof does not
select the finite tensor inventory, certificate, catalyst, or near-optimal
matrix block needed to turn that route into a practical generator.

There is a stronger derived constructive route: an explicit catalytic
restriction can be amplified while retaining one fixed auxiliary tensor, then
combined repeatedly with any supplied matrix multiplication scheme. This gives
a coefficient-producing recurrence with provable convergence below a selected
slack exponent. The remaining central problem is acquiring a useful finite
catalyst and its explicit restriction matrices. This derivation has not been
formalized in Lean or implemented in this research pass.

## Where the finite constructions cease to select witnesses

1. `AuxiliarySeparation/Character/FiniteSeparation.lean:57` defines
   `separationPolynomialApproximation` with degree `(M-1)^2`, leading degree
   zero, and a coordinate rank bound. Its polynomial and finite local maps
   are genuine construction data. The finite Fourier period satisfies
   `5M <= L <= 6M`; the scalar image of `L` is nonzero.
2. `Character/FiniteSeparation.lean:77–99`,
   `Character.value_separationTarget_le`, passes from finite evaluations to a
   character inequality. Its decisive input is
   `Character.value_le_of_polynomialApproximation`
   (`Character/Degeneration.lean:103–139`). The latter takes all powers and
   removes the interpolation factor `D*n+1` using an asymptotic inequality.
   A finite implementation must retain that factor. It cannot treat the
   character inequality as a single finite restriction of identical cost.
3. `Spectrum/StateObstruction.lean:283–308`,
   `normalizedStates_nonempty_of_no_catalyst`, combines every finite state
   system by infinite-product compactness. It does not choose one finite
   system or give a bound on its size.
4. `Spectrum/MultiplicativeStates.lean:117–156`,
   `exists_finite_multiplicative_slice`, uses a compact convex fixed-point
   theorem at each finite step. Although the set of distinguished
   multiplicative elements is finite there, the state space still has all
   tensor-class coordinates. This is not a finite polynomial solver.
   `exists_multiplicative_of_invariant_states` at lines 161–190 then applies
   compactness again to all distinguished elements.
5. `Character/Existence.lean:39–57`, `exists_detecting_character`, puts those
   two existence arguments together. It produces a real-valued character;
   it does not return decomposition coefficients.
6. `Arithmetic/CharacterRounding.lean:93–105`,
   `exponent_le_of_detecting_characters`, turns the character bound into an
   infimum bound. `Arithmetic/RankExponent.lean:219–234`,
   `exists_exactMatrixRank_lt_rpow` and
   `exists_rankAtMost_of_exponent_slack`, then extract an unspecified finite
   block from that infimum. The block size and coefficient arrays are not
   selected effectively in these definitions.

Thus removing `Classical.choice` from Fourier roots or interpolation nodes
would not resolve the central coefficient-generation question. Those finite
choices are already straightforward over suitably supplied finite extensions.

## A finite certificate that could actually be checked

Let `S` denote finite tensors modulo mutual linear restriction, with order
`x <= y` meaning that `x` is a restriction of `y`.
`Tensor/Semiring.lean:46–48` defines the required concrete order witness:
three finite matrices and the exact tensor equality.

For a finite selected tensor inventory, assign one variable `f(x)` per
selected class. `Spectrum/StateObstruction.lean:34–52` gives the integer
constraint vectors:

| Constraint | Required scalar inequality |
| --- | --- |
| lower | `f(x) >= 0` |
| upper | `R(x)*f(1) - f(x) >= 0` |
| addPos / addNeg | both signs of `f(x+y)-f(x)-f(y)` |
| order | `f(y)-f(x) >= 0` for a witnessed `x <= y` |
| detector | `f(M_d*x)-k*f(x) >= 0` |

Normalize `f(1)=1`. `R(x)` can be a supplied rank upper bound with an explicit
decomposition. It need not be the unknown minimum rank: the state-space
lemmas take arbitrary `R` satisfying their domination assumptions. Using
naive decompositions is sound, though likely weak.

For integer columns `A_i`, a checkable dual certificate is

```
m > 0, n_i >= 0 integers,
sum_i n_i * A_i = -m * e_1.
```

This is the exact certificate in
`Spectrum/StateObstruction.lean:162–181`,
`no_stateConstraint_certificate`.
`Convex/Farkas.lean:65–74` is the finite alternative;
`Convex/RationalCone.lean:162–205` rationalizes a real solution and clears
denominators. An executable implementation can use rational linear programming
or rational elimination and verify the resulting integer vector directly.

The certificate must also carry:

- every tensor's field, shape, and exact coefficients;
- decomposition arrays establishing each upper bound;
- three local matrices establishing each order constraint;
- explicit identifications/restrictions when two representatives share a
  class; quotient-class equality cannot be guessed from matching dimensions;
- addition/product coordinate conventions and the associated reindexings.

An infeasible numerical LP or a list of scalar inequalities without these
tensor witnesses does not extract decomposition coefficients.

### Existence of a finite dual certificate: conditional argument

Choose `d >= 2` and integer `k > d^(9/4)` over an algebraically closed scalar
field where the proof's uniform character bound applies. If the normalized
state space were nonempty, `normalizedStates_exists_multiplicative`
(`Spectrum/MultiplicativeStates.lean:194–209`) would produce a character
detecting at least `k` on `M_d`, contradicting that bound. Therefore the state
space is empty. The compactness argument used in
`StateObstruction.lean:289–306`, read contrapositively, entails an inconsistent
finite constraint family; the finite rational cone alternative yields a dual
integer certificate.

This is an existence argument for a search target, not a bound on the search.
Over a finite field, tensor arrays and finite restriction matrices can be
enumerated, and fixed-pair restrictions can be decided by finite enumeration.
Over the algebraic closure of a finite field, witnesses can be enumerated over
finite extensions. One must keep the field/extension embeddings consistent.
The above character argument must not simply be applied to the finite base
field with an assumption of infinite interpolation points.

## From a dual certificate to a catalyst

`Spectrum/Catalyst.lean:28–56`, `exists_catalyst_of_completion_eq`, converts an
additive-completion equality to

```
D + m + M_d*s <= D + k*s.
```

`exists_catalyst_of_completion_sum_eq` at lines 62–86 aggregates the selected
relations. The current implementation selects an additive-localization witness
`C` at line 49. An executable certificate needs the corresponding finite
catalyst and actual local maps; a bare equality in the Grothendieck group is
insufficient data for applying a scheme.

Proposed implementation strategies are either to retain finite relation
rewrites while building the certificate and construct their cancellation
catalyst, or to enumerate candidate catalysts/restrictions and check them. The
latter is exhaustive and potentially enormous. Neither strategy has been
implemented or proved quantitatively in this research pass.

## A conditional finite bootstrapping formula

Suppose the following concrete data have been supplied and certified over one
field:

- a restriction `D + M_d*s <= D + k*s`;
- a scalar restriction `1 <= s`;
- a restriction `D <= C*s`, and an `R`-term decomposition of `s`;
- an `r`-term decomposition of the square matrix multiplication tensor `M_u`.

For natural `a` and positive `j`, these data construct a scheme for

```
block size b = (d*u^a)^j,
term bound  R * (r^a*k + C)^j.
```

Derivation: the seed gives `M_(u^(a*j)) <= (r^a)^j`. Tensor with `M_d^j`,
insert the scalar restriction `1 <= s`, and apply
`Spectrum/Obstruction.lean:186–225`, `catalytic_power_comparison`, with
`n=r^a`. Finally decompose `(r^a*k+C)^j` copies of `s`. This is exactly the
finite comparison used by `catalytic_semiring_obstruction` at lines 248–261,
with a supplied seed in place of its minimum-rank budget witness. Tensor
products and composition of the local matrices return actual coefficients.

If `tau=log_u(r)` and

```
C < r^a * (d^tau - k),
```

then `log_(d*u^a)(r^a*k+C) < tau`; sufficiently large `j` absorbs the fixed
factor `R` and yields an improvement of this seed's certified exponent.
This inequality is a derived proposal, not a repository theorem or a produced
useful catalyst. All rank upper bounds may be conservative.

The source's obstruction proof chooses an already near-optimal block from the
infimum at `Spectrum/Obstruction.lean:113–114`; that is the witness we are
trying to generate. Naively iterating the displayed formula with a large `j`
does not establish convergence. The following amplification removes that
problem.

### Powered catalyst and a convergent constructive recurrence

The main agent proposed the following fix; this agent checked its ordered-
semiring derivation. It is a derived proposal, not a repository theorem.

Write `T=M_d`. For a fixed positive `b`, define

```
D_b = sum_(i=0)^(b-1) k^(b-1-i) * D * T^i.
```

The original catalyst restriction entails

```
D_b + T^b*s <= D_b + k^b*s.
```

Proof by induction uses `D_(b+1)=D*T^b+k*D_b`. Multiply the original
inequality by `T^b`, add `k*D_b`, then compose with the `b`-stage inequality
multiplied by `k` and augmented by `D*T^b`. These operations use only sums,
products, and composition of restrictions. They never cancel a catalyst and
they keep **the same tensor `s`** at every stage. Coordinate maps can therefore
be constructed from the original maps using the existing finite operations.

With a `C_D`-term scheme for `D` and the cubic scheme for `M_d`, a conservative
explicit bound is

```
C_b = C_D * sum_(i=0)^(b-1) k^(b-1-i) * (d^3)^i.
```

Consequently `D_b <= C_b <= C_b*s`. Apply
`catalytic_power_comparison` with `n=r_t` and `j=1` to any supplied
`r_t`-term scheme for `M_(u_t)`. This produces a new exact scheme with

```
u_(t+1) = u_t * d^b,
r_(t+1) = R * (k^b*r_t + C_b).
```

`r_t` here is a certified term upper bound, never an unknown minimum rank.
Padding a scheme by zero terms makes this recurrence exact as a term count if
desired. The integer `k` denotes tensor copies, not scalar reduction modulo
the field characteristic.

Choose a target `tau_goal > log_d(k)`, and then choose `b` such that

```
Q = R*k^b < d^(b*tau_goal).
```

Such a `b` exists because `R` is fixed. For `Q>1`, setting `c=R*C_b` gives

```
r_t = Q^t*r_0 + c*(Q^t-1)/(Q-1)
    <= Q^t * (r_0 + c/(Q-1)),
u_t = u_0*d^(b*t).
```

Let `eta=d^(b*tau_goal)/Q > 1`. A sufficient stopping condition is

```
eta^t > (r_0 + c/(Q-1)) / u_0^tau_goal.
```

This can be found by increasing the integer `t` and making exact rational-
exponent comparisons. Thus an explicit catalyst with `log_d(k)<tau_goal`
gives a finite terminating coefficient generator, with a conditional parameter
bound. It does not require searching for a near-optimal seed.

For `tau_goal>9/4`, sufficiently large `d` admit an integer
`d^(9/4)<k<d^tau_goal`. The finite-certificate existence argument above then
provides a catalyst in principle over the algebraic closure. Its finitely
many coefficients, restriction matrices, and scalar-unit witness belong to
one fixed finite extension of a finite base field. All powered catalysts and
recurrence steps remain in that extension. If its degree is `e`, descend the
**final** `r_t`-term scheme at cost `e^2*r_t`, replacing the stopping condition
by

```
eta^t > e^2 * (r_0 + c/(Q-1)) / u_0^tau_goal.
```

The extension factor is fixed across the recurrence. Descending every seed
and charging that factor again at each iteration would invalidate this bound.

This improves the theoretical extraction picture substantially: a useful
catalyst is now sufficient data for convergence, not merely a scalar exponent
certificate. The source still supplies no size bound or efficient algorithm
for finding that catalyst. Formalization and a small finite example with all
maps verified would be appropriate next checks before treating this route as
implemented.

## Bounds available and bounds still missing

The source gives Fourier period at most `6M`, separation degree `(M-1)^2`,
linear interpolation overhead in the power, and fixed coefficient-algebra
descent overhead. These bound costs **after choosing the finite data**.

For a fixed integer LP, standard basic-feasible-solution extraction gives a
useful independent proposal for bounding denominators: a nonnegative solution
of `A*x=-e_1` may be supported on `s <= number_of_coordinates` independent
columns. If all relevant integer entries have absolute value at most `H >= 1`,
Cramer's rule and Hadamard's bound give a common denominator and cleared
coefficients bounded by `s^(s/2)*H^s`, using a nonsingular `s`-row minor.
This elementary bound is not established as a quantitative theorem in the
repository. It controls only one already selected finite system.

No inspected source selects or bounds:

- tensor dimensions and number of tensors sufficient for an inconsistent LP;
- finite-extension degree containing the useful restriction witnesses;
- catalyst dimensions or ranks after completion;
- the block selected by `exists_exactMatrixRank_lt_rpow`;
- the amount of work required to obtain a practically useful scheme.

## Guaranteed extraction versus useful extraction

For a fixed finite field and rational `epsilon=p/q > 0`, enumerate `n >= 2`,
term counts, and all arrays defining decompositions of `M_n`. Verify all
coefficients and the exact integer inequality

```
R^(4*q) <= n^(9*q + 4*p).
```

`Arithmetic/RankExponent.lean:230–234` together with the all-fields rank
bound supplies existence, so the exhaustive procedure eventually terminates.
Each fixed `(n,R)` search is finite. This assertion depends on the theorem's
validity; it is not an observed runtime result. For arbitrary fields without
effective coefficient representations and equality, no analogous executable
algorithm follows merely from the arbitrary-field theorem.

The next meaningful practical research artifact would be a **small explicit
catalytic restriction plus a coefficient-producing compiler for the powered-
catalyst recurrence**, or a finite witnessed LP that produces one. A tiny
demonstration of rational LP dual certificates is useful infrastructure, but
cannot be reported as solving the `9/4` coefficient-extraction problem.
