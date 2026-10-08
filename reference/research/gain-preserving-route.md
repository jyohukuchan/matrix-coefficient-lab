# Keeping positive scalar gains improves the coefficient recurrence

This is an independently audited mathematical derivation proposed by main.
It strengthens the earlier conditional catalyst compiler by keeping scalar
gains throughout powering and absorption. It is not a newly Lean-formalized
theorem, and it does not find the useful original catalyst.

All copy counts in this report are ordinary nonnegative integers. A scalar
gain of two denotes two independent tensor blocks even in characteristic two.

## Preserve gains while powering the same auxiliary

Suppose actual local matrices certify

```
D ⊕ m Unit ⊕ (M_d tensor S) <= D ⊕ k*S,
```

with `d>=2`, `m>=1`, `S!=0`. Set

```
D_1=D, g_1=m,
D_(b+1)=(D tensor M_(d^b)) ⊕ k*D_b,
g_(b+1)=k*g_b+m*d^b.
```

These data construct actual maps for

```
D_b ⊕ g_b Unit ⊕ (M_(d^b) tensor S) <= D_b ⊕ k^b*S.       (1)
```

Induction has two steps. First apply `k` copies of the previous comparison
using the same fixed auxiliary S, retaining all `k*g_b` unit blocks. Then
apply the **original positive-gain** comparison tensored with `M_(d^b)` to
the remaining `D tensor M_(d^b)` and `k*(M_(d^b) tensor S)` blocks. This yields
an additional `m` independent copies of `M_(d^b)`, not just a discarded gain.

Each matrix tensor `M_n` restricts to `n` independent scalar units by
selecting coordinates `(i,i)` on all three matrix legs. The coefficient at
three selected diagonal coordinates is one precisely when their indices
agree. Thus the extra `m` matrix blocks supply `m*d^b` units. The other matrix
product reindexes to `M_(d^(b+1)) tensor S` as in the existing core compiler.

No catalyst cancellation occurs, no scalar `m` is divided out, and S is not
replaced by a new tensor. Previously produced scalar blocks are spectators
during the second step. Their order relative to the new gains can be fixed
by direct-sum coordinate permutations.

The explicit count is

```
g_b=m*sum_(i=0)^(b-1) k^(b-1-i)*d^i
   =m*(k^b-d^b)/(k-d).
```

The quotient formula is a numerical identity. Implementations can use the
integer recurrence without division; gain blocks remain separate objects.

## Use the catalyst's decomposition directly

Let a supplied exact scheme for `D_b` have `C_b` terms, and a supplied exact
scheme for S have `R_s` terms. Let the current square matrix seed have size u
and r terms. Set `K=k^b` and `v=d^b`.

The stronger construction starts from

```
C_b scalar units ⊕ r*K copies of S.
```

Its scheme has exactly `C_b+r*K*R_s` terms. The catalyst scheme's three
families restrict the first `C_b` units to `D_b`. This differs from the old
promotion through `C_b` copies of S, which unnecessarily multiplied that
catalyst cost by `R_s`.

Apply (1) successively to the same retained catalyst and each group of K
auxiliary copies. Preserve the accumulated gain units as spectators. The
result contains

```
D_b ⊕ r*g_b Unit ⊕ r*(M_v tensor S).
```

Discard the final catalyst, apply the r-term seed's families to the r matrix
slots, and extract a scalar unit from nonzero S. The exact output becomes

```
r*g_b Unit ⊕ M_(u*v).
```

This transformation must use the full cross-slot catalyst maps, not assume
that `D_b` is passed through unchanged. It produces actual coefficient arrays
with `C_b+r*K*R_s` terms before the final scalar elimination.

## Exact scalar-output elimination

An exact r-term decomposition of `m Unit ⊕ T` yields an `(r-m)`-term
decomposition of T. This special property needs no general direct-sum
additivity theorem.

For one scalar block at output coordinate s, choose a term j with
`c_j[s]!=0`. Such a term exists because the scalar tensor's coefficient is
one. On the remaining input axes, the scalar output identity says

```
sum_i c_i[s]*a_i_old tensor b_i_old=0.
```

Eliminate the j-th old input product from each remaining output. For `i!=j`
set

```
a_i'=a_i_old,
b_i'=b_i_old,
c_i'=c_i_old-(c_i[s]/c_j[s])*c_j_old.
```

Remove that scalar coordinate on all three axes and drop term j. The new
families agree with the remaining target at every coefficient. Repeat for all
scalar blocks. Divisions are by verified nonzero field elements, never by
the integer number of blocks. The full direct-sum target identity, including
zero cross-block coefficients, is essential to this argument.

Applying this elimination to the retained `r*g_b` units gives a matrix scheme
with the improved term upper bound

```
r_next=C_b+(k^b*R_s-g_b)*r,
u_next=d^b*u.                                               (2)
```

Zero or duplicate terms may be simplified afterward. Such simplification is
not needed for (2) and must not be confused with finding the minimum rank.

## Positivity and asymptotic growth

Every flattening of the original certificate gives the necessary inequality

```
m <= (k-d^2)*f_i(S) <= (k-d^2)*R_s.
```

Since `m>0` and `S!=0`, this entails `k>d^2>d`. In particular

```
R_s-m/(k-d)>0,
Q_b=k^b*R_s-g_b
   =k^b*(R_s-m/(k-d))+m*d^b/(k-d)>0.
```

The powered flattening also gives
`g_b <= (k^b-d^(2b))*f_i(S)`, hence
`Q_b >= d^(2b)*f_i(S)>1`. Thus the usual geometric recurrence formula applies.

For a target exponent `tau>log_d(k)`, choose a fixed b satisfying
`Q_b<d^(b*tau)`. Such a b exists: `log_d(Q_b)/b` tends to `log_d(k)`.
For a rational target `tau=A/B`, the exact integer test is

```
Q_b^B < d^(b*A).
```

With that b fixed, (2) gives

```
r_t <= Q_b^t*(r_0+C_b/(Q_b-1)),
u_t=u_0*d^(b*t).
```

It reaches the requested strict rank inequality in finite time. Explicit
integer comparisons, rather than floating-point logarithms, should decide
when the bound is reached.

All witness matrices and the initial schemes lie in one fixed supplied
coefficient field. Products, sums, and scalar-output elimination stay in that
field. If its extension degree is e, descend the final matrix scheme once;
the condition becomes `e^2*r_t < u_t^tau`. Do not charge e squared again at
each iteration.

## Finite control values

For the conservative coordinate certificate `d=2,k=9,m=1`:

```
g_2=9+2=11, K=81.
```

With S a scalar unit, a scalar matrix seed, and an empty D, the actual
`compile_gain_powered_catalyst` and `compile_gain_absorption` arrays produce
`81-11=70` matrix terms. With a cross-mixed scalar D, the generated powered
catalyst scheme has `C_2=8+9=17` terms; absorption produces `17+70=87` terms.
Both final size-four matrix coefficient identities passed exact checks over
F_2. The gain-free powered path generated 81 terms in the empty-D case.

These are verified compiler control values for conservative witnesses. The
70-term scheme does not beat the naive 64-term size-four matrix scheme and
does not exhibit a useful near-9/4 catalyst. No Lean theorem was added or
built for this derived route.

The tests also retain 24 separate gains for `d=2,k=10,m=2` over F_4, even
though 24 embeds as zero in that field. The resulting 94-term F_4 scheme
descends once to 376 prime-field terms. An additional independent check used
D of shape `(1,2,1)` and S of shape `(2,1,3)`, both with non-prime-field
coefficients. The powered source and target shapes were `(175,107,256)` and
`(56,53,72)`. This check explicitly raised the dense tensor limit using
`CompilerLimits(max_tensor_coefficients=5_000_000)`; the source exceeds the
default two-million-coefficient cap. Exact maps and the final 185-term
size-four scheme verified,
matching `34+(81*2-11)`.

## Reviewed automatic planning interface

`reference/gain_planning.py` implements the integer criterion above, the
recurrence `r_next=C_b+Q_b*r`, and one final degree-squared field descent.
Its numeric plan is separate from coefficient generation: a plan may exceed
the dense compiler caps, in which case construction is reported incomplete.
The final tensor and scheme sizes are checked before powered compilation.

For a fixed b, `Q_b^B<d^(b*A)` is sufficient even when `k^B>=d^A`; the
actual retained gains improve that finite block. The implementation correctly
accepts this case. In particular, the conservative `d=2,k=9,m=1` scalar
example has `Q_1=8`, so the target `31/10` is possible despite the old k-based
criterion failing. This control test does not establish exponent 9/4.

Independent review ran the compiler, gain-planner, and finite-type test
modules together: 34 tests passed. A subsequently added empty-axis metadata
regression also passed separately. All verification here concerns actual
finite coefficients and Python behavior, not formal verification of the
asymptotic argument.
