# Same-field geometric positive-gain certificates

This note independently audits the main agent's geometric direct-sum
construction. It confirms that, assuming the repository's exact-rank
exponent theorem, positive-gain certificates at `d=2,k=5` exist already over
`F₂`, with `D=0` and a specified direct-sum family of auxiliaries. This is an
existence and certificate-compilation reduction; no useful matrix scheme or
new catalyst coefficients have been extracted here. No Lean files were
changed or built.

## Conditional construction from one supplied exact scheme

Let `K` be any field, `d,k,b` positive integers, and let a supplied exact
scheme for `M_(d^b)` have `r` terms with `r<k^b`. The supplied term count may
be an upper bound; it need not be the least rank. Define

```
S_b = ⊕_{i=0}^{b−1} k^(b−1−i) copies of M_(d^i),
m   = k^b−r > 0,
D   = 0.
```

Here `M₁` is one scalar tensor and every displayed integer counts independent
direct-sum copies. No integer is interpreted as a possibly vanishing field
scalar. Multiplying by `M_d` and using the explicit coordinate identification
`M_d⊗M_(d^i) ≅ M_(d^(i+1))` gives

```
k S_b     = k^b Unit ⊕ ⊕_{i=1}^{b−1} k^(b−i) M_(d^i),
M_d⊗S_b   = M_(d^b)  ⊕ ⊕_{i=1}^{b−1} k^(b−i) M_(d^i).
```

The middle blocks match exactly. Among the `k^b` source scalar blocks, use
`r` blocks to generate `M_(d^b)` by the supplied rank scheme, and retain the
remaining `m` blocks as output scalar gains. Match the middle blocks by
identity maps and coordinate permutations. Therefore actual local linear
maps over `K` give

```
m Unit ⊕ (M_d⊗S_b) ≤ k S_b.
```

For a rank scheme `(a_j,b_j,c_j)`, the scalar-to-matrix part of the three
maps has precisely those vectors as its columns. Its remaining gain columns
are coordinate vectors, and the middle-block columns are identity columns.
Every source block maps into one target block on all three legs, so the
coefficient identity follows block by block. There are no unaccounted cross
terms. All coefficients remain in the supplied scheme's field.

This proof makes no tensor-rank direct-sum additivity assumption. In the
boundary case `b=1`, the middle sum is empty and it reduces to the usual
rank-scheme restriction plus unused scalar terms. It does not require an
algebraically closed field, interpolation nodes, or a field extension.

## Why the exact exponent supplies a power-of-d block

Let `d≥2` and put `c=log_d k`. Assume `ν_K<c`, where `ν_K` is the infimum
of the finite exact-rank exponents. Choose `τ` with `ν_K<τ<c`. The infimum
witness theorem supplies an integer `n≥2` and an exact scheme over **K**
with term count `r<n^τ<n^c`. This gives a strictly positive gap

```
ρ = r/n^c < 1.
```

For every integer `t≥1`, tensoring this scheme gives an exact scheme for
`M_(n^t)` with at most `r^t` terms. Define `b_t` by integer comparisons as

the largest integer satisfying `d^b_t≤n^t`. Equivalently,
`b_t=floor(t log_d n)`. Restrict the coefficient vectors to the upper-left
`d^b_t` matrix coordinates on each leg; zero padding of the two input
matrices and selecting the output submatrix shows that the result is an
exact scheme for `M_(d^b_t)` over the same field. Its rank upper bound is
`r^t`.

Since `b_t>t log_d n−1`,

```
r^t / k^b_t < k (r/n^c)^t = k ρ^t → 0.
```

For sufficiently large finite `t`, we have `b_t≥1` and `r^t<k^b_t`. Thus a
power-of-d block satisfying the conditional construction exists. Once a
seed `(n,r)` is actually supplied, the correct `t,b_t` can be found using
only integer powers and strict comparisons; floating-point logarithms are
not needed by a certificate generator. There is no uniform effective bound
on the initial seed provided by the infimum theorem.

The submatrix restriction direction is important: it uses `d^b_t≤n^t`.
Rounding upward and attempting to restrict a smaller matrix algorithm to a
larger matrix would not justify this argument.

## Repository theorem route

The following existing sources supply exact coefficient statements, rather
than identities only on finite-field input points:

- `lean/OAI/LinearAlgebra/MatrixMultiplication/AuxiliarySeparation/Main.lean:36`
  contains `exactRankExponent_le_nine_quarters_allFields`; its statement is
  `exactRankExponent F≤9/4` for every field `F`.
- `Arithmetic/RankExponent.lean:51` defines the least exact rank by `Nat.find`,
  and line 57, `exactRank_spec`, supplies an actual `Tensor.RankAtMost`
  decomposition at that rank. Here and below the arithmetic paths are
  relative to the same `AuxiliarySeparation` directory.
- `Arithmetic/RankExponent.lean:219`,
  `exists_exactMatrixRank_lt_rpow`, obtains the strictly smaller finite
  exact-rank block whenever the infimum is strictly below the requested
  exponent. This is the existential seed used above.
- `Arithmetic/RankExponent.lean:69`, `exactRank_restrict_le`, and line 91,
  `exactRank_product_le`, supply restriction monotonicity and exact
  product-decomposition upper bounds.
- `Arithmetic/FieldExtension.lean:20`,
  `rankAtMost_matrixMultiplication_of_power`, identifies a tensor power of a
  matrix tensor with the matrix tensor whose dimension is the corresponding
  integer power.

The geometric direct-sum construction and floor-size specialization in this
note have been mathematically audited; they are not claimed as newly
kernel-checked Lean lemmas. The all-field exponent theorem internally uses
algebraic-closure descent, but its conclusion concerns ranks **over the base
field itself**. Using that conclusion and its infimum witness therefore
requires no extension in the resulting existential scheme.

For `d=2,k=5`, the strict inequality needed here is `9/4<log₂5`, witnessed
without numerical approximation by `2^9=512<625=5^4`. Consequently the
repository theorem implies some finite exact binary scheme

```
M_(2^b) with r<5^b,
```

and hence a positive-gain `D=0` certificate over `F₂` in the displayed
`S_b` family. The same statement holds over each other field separately.
This is stronger than the earlier algebraically-closed catalyst-existence
argument for purposes of a same-field search space.

## Complete in principle, but not a practical size bound

Over `F₂`, for each fixed `b`, exhaustive enumeration of the three arrays of
a rank-`(5^b−1)` scheme is finite. Any scheme with fewer terms can be padded
with zero terms. Coefficientwise verification is decidable. Enumeration for
unbounded `b` is therefore complete in principle, and every success compiles
into the geometric catalyst. No extension degree, arbitrary `D`, or arbitrary
auxiliary support needs to be included in this theoretical search class.
The rank-block search is already the unresolved coefficient-extraction
problem: converting a successful block to a catalyst does not discover that
block or explain how to find it cheaply.

The directly compiled certificates are immediately large. Each of the
auxiliary's three axes has dimension

```
s_b = Σ_{i=0}^{b−1} 5^(b−1−i) 4^i = 5^b−4^b.
```

The previously source-audited Alder–Strassen theorem gives
`rank(M_n)≥2n²−1` over every field: matrix multiplication is multiplication
in the associative unital simple algebra `Mat_n(K)`, of dimension `n²` and
one maximal two-sided ideal. Thus the supplied-scheme condition `r<5^b` for this direct compilation is
ruled out for the first three blocks:

| b | Matrix dimension | External lower bound `2·4^b−1` | Strict target `r<5^b` | Auxiliary axis size |
| --- | --- | --- | --- | --- |
| 1 | 2 | 7 | r<5 | 1 |
| 2 | 4 | 31 | r<25 | 9 |
| 3 | 8 | 127 | r<125 | 61 |
| 4 | 16 | 511 | r<625 | 369 |

The row for `b=4` only survives this lower bound; it is not a supplied scheme
or an existence claim specifically at `b=4`. The lower bound is the external
classical theorem documented in `catalyst-obstructions.md`, not a newly
formalized Lean result. Its algebra assumptions hold in positive
characteristic as well.

Already at the smallest block not excluded by this supplied-scheme test,
the source `5S_b` has axis size
1,845 and a dense coefficient cube of `1,845³=6,280,426,125` entries. It
exceeds the current two-million-entry dense verifier budget by a large
margin. The direct-sum support is sparse, so this count concerns that dense
backend, not an intrinsic need to materialize every zero coefficient. A
sparse or symbolic block compiler could verify a *supplied* scheme without
forming that cube; it would still need the missing good matrix scheme.
Larger blocks only increase these dimensions. This result narrows a complete
theoretical search class and clarifies field requirements; it does not make
an unbounded exhaustive coefficient search practical.


### Separate bounds for arbitrary mixing on the same auxiliaries

The `b≥4` conclusion above concerns **direct compilation from a supplied
rank-`r<5^b` scheme**. It does not rule out every possible mixed restriction
map on `S_b` when `b≤3`: the stated construction is a sufficient method, not
a necessary normal form for all maps on that support.

Two small members are independently excluded by the earlier obstructions.
`S₁=Unit` is the scalar case. For `S₂=5 Unit⊕M₂`, use its explicit rank upper
bound twelve and identify it with multiplication in
`A=K⁵×Mat₂(K)`. This algebra has dimension nine and six maximal two-sided
ideals. The algebra filter from `catalyst-obstructions.md` requires

```
5 rank_upper(S₂) ≥ 1+8 dim(A)−t(A),
60 ≥ 67,
```

which fails. This excludes arbitrary catalyst maps on `S₂`, including an
extra fixed `D`, by the repeated-target lower-bound argument.

For `S₃=25 Unit⊕5M₂⊕M₄`, Strassen's supplied seven-term scheme and its square
give an upper certificate of `25+5·7+49=109`. The corresponding algebra has
dimension sixty-one and thirty-one maximal two-sided ideals. The same filter
only gives `545≥458`, which passes; this note makes no impossibility claim
for mixed maps on `S₃`. Its dense source cube already has
`305³=28,372,625` entries, beyond the current backend cap. The direct
scheme-to-catalyst compiler still cannot use `b=3`, because that would
require the independently impossible `rank(M₈)<125`.
