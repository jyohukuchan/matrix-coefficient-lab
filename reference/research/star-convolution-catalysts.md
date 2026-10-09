# Independent exclusion of convolution catalysts for the singular star

For every field, neither full convolution D of concise shape `(3,4,6)` nor
its highest-output-degree truncation of shape `(3,4,5)` can satisfy an exact
restriction

```
D direct_sum m Unit direct_sum (M_2 tensor S_star)
    <= D direct_sum 5 S_star,                    m >= 1.
```

Here S_star is the `(2,3,3)` singular star pencil. These two exclusions are
mathematical derivations proposed by the finite-search/main agents and
independently audited by the constructive-research agent. They are not Lean
theorems, and no Lean build or proof-source edit was made for this report.
They are independent of SAT timeout outcomes. No catalyst maps were found.

A subsequent independently audited combination of these arguments gives a
stronger necessary condition: if the first flattening rank of D is three,
its generic first-slice rank must be at least five. Thus both other
flattening ranks must also be at least five. The componentwise minimal
remaining concise-shape boundaries are now `(3,5,5)`, `(4,4,5)`, and
`(4,5,4)`. These boundaries describe candidates still needing analysis, not
existing catalyst witnesses. The section on first dimension three proves
this strengthening; later sections give additional rank-drop conditions
and explain why a fixed-star witness would improve the original 9/4 bound.

The full convolution definition agrees with `convolution_tensor(field,3,4)`
in `reference/sectors.py:29`: coefficient `(i,j,k)` is one exactly when
`i+j=k`. The truncation deletes output coordinate k=5. The first-concision,
slice-rank, and scalar-extension reductions used below are also detailed in
`reference/research/star-nonzero-catalyst.md`.

## Exact slice rank-drop loci

Write the first polynomial as `delta_0+delta_1*t+delta_2*t^2`. Acting on
the four coefficients of a polynomial of degree at most three gives the
following six-by-four slice matrix, with output coordinates as rows:

```
delta_0       0        0        0
delta_1 delta_0        0        0
delta_2 delta_1  delta_0        0
      0 delta_2  delta_1  delta_0
      0       0  delta_2  delta_1
      0       0        0  delta_2
```

A nonzero polynomial multiplies injectively in the polynomial ring over
any field. Therefore the full matrix has rank four for every nonzero
delta, and its rank-drop locus is exactly `{delta=0}`, of codimension three.

For the truncated five-by-four matrix, if delta_0 is nonzero, rows 0 through
3 have determinant delta_0^4. If delta_0=0 but delta_1 is nonzero, rows 1
through 4 have determinant delta_1^4. If both vanish, the remaining
delta_2-only matrix has rank three for delta_2 nonzero, and rank zero at
the origin. Consequently its exact rank-drop locus is

```
{delta_0=delta_1=0},
```

a single linear subspace of dimension one, not a larger exceptional variety.
Both tensors have generic first-slice rank four and first flattening rank
three. All statements remain true in every characteristic.

An independent bounded check enumerated every delta over F_2,F_3,F_4,F_5.
The full rank distributions were respectively `{0:1,4:7}`, `{0:1,4:26}`,
`{0:1,4:63}`, `{0:1,4:124}`. Truncated distributions were `{0:1,3:1,4:6}`,
`{0:1,3:2,4:24}`, `{0:1,3:3,4:60}`, `{0:1,3:4,4:120}`. This corroborates
the explicit minors; the all-field conclusion follows from the algebraic
argument, not this finite enumeration.

## Restriction and low-rank geometry

Discard surplus target scalar units to reduce to m=1. Extend any purported
exact coefficient restriction to an algebraic closure. Linear ranks and the
coefficient identity are preserved. The following dimension arguments take
place over this algebraically closed field.

The target first space V has coordinates

```
(delta_0,delta_1,delta_2, epsilon, A, B),
dimension(V)=3+1+8=12.
```

The source first space has three D coordinates and five star pairs, hence
dimension thirteen. Target first conciseness forces the first covector
pullback F from V to that source space to be injective. Let Q_i be its
two-row projection to source star pair i, and let `H_i=kernel(Q_i)`.
Then `dimension(H_i)>=10`.

On H_i, the source slice rank is at most `4+4*2=12`. Thus every target slice
on H_i has rank at most twelve, independently of any mixing in the other
two restriction maps. The target slice has rank

```
rank(D(delta)) + indicator(epsilon != 0)
 + 2 rank(vertical_stack(A,B)) + 2 rank(horizontal_concat(A,B)).
```

The core has maximum rank eight. Its lower-rank locus is the union of the
two rectangular rank-at-most-one loci, each of dimension five in the eight
core coordinates. Let B_bad denote this union. The target low-rank locus is
therefore exactly the union of

```
E       = {epsilon=0},                         dimension 11,
Z_D     = {D(delta) has rank < 4} x K x K^8,
Z_core  = K^3 x K x B_bad,                    dimension 9.
```

For full convolution, dimension(Z_D)=9; for the truncation it is ten.
Every H_i is an irreducible linear space. If it lies in a finite union of
closed varieties, it lies in one of them. It cannot lie in a variety whose
dimension is smaller than its own. In particular it cannot lie in Z_core.

## Full convolution

All five H_i must lie in E, since the other two components have dimension
nine. Equivalently, the epsilon form belongs to the row space of every Q_i.
Each Q_i has at most two independent rows, so their combined row space has
dimension at most `1+5=6`: one common epsilon form and at most one additional
form per pair. The three source-D rows add at most three dimensions.

The row rank of F is thus at most nine, contradicting its injective rank
twelve. Full convolution D of shape `(3,4,6)` is excluded.

## Truncated convolution

Each H_i is either contained in E (type A), or contained in Z_D (type B).
For type B, both H_i and Z_D have dimension ten, so they are equal:

```
H_i = {delta_0=delta_1=0},
rowspace(Q_i) = span(delta_0,delta_1).
```

If every pair is type A, the previous row-rank bound of nine applies. If
some pair is type B, there are at most four type-A pairs. All type-B rows
lie in one fixed two-dimensional space; the type-A rows lie in the span of
epsilon and at most four additional forms. The five pairs together have
row rank at most `2+1+4=7`. Including the three source-D rows gives rank
at most ten, still smaller than twelve. The truncated `(3,4,5)` tensor is
also excluded over every field.

## Audited conditional extensions

The same argument applies to a first-concise arbitrary D with first
dimension h and generic slice rank rho. Let Delta be its generic rank-drop
locus in K^h, defined by all rho-by-rho minors. Its definition and dimension
are over an algebraic closure; finite-field point samples alone do not
determine its algebraic codimension.

The target first dimension is `h+9`. Each H_i has dimension at least `h+7`,
and the rank-at-most-`rho+8` target locus is the union of the gain-zero
hyperplane, `Delta x K x K^8`, and `K^h x K x B_bad`. The last component
has dimension `h+6`.

If Delta has codimension at least three, its component also has dimension
at most `h+6`. Every H_i must then have type A. The total first pullback
row rank is at most `h+6`, less than `h+9`, a contradiction.

If Delta is contained in a single linear codimension-two subspace, any
type-B H_i must equal that subspace times the remaining coordinates, by
dimension. Its two pair rows therefore lie in one fixed two-dimensional
space of target-D forms. With some type B, the row rank is at most `h+7`,
again insufficient. If there are no type-B spaces, the type-A contradiction
applies. The rank-drop locus need not fill the bounding subspace for this
sufficient exclusion; if its dimension is smaller, type B is impossible.

There is a further valid strengthening specifically for h=3: any rank-drop
locus of codimension at least two is excluded, even if it has multiple or
nonlinear components. A type-B H_i has dimension at least ten and its
projection to target-D coordinates has dimension at most one. Its kernel
under that projection has dimension at least nine, which is exactly the
full kernel of the target-D projection. It follows that H_i has dimension
ten and contains that full nine-dimensional kernel. Thus every type-B pair
row lies in the three-dimensional space of target-D covectors.

If some pair is type B, at most four are type A. All pair rows together
then span at most `3+1+4=8` dimensions; the source-D rows add at most three,
giving at most eleven rather than the required twelve. This h=3 refinement
does not extend by the same bound to arbitrary h: the resulting estimate
would be `2h+5`, which need not be smaller than `h+9` for h>=4.

These conditional lemmas exclude families of D, not all possible catalysts.
For example, a generic-rank-drop hypersurface has codimension one, so none
of these rank-drop criteria alone excludes it. The additional necessary
slice plane does, however, rule out all generic-rank-four cases when h=3,
as follows.

## Further strengthening: first dimension three requires generic rank five

The already audited nonzero-catalyst argument supplies a K-defined
two-dimensional subspace L of source-D first covectors such that every
nonzero slice in L, including after algebraic-closure extension, has matrix
rank at least four. This is the all-h slice-plane condition proved in
`reference/research/star-nonzero-catalyst.md`; it follows by taking the
kernel of the target-D/gain projection on the source-D part of the first
pullback image. In particular a purported catalyst cannot have generic
slice rank smaller than four.

Suppose now h=3 and the generic slice rank rho is exactly four. Work over
the algebraic closure. Its rank-drop set

```
R = {delta in Kbar^3 : rank D(delta) < 4}
```

is a closed homogeneous cone, defined by all four-by-four minors of the
linear slice matrix. At least one minor is a nonzero polynomial because
rho=4. Every slice has rank at most four over this field. Therefore the
necessary plane L has

```
R intersect L = {0}.
```

This forces dimension(R)<=1. To check the geometric step, any irreducible
cone component of affine dimension at least two contains the origin and
projectivizes to a positive-dimensional closed variety in projective
two-space. Such a variety meets every projective line: the projective
dimension inequality gives intersection dimension at least
`dimension(variety)+1-2>=0`. The projective line associated to L would
therefore contain a point of that component, contradicting the displayed
intersection. Equivalently, intersecting an affine cone component with
the linear hyperplane L lowers its dimension by at most one; a component
of dimension at least two cannot meet L only at the origin.

Thus the generic rank-drop locus has codimension at least two. The h=3
rank-drop exclusion proved above applies and contradicts existence of the
catalyst. Together with the plane's rank-at-least-four condition, this
excludes every h=3 catalyst with generic slice rank at most four, leaving
the necessary bound

```
h=3  ==>  generic first-slice rank rho>=5.
```

This uses algebraic dimension over an algebraic closure and the plane's
rank bound for all nonzero extended-field vectors. A rank bound tested
only on the original finite-field points would not suffice for the
projective intersection step. Exact coefficient restrictions supply the
needed extension-stable condition.

Generic slice rank is at most each of the other two flattening ranks after
concision. Hence with h=3 both are at least five. For h>=4 the existing
all-h condition still requires both other ranks at least four and excludes
the pair `(4,4)`. The revised minimal componentwise boundaries are exactly

```
(3,5,5), (4,4,5), (4,5,4).
```

In particular, every tensor with first flattening rank three and either
other flattening rank four is excluded, including every concise `(3,4,5)`
tensor, not only the truncated convolution example. A rank-six upper
decomposition for a `(3,4,5)` tensor does not circumvent this obstruction:
tensor rank and generic slice matrix rank are different quantities, and
the latter is bounded by four for that shape. A generic-rank-five tensor
may drop to rank four at intersections with L, so this particular argument
does not automatically extend the lower bound to six.

## Further obstruction: rank-drop must span the whole first space

There is another all-h necessary condition for an exact star catalyst:
the generic first-slice rank-drop locus Delta must not be contained in any
proper linear hyperplane of the first-coordinate space. Equivalently, its
linear span over the algebraic closure must be the whole first space.

Suppose instead Delta is contained in `{b=0}` for a nonzero linear form b
on the target-D first coordinates. The notation b in this section means
this fixed linear form, not a matrix dimension. Continue with h first
coordinates, generic slice rank rho, one scalar gain epsilon, and an
injective first pullback from the target's `h+9` coordinates. All geometric
sets and dimensions below are over an algebraic closure.

For each source star pair, its pullback Q_i has at most two independent
rows, and H_i=kernel(Q_i) has dimension at least `h+7`. On H_i the source
slice rank is at most `rho+8`. Consequently H_i lies in the target
rank-at-most-`rho+8` set, the union of

```
{epsilon=0},
Delta x K x K^8,
K^h x K x B_bad.
```

The last component has dimension `h+6`, so it cannot contain H_i. Since
H_i is irreducible, it lies in one of the first two closed sets. In the
first case epsilon annihilates H_i; in the second case b does. Thus the
row space of Q_i contains epsilon or b. The second assertion needs no
codimension bound on Delta and no assumption that its components are
linear. It only uses its containment in the hyperplane.

Each pair row space therefore lies in the span of `{epsilon,b}` and at
most one additional form. Their five row spaces together span at most
`2+5=7` dimensions. The h source-D rows of the pullback add at most h,
so its rank is at most `h+7`, strictly smaller than the required `h+9`.
This proves the obstruction. If all pair spaces share the same one of
epsilon and b, the bound is even smaller, but this is not needed.

The scalar-extension argument is essential when the base field is finite:
Delta here is the homogeneous determinantal variety defined by the generic
maximal minors, not a set inferred from finitely many evaluated slices.
The generic rank and all exact restriction maps are unchanged by extension.

## Application to the `(4,4,5)` band-pencil candidate

Let `L_4(a,b)` denote the four-by-five matrix with a on positions `(j,j)`
and b on `(j,j+1)`, for j=0,1,2,3. The candidate

```
D(a,b,c,d) = L_4(a,b) + c E_00 + d E_11
```

has the slice matrix

```
a+c     b     0     0     0
  0   a+d     b     0     0
  0     0     a     b     0
  0     0     0     a     b
```

Its five maximal minors, in the order obtained by deleting columns
0,1,2,3,4, are

```
b^4,
(a+c)*b^3,
(a+c)*(a+d)*b^2,
(a+c)*(a+d)*a*b,
(a+c)*(a+d)*a^2.
```

These formulas are valid in every characteristic. If b is nonzero, the
first minor gives rank four. If b=0, the remaining matrix is diagonal on
its first four columns, with entries `(a+c,a+d,a,a)`. Its exact rank-drop
locus is therefore

```
{b=0,a=0} union {b=0,a+c=0} union {b=0,a+d=0}.
```

It consists of three distinct codimension-two planes in the four first
coordinates. Their union is not contained in one codimension-two plane,
so the earlier fixed-plane criterion does not directly apply. However all
three lie in the single hyperplane `{b=0}`. The newly audited obstruction
excludes this D over every field, for every positive gain count.

The first pullback would need rank thirteen. In this example the same
argument bounds the five star-pair rows by seven dimensions and the
source-D rows by four, giving rank at most eleven. The candidate does have
a two-dimensional full-rank-four plane `(c,d)=(0,0)`, so passing the earlier
minimum-rank plane filter is insufficient.

Independent enumeration of every `(a,b,c,d)` over F_2,F_3,F_4,F_5 checked
the displayed exact rank-drop condition. Slice rank histograms were
`{0:1,1:2,2:2,3:2,4:9}`, `{0:1,1:4,2:6,3:8,4:62}`,
`{0:1,1:6,2:12,3:18,4:219}`, and
`{0:1,1:8,2:20,3:32,4:564}`. This computation confirms the minors at
these points; the exclusion is the all-field linear-algebra and dimension
proof above, independently of bounded SAT outcomes.

## Diagonal first-dimension-three catalysts are also excluded

For

```
D(x,y,z) = diag(x, x+z, y, y+z, x+y, x+y+z),
```

the generic slice rank is six. The plane z=0 has slice rank at least four
for every nonzero `(x,y)` over the algebraic closure: its diagonal entries
are two copies each of x, y, x+y; at most one of these three values can
vanish for a nonzero pair. The generic rank-drop set is the union of the
six displayed linear hyperplanes in the first-coordinate space. It has
codimension one and spans that whole three-dimensional space.

The plane requirement, the h=3 generic-rank-five lower bound, the
codimension-based exclusions, and the hyperplane-containment exclusion
alone do not exclude this candidate. A further argument proposed by the
finite-search agent, then independently audited by the constructive and
existence-audit agents, does exclude it. In fact it excludes every diagonal
D whose first flattening rank is three, regardless of its generic slice
rank. No regular-pencil commuting-source argument is needed.

Remove identically zero diagonal entries first. Write a general such D as

```
D(delta) = diag(ell_1(delta), ..., ell_N(delta)),
```

where the nonzero linear forms ell_j span the three-dimensional source-D
dual space. Its generic slice rank is N; its target-D rank-drop locus is
the union of the hyperplanes ell_j=0. Repeated forms are allowed. Work over
an algebraic closure as before, reducing surplus gains to one.

Let V be the twelve-dimensional target first space. Let Y_star be the
four-dimensional subspace of V dual forms spanned by the three target-D
coordinate forms and the scalar-gain form epsilon. For each source star
pair, let W_i be its pullback row space, of dimension at most two. Its kernel
H_i has dimension at least ten. On H_i the source slice rank is at most
`N+8`, so the target slice lies in the union of epsilon=0, the target-D
rank-drop hyperplanes, and the core bad locus times D/gain coordinates.
The last component has dimension nine and cannot contain H_i.

Irreducibility therefore implies that W_i contains a nonzero form u_i in
Y_star, either epsilon or one of the target-D ell_j. Define

```
P = sum_i W_i,
Z = rowspace of the three source-D rows of the first pullback.
```

Each W_i contributes at most one additional form beyond u_i. Thus
`dimension(P)<=4+5=9`, while `dimension(Z)<=3`. Target first conciseness
forces the full pullback row rank to be twelve. Consequently all bounds
are equalities:

```
dimension(P)=9, dimension(Z)=3, P intersect Z={0}.
```

Moreover the chosen u_i must span all of Y_star. Otherwise their span has
dimension at most three and P would have dimension at most eight. Hence
`Y_star subset P`, and in particular `Z intersect Y_star={0}`.

Choose two linearly independent source diagonal-entry forms ell_j and
ell_k, which exist because all entry forms span the source-D dual space
of dimension three. The source-D pullback has row rank three, so their
pullbacks are still independent elements of Z. Their simultaneous kernel
K in V has dimension ten. On K those two distinct source diagonal entries
vanish, so the source-D slice rank is at most `N-2`; adding all five star
blocks gives source rank at most `N+8`. The target slice again lies in
the same low-rank union.

The ten-dimensional irreducible K cannot lie in the nine-dimensional core
bad component. It must therefore lie in epsilon=0 or in one of the target-D
rank-drop hyperplanes. The corresponding nonzero normal belongs to Y_star
and annihilates K. But the annihilator of K is exactly the span of its two
chosen source-D pullback forms, a subspace of Z. This gives a nonzero
element of `Z intersect Y_star`, contradicting its being zero.

The proof works in every characteristic. For the concrete `(3,6,6)` D,
one can choose the source forms x and y; their vanishing actually leaves
only three nonzero diagonal entries, a stronger bound than needed. This
candidate is thus excluded over all fields, independently of SAT results.
The deduction requires D to have the stated common diagonal slice form;
commuting or upper-triangular slices alone do not justify the two-entry
rank reduction. Arbitrary non-diagonal D remains subject to the earlier
necessary conditions and is not globally excluded by this argument.

## General first-dimension-three refinement

The main agent subsequently derived the following stronger necessary
conditions for arbitrary, not necessarily diagonal, first-concise D with
h=3. The constructive and existence-audit agents independently checked the
linear-algebra, dimension, and coefficient-field steps.

Every nonzero D first slice over the algebraic closure must have matrix
rank at least four. Moreover its generic rank-drop cone must contain at
least three linear plane components whose normal forms are independent.
Those planes are defined over the coefficient field of an actual assumed
restriction witness. For a `(3,5,5)` D, its determinant must therefore have
at least three independent linear factors over that witness field.

To prove this, let Y_star again be the four-dimensional space of target-D
and gain forms, let W_i be the five source-star pair row spaces, and let
P be their sum. For arbitrary D the generic rank-drop cone Delta is a
proper closed subset of the three target-D coordinates. Each H_i=ker W_i
has dimension at least ten and lies in the usual target low-rank union.
The core-bad component has dimension nine and is impossible.

If H_i lies in the gain-zero hyperplane (type A), W_i contains epsilon.
Otherwise its projection to the target-D coordinates is a linear subspace
contained in Delta (type B). That projection cannot be the full space,
since Delta is proper, and hence has dimension at most two. A nonzero
target-D form annihilates it and lies in W_i. Thus in both cases W_i
contains a nonzero Y_star form u_i, and the same bounds as above give

```
dimension(P)<=4+5=9, dimension(Z)<=3,
dimension(P+Z)=12.
```

It follows that P has dimension nine, Z dimension three, their intersection
is zero, and the chosen u_i span all of Y_star. Consequently Y_star is a
subspace of P.

The source first space has dimension thirteen and the image H of the
injective pullback has dimension twelve. Its single normal eta has no
source-D component: the ten source-star forms restrict to a nine-dimensional
P, so their one-dimensional restriction kernel already supplies the entire
normal line of H. Therefore the whole three-dimensional source-D coordinate
space, with all star pairs zero, lies in H. Its target preimage annihilates
P and hence Y_star, so it has only core coordinates. Every nonzero such
vector yields a nonzero core slice of rank at least four. Slice-rank
monotonicity gives the asserted rank-at-least-four condition for every
nonzero source-D slice, including after scalar extension.

Since Y_star is contained in P, the quotient P/Y_star has dimension five.
Each W_i has dimension at most two and meets Y_star in at least one
dimension; its quotient contribution is therefore at most one. Five
contributions must produce dimension five, so all five W_i have dimension
two, their intersections with Y_star have dimension exactly one, and
their quotient contributions are independent.

For type B, a target-D projection of dimension at most one would yield
at least two independent target-D forms annihilating H_i, contradicting
the just-proved one-dimensional intersection. Its projection is therefore
exactly a two-dimensional plane contained in Delta. Because Delta is
proper in three dimensions, such a plane is an irreducible component of
Delta. Its normal spans W_i intersect Y_star and is a pure target-D form.
For type A the intersection instead is exactly the epsilon line; the two
types cannot coincide.

The five one-dimensional intersections span all four dimensions of
Y_star. Hence there is at least one type-A pair and at least three type-B
pairs whose plane normals span the three target-D dual coordinates. The
number of type-B pairs is three or four, so this assertion does not assume
all five normals are independent.

All these spaces, kernels, projections, and annihilators can be chosen
over the field of the actual restriction maps: linear ranks are invariant
under extension. This yields three independent planes over that coefficient
field. If an assumed witness uses an extension E of an original field K,
the argument asserts E-defined planes, not necessarily K-defined planes.
For exclusion without fixing the witness field, inspect plane components
over the algebraic closure; factorization only over a smaller base field
may miss factors that appear after extension.

For a `(3,5,5)` D, the preceding generic-rank lower bound forces generic
rank five. Its nonzero homogeneous determinant vanishes on each necessary
plane, so the corresponding independent linear forms divide it. Thus a
verified geometrically irreducible determinant, or one without three
independent geometric linear factors, rules out that D. A random matrix
net is not automatically excluded merely because it is called random;
the factor condition must actually be checked. Conversely, satisfying
these necessary conditions does not construct a catalyst.

## Calibration of the fixed star search target

The fixed auxiliary S_star has an explicit rank-three degeneration over
every original coefficient field, including characteristic two. With
first-leg basis a_1,a_2 and other-leg basis e_0,e_1,e_2, set

```
P(t) = a_1 tensor (e_0+t e_1) tensor (e_0+t e_1)
     + a_2 tensor (e_0+t e_2) tensor (e_0+t e_2)
     - (a_1+a_2) tensor e_0 tensor e_0
     = t S_star + t^2 E,

E = a_1 tensor e_1 tensor e_1 + a_2 tensor e_2 tensor e_2.
```

Constant terms cancel in every characteristic. This also gives an exact
coefficient-extraction bound for every n without interpolation nodes or
any field extension. Expand P(t)^n as a sum of 3^n rank-one polynomial
word terms. Each word has a constant first factor A_w and two vector
polynomials B_w(t), C_w(t), each of degree at most n. Since the coefficient
of t^n in P(t)^n is exactly S_star^n, the termwise formula is

```
[t^n] (A_w tensor B_w(t) tensor C_w(t))
  = sum_(i=0..n) A_w tensor [t^i]B_w(t) tensor [t^(n-i)]C_w(t).
```

It supplies at most n+1 rank-one terms per word, all over the original
field. Therefore

```
rank(S_star^n) <= (n+1)*3^n.
```

Independent checks constructed these actual coefficient schemes for
n=0,1,2,3 over F_2,F_3,F_4,F_5 and verified every tensor coefficient. For
n=3 the unpruned scheme supplied 108 terms, as the formula predicts. This
is a decomposition of powers of the auxiliary, not a new low-rank matrix
multiplication scheme.

For any normalized nonnegative, additive, multiplicative, restriction-
monotone character phi, the rank bound yields
`phi(S_star)^n <= (n+1)*3^n`, hence `phi(S_star)<=3`. A scalar unit restricts
from S_star, so its character value is positive. A real exact star witness
with k=5 and positive gain m would imply

```
phi(D) + m + phi(M_2)*phi(S_star)
    <= phi(D) + 5*phi(S_star),
phi(M_2) <= 5 - m/phi(S_star) <= 5-m/3 <= 14/3.
```

Subtraction here is ordinary subtraction of finite real character values;
it does not cancel tensors in the restriction preorder.

The character inequality alone should not be presented as an arithmetic
algorithm. To translate it to an exponent bound, use the detecting-character
bridge `exists_detecting_character` in
`Character/Existence.lean:39`, the product identity
`M_(2^j) isomorphic to M_2^j`, and the exact-rank-to-arithmetic bridge
`omega_le_exactRankExponent` in `Arithmetic/Exponent.lean:117`.
Explicitly, if the exact-rank exponent nu were greater than
`log_2(14/3)`, for large j there would be an integer k strictly between
`(14/3)^j` and `2^(j*nu)`. Detection for dimension 2^j gives a character
with value at least k, whereas the witness gives every character value
on that matrix tensor at most `(14/3)^j`, a contradiction. Under these
existing bridge inputs the witness would imply

```
omega <= log_2(14/3) = approximately 2.222392421 < 9/4.
```

This implication is conditional on an actual fixed-star witness. If its
maps are only supplied over a finite extension of an original base field,
passing the arithmetic bound back also uses the fixed finite-coefficient-
algebra descent bridge; its overhead must not be charged repeatedly.
The auxiliary coefficient bound above itself requires no such extension.
It is separate from the current compiler's supplied four-term S scheme
and its corresponding finite recurrence costs.

Consequently the original 9/4 theorem does not guarantee a gain-one,
`d=2,k=5` witness with this particular S_star. Its existence would establish
a strictly stronger exponent bound. This calibrates the search target; it
does not prove that such a witness is impossible. The existing
`normalizedStates_nonempty_of_no_catalyst` in
`Spectrum/StateObstruction.lean:283` assumes absence of catalysts for all
auxiliary s. Its contrapositive produces some auxiliary and cannot by
itself promise this fixed one, or a chosen direct sum of copies of it.
