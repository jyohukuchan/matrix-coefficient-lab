> This earlier slice-space audit is superseded by the
> [replicated rank proof](mixed-dot-rank-catalysts.md), which excludes
> every finite catalyst. Its narrower claims and former unresolved boundary
> are retained below as intermediate results.

# Mixed-orientation dot catalyst audit

This is an independent mathematical audit, not a Lean theorem or an exact
certificate. It concerns the supplied auxiliary

\[
S=\operatorname{DotX}_2\oplus\operatorname{DotY}_2\oplus
\operatorname{DotZ}_2,
\]

of shape \((5,5,5)\), where the named leg is the singleton leg of each
length-two dot tensor. Its supplied decomposition has six terms. The
question is whether the exact restriction

\[
D\oplus\mathrm{Unit}\oplus(M_2\otimes S)
\ \leq\ D\oplus5S                                      \tag{1}
\]

can exist. Any larger positive gain restricts to gain one, so excluding
(1) excludes all positive gains.

## Audited result: two differently oriented dots

For

\[
D=\operatorname{DotX}_2\oplus\operatorname{DotY}_2
\quad\text{(shape \((3,3,4)\))},
\]

(1) is impossible over every field. By leg permutation, this also excludes
the direct sum of any two distinct singleton-leg orientations. This result
does not assert an exclusion for arbitrary mixtures of dot tensors.

We may extend scalars to the algebraic closure: an exact restriction over
the original field remains exact there. The following linear-space and
irreducibility arguments therefore take place over an infinite
algebraically closed field, including characteristic two.

### Slice geometry

A first-leg slice of one source copy of \(S\), with scalar coordinate
\(s\) and two length-two vector coordinates \(u,v\), has rank

\[
2\mathbf1_{s\ne0}+\mathbf1_{u\ne0}+\mathbf1_{v\ne0}.
\]

The first-slice rank of \(D\) is
\(2\mathbf1_{d\ne0}+\mathbf1_{w\ne0}\), with scalar \(d\) and vector
\(w\in K^2\), and has generic rank three.

The three disjoint blocks in a first slice of the matrix-product core
\(M_2\otimes S\) are parametrized by

\[
A\in\operatorname{Mat}_{2\times2},\qquad
B\in\operatorname{Mat}_{4\times2},\qquad
C\in\operatorname{Mat}_{2\times4}.
\]

Their total slice rank is

\[
4\operatorname{rank}A+2\operatorname{rank}B+
2\operatorname{rank}C,                                   \tag{2}
\]

with generic value 16. Together with \(D\) and the gain scalar \(g\),
the target first-slice rank is generically 20.

The generic rank-drop cone of the target is the finite union of the
components

* \(g=0\);
* \(d=0\);
* \(w=0\), a codimension-two linear space;
* \(\det A=0\);
* \(\operatorname{rank}B\le1\);
* \(\operatorname{rank}C\le1\).

For a matrix linear space of rank at most one, all matrices have either a
common one-dimensional image or a common one-dimensional row space.
One elementary proof starts with a nonzero rank-one matrix \(xy^T\):
every other \(uv^T\) must share its image or row direction, since otherwise
a generic linear combination has rank two. Mixing the two possibilities
again gives a rank-two combination. Thus the maximal dimensions of such
spaces in the three matrix blocks are respectively 2, 4, and 4. Their
codimension losses in the full matrix spaces are respectively 2, 4, and 4.

### The pair-normal lemma

An exact restriction induces a first-leg pullback. Let \(c_1,\ldots,c_5\)
be the five source \(S\)-scalar coordinates pulled back to the target
first-input space. Killing any two of them bounds the source slice rank
by

\[
3+3\cdot4+2\cdot2=19.
\]

The same target slices therefore have rank at most 19. Their linear
space, of codimension at most two, lies in the target rank-drop cone.
An irreducible linear space contained in a finite union of closed sets
lies in one member. The \(B,C\) components cannot contain such a space
because their maximal linear spaces have codimension four. Consequently
every pair-normal span \(\langle c_i,c_j\rangle\) either contains
\(g\), contains \(d\), is the entire pure \(w\)-normal plane, or is a
pure \(A\)-normal plane whose annihilator is a rank-one matrix plane.

Here and below, "pure" means supported on the indicated direct summand.
This conclusion also covers dependent pairs: their kernels have still
smaller codimension, so only a scalar-hyperplane component can contain
them.

### Classifying the five normals

Put \(Y=\langle g,d\rangle\), and let
\(W=\langle c_1,\ldots,c_5\rangle\). If \(W\) has quotient dimension
at most one modulo \(Y\), then the projection of \(\ker W\) onto the
20-dimensional core first space has codimension at most one. Its
projections onto the \(A,B,C\) matrix blocks have dimensions at least
3, 7, 7. None is contained in the corresponding rank-drop cone, so the
generic core slice rank remains 16. But killing all five source scalars
bounds the source slice rank by \(3+10=13\), a contradiction.

Suppose instead that two \(c_i\) have independent classes modulo \(Y\).
Their pair cannot contain a scalar axis. By the pair-normal lemma, they
are pure normals in one common two-dimensional block: either the
\(w\) block or one \(A\)-normal rank-one plane. Pairing each remaining
normal with these two shows that every remaining nonzero quotient class
is pure in this same block. Any remaining zero quotient class is pure
\(g\) or pure \(d\), since pairing it with a nonzero pure-block normal
must contain one of those scalar axes.

If the common block is \(w\), then \(\ker W\) leaves the entire matrix
core intact. Its generic slice rank is at least 16, again exceeding 13.

In the remaining case, \(\ker W\) leaves a rank-one \(A\)-plane,
unrestricted \(B,C\), and unrestricted \(w\). Its generic contribution
from (2) and \(w\) is \(12+1=13\). The source bound of 13 forces both
scalar axes \(g,d\) to vanish on this common kernel. Thus the restricted
target pencil has generic rank exactly 13.

Let \(t\) be the source \(D\)-summand's singleton-first scalar pulled
back to this kernel. It cannot be zero identically: that would reduce
the source bound to \(1+10=11\). Killing this nonzero \(t\) therefore
gives a hyperplane on which source slice rank is at most 11.

The restricted target pencil has no rank-drop hyperplane. Its rank-drop
components are \(A=0\) inside the two-dimensional rank-one plane,
\(w=0\) inside the two-dimensional vector block, and the unrestricted
\(B,C\) rank-drop components. All have maximal linear codimension at
least two. A hyperplane cannot lie in their finite union. Its generic
target slice rank remains 13, contradicting the source bound of 11.

This completes the all-field exclusion.

## Stronger audited result: an orientation count at most two

Let

\[
D=x\operatorname{DotX}_2\oplus y\operatorname{DotY}_2\oplus
z\operatorname{DotZ}_2,
\qquad x,y,z\in\mathbb N.
\]

Then (1) is impossible over every field whenever
\(\min(x,y,z)\le2\). In particular, arbitrary direct sums using at most
two singleton-leg orientations are excluded. The two-dot result above
is a special case. This theorem does not claim an exclusion when all
three counts are at least three.

Permute the legs so that \(x\le2\), and put \(v=y+z\). The generic
first-slice rank of \(D\) is \(2x+v\). Its first pencil has \(x\)
scalar blocks of weight two and \(v\) length-two vector blocks of
weight one. The target generic first-slice rank is \(2x+v+17\).
Killing two source auxiliary scalar forms bounds source rank by
\(2x+v+16\), so the same pair-normal argument applies.

Put \(Y=\langle g,d_1,\ldots,d_x\rangle\). Every pair-normal span
either contains one actual scalar coordinate \(g,d_j\), is a pure
normal plane for one of the \(v\) target vector blocks, or is a pure
\(A\)-normal plane whose annihilator is rank one. The \(B,C\)
components still have maximal linear codimension four and cannot occur.
If two normals have independent classes modulo \(Y\), the pair lemma
forces all nonzero quotient normals into one common nonscalar plane.
The remaining normals, if any, are pure individual scalar axes. This
follows by pairing each with the two independent normals, just as in the
two-dot proof; it does not assume coordinate source maps.

Let \(J\) be the common kernel of all five auxiliary scalar forms.
The source slice rank on \(J\) is at most

\[
2x+v+10.                                                \tag{3}
\]

If the quotient-normal span modulo \(Y\) has dimension at most one,
the projection of \(J\) onto the joint matrix-core and target-vector
space has codimension at most one. All three matrix projections retain
full generic ranks, and each length-two target vector projection remains
nonzero. Thus the target rank is at least \(16+v\), greater than (3)
for \(x\le2\).

If the common nonscalar plane is a target vector block, the entire
matrix core survives and at most that one vector block is killed. The
target generic rank is at least \(16+v-1\), also greater than (3).

It remains to consider a common \(A\)-normal rank-one plane. Then
\(J\) leaves unrestricted target vector blocks and unrestricted
\(B,C\). Write \(p\) for the number of target \(D\)-scalar axes
vanishing on \(J\), and \(q\in\{0,1\}\) for whether the gain axis
vanishes. Since all scalar normals in this branch are individual axes,
the surviving scalar coordinates are independent. The target generic
rank on \(J\) is exactly

\[
12+v+2(x-p)+(1-q).
\]

Comparison with (3) requires \(2p+q\ge3\).

* For \(x=0\), this is impossible.
* For \(x=1\), it forces \((p,q)=(1,1)\). The target generic rank
  is \(12+v\). The source singleton-first \(D\)-scalar pullback
  cannot vanish identically on \(J\), since that would reduce (3) by
  two. Its kernel is a hyperplane costing two source rank. But the
  target pencil on \(J\) has no rank-drop hyperplane: all surviving
  rank-drop components have maximal linear codimension at least two.
  This is a contradiction.
* For \(x=2\), the possible pairs are
  \((p,q)=(1,1),(2,0),(2,1)\). The first two are handled by the
  single remaining scalar coordinate, as follows. The final pair is
  handled by concision.

For \((p,q)=(1,1)\), the target generic rank is \(14+v\), equal to
the source bound. For \((p,q)=(2,0)\), it is \(13+v\), one below
the source bound. Each of the two source singleton-first \(D\)-scalar
pullbacks is nonzero on \(J\); if one were zero the source bound
would already be \(12+v\). Killing either nonzero form costs two
source rank, and therefore forces a target rank-drop hyperplane. The
only such component is the one surviving target scalar axis. Both
source forms must therefore be proportional to that same scalar
coordinate. Killing them both makes the source bound \(10+v\),
while the surviving target matrix and vector blocks have generic rank
\(12+v\). This contradiction covers both pairs.

For \((p,q)=(2,1)\), all target scalar axes vanish. Restricting the
\(A\) input to a rank-one two-plane leaves the core's other two
flattening dimensions either \((20,16)\) or \((16,20)\). Here the
\(A\)-block contributes \((8,4)\) or \((4,8)\), while the two
unrestricted other core blocks contribute \((4,8)\) and \((8,4)\).
These contributions occupy direct summands.

On the source side, each of the five auxiliary copies with its first
scalar killed has other flattening dimensions at most \((3,3)\).
The two source \(\operatorname{DotX}_2\) summands together contribute
\((4,4)\). Thus, after adding the common contributions \((b,c)\)
of all \(\operatorname{DotY}_2,\operatorname{DotZ}_2\) summands, where
\(b=y+2z\) and \(c=2y+z\),
the source other flattenings are bounded by \((b+19,c+19)\).
The target has either \((b+20,c+16)\) or \((b+16,c+20)\).
One target flattening exceeds the corresponding source flattening,
which is impossible under a restriction. This completes the proof.

The source inequalities here do not presume that its first pullback is
onto any direct summand. All its first slices on \(J\) have auxiliary
singleton-first coordinate zero, so the restricted source tensor factors
through the corresponding \((4,3,3)\) auxiliary tensor in each copy.
Every source \(D\) block contributes at most its original flattening
dimension. These yield the stated upper bounds. By contrast, the target
blocks on \(J\) occupy independent second- and third-coordinate
summands and have the stated actual flattening dimensions. This is the
concision comparison needed by an arbitrary exact restriction.

### Canonical coefficients and executable controls

A canonical \(D\) is constructed by successive direct sums of blocks
with supports \((0,i,i)\), \((i,0,i)\), or \((i,i,0)\), for
\(i=0,1\), translating each block into disjoint coordinates. With
\(n=x+y+z\), its concise shape is

\[
(x+2y+2z,\ 2x+y+2z,\ 2x+2y+z)
=(2n-x,\ 2n-y,\ 2n-z).
\]

Thus counts inferred from a canonical tensor's shape are
\(x=2n-\dim_1\), \(y=2n-\dim_2\), \(z=2n-\dim_3\), with
\(n=(\dim_1+\dim_2+\dim_3)/5\). Shape alone does not identify an
arbitrary tensor as this direct sum; recognition must also check its
actual coefficients or a certified isomorphism.

Independent one-off exact coefficient controls checked every canonical
count triple \(0\le x,y,z\le3\) over characteristics two and three:
128 flattening-concision checks agreed with the formula. The matrix
core was separately constructed from the eight \(M_2\) coefficient
triples and six \(S\) triples. Restricting its \(A\) coordinates to
\(\langle E_{00},E_{01}\rangle\) gave flattenings \((18,20,16)\);
restricting to \(\langle E_{00},E_{10}\rangle\) gave
\((18,16,20)\), in both characteristics. These four additional checks
confirm the two concision normal forms. They are controls, not a
finite-characteristic substitute for the all-field proof above.

## Remaining scope

Only the region \(x,y,z\ge3\) remains outside the theorem above.
There, a common first-normal plane can coexist with three removed
singleton-first \(D\)-scalar axes and the retained gain. Alternating
slice arguments in the other legs introduces further normal-plane
branches. No all-mixture conclusion is claimed, and no numerical or
finite-field sample is used to replace those missing arguments.
