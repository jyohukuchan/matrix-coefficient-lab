# Rank obstruction for separate oriented dots

**Status:** independently audited mathematical proof over every field;
not formalized in Lean. This supersedes the earlier slice-space exclusions
for dot-sum catalysts: the exclusion now covers **every finite tensor
catalyst**. It supplies a rank obstruction, not a normalized-state lower
row or a coefficient witness.

Let \(\operatorname{DotX}_a\), \(\operatorname{DotY}_b\), and
\(\operatorname{DotZ}_c\) be matrix tensors of ranks \(a,b,c\), with
their named leg singleton. For positive integers \(a,b,c\), put

\[
S=\operatorname{DotX}_a\oplus\operatorname{DotY}_b
  \oplus\operatorname{DotZ}_c.
\]

If \(D\) is any finite tensor, then no exact
restriction

\[
D\oplus m\mathrm{Unit}\oplus(M_2\otimes S)
\ \leq\ D\oplus5S                                      \tag{1}
\]

exists for any \(m\ge1\), over any field. In particular, this excludes
every finite catalyst for the three separate length-two dots of shape
\((5,5,5)\).

The proof uses ordinary exact tensor rank. Its ingredients are elementary
matrix-summand additivity, first-leg substitution, and the all-field
rectangular bound \(R(M_{2,2,k})\ge7k/2\). A self-contained proof of the
last bound and its direct-sum replication appear below, following the
argument restated in
[Jason Yang, *New lower bounds on tensor rank of (2,n,m) matrix
multiplication with GPT-6*, §2](https://arxiv.org/html/2609.14393v1#S2).
That section attributes the original argument to Yaroslav Shitov.

## Matrix summands add their exact rank

Let \(T\in X_T\otimes Y_T\otimes Z_T\) be arbitrary and let
\(E\in K\otimes Y_E\otimes Z_E\) be a matrix tensor of matrix rank
\(r\). Then

\[
R(T\oplus E)=R(T)+r.                                    \tag{2}
\]

For the lower bound, take any rank-\(R\) decomposition of \(T\oplus E\).
The \(Y_E\) components of its second factors span a space of dimension
at least \(r\): contracting its first leg to the scalar \(E\)-coordinate
and projecting its last two legs retrieves the rank-\(r\) matrix \(E\).
Select \(r\) terms with independent \(Y_E\) components. There is a
linear map \(L:Y_E\to Y_T\) sending each selected component to the
negative of that term's \(Y_T\) component. Extend the map linearly to
the remaining \(Y_E\) directions.

Apply the first-leg projection onto \(X_T\), the second-leg map
\((y_T,y_E)\mapsto y_T+L(y_E)\), and the third-leg projection onto
\(Z_T\). These maps retrieve exactly \(T\), and the selected \(r\)
terms vanish. Thus \(R(T)\le R-r\). The reverse inequality follows by
concatenating decompositions. This proof uses no characteristic
assumption and applies after any permutation of the three legs.

Consequently, if \(D\) is a finite direct sum of matrix tensors with
ranks \(r_1,\ldots,r_s\), then

\[
R(T\oplus D)=R(T)+\sum_i r_i.                           \tag{3}
\]

## First-leg substitution with a protected summand

Suppose \(T=T_0\oplus T_1\), and the first flattening of \(T_1\) has
rank \(b\). Then

\[
R(T)\ge R(T_0)+b.                                      \tag{4}
\]

First compress each summand to its actual first support. Write the
first support of \(T_0\) as \(X_0\), of dimension \(a\). The full
first support of \(T\) has dimension \(a+b\).

In a rank decomposition, its first factors span this support. If the
current first support has dimension greater than \(a\), some first
factor lies outside the protected copy of \(X_0\). Quotient by that
factor. The quotient is injective on \(X_0\), and can be identified with
a smaller first space while keeping \(X_0\) fixed. It kills one
decomposition term and reduces the first support dimension by exactly
one. The projections of the other two legs onto the \(T_0\) blocks
continue to retrieve unchanged \(T_0\): all original \(T_1\)
contributions vanish under those projections, including any first-leg
components introduced by the quotient.

Repeat \(b\) times. The resulting tensor has rank at most \(R(T)-b\)
and still restricts to \(T_0\), proving (4). This is exact substitution;
it does not assert direct-sum rank additivity for arbitrary summands.

## The rectangular bound over every field

For every positive integer \(k\),

\[
R(M_{2,2,k})\ge\frac72 k.                              \tag{5}
\]

Extend the field to its algebraic closure. Any rank scheme over the
original field remains a scheme of the same size there, so a lower bound
over the closure implies the original-field lower bound. We now work
over this infinite field.

Let a rank-\(R\) scheme for multiplication of a \(2\times2\) matrix
\(A\) by a \(2\times k\) matrix \(B\) be written

\[
AB=W\Delta(A)\beta(B),
\]

where \(\Delta(A)\) is the diagonal matrix of the \(R\) first linear
forms, \(\beta:K^{2\times k}\to K^R\) contains the second linear
forms, and \(W:K^R\to K^{2\times k}\) contains the output factors.
Delete any term with identically zero first form.

Every remaining first form is nonzero on the rank-one variety: rank-one
matrices span the \(2\times2\) matrix space. Parametrize that variety
by \(uv^T\). The product of finitely many nonzero first-form
polynomials is nonzero, so over the infinite field there is a rank-one
matrix \(M\) on which all first forms are nonzero.

Choose invertible matrices \(P,Q\) with \(M=PE_{00}Q^{-1}\).
Replace the first forms by \(\alpha_i(PAQ^{-1})\), second forms by
\(\beta_i(QB)\), and output factors by \(P^{-1}C_i\). This is again
a scheme for \(AB\), because
\((PAQ^{-1})(QB)=PAB\). Thus we may assume
\(\Delta_{00}:=\Delta(E_{00})\) is invertible. This normalization
states the compatible three changes explicitly.

The map \(\beta\) is injective: if \(\beta(B)=0\), taking \(A=I\)
gives \(B=0\). The map \(W\) is surjective, again because products
with \(A=I\) cover all outputs. Hence \(\dim\ker W=R-2k\).

Let \(U\) be the \(k\)-dimensional space of matrices whose top row is
zero. Define

\[
\kappa_0(u)=\Delta_{00}\beta(u),\qquad
\kappa_1(u)=\Delta_{10}\beta(u),
\]

\[
\mu(u)=\Delta_{01}\beta(u)-\Delta_{00}\beta(E_{01}u).
\]

All three maps take values in \(\ker W\). Indeed,
\(E_{00}u=E_{10}u=0\) and
\(E_{00}E_{01}u=E_{01}u\). Also \(\kappa_0\) is injective. Put
\(J=\kappa_0(U)\), so \(\dim J=k\).

The map

\[
\Phi:U\longrightarrow(\ker W/J)^2,\qquad
u\longmapsto(\kappa_1(u)+J,\mu(u)+J)
\]

is injective. To check this, suppose both components lie in \(J\).
The diagonal matrices commute, giving

\[
\Delta_{01}\kappa_1(u)-\Delta_{10}\mu(u)
=\Delta_{10}\Delta_{00}\beta(E_{01}u)
\in\Delta_{01}J-\Delta_{10}J.
\]

Cancel the invertible \(\Delta_{00}\), using its commutation with
the other diagonal matrices, to obtain

\[
\Delta_{10}\beta(E_{01}u)
\in\Delta_{01}\beta(U)-\Delta_{10}\beta(U).
\]

Applying \(W\) gives
\(E_{11}u\in E_{01}U-E_{10}U=E_{01}U\). The left side equals
\(u\) and has only its bottom row possibly nonzero; the right side
has only its top row possibly nonzero. Thus \(u=0\).

Dimension comparison now yields

\[
k\le2\bigl((R-2k)-k\bigr),
\]

which is (5). Nothing in this argument divides by two as a field
element. The inequalities concern integer dimensions, and the proof
includes characteristic two.

## Applying the bounds to separate dots

Permute the legs so that \(a=\max(a,b,c)\). The square matrix
multiplication tensor is isomorphic to each such leg permutation (cyclic
index relabeling and transposition), so this does not change its rank.
The first block
\(M_2\otimes\operatorname{DotX}_a\) is the matrix multiplication
tensor \(M_{2,2,2a}\): the length-\(a\) index joins the output-column
index of the \(2\times2\) multiplication. Its first flattening has
dimension four. The other two blocks have first flattening dimensions
\(4b\) and \(4c\), and occupy disjoint coordinate summands.
Applying (4) and then (5) gives

\[
R(M_2\otimes S)
\ge 4b+4c+R(M_{2,2,2a})
\ge7a+4b+4c
\ge5(a+b+c),                                           \tag{6}
\]

where the last inequality is \(2a\ge b+c\).

As an immediate first consequence, suppose \(D\) is a direct sum of
matrix tensors and write \(r_D=\sum_i r_i\) for their matrix ranks.
By (3) and (6), the left tensor in (1) has exact rank at least

\[
r_D+m+5(a+b+c).
\]

The right tensor has an explicit decomposition with
\(r_D+5(a+b+c)\) terms. Exact rank cannot increase under a restriction,
so (1) would require \(m\le0\). This contradicts positive gain.

For equal lengths \(a=b=c=t\), the core bound reads
\(R(M_2\otimes S_t)\ge15t\). For the literal length-two auxiliary,
it reads \(R(M_2\otimes S_2)\ge30\), against the supplied source
rank upper bound 30, leaving no positive-gain margin even with any
finite matrix-sum catalyst. The following replication argument removes
this assumption on \(D\).

## Replicating the rectangular bound

For every positive integer \(N\), the same proof gives

\[
R\bigl(\bigoplus_{i=1}^N M_{2,2,k}\bigr)
\ge\frac72 Nk.                                        \tag{7}
\]

Here is the required adaptation, without assuming arbitrary direct-sum
rank additivity. The first-input space is a direct sum of \(N\) copies
of \(\operatorname{Mat}_{2\times2}\), and second and output spaces
are direct sums of \(N\) copies of \(\operatorname{Mat}_{2\times k}\).
Multiplication acts separately in the blocks.

The product of the \(N\) rank-one first-matrix varieties is irreducible
and spans the full first-input direct sum. To verify the spanning point
explicitly, a nonzero linear form has a nonzero component on some
matrix block; independently varying and scaling a rank-one matrix in
that block prevents its restriction to the product from vanishing
identically. Thus a tuple of nonzero rank-one matrices can avoid all
finitely many nonzero first forms of a rank scheme simultaneously.

Normalize each tuple component to \(E_{00}\) by its own compatible
\(P_i,Q_i\) changes of first, second, and output coordinates. The
diagonal scheme matrix \(\Delta(\mathbf E_{00})\) for the tuple
\(\mathbf E_{00}=(E_{00},\ldots,E_{00})\) is invertible.

Now use the global tuples \(\mathbf E_{10},\mathbf E_{01},
\mathbf E_{11}\), and the space \(U\) of tuples whose top row is
zero in every block. The maps \(\beta,W\) are respectively injective
and surjective by taking the tuple of identity matrices. Their second
and output dimensions are \(2Nk\), while \(\dim U=Nk\).
The same definitions of \(\kappa_0,\kappa_1,\mu,J\) give
\(\dim J=Nk\) and an injective map

\[
U\longrightarrow(\ker W/J)^2.
\]

The preceding proof applies word for word with the global tuples:
diagonal matrices commute, \(\mathbf E_{10}U=0\), and
\(\mathbf E_{11}u\in\mathbf E_{01}U\) forces the bottom row of
every block to vanish. Comparing dimensions gives
\(Nk\le2(R-3Nk)\), proving (7) in every characteristic.

## Eliminating an arbitrary finite catalyst

Put \(C=M_2\otimes S\) and \(L=a+b+c\), again rotating the legs so
that \(a\) is the largest length. In \(N C\), protect the \(N\)
copies of its first component \(M_{2,2,2a}\). The remaining first
support dimension is \(4N(b+c)\). Substitution (4), applied to this
whole direct sum, and replication (7) yield

\[
R(NC)\ge4N(b+c)+7Na\ge5NL.                            \tag{8}
\]

The \(Nm\) separate unit summands each add one to exact rank by (2).
This is an iteration of rank-one matrix-summand additivity; the diagonal
tensor \(Nm\mathrm{Unit}\) is not itself a singleton-leg matrix.
Consequently,

\[
R\bigl(N(m\mathrm{Unit}\oplus C)\bigr)
\ge Nm+5NL.                                           \tag{9}
\]

Suppose (1) held for some arbitrary finite tensor \(D\). Let
\(T=m\mathrm{Unit}\oplus C\) and \(B=5S\), so the supplied maps
would give \(D\oplus T\le D\oplus B\). Applying the same maps
successively, with the other summands as spectators, gives for every
\(N\ge1\)

\[
D\oplus NT\le D\oplus NB.
\]

Dropping \(D\) from the left tensor gives
\(NT\le D\oplus NB\). The right tensor has rank at most
\(R(D)+5NL\), because the supplied decomposition of \(S\) has
\(L\) terms. On the other hand, (9) gives

\[
Nm+5NL\le R(D)+5NL,
\qquad\text{hence}\qquad Nm\le R(D).
\]

Choosing \(N>R(D)\) contradicts \(m\ge1\). No rank additivity
assumption on \(D\) is used. This proves the announced exclusion for
every finite tensor catalyst, equal or unequal positive dot lengths,
and every scalar field.

## Scope and verification

The proof concerns ordinary exact rank and formal tensor restrictions,
so it applies over finite fields as well as infinite fields. It uses
neither a finite search timeout nor a spectral assignment. In
particular, a dot's rank upper bound does not force every normalized
spectral character to take that value; no such assertion is needed.

The mathematical theorem covers arbitrary positive lengths and arbitrary
finite tensor catalysts. Executable recognition may intentionally
cover the smaller class of literal \(S_2\);
recognition must compare coefficients or a certified isomorphism, not
infer membership from shape alone. No Lean theorem or complete
coefficient-extraction certificate is claimed here.

[`mixed_dot_rank.py`](../mixed_dot_rank.py) checks the actual rectangular
projection and first concision for equal lengths and repeated blocks.
[`mixed_dot_catalyst_screen.py`](../mixed_dot_catalyst_screen.py) recognizes
literal S2 and supports every finite D, using a coordinate rank upper bound
to compute a sufficient contradictory repetition count. It does not
materialize that iteration or treat a tensor-rank lower bound as an ordinary
normalized-state row.
