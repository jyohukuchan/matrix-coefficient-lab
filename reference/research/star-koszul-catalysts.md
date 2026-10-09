# A Koszul obstruction to every positive-gain star catalyst

This note proves an exclusion for one fixed auxiliary tensor over **every
field**, including finite fields and characteristic two. Let

\[
 S=a_1\otimes(e_0\otimes e_1+e_1\otimes e_0)
   +a_2\otimes(e_0\otimes e_2+e_2\otimes e_0).
\]

For any finite three-leg tensor \(D\), any integer \(m\geq1\), and any field
\(K\), there cannot be three linear local restriction maps realizing

\[
 D\oplus m\langle1\rangle\oplus(M_2\otimes S)
       \ \leq\ D\oplus5S.                                      \tag{1}
\]

The result concerns exact local restrictions and a positive number of
independent scalar summands. It supplies no matrix multiplication scheme,
optimal rank statement, or exclusion for a different auxiliary tensor. It
also excludes fewer than five source copies, since extra copies may be
projected away. The algebraic argument below is a written proof backed by a
replayable integer certificate; it has not been formalized in Lean.

## The exact integer target and its flattening

Put \(T=\langle1\rangle\oplus(M_2\otimes S)\), of shape \((9,13,13)\).
The matrix tensor uses row-major triples
\((2i+j,2j+h,2i+h)\), for \(i,j,h\in\{0,1\}\).
The integer support of \(T\) consists of \((0,0,0)\) and the 32 triples

\[
 (1+2(2i+j)+a,\;1+3(2j+h)+b,\;1+3(2i+h)+c),
\]

where \((a,b,c)\) runs through
\((0,0,1),(0,1,0),(1,0,2),(1,2,0)\). Every supported coefficient is one.
This literal construction was independently compared, at every coordinate,
with the reference tensor constructors.

For \(A=K^9\), form the Koszul flattening

\[
 \mathcal K_T:\Lambda^4 A\otimes B^*
                      \longrightarrow\Lambda^5 A\otimes C.
\]

Its column indexed by \((I,b)\), where \(|I|=4\), has the entry
\((-1)^{|\{i\in I:i<a\}|}T_{abc}\) in row
\((I\cup\{a\},c)\), whenever \(a\notin I\). Thus its integer matrix has
1638 rows, 1638 columns, and 2310 nonzero entries. A rank-one tensor contributes
matrix rank at most

\[
 \binom{8}{4}=70.                                                \tag{2}
\]

Indeed, for nonzero first factor \(a\), exterior multiplication by \(a\)
has this rank after choosing a basis with \(a\) first; the other two factors
contribute a rank-one map. This argument works in every characteristic.

The bipartite graph of the nonzero matrix entries has 550 nonempty connected
components, each with at most 12 rows and 12 columns. Components have disjoint
row and column sets and exhaust every nonzero entry. Recorded integer row
and column swaps and additions transform each component to expose a
square minor of determinant \(+1\) or \(-1\). The sum of their minor sizes
is **1054**. All operations are unimodular, so the matrices remain equivalent
after reduction to any field. Every selected determinant remains nonzero.
Consequently

\[
 \operatorname{rank}\mathcal K_T\geq1054
                 \quad\text{over every field}.                 \tag{3}
\]

The certificate was independently checked by reconstructing the signed
matrix from the support above, replaying all 1244 elementary operations, and
computing the selected determinants using rational Gaussian elimination
rather than the certificate generator's elimination routine. An independent
bitset elimination gives rank exactly 1054 over \(\mathbb F_2\).
The all-field statement follows from the unimodular operations and unit
minors, not from that binary rank computation alone.

## Repetition keeps the exterior dimension fixed

Assume (1). Applying the same restriction successively to a single retained
copy of \(D\), each time with five fresh copies of \(S\), gives

\[
 D\oplus5NS\ \geq\ D\oplus Nm\langle1\rangle
                          \oplus N(M_2\otimes S)
\]

for every positive integer \(N\). Project away \(D\) and select one scalar
summand from each of the \(N\) gain groups. This gives a restriction to
\(N\) direct copies of \(T\).

Now identify the first coordinate spaces of those copies by the linear map
\(\bigoplus_{n=1}^N K^9\to K^9\) that sends each copy's basis vector \(e_a\)
to the same \(e_a\). **Keep the second and third coordinate spaces separate.**
The resulting tensor \(\widetilde T_N\) has shape \((9,13N,13N)\).
Its Koszul matrix, using the same exterior degree four, is a block diagonal
sum of \(N\) copies of \(\mathcal K_T\), after permuting rows and columns.
In particular,

\[
 \operatorname{rank}\mathcal K_{\widetilde T_N}\geq1054N.         \tag{4}
\]

This is additivity of the rank of a block diagonal **matrix**. It does not
use additivity of tensor rank or border rank. The rank-one contribution
remains 70 because the first coordinate space still has dimension nine.
As a finite control, independently assembling the two-copy merged matrix
and eliminating over \(\mathbb F_2\) gives rank 2108.

## Formal polynomial transfer and the contradiction

There is an integer, division-free three-term polynomial identity

\[
\begin{aligned}
 P(t)={}&a_1\otimes(e_0+te_1)\otimes(e_0+te_1)\\
       &+a_2\otimes(e_0+te_2)\otimes(e_0+te_2)\\
       &-(a_1+a_2)\otimes e_0\otimes e_0\\
     ={}&tS+t^2E,
\end{aligned}
\]

where \(E=a_1\otimes e_1\otimes e_1+a_2\otimes e_2\otimes e_2\).
This identity holds formally over every field, including characteristic
two. No evaluation nodes or field extension are needed.

Choose any finite exact rank-one decomposition of \(D\) with \(C_D\) terms;
the coefficientwise decomposition always supplies one. Multiplying that
decomposition by \(t\) and taking \(5N\) copies of \(P(t)\) yields a
polynomial tensor with \(C_D+15N\) rank-one terms and leading coefficient
\(D\oplus5NS\) at degree one. Apply the iterated restriction, the scalar
projections, and the first-space identification to every term. The result
has the form

\[
 Q(t)=t\widetilde T_N+t^2 F.
\]

Over the rational function field \(K(t)\), the tensor \(Q(t)/t\) still has
at most \(C_D+15N\) rank-one terms. Its coordinates are regular at \(t=0\)
and specialize to \(\widetilde T_N\). A nonzero minor of the specialized
Koszul matrix remains a nonzero polynomial minor; hence (4) is a lower bound
for its rank over \(K(t)\). By (2) and subadditivity of matrix rank,

\[
 1054N\leq70(C_D+15N)=70C_D+1050N,
 \qquad 4N\leq70C_D.                                           \tag{5}
\]

The right side is fixed while \(N\) is arbitrary. Taking
\(N>70C_D/4\) contradicts (5), proving the exclusion (1).
The use of \(K(t)\) is a proof step; it does not claim an executable
decomposition over the original field or silently extend any finite
implementation's coefficient field.

## Certificate artifacts

The compact [integer certificate](../star_koszul.json) records the support
component indices, every elementary operation and each unit minor determinant.
[`star_koszul.py`](../star_koszul.py) reconstructs the integer target and signed
sparse matrix from the explicit formulas above and verifies all component
partitions and operations, with explicit allocation and integer caps.
[`test_star_koszul.py`](../test_star_koszul.py) separately constructs the
exterior matrix by deleting an axis from its output wedge, compares every
target coefficient and checks corruption controls. The reference
[`star_degeneration.py`](../star_degeneration.py) verifies the complete
formal polynomial identity at every tensor coordinate. Large generation and
independent audit logs remain outside the repository.
