# Border rank of the three-oriented dot tensor

The mixed-dot tensor has border rank exactly six over every field.

Use the literal support supplied in the task:

  (0,0,0), (0,1,1), (1,2,2), (2,2,3), (3,3,4), (4,4,4).

Its first-leg slice at v=(s,u0,u1,w0,w1), with second-leg coordinates as
rows and third-leg coordinates as columns, is

  [ s 0 0  0  0  ]
  [ 0 s 0  0  0  ]
  [ 0 0 u0 u1 0  ]
  [ 0 0 0  0  w0 ]
  [ 0 0 0  0  w1 ].

For every tensor of rank at most five with second and third dimensions
five, its slice matrices satisfy the degree-six polynomial identity

  A_x adj(A_y) A_z - A_z adj(A_y) A_x = 0.

Indeed, pad a rank decomposition to five terms and put their second and
third vectors into the columns of 5-by-5 matrices B,C. Write D_s for the
diagonal matrix of first-leg evaluations. Then A_s=B D_s C^T. The
adjugate product identity and the identities M adj(M)=adj(M) M=det(M)I
give

  A_x adj(A_y) A_z
    = det(B)det(C) B D_x adj(D_y) D_z C^T.

All three middle matrices are diagonal, so the expression is symmetric
in x,z. This proof uses polynomial identities over the integers and
requires no inverses, division, generic invertible slice, or
characteristic restriction. Thus the identity also vanishes on the
Zariski closure of tensors of rank at most five, including after extension
to the algebraic closure of any field.

For the displayed tensor take

  x=e1, y=e0+e2+e4, z=e3.

The only nonzero entry of adj(A_y) is (2,3), equal to -1. The cofactor
minor is the 4-by-4 identity matrix; its sign is (-1)^(2+3)=-1. Therefore
the displayed commutator has its only nonzero entry at (2,4), also -1.
That value remains nonzero in every characteristic. Border rank is
therefore at least six. The six literal support terms give an exact
rank-one decomposition, proving the matching border upper bound six.

The adjacent Python script independently constructs all integer slices,
cofactors and products, checks the identity on twelve explicitly exported
integer rank-five tensors, and exports the actual obstruction. The tests
support the symbolic proof; the symbolic argument supplies its universal
scope. The conclusion rules out border rank five and does not itself
exclude arbitrary catalytic D for the six-term source budget.
