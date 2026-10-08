"""Exact convolution using finite evaluation nodes and the leading coefficient.

The extra projective node is infinity: it selects the top coefficient of each
input polynomial. This is an executable finite construction, not a new Lean
formalization or a near-9/4 matrix multiplication scheme.
"""

from .convolution import lagrange_basis
from .fields import FiniteField
from .polynomials import Polynomial
from .sectors import convolution_tensor
from .tensor_schemes import TensorScheme


def projective_convolution_scheme(field, a, b, nodes=None):
    """Generate C(a,b) with a+b-2 finite nodes and one infinity term.

    For positive axes, exactly degree=a+b-2 distinct finite nodes are required;
    they may include zero. Infinity selects input coefficients a-1 and b-1.
    Its output polynomial is t^degree-sum_i nodes[i]^degree*L_i(t).

    Canonical field encodings enumerate default nodes; they are not integer
    scalar casts. The construction requires field.order >= degree. Failure of
    this node condition does not imply nonexistence of a different exact scheme.
    Empty input axes produce a zero-term scheme and require an empty node list.
    a=b=1 requires no finite nodes and returns the one-term scalar scheme.
    """
    if not isinstance(field, FiniteField):
        raise ValueError("projective interpolation requires an explicitly represented finite field")
    target = convolution_tensor(field, a, b)
    degree = a+b-2 if a and b else 0
    if nodes is None:
        if degree > field.order:
            raise ValueError("not enough finite field elements for projective convolution interpolation")
        nodes = tuple(range(degree))
    else:
        nodes = tuple(nodes)
    if len(nodes) != degree:
        raise ValueError(f"projective convolution requires exactly {degree} finite nodes")
    basis = lagrange_basis(field, nodes)
    if not a or not b:
        return TensorScheme(target, (), (), ())
    infinite_output = Polynomial.monomial(field, degree)
    for node, polynomial in zip(nodes, basis):
        infinite_output = infinite_output-polynomial*Polynomial.constant(field, field.pow(node, degree))
    left = tuple(tuple(field.pow(node, i) for i in range(a)) for node in nodes)
    right = tuple(tuple(field.pow(node, j) for j in range(b)) for node in nodes)
    left += (tuple(int(i == a-1) for i in range(a)),)
    right += (tuple(int(j == b-1) for j in range(b)),)
    output = tuple(tuple(polynomial.coefficient(k) for k in range(degree+1)) for polynomial in basis)
    output += (tuple(infinite_output.coefficient(k) for k in range(degree+1)),)
    scheme = TensorScheme(target, left, right, output)
    scheme.require_exact()
    return scheme
