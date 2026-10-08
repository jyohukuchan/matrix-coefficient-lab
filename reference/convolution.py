"""Evaluation/interpolation decomposition from Convolution/Rank.lean."""

from .polynomials import Polynomial
from .sectors import convolution_tensor
from .tensor_schemes import TensorScheme


def lagrange_basis(field, nodes):
    """Full Lagrange basis polynomials; zero is a valid evaluation node."""
    nodes = tuple(nodes)
    for node in nodes:
        field.check(node)
    if len(set(nodes)) != len(nodes):
        raise ValueError("interpolation nodes must be distinct")
    basis = []
    for i, node in enumerate(nodes):
        polynomial = Polynomial.constant(field, 1)
        denominator = 1
        for j, other in enumerate(nodes):
            if i != j:
                polynomial = polynomial*Polynomial(field, (field.neg(other), 1))
                denominator = field.mul(denominator, field.sub(node, other))
        basis.append(polynomial*Polynomial.constant(field, field.inv(denominator)))
    return tuple(basis)


def convolution_scheme(field, a, b, nodes=None):
    """Decompose C(a,b) using a+b-1 evaluations and full interpolation.

    No leading power is removed here, so the evaluation nodes may include 0.
    Canonical encodings enumerate the supplied field, not integer scalar casts.
    Extra distinct nodes are allowed; empty input axes need no evaluations.
    """
    target = convolution_tensor(field, a, b)
    count = target.shape[2] if a and b else 0
    if nodes is None:
        if count > field.order:
            raise ValueError("not enough field elements for convolution interpolation")
        nodes = tuple(range(count))
    else:
        nodes = tuple(nodes)
    basis = lagrange_basis(field, nodes)
    if len(nodes) < count:
        raise ValueError("not enough distinct nodes for the convolution degree bound")
    if not a or not b:
        return TensorScheme(target, (), (), ())
    scheme = TensorScheme(target,
                          tuple(tuple(field.pow(node, i) for i in range(a)) for node in nodes),
                          tuple(tuple(field.pow(node, j) for j in range(b)) for node in nodes),
                          tuple(tuple(p.coefficient(k) for k in range(count)) for p in basis))
    scheme.require_exact()
    return scheme
