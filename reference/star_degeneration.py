"""A three-term star-pencil degeneration and direct coefficient extraction.

This is a finite mathematical calibration, not a matrix multiplication scheme.
The constant first factors let [t^e] of the powered polynomial be decomposed
into (e+1)*3^e actual rank-one terms over the supplied field itself. No
interpolation nodes or extension of that field are needed.
"""

from dataclasses import dataclass
from itertools import product

from .catalyst_compiler import CompilerLimits, CompilerResourceLimit
from .fields import FiniteField
from .polynomials import Polynomial
from .tensor_degenerations import TensorDegeneration
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor


DEFAULT_LIMITS = CompilerLimits()


def star_pencil_tensor(field):
    """S=a1 tensor (e0 e1+e1 e0)+a2 tensor (e0 e2+e2 e0)."""
    return FiniteTensor.from_function(field, (2, 3, 3),
                                     lambda i, j, k: int((j, k) in ((0, i+1), (i+1, 0))))


def star_error_tensor(field):
    """E=a1 tensor e1 tensor e1 + a2 tensor e2 tensor e2."""
    return FiniteTensor.from_function(field, (2, 3, 3),
                                     lambda i, j, k: int(j == k == i+1))


def require_star_identity(degeneration):
    """Check the complete formal identity P(t)=t*S+t^2*E at all 18 coordinates."""
    if not isinstance(degeneration, TensorDegeneration):
        raise ValueError("supply an actual polynomial tensor degeneration")
    field = degeneration.field
    if degeneration.target != star_pencil_tensor(field) or degeneration.leading != 1:
        raise ValueError("the pencil must target the star tensor at leading degree one")
    error = star_error_tensor(field)
    for coordinate in product(range(2), range(3), range(3)):
        expected = Polynomial(field, (0, degeneration.target.coefficient(*coordinate),
                                      error.coefficient(*coordinate)))
        if degeneration.tensor_polynomial(*coordinate) != expected:
            raise ValueError(f"star pencil formal identity fails at {coordinate}")
    degeneration.require_valid()


def star_polynomial_degeneration(field):
    """The actual three terms, including -(a1+a2) tensor e0 tensor e0."""
    if not isinstance(field, FiniteField):
        raise ValueError("supply a represented finite field")
    zero, one = Polynomial(field), Polynomial.constant(field, 1)
    t, minus_one = Polynomial.monomial(field, 1), Polynomial.constant(field, field.neg(1))
    degeneration = TensorDegeneration(
        star_pencil_tensor(field), 1,
        ((one, zero), (zero, one), (minus_one, minus_one)),
        ((one, t, zero), (one, zero, t), (one, zero, zero)),
        ((one, t, zero), (one, zero, t), (one, zero, zero)))
    require_star_identity(degeneration)
    return degeneration


@dataclass(frozen=True)
class StarPowerExtraction:
    """Formal powered input and its independently verified exact decomposition."""

    exponent: int
    powered: TensorDegeneration
    scheme: TensorScheme

    @property
    def term_bound(self):
        return (self.exponent+1)*3**self.exponent


def extract_star_power(field, exponent=1, *, limits=DEFAULT_LIMITS):
    """Extract [t^e]P(t)^tensor e by summing b_i tensor c_(e-i).

    Every first-leg polynomial is constant. Thus one term per split suffices,
    without a third degree index or characteristic-dependent interpolation.
    Zero rank-one terms are retained so the displayed length is exactly the
    elementary upper bound; no rank optimality claim is made.
    """
    if not isinstance(field, FiniteField):
        raise ValueError("supply a represented finite field")
    if type(exponent) is not int or exponent < 0:
        raise ValueError("the exponent must be a nonnegative integer")
    if not isinstance(limits, CompilerLimits):
        raise ValueError("supply explicit compiler limits")
    # The target has 18^e coefficients, hence at least 2^(4e). Reject huge e
    # before any exponentiation, polynomial construction, or power loop.
    if 4*exponent >= limits.max_tensor_coefficients.bit_length():
        raise CompilerResourceLimit("star power tensor cap reached; construction incomplete")
    shape = (2**exponent, 3**exponent, 3**exponent)
    terms = (exponent+1)*3**exponent
    limits.scheme(shape, terms)
    if exponent == 0:
        scheme = TensorScheme.from_tensor(FiniteTensor.unit(field))
        return StarPowerExtraction(0, TensorDegeneration.from_scheme(scheme), scheme)
    powered = star_polynomial_degeneration(field).tensor_power(exponent)
    # This check is deliberately formal, including coefficients below e.
    powered.require_valid()
    if any(polynomial.degree > 0 for row in powered.a for polynomial in row):
        raise ValueError("direct two-leg extraction requires constant first factors")
    aa, bb, cc = [], [], []
    for a, b, c in zip(powered.a, powered.b, powered.c):
        first = tuple(polynomial.coefficient(0) for polynomial in a)
        for split in range(exponent+1):
            aa.append(first)
            bb.append(tuple(polynomial.coefficient(split) for polynomial in b))
            cc.append(tuple(polynomial.coefficient(exponent-split) for polynomial in c))
    scheme = TensorScheme(powered.target, aa, bb, cc)
    scheme.require_exact()
    return StarPowerExtraction(exponent, powered, scheme)
