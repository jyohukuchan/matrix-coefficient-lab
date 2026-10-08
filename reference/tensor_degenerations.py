"""Polynomial approximations to arbitrary explicit finite tensor targets."""

from dataclasses import dataclass
from functools import cached_property
from itertools import product

from .constructions import constant_weights, nonzero_nodes
from .polynomials import Polynomial
from .tensor_schemes import TensorScheme, _families, _product_families
from .tensors import FiniteTensor


@dataclass(frozen=True)
class TensorDegeneration:
    """P(t)=sum_q a_q(t) tensor b_q(t) tensor c_q(t), targeting an explicit T.

    All coefficients below leading must vanish; [t^leading]P must equal T at
    every coordinate. The target need not be a matrix multiplication tensor.
    Higher coefficients are unrestricted, including non-proportional noise.
    """

    target: FiniteTensor
    leading: int
    a: tuple[tuple[Polynomial, ...], ...]
    b: tuple[tuple[Polynomial, ...], ...]
    c: tuple[tuple[Polynomial, ...], ...]

    def __post_init__(self):
        if not isinstance(self.target, FiniteTensor):
            raise ValueError("target must be a FiniteTensor")
        if type(self.leading) is not int or self.leading < 0:
            raise ValueError("leading degree must be a nonnegative integer")

        def check(polynomial):
            if not isinstance(polynomial, Polynomial) or polynomial.field != self.field:
                raise ValueError("every coefficient must be a Polynomial over the same field")

        for name, family in zip(("a", "b", "c"),
                                _families(self.shape, self.a, self.b, self.c, check)):
            object.__setattr__(self, name, family)

    @property
    def field(self):
        return self.target.field

    @property
    def shape(self):
        return self.target.shape

    @property
    def terms(self):
        return len(self.a)

    @property
    def degree_bound(self):
        return sum(max((max(0, p.degree) for row in family for p in row), default=0)
                   for family in (self.a, self.b, self.c))

    @property
    def recovery_node_count(self):
        return max(0, self.degree_bound-self.leading) + 1

    @cached_property
    def _tensor_coefficients(self):
        zero = Polynomial(self.field)
        coefficients = {}
        for a, b, c in zip(self.a, self.b, self.c):
            aa = [(i, x) for i, x in enumerate(a) if x.coefficients]
            bb = [(j, y) for j, y in enumerate(b) if y.coefficients]
            cc = [(k, z) for k, z in enumerate(c) if z.coefficients]
            for (i, x), (j, y), (k, z) in product(aa, bb, cc):
                key = i, j, k
                coefficients[key] = coefficients.get(key, zero) + x*y*z
        return coefficients

    def tensor_polynomial(self, i, j, k):
        # The target's accessor validates all three coordinates.
        self.target.coefficient(i, j, k)
        return self._tensor_coefficients.get((i, j, k), Polynomial(self.field))

    @cached_property
    def _first_mismatch(self):
        for i, j, k in product(*(range(d) for d in self.shape)):
            polynomial = self.tensor_polynomial(i, j, k)
            for degree, coefficient in enumerate(polynomial.coefficients[:self.leading]):
                if coefficient:
                    return ((i, j, k), degree, 0, coefficient)
            expected = self.target.coefficient(i, j, k)
            actual = polynomial.coefficient(self.leading)
            if actual != expected:
                return ((i, j, k), self.leading, expected, actual)
        return None

    def verify(self):
        return self._first_mismatch is None

    def require_valid(self):
        if self._first_mismatch is not None:
            coordinate, degree, expected, actual = self._first_mismatch
            raise ValueError(f"invalid degeneration at {coordinate}, degree {degree}:"
                             f" expected {expected}, got {actual}")

    @classmethod
    def from_scheme(cls, scheme, leading=0):
        scheme.require_exact()
        if type(leading) is not int or leading < 0:
            raise ValueError("leading degree must be a nonnegative integer")
        f = scheme.field
        return cls(scheme.target, leading,
                   [[Polynomial.monomial(f, leading, x) for x in row] for row in scheme.a],
                   [[Polynomial.constant(f, x) for x in row] for row in scheme.b],
                   [[Polynomial.constant(f, x) for x in row] for row in scheme.c])

    def evaluate_at(self, node):
        """Intermediate coefficient families for P(node), usually not equal to T.

        The returned scheme retains T as its declared target. Its apply() method
        rejects it unless its own exact coefficient certificate passes.
        """
        self.field.check(node)
        evaluate = lambda family: [[p.evaluate(node) for p in row] for row in family]
        return TensorScheme(self.target, evaluate(self.a), evaluate(self.b), evaluate(self.c))

    def recover_scheme(self, nodes=None):
        self.require_valid()
        f = self.field
        nodes = nonzero_nodes(f, self.recovery_node_count) if nodes is None else tuple(nodes)
        weights = constant_weights(f, nodes)
        if len(nodes) < self.recovery_node_count:
            raise ValueError("not enough nodes for the normalized tensor degree bound")
        aa, bb, cc = [], [], []
        for node, weight in zip(nodes, weights):
            evaluated = self.evaluate_at(node)
            scale = f.div(weight, f.pow(node, self.leading))
            aa.extend(evaluated.a)
            bb.extend(evaluated.b)
            cc.extend(tuple(f.mul(scale, x) for x in row) for row in evaluated.c)
        scheme = TensorScheme(self.target, aa, bb, cc)
        scheme.require_exact()
        return scheme

    def tensor_product(self, other):
        self.require_valid()
        other.require_valid()
        target = self.target.tensor_product(other.target)
        families = _product_families(self.a, self.b, self.c, other.a, other.b, other.c,
                                     lambda a, b: a*b)
        return TensorDegeneration(target, self.leading+other.leading, *families)

    def tensor_power(self, exponent):
        if type(exponent) is not int or exponent < 0:
            raise ValueError("tensor exponent must be a nonnegative integer")
        self.require_valid()
        result = TensorDegeneration.from_scheme(TensorScheme.from_tensor(FiniteTensor.unit(self.field)))
        for _ in range(exponent):
            result = result.tensor_product(self)
        return result
