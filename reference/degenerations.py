"""Generate exact bilinear schemes from certified polynomial degenerations.

Implements the coefficient-level PolynomialApproximation interface and
normalized nonzero-node recovery from Polynomial/Interpolation.lean. The target
is a matrix multiplication tensor. General source tensors and restriction maps,
and the existence proof's 9/4 decomposition generator, are not implemented.
"""

from dataclasses import dataclass
from functools import cached_property
from itertools import product

from .constructions import constant_weights, nonzero_nodes
from .fields import FiniteField
from .polynomials import Polynomial
from .schemes import BilinearScheme, naive_scheme


@dataclass(frozen=True)
class PolynomialDegeneration:
    """P(t) = sum_q a_q(t) tensor b_q(t) tensor c_q(t).

    Each entry of a,b,c is a Polynomial. Vectors use BilinearScheme's row-major
    coordinates. A valid certificate has [t^j]P=0 for j<leading and
    [t^leading]P equal to the (n,m,k) matrix multiplication tensor. Higher
    coefficients need not vanish or be proportional to the target tensor.
    """

    field: FiniteField
    n: int
    m: int
    k: int
    leading: int
    a: tuple[tuple[Polynomial, ...], ...]
    b: tuple[tuple[Polynomial, ...], ...]
    c: tuple[tuple[Polynomial, ...], ...]

    def __post_init__(self):
        if any(type(d) is not int or d < 1 for d in (self.n, self.m, self.k)):
            raise ValueError("matrix dimensions must be positive integers")
        if type(self.leading) is not int or self.leading < 0:
            raise ValueError("leading degree must be a nonnegative integer")
        for name in ("a", "b", "c"):
            object.__setattr__(self, name, tuple(tuple(row) for row in getattr(self, name)))
        if not len(self.a) == len(self.b) == len(self.c):
            raise ValueError("coefficient families must have equal term counts")
        for family, size in ((self.a, self.n * self.m),
                             (self.b, self.m * self.k), (self.c, self.n * self.k)):
            if any(len(row) != size for row in family):
                raise ValueError("polynomial coefficient vector has the wrong length")
            for row in family:
                for polynomial in row:
                    if not isinstance(polynomial, Polynomial) or polynomial.field != self.field:
                        raise ValueError("every coefficient must be a Polynomial over the same field")

    @property
    def terms(self):
        return len(self.a)

    @property
    def degree_bound(self):
        """Sum of maximum local degrees, a safe bound even with cancellation."""
        return sum(max((max(0, p.degree) for row in family for p in row), default=0)
                   for family in (self.a, self.b, self.c))

    @property
    def recovery_node_count(self):
        """Enough nodes for P(t)/t^leading, assuming the certificate passes."""
        return max(0, self.degree_bound - self.leading) + 1

    @cached_property
    def _tensor_coefficients(self):
        # Work in the formal polynomial ring; cancellation is across all terms.
        zero = Polynomial(self.field)
        coefficients = {}
        for a, b, c in zip(self.a, self.b, self.c):
            aa = [(i, x) for i, x in enumerate(a) if x.coefficients]
            bb = [(j, y) for j, y in enumerate(b) if y.coefficients]
            cc = [(h, z) for h, z in enumerate(c) if z.coefficients]
            for (i, x), (j, y), (h, z) in product(aa, bb, cc):
                key = i, j, h
                coefficients[key] = coefficients.get(key, zero) + x * y * z
        return coefficients

    def tensor_polynomial(self, left_index, right_index, output_index):
        """Inspect one formal tensor coordinate, even for an invalid certificate."""
        for index, size in ((left_index, self.n * self.m), (right_index, self.m * self.k),
                            (output_index, self.n * self.k)):
            if type(index) is not int or not 0 <= index < size:
                raise ValueError("tensor coordinate is out of range")
        return self._tensor_coefficients.get((left_index, right_index, output_index),
                                             Polynomial(self.field))

    @cached_property
    def _first_mismatch(self):
        for i, j, h in product(range(self.n * self.m), range(self.m * self.k),
                               range(self.n * self.k)):
            polynomial = self.tensor_polynomial(i, j, h)
            for degree, coefficient in enumerate(polynomial.coefficients[:self.leading]):
                if coefficient != 0:
                    return ((i, j, h), degree, 0, coefficient)
            row, inner = divmod(i, self.m)
            inner2, column = divmod(j, self.k)
            outrow, outcolumn = divmod(h, self.k)
            expected = int(row == outrow and inner == inner2 and column == outcolumn)
            actual = polynomial.coefficient(self.leading)
            if actual != expected:
                return ((i, j, h), self.leading, expected, actual)
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
        """Constant/monomial lift of a known scheme, useful as a control."""
        scheme.require_exact()
        if type(leading) is not int or leading < 0:
            raise ValueError("leading degree must be a nonnegative integer")
        f = scheme.field
        return cls(f, scheme.n, scheme.m, scheme.k, leading,
                   [[Polynomial.monomial(f, leading, x) for x in row] for row in scheme.a],
                   [[Polynomial.constant(f, x) for x in row] for row in scheme.b],
                   [[Polynomial.constant(f, x) for x in row] for row in scheme.c])

    def evaluate_at(self, node):
        """Return the evaluated coefficient families, NOT necessarily an exact scheme.

        This is P(node), not its leading coefficient. The returned BilinearScheme
        permits inspection/verification, but rejects multiplication unless its
        own exact matrix-tensor certificate happens to pass.
        """
        self.field.check(node)
        evaluate = lambda family: [[p.evaluate(node) for p in row] for row in family]
        return BilinearScheme(self.field, self.n, self.m, self.k,
                              evaluate(self.a), evaluate(self.b), evaluate(self.c))

    def recover_scheme(self, nodes=None):
        """Generate a certified exact scheme using nonzero-node interpolation.

        At each node t, append r evaluated terms, multiplying c by w_t/t^leading.
        With D=self.degree_bound, D-leading+1 nodes suffice. This uses the
        degree of the normalized polynomial, tighter than the Lean power-rank
        theorem's D+1 bound; it follows its normalized recovery identity.
        Extra supplied nodes are allowed. Zero terms are retained.
        """
        self.require_valid()
        f = self.field
        if nodes is None:
            nodes = nonzero_nodes(f, self.recovery_node_count)
        else:
            nodes = tuple(nodes)
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
        scheme = BilinearScheme(f, self.n, self.m, self.k, aa, bb, cc)
        scheme.require_exact()
        return scheme

    def tensor_product(self, other):
        """Tensor polynomial families BEFORE recovery, using the same parameter t."""
        if self.field != other.field:
            raise ValueError("tensor factors must use the same field representation")
        self.require_valid()
        other.require_valid()

        def combine(a, b, rows1, columns1, rows2, columns2):
            return tuple(a[i1 * columns1 + j1] * b[i2 * columns2 + j2]
                         for i1 in range(rows1) for i2 in range(rows2)
                         for j1 in range(columns1) for j2 in range(columns2))

        aa, bb, cc = [], [], []
        for q, r in product(range(self.terms), range(other.terms)):
            aa.append(combine(self.a[q], other.a[r], self.n, self.m, other.n, other.m))
            bb.append(combine(self.b[q], other.b[r], self.m, self.k, other.m, other.k))
            cc.append(combine(self.c[q], other.c[r], self.n, self.k, other.n, other.k))
        return PolynomialDegeneration(self.field, self.n * other.n, self.m * other.m,
                                      self.k * other.k, self.leading + other.leading, aa, bb, cc)

    def tensor_power(self, exponent):
        if type(exponent) is not int or exponent < 0:
            raise ValueError("tensor exponent must be a nonnegative integer")
        self.require_valid()
        result = PolynomialDegeneration.from_scheme(naive_scheme(self.field, 1, 1, 1))
        for _ in range(exponent):
            result = result.tensor_product(self)
        return result


def left_perturbation_fixture(scheme):
    """Synthetic TEST INPUT built from a known exact scheme, not a 9/4 construction.

    Replace each a_q[j] by t*a_q[j] + t^2*a_q[(j+1) mod (n*m)].
    Thus P(t)=t*T+t^2*S with a cyclic permutation on the left coordinate.
    The higher tensor S need not be proportional to T. This exercises recovery
    from genuine tensor-valued noise, while avoiding any claim of a new scheme.
    """
    scheme.require_exact()
    f = scheme.field
    aa = [[Polynomial(f, (0, x, row[(i + 1) % len(row)])) for i, x in enumerate(row)]
          for row in scheme.a]
    constant = lambda family: [[Polynomial.constant(f, x) for x in row] for row in family]
    return PolynomialDegeneration(f, scheme.n, scheme.m, scheme.k, 1,
                                  aa, constant(scheme.b), constant(scheme.c))
