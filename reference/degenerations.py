"""Generate exact bilinear schemes from certified polynomial degenerations.

Implements the coefficient-level PolynomialApproximation interface and
normalized nonzero-node recovery from Polynomial/Interpolation.lean. This is the
matrix-specific adapter for the generic TensorDegeneration interface. General
polynomial restriction maps and the 9/4 decomposition generator are absent.
"""

from dataclasses import dataclass
from functools import cached_property
from itertools import product

from .fields import FiniteField
from .polynomials import Polynomial
from .schemes import BilinearScheme, naive_scheme
from .tensor_degenerations import TensorDegeneration
from .tensors import matrix_multiplication_tensor


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
        degeneration = self.as_tensor_degeneration()
        for name in ("a", "b", "c"):
            object.__setattr__(self, name, getattr(degeneration, name))

    @cached_property
    def _tensor_degeneration(self):
        return TensorDegeneration(matrix_multiplication_tensor(self.field, self.n, self.m, self.k),
                                  self.leading, self.a, self.b, self.c)

    def as_tensor_degeneration(self):
        """Expose the same certificate with an explicit finite-tensor target."""
        return self._tensor_degeneration

    @property
    def terms(self):
        return len(self.a)

    @property
    def degree_bound(self):
        """Sum of maximum local degrees, a safe bound even with cancellation."""
        return self.as_tensor_degeneration().degree_bound

    @property
    def recovery_node_count(self):
        """Enough nodes for P(t)/t^leading, assuming the certificate passes."""
        return self.as_tensor_degeneration().recovery_node_count

    @cached_property
    def _tensor_coefficients(self):
        return self.as_tensor_degeneration()._tensor_coefficients

    def tensor_polynomial(self, left_index, right_index, output_index):
        """Inspect one formal tensor coordinate, even for an invalid certificate."""
        return self.as_tensor_degeneration().tensor_polynomial(left_index, right_index, output_index)

    @cached_property
    def _first_mismatch(self):
        return self.as_tensor_degeneration()._first_mismatch

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
        evaluated = self.as_tensor_degeneration().evaluate_at(node)
        return BilinearScheme(self.field, self.n, self.m, self.k,
                              evaluated.a, evaluated.b, evaluated.c)

    def recover_scheme(self, nodes=None):
        """Generate a certified exact scheme using nonzero-node interpolation.

        At each node t, append r evaluated terms, multiplying c by w_t/t^leading.
        With D=self.degree_bound, D-leading+1 nodes suffice. This uses the
        degree of the normalized polynomial, tighter than the Lean power-rank
        theorem's D+1 bound; it follows its normalized recovery identity.
        Extra supplied nodes are allowed. Zero terms are retained.
        """
        return self.as_tensor_degeneration().recover_scheme(nodes).to_matrix_scheme(self.n, self.m, self.k)

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
