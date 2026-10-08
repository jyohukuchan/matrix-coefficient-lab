"""Exact rank-one decompositions for arbitrary finite three-leg tensors."""

from dataclasses import dataclass
from functools import cached_property
from itertools import product

from .tensors import FiniteTensor, _vector, matrix_multiplication_tensor


def _families(shape, a, b, c, check):
    families = tuple(tuple(tuple(row) for row in family) for family in (a, b, c))
    if not len(families[0]) == len(families[1]) == len(families[2]):
        raise ValueError("coefficient families must have equal term counts")
    for family, size in zip(families, shape):
        if any(len(row) != size for row in family):
            raise ValueError("coefficient vector has the wrong length")
        for row in family:
            for value in row:
                check(value)
    return families


def _product_families(a, b, c, other_a, other_b, other_c, multiply):
    # Each leg uses the same lexicographic pair indexing as FiniteTensor.
    result = [[], [], []]
    for q, r in product(range(len(a)), range(len(other_a))):
        for out, first, second in zip(result, (a, b, c), (other_a, other_b, other_c)):
            out.append(tuple(multiply(x, y) for x in first[q] for y in second[r]))
    return tuple(tuple(rows) for rows in result)


@dataclass(frozen=True)
class TensorScheme:
    """T = sum_q a[q] tensor b[q] tensor c[q], for an explicit target T.

    terms is the supplied decomposition length, a rank upper bound if verified.
    It is not minimum rank, scalar-operation count, or execution time.
    """

    target: FiniteTensor
    a: tuple[tuple[int, ...], ...]
    b: tuple[tuple[int, ...], ...]
    c: tuple[tuple[int, ...], ...]

    def __post_init__(self):
        if not isinstance(self.target, FiniteTensor):
            raise ValueError("target must be a FiniteTensor")
        for name, family in zip(("a", "b", "c"),
                                _families(self.shape, self.a, self.b, self.c, self.field.check)):
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

    @cached_property
    def _first_mismatch(self):
        f = self.field
        coefficients = {}
        for a, b, c in zip(self.a, self.b, self.c):
            aa = [(i, x) for i, x in enumerate(a) if x]
            bb = [(j, y) for j, y in enumerate(b) if y]
            cc = [(k, z) for k, z in enumerate(c) if z]
            for (i, x), (j, y), (k, z) in product(aa, bb, cc):
                key = i, j, k
                coefficients[key] = f.add(coefficients.get(key, 0), f.mul(f.mul(x, y), z))
        for i, j, k in product(*(range(d) for d in self.shape)):
            expected = self.target.coefficient(i, j, k)
            actual = coefficients.get((i, j, k), 0)
            if actual != expected:
                return ((i, j, k), expected, actual)
        return None

    def verify(self):
        return self._first_mismatch is None

    def require_exact(self):
        if self._first_mismatch is not None:
            coordinate, expected, actual = self._first_mismatch
            raise ValueError(f"invalid tensor coefficient {coordinate}: expected {expected}, got {actual}")

    def apply(self, left, right):
        """Evaluate the certified bilinear map on two input vectors."""
        self.require_exact()
        f = self.field
        left, right = _vector(f, left, self.shape[0]), _vector(f, right, self.shape[1])
        values = [f.mul(f.dot(a, left), f.dot(b, right)) for a, b in zip(self.a, self.b)]
        return tuple(f.sum(f.mul(c[k], value) for c, value in zip(self.c, values))
                     for k in range(self.shape[2]))

    @classmethod
    def from_tensor(cls, target):
        """A trivial exact decomposition, at most dim(X)*dim(Y) terms.

        One term per input-coordinate pair, with its output coefficient vector;
        pairs with a zero output vector are omitted. No low-rank search.
        """
        if not isinstance(target, FiniteTensor):
            raise ValueError("target must be a FiniteTensor")
        aa, bb, cc = [], [], []
        x, y, z = target.shape
        for i, j in product(range(x), range(y)):
            row = tuple(target.coefficient(i, j, k) for k in range(z))
            if any(row):
                aa.append(tuple(int(h == i) for h in range(x)))
                bb.append(tuple(int(h == j) for h in range(y)))
                cc.append(row)
        return cls(target, aa, bb, cc)

    def tensor_product(self, other):
        self.require_exact()
        other.require_exact()
        target = self.target.tensor_product(other.target)
        return TensorScheme(target, *_product_families(self.a, self.b, self.c,
                                                      other.a, other.b, other.c, self.field.mul))

    def tensor_power(self, exponent):
        if type(exponent) is not int or exponent < 0:
            raise ValueError("tensor exponent must be a nonnegative integer")
        self.require_exact()
        result = TensorScheme.from_tensor(FiniteTensor.unit(self.field))
        for _ in range(exponent):
            result = result.tensor_product(self)
        return result

    def rescale_terms(self, left_scale, right_scale):
        f = self.field
        inverse = f.inv(f.mul(left_scale, right_scale))
        return TensorScheme(self.target,
                            [[f.mul(left_scale, x) for x in row] for row in self.a],
                            [[f.mul(right_scale, x) for x in row] for row in self.b],
                            [[f.mul(inverse, x) for x in row] for row in self.c])

    def to_matrix_scheme(self, n, m, k):
        """Convert only if the target is exactly the requested matrix tensor."""
        from .schemes import BilinearScheme

        self.require_exact()
        if self.target != matrix_multiplication_tensor(self.field, n, m, k):
            raise ValueError("target is not the requested matrix multiplication tensor")
        return BilinearScheme(self.field, n, m, k, self.a, self.b, self.c)

    def project_to_prime_field(self):
        """The project_basis formula: certify a decomposition of pi(target).

        pi is the constant power-basis coordinate, even for non-base-valued
        targets. This is linear projection, not descent of an unchanged tensor.
        Fix the extension and take powers before projection to avoid repeating
        its dimension-squared overhead.
        """
        self.require_exact()
        extension = self.field
        d = extension.degree
        aa, bb, cc = [], [], []
        for q, i, j in product(range(self.terms), range(d), range(d)):
            aa.append(tuple(extension.coordinates(x)[i] for x in self.a[q]))
            bb.append(tuple(extension.coordinates(x)[j] for x in self.b[q]))
            basis_product = extension.mul(extension.p**i, extension.p**j)
            cc.append(tuple(extension.coordinates(extension.mul(basis_product, x))[0]
                            for x in self.c[q]))
        return TensorScheme(self.target.project_to_prime_field(), aa, bb, cc)

    def descend(self):
        """Descend an unchanged tensor only when all its coefficients lie in F_p."""
        self.require_exact()
        if any(c >= self.field.p for c in self.target.coefficients):
            raise ValueError("target has non-base-field coefficients; use explicit linear projection instead")
        return self.project_to_prime_field()
