"""Coefficient-level schemes, block composition, and fixed-algebra descent.

These implement formulas from RecursiveBlockPrograms.lean and FieldDescent.lean.
Output coefficients use (row, column), transposing Lean's (column, row) index.
"""

from dataclasses import dataclass
from functools import cached_property
from itertools import product

from .fields import FiniteField


def _matrix(field, matrix, rows, columns):
    if len(matrix) != rows or any(len(row) != columns for row in matrix):
        raise ValueError(f"expected a {rows} by {columns} matrix")
    return tuple(field.check(x) for row in matrix for x in row)


def naive_multiply(field, left, right):
    """Independent ordinary matrix multiplication, also used by the tests."""
    if not left or not right or not right[0]:
        raise ValueError("matrix dimensions must be positive")
    n, m, k = len(left), len(right), len(right[0])
    _matrix(field, left, n, m)
    _matrix(field, right, m, k)
    return [[field.dot(left[i], [right[j][h] for j in range(m)])
             for h in range(k)] for i in range(n)]


@dataclass(frozen=True)
class BilinearScheme:
    """AB = sum_q c[q] * dot(a[q], A) * dot(b[q], B).

    a, b, c are row-major coefficient vectors. terms is a decomposition
    length / rank upper bound, not the minimal tensor rank or operation count.
    """

    field: FiniteField
    n: int
    m: int
    k: int
    a: tuple[tuple[int, ...], ...]
    b: tuple[tuple[int, ...], ...]
    c: tuple[tuple[int, ...], ...]

    def __post_init__(self):
        if any(type(d) is not int or d < 1 for d in (self.n, self.m, self.k)):
            raise ValueError("matrix dimensions must be positive integers")
        for name in ("a", "b", "c"):
            object.__setattr__(self, name, tuple(tuple(row) for row in getattr(self, name)))
        if not len(self.a) == len(self.b) == len(self.c):
            raise ValueError("coefficient families must have equal term counts")
        for family, size in ((self.a, self.n * self.m),
                             (self.b, self.m * self.k), (self.c, self.n * self.k)):
            if any(len(row) != size for row in family):
                raise ValueError("coefficient vector has the wrong length")
            for row in family:
                for x in row:
                    self.field.check(x)

    @property
    def terms(self):
        return len(self.a)

    @cached_property
    def _first_mismatch(self):
        # Accumulate only nonzero contributions, then check EVERY tensor
        # coordinate. This certifies a formal bilinear identity even over F_2.
        f = self.field
        coefficients = {}
        for a, b, c in zip(self.a, self.b, self.c):
            aa = [(i, x) for i, x in enumerate(a) if x]
            bb = [(j, y) for j, y in enumerate(b) if y]
            cc = [(h, z) for h, z in enumerate(c) if z]
            for (i, x), (j, y), (h, z) in product(aa, bb, cc):
                key = i, j, h
                coefficients[key] = f.add(coefficients.get(key, 0), f.mul(f.mul(x, y), z))
        for i, j, h in product(range(self.n * self.m), range(self.m * self.k),
                               range(self.n * self.k)):
            row, inner = divmod(i, self.m)
            inner2, column = divmod(j, self.k)
            outrow, outcolumn = divmod(h, self.k)
            expected = int(row == outrow and inner == inner2 and column == outcolumn)
            actual = coefficients.get((i, j, h), 0)
            if actual != expected:
                return ((i, j, h), expected, actual)
        return None

    def verify(self):
        return self._first_mismatch is None

    def require_exact(self):
        if self._first_mismatch is not None:
            coordinate, expected, actual = self._first_mismatch
            raise ValueError(f"invalid tensor coefficient {coordinate}: expected {expected}, got {actual}")

    def multiply(self, left, right):
        """Evaluate a certified scheme on scalar matrices."""
        self.require_exact()
        f = self.field
        left = _matrix(f, left, self.n, self.m)
        right = _matrix(f, right, self.m, self.k)
        values = [f.mul(f.dot(a, left), f.dot(b, right)) for a, b in zip(self.a, self.b)]
        return [[f.sum(f.mul(c[i * self.k + j], value)
                       for c, value in zip(self.c, values))
                 for j in range(self.k)] for i in range(self.n)]

    def multiply_recursive(self, left, right):
        """The square block identity of RecursiveBlockPrograms.lean.

        Both inputs must be square of size n^depth, for n=m=k>=2.
        Each level forms linear combinations of blocks, recurses, then combines
        the products. No performance claim or arithmetic-cost instrumentation.
        """
        self.require_exact()
        if not self.n == self.m == self.k or self.n < 2:
            raise ValueError("recursion requires a square block scheme of size >= 2")
        f, size = self.field, len(left)
        _matrix(f, left, size, size)
        _matrix(f, right, size, size)
        power = size
        while power > 1 and power % self.n == 0:
            power //= self.n
        if power != 1:
            raise ValueError("input size must be a nonnegative power of the block size")

        def recurse(a, b):
            size = len(a)
            if size == 1:
                return [[f.mul(a[0][0], b[0][0])]]
            width = size // self.n

            def combine(matrix, weights):
                return [[f.sum(f.mul(weights[i * self.n + j], matrix[i * width + s][j * width + t])
                                for i in range(self.n) for j in range(self.n))
                         for t in range(width)] for s in range(width)]

            products = [recurse(combine(a, x), combine(b, y)) for x, y in zip(self.a, self.b)]
            return [[f.sum(f.mul(c[(i // width) * self.n + j // width], value[i % width][j % width])
                           for c, value in zip(self.c, products))
                     for j in range(size)] for i in range(size)]

        return recurse(left, right)

    def tensor_product(self, other):
        """Tensor two schemes and identify product coordinates with matrices."""
        if self.field != other.field:
            raise ValueError("tensor factors must use the same field representation")
        self.require_exact()
        other.require_exact()
        f = self.field

        def combine(a, b, rows1, columns1, rows2, columns2):
            return tuple(f.mul(a[i1 * columns1 + j1], b[i2 * columns2 + j2])
                         for i1 in range(rows1) for i2 in range(rows2)
                         for j1 in range(columns1) for j2 in range(columns2))

        aa, bb, cc = [], [], []
        for q, r in product(range(self.terms), range(other.terms)):
            aa.append(combine(self.a[q], other.a[r], self.n, self.m, other.n, other.m))
            bb.append(combine(self.b[q], other.b[r], self.m, self.k, other.m, other.k))
            cc.append(combine(self.c[q], other.c[r], self.n, self.k, other.n, other.k))
        return BilinearScheme(f, self.n * other.n, self.m * other.m, self.k * other.k, aa, bb, cc)

    def tensor_power(self, exponent):
        if type(exponent) is not int or exponent < 0:
            raise ValueError("tensor exponent must be a nonnegative integer")
        self.require_exact()
        result = naive_scheme(self.field, 1, 1, 1)
        for _ in range(exponent):
            result = result.tensor_product(self)
        return result

    def rescale_terms(self, left_scale, right_scale):
        """Change a decomposition's coefficients without changing its tensor."""
        f = self.field
        inverse = f.inv(f.mul(left_scale, right_scale))
        return BilinearScheme(f, self.n, self.m, self.k,
                              [[f.mul(left_scale, x) for x in a] for a in self.a],
                              [[f.mul(right_scale, x) for x in b] for b in self.b],
                              [[f.mul(inverse, x) for x in c] for c in self.c])


def naive_scheme(field, n, m, k):
    if any(type(d) is not int or d < 1 for d in (n, m, k)):
        raise ValueError("matrix dimensions must be positive integers")
    aa, bb, cc = [], [], []
    for i, j, h in product(range(n), range(m), range(k)):
        aa.append(tuple(int(x == i * m + j) for x in range(n * m)))
        bb.append(tuple(int(x == j * k + h) for x in range(m * k)))
        cc.append(tuple(int(x == i * k + h) for x in range(n * k)))
    return BilinearScheme(field, n, m, k, aa, bb, cc)


def strassen_scheme(field):
    """Known 2x2, seven-term TEST SEED; not a scheme extracted from the 9/4 proof."""
    a = ((1, 0, 0, 1), (0, 0, 1, 1), (1, 0, 0, 0), (0, 0, 0, 1),
         (1, 1, 0, 0), (-1, 0, 1, 0), (0, 1, 0, -1))
    b = ((1, 0, 0, 1), (1, 0, 0, 0), (0, 1, 0, -1), (-1, 0, 1, 0),
         (0, 0, 0, 1), (1, 1, 0, 0), (0, 0, 1, 1))
    c = ((1, 0, 0, 1), (0, 0, 1, -1), (0, 1, 0, 1), (1, 0, 1, 0),
         (-1, 1, 0, 0), (0, 0, 0, 1), (1, 0, 0, 0))
    convert = lambda family: tuple(tuple(field.embed(x) for x in row) for row in family)
    return BilinearScheme(field, 2, 2, 2, convert(a), convert(b), convert(c))


def descend(scheme):
    """Project to F_p using pi(v)=the constant power-basis coordinate.

    For each q,i,j: a'=coord_i(a_q), b'=coord_j(b_q),
    c'=pi(X^i * X^j * c_q). Thus r terms become r*d^2 terms.
    To retain a FIXED overhead across powers, call descend(scheme.tensor_power(t)),
    rather than descend(scheme).tensor_power(t). The supplied extension stays fixed.
    No minimal coefficient algebra is computed; the whole supplied field is used.
    """
    scheme.require_exact()
    extension = scheme.field
    base = FiniteField(extension.p)
    d = extension.degree
    aa, bb, cc = [], [], []
    for q, i, j in product(range(scheme.terms), range(d), range(d)):
        aa.append(tuple(extension.coordinates(x)[i] for x in scheme.a[q]))
        bb.append(tuple(extension.coordinates(x)[j] for x in scheme.b[q]))
        basis_product = extension.mul(extension.p ** i, extension.p ** j)
        cc.append(tuple(extension.coordinates(extension.mul(basis_product, x))[0]
                        for x in scheme.c[q]))
    return BilinearScheme(base, scheme.n, scheme.m, scheme.k, aa, bb, cc)
