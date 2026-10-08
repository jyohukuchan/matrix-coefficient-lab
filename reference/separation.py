"""Finite Fourier separation and square filtering from Separation/Basic.lean.

The output is a full direct sum of branch tensors with dot-product factors.
This implements finite coefficient identities, not spectral existence or a
generator for the matrix multiplication schemes promised by the 9/4 theorem.
"""

from dataclasses import dataclass
from functools import cached_property
from itertools import product

from .constructions import (constant_weights, fourier_filter, is_primitive_root,
                            nonzero_nodes, primitive_root, separation_period)
from .maps import LinearMap
from .polynomials import Polynomial
from .tensor_degenerations import TensorDegeneration
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, _tensor_family, shared_first_tensor


def phase(g, h, u, v):
    return u-v+2*(g-h)


def square_weights(g, h, u, v):
    """Signed INTEGER weights; their sum on phase zero is (g-h)^2."""
    return g*g, h*u-h*h, -h*v


@dataclass(frozen=True)
class FiniteSeparation:
    branches: tuple[FiniteTensor, ...]
    period: int | None = None
    root: int | None = None

    def __post_init__(self):
        branches = _tensor_family(self.branches, uniform=True)
        object.__setattr__(self, "branches", branches)
        period = separation_period(self.field, self.count) if self.period is None else self.period
        if type(period) is not int or period < 5*self.count or self.field.embed(period) == 0:
            raise ValueError("separation period must be at least 5*M and invertible in the field")
        root = primitive_root(self.field, period) if self.root is None else self.root
        if not is_primitive_root(self.field, root, period):
            raise ValueError("separation root must have the required exact order")
        object.__setattr__(self, "period", period)
        object.__setattr__(self, "root", root)

    @property
    def field(self):
        return self.branches[0].field

    @property
    def count(self):
        return len(self.branches)

    @property
    def branch_shape(self):
        return self.branches[0].shape

    @property
    def shape(self):
        x, y, z = self.branch_shape
        m = self.count
        return x*m, m*y*m, m*z*m

    @cached_property
    def source(self):
        return shared_first_tensor(self.branches)

    def _decode(self, i, j, k):
        for index, size in zip((i, j, k), self.shape):
            if type(index) is not int or not 0 <= index < size:
                raise ValueError("separation coordinate is out of range")
        x, g = divmod(i, self.count)
        hy, u = divmod(j, self.count)
        hz, v = divmod(k, self.count)
        h, y = divmod(hy, self.branch_shape[1])
        h2, z = divmod(hz, self.branch_shape[2])
        return x, y, z, g+1, h+1, h2+1, u+1, v+1

    @cached_property
    def target(self):
        def coefficient(i, j, k):
            x, y, z, g, h, h2, u, v = self._decode(i, j, k)
            return self.branches[h-1].coefficient(x, y, z) if g == h == h2 and u == v else 0
        return FiniteTensor.from_function(self.field, self.shape, coefficient)

    @cached_property
    def projected(self):
        """Independent expected support after Fourier averaging, before weights."""
        def coefficient(i, j, k):
            x, y, z, g, h, h2, u, v = self._decode(i, j, k)
            return self.branches[h-1].coefficient(x, y, z) if h == h2 and phase(g, h, u, v) == 0 else 0
        return FiniteTensor.from_function(self.field, self.shape, coefficient)

    @cached_property
    def projection_maps(self):
        """The proof's three maps from period complete independent source copies."""
        f, m, length = self.field, self.count, self.period
        x, y, z = self.source.shape
        inverse = f.inv(f.embed(length))
        first = LinearMap(f, length*x, tuple(tuple(f.pow(self.root, r*(2*(g+1))) if old == i else 0
                                                  for r in range(length) for old in range(x))
                                             for i in range(x) for g in range(m)))
        second = LinearMap(f, length*y, tuple(tuple(f.pow(self.root, r*((u+1)-(j//self.branch_shape[1]+1))) if old == j else 0
                                                   for r in range(length) for old in range(y))
                                              for j in range(y) for u in range(m)))
        third = LinearMap(f, length*z, tuple(tuple(f.mul(inverse, f.pow(self.root, r*(-(v+1)-(k//self.branch_shape[2]+1)))) if old == k else 0
                                                  for r in range(length) for old in range(z))
                                             for k in range(z) for v in range(m)))
        return first, second, third

    @cached_property
    def weights(self):
        x, y, z = self.source.shape
        m = self.count
        return (tuple((g+1)**2 for _ in range(x) for g in range(m)),
                tuple((j//self.branch_shape[1]+1)*(u+1)-(j//self.branch_shape[1]+1)**2
                      for j in range(y) for u in range(m)),
                tuple(-(k//self.branch_shape[2]+1)*(v+1) for k in range(z) for v in range(m)))

    @property
    def shifts(self):
        # Uniform shifts work even when an ambient axis is empty.
        m = self.count
        return 0, m*(m-1), m*m

    @property
    def leading(self):
        return sum(self.shifts)

    @property
    def normalized_degree_bound(self):
        return (self.count-1)**2

    @cached_property
    def _fourier_values(self):
        bound = 3*(self.count-1)
        return {e: fourier_filter(self.field, self.period, self.root, e) for e in range(-bound, bound+1)}

    def polynomial_coefficient(self, i, j, k):
        """Actual Fourier average followed by shifted integer square weights."""
        x, y, z, g, h, h2, u, v = self._decode(i, j, k)
        coefficient = 0 if h != h2 else self.field.mul(self.branches[h-1].coefficient(x, y, z),
                                                       self._fourier_values[phase(g, h, u, v)])
        degree = self.leading+sum(weights[index] for weights, index in zip(self.weights, (i, j, k)))
        return Polynomial.monomial(self.field, degree, coefficient)

    @cached_property
    def _first_mismatch(self):
        for g, h, u, v in product(range(1, self.count+1), repeat=4):
            e = phase(g, h, u, v)
            if abs(e) >= self.period or self._fourier_values[e] != int(e == 0):
                return "Fourier support is not the integer phase-zero equation"
            if e == 0 and sum(square_weights(g, h, u, v)) != (g-h)**2:
                return "surviving weight is not the stated square"
        for coordinate in product(*(range(d) for d in self.shape)):
            p = self.polynomial_coefficient(*coordinate)
            if p.coefficient(self.leading) != self.target.coefficient(*coordinate):
                return f"leading coefficient is not the separated target at {coordinate}"
            if p.degree > self.leading+self.normalized_degree_bound:
                return f"square-filtered tensor exceeds its degree bound at {coordinate}"
            if any(p.coefficients[:self.leading]):
                return f"polynomial has a lower coefficient at {coordinate}"
        return None

    def verify(self):
        return self._first_mismatch is None

    def require_valid(self):
        if self._first_mismatch is not None:
            raise ValueError(f"invalid finite separation: {self._first_mismatch}")

    def projection_scheme(self, source_scheme=None):
        self.require_valid()
        if source_scheme is None:
            source_scheme = TensorScheme.from_tensor(self.source)
        if not isinstance(source_scheme, TensorScheme) or source_scheme.target != self.source:
            raise ValueError("source scheme must target the exact shared-input branch tensor")
        source_scheme.require_exact()
        copies = TensorScheme.direct_sum((source_scheme,)*self.period)
        result = copies.restrict(*self.projection_maps)
        if result.target != self.projected:
            raise ValueError("local Fourier maps do not realize the expected projection")
        return result

    def weighted_maps(self, node):
        """Unshifted proof maps at a NONZERO parameter, permitting negative powers."""
        self.field.check(node)
        if node == 0:
            raise ValueError("unshifted square weights require a nonzero parameter")
        return tuple(LinearMap(self.field, m.input_size,
                               tuple(tuple(self.field.mul(self.field.pow(node, weight), value) for value in row)
                                     for row, weight in zip(m.rows, weights)))
                     for m, weights in zip(self.projection_maps, self.weights))

    def generate_degeneration(self, source_scheme=None):
        projected = self.projection_scheme(source_scheme)
        families = tuple(tuple(tuple(Polynomial.monomial(self.field, weight+shift, coefficient)
                                     for weight, coefficient in zip(weights, row)) for row in family)
                         for weights, shift, family in zip(self.weights, self.shifts, (projected.a, projected.b, projected.c)))
        result = TensorDegeneration(self.target, self.leading, *families)
        result.require_valid()
        for coordinate in product(*(range(d) for d in self.shape)):
            if result.tensor_polynomial(*coordinate) != self.polynomial_coefficient(*coordinate):
                raise ValueError(f"weighted families do not realize the formal separation polynomial at {coordinate}")
        return result

    def recover_power(self, exponent, nodes=None, source_scheme=None):
        if type(exponent) is not int or exponent < 0:
            raise ValueError("tensor exponent must be a nonnegative integer")
        count = exponent*self.normalized_degree_bound+1
        nodes = nonzero_nodes(self.field, count) if nodes is None else tuple(nodes)
        constant_weights(self.field, nodes)
        if len(nodes) < count:
            raise ValueError("not enough nodes for the square-filtered tensor degree bound")
        degeneration = self.generate_degeneration(source_scheme)
        powered = degeneration if exponent == 1 else degeneration.tensor_power(exponent)
        bound = exponent*(self.leading+self.normalized_degree_bound)
        return powered.recover_scheme(nodes, degree_bound=bound)

    def recover_scheme(self, nodes=None, source_scheme=None):
        return self.recover_power(1, nodes, source_scheme)

    @cached_property
    def direct_sum_target(self):
        dot = FiniteTensor.from_function(self.field, (1, self.count, self.count), lambda i, j, k: int(j == k))
        return FiniteTensor.direct_sum(b.tensor_product(dot) for b in self.branches)

    def as_direct_sum(self, scheme):
        """Reorder X so the separated tensor is literally sum_h B_h tensor dot_M."""
        if not isinstance(scheme, TensorScheme) or scheme.target != self.target:
            raise ValueError("scheme must target this construction's separated tensor")
        x = self.branch_shape[0]
        result = scheme.restrict(LinearMap.selection(self.field, self.shape[0],
                                                     (i*self.count+h for h in range(self.count) for i in range(x))),
                                 LinearMap.identity(self.field, self.shape[1]),
                                 LinearMap.identity(self.field, self.shape[2]))
        if result.target != self.direct_sum_target:
            raise ValueError("separated target does not match the branch/dot-product direct sum")
        return result
