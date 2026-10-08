"""Finite determinant filtration from Determinant/Filtration.lean.

Keep X fixed; change Y and Z from polynomial-pair monomials to the actual
quotient/kernel bases. The triangular cross block has polynomial order two,
while the two matched convolution blocks have order one. All identities are
checked by exact tensor and formal polynomial coefficients in every field.
"""

from dataclasses import dataclass
from functools import cached_property
from itertools import product

from .constructions import constant_weights, nonzero_nodes
from .constraint_generators import WitnessedOrderExport
from .fields import FiniteField
from .maps import LinearMap
from .polynomials import Polynomial
from .sectors import convolution_tensor
from .tensor_degenerations import TensorDegeneration
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, shared_first_tensor
from .witnessed_constraints import DEFAULT_LIMITS, StateAssemblyLimit, WitnessedOrder


def _natural(value):
    if type(value) is not int or value < 0:
        raise ValueError("filtration degrees must be nonnegative integers")
    return value


def adapted_vectors(field, degree, *, limits=DEFAULT_LIMITS):
    """Rows q_0,...,q_(n+1),k_0,...,k_(n-1) in old (A,B) monomials.

    q_j=(X^j,0) for j<=n; q_(n+1)=(0,X^n);
    k_j=(-X^(j+1),X^j). The same basis works in characteristic two.
    """
    _natural(degree)
    size = 2*(degree+1)
    limits.maps((1, size, size), (1, size, size))
    rows = []
    for j in range(degree+1):
        rows.append(tuple(int(old == j) for old in range(size)))
    rows.append(tuple(int(old == size-1) for old in range(size)))
    for j in range(degree):
        rows.append(tuple(field.neg(1) if old == j+1 else 1 if old == degree+1+j else 0
                          for old in range(size)))
    return LinearMap(field, size, rows)


def adapted_coordinates(field, degree, *, limits=DEFAULT_LIMITS):
    """Output-dual coordinates q_i=A_i+B_(i-1), q_top=B_n, k_j=B_j."""
    _natural(degree)
    size = 2*(degree+1)
    limits.maps((1, size, size), (1, size, size))
    rows = []
    for i in range(degree+1):
        rows.append(tuple(int(old == i or (i > 0 and old == degree+i)) for old in range(size)))
    rows.append(tuple(int(old == size-1) for old in range(size)))
    for j in range(degree):
        rows.append(tuple(int(old == degree+1+j) for old in range(size)))
    return LinearMap(field, size, rows)


def _transpose(mapping):
    return LinearMap(mapping.field, mapping.output_size,
                     tuple(tuple(row[column] for row in mapping.rows) for column in range(mapping.input_size)))


@dataclass(frozen=True)
class DeterminantFiltration:
    field: FiniteField
    d: int
    e: int
    limits: object = DEFAULT_LIMITS

    def __post_init__(self):
        if not isinstance(self.field, FiniteField):
            raise ValueError("filtration requires a represented finite field")
        _natural(self.d)
        _natural(self.e)
        self.limits.tensor(self.shape)
        self.limits.maps(self.shape, self.shape)

    @property
    def shape(self):
        return self.d+1, 2*(self.e+1), 2*(self.d+self.e+1)

    @cached_property
    def source(self):
        convolution = convolution_tensor(self.field, self.d+1, self.e+1)
        return shared_first_tensor((convolution, convolution))

    def _coefficient(self, i, j, k, include_cross):
        n = self.d+self.e
        if j < self.e+2:
            if k < n+2:
                return int(i+j == k)
            return int(include_cross and j == self.e+1 and i < self.d and self.e+i == k-(n+2))
        if k < n+2:
            return 0
        return int(i+j-(self.e+2) == k-(n+2))

    @cached_property
    def adapted(self):
        return FiniteTensor.from_function(self.field, self.shape,
                                          lambda i, j, k: self._coefficient(i, j, k, True))

    @cached_property
    def graded(self):
        return FiniteTensor.from_function(self.field, self.shape,
                                          lambda i, j, k: self._coefficient(i, j, k, False))

    @cached_property
    def cross(self):
        return FiniteTensor.from_function(self.field, self.shape,
                                          lambda i, j, k: self.field.sub(self.adapted.coefficient(i, j, k),
                                                                        self.graded.coefficient(i, j, k)))

    @cached_property
    def basis_maps(self):
        return (LinearMap.identity(self.field, self.d+1),
                adapted_vectors(self.field, self.e, limits=self.limits),
                adapted_coordinates(self.field, self.d+self.e, limits=self.limits))

    @cached_property
    def diagonal_maps(self):
        f = self.field
        return ((Polynomial.constant(f, 1),)*(self.d+1),
                (Polynomial.monomial(f, 1),)*(self.e+2)+(Polynomial.constant(f, 1),)*self.e,
                (Polynomial.constant(f, 1),)*(self.d+self.e+2)+(Polynomial.monomial(f, 1),)*(self.d+self.e))

    @cached_property
    def polynomial_maps(self):
        """Actual original-source-to-graded maps, including the basis change."""
        return tuple(tuple(tuple(scale*Polynomial.constant(self.field, coefficient) for coefficient in row)
                           for scale, row in zip(diagonal, mapping.rows))
                     for diagonal, mapping in zip(self.diagonal_maps, self.basis_maps))

    def polynomial_coefficient(self, i, j, k):
        """Contract actual source coefficients with the polynomial Y,Z maps."""
        self.source.coefficient(i, j, k)  # bounds check; source/adapted axis widths agree
        result = Polynomial(self.field)
        middle, last = self.polynomial_maps[1][j], self.polynomial_maps[2][k]
        for y, by in enumerate(middle):
            if by.degree < 0:
                continue
            for z, cz in enumerate(last):
                coefficient = self.source.coefficient(i, y, z)
                if coefficient and cz.degree >= 0:
                    result = result+by*cz*Polynomial.constant(self.field, coefficient)
        return result

    @cached_property
    def _first_mismatch(self):
        for degree in (self.e, self.d+self.e):
            vectors = adapted_vectors(self.field, degree, limits=self.limits)
            coordinates = adapted_coordinates(self.field, degree, limits=self.limits)
            basis = _transpose(vectors)
            identity = LinearMap.identity(self.field, basis.input_size)
            if coordinates.compose(basis) != identity or basis.compose(coordinates) != identity:
                return "quotient/kernel basis and dual coordinate maps are not inverse"
        if self.source.restrict(*self.basis_maps) != self.adapted:
            return "original polynomial multiplication does not reconstruct the adapted tensor"
        for matrix, bound in zip(self.polynomial_maps, (0, 1, 1)):
            if any(p.degree > bound for row in matrix for p in row):
                return "local polynomial map exceeds its formal degree bound"
        for coordinate in product(*(range(size) for size in self.shape)):
            expected = Polynomial(self.field, (0, self.graded.coefficient(*coordinate), self.cross.coefficient(*coordinate)))
            if self.polynomial_coefficient(*coordinate) != expected:
                return f"formal polynomial is not t*graded+t^2*cross at {coordinate}"
        return None

    def verify(self):
        return self._first_mismatch is None

    def require_valid(self):
        if self._first_mismatch is not None:
            raise ValueError(f"invalid determinant filtration: {self._first_mismatch}")

    def generate_degeneration(self, source_scheme=None):
        """Apply the polynomial local maps to a supplied exact source scheme."""
        self.require_valid()
        if source_scheme is None:
            upper = self.shape[0]*self.shape[1]
            self.limits.maps((upper,)*3, self.shape)
            source_scheme = TensorScheme.from_tensor(self.source)
        if not isinstance(source_scheme, TensorScheme) or source_scheme.target != self.source:
            raise ValueError("source scheme must target the exact original two-copy convolution")
        source_scheme.require_exact()
        self.limits.maps((source_scheme.terms,)*3, self.shape)
        families = []
        for matrix, family in zip(self.polynomial_maps, (source_scheme.a, source_scheme.b, source_scheme.c)):
            transformed = []
            for vector in family:
                output = []
                for row in matrix:
                    value = Polynomial(self.field)
                    for polynomial, scalar in zip(row, vector):
                        if scalar:
                            value = value+polynomial*Polynomial.constant(self.field, scalar)
                    output.append(value)
                transformed.append(tuple(output))
            families.append(tuple(transformed))
        result = TensorDegeneration(self.graded, 1, *families)
        result.require_valid()
        for coordinate in product(*(range(size) for size in self.shape)):
            if result.tensor_polynomial(*coordinate) != self.polynomial_coefficient(*coordinate):
                raise ValueError("generated coefficient families do not realize the actual polynomial maps")
        return result

    def export_order(self, source_key, graded_key, *, nodes=None, limits=None):
        """An actual ordinary row: len(nodes) original copies -> one graded tensor.

        At least two distinct nonzero nodes are required. This finite row
        retains interpolation overhead and preserves the shared first leg.
        """
        limits = self.limits if limits is None else limits
        if any(not isinstance(key, str) or not key for key in (source_key, graded_key)) or source_key == graded_key:
            raise ValueError("source and graded IDs must be fresh distinct nonempty strings")
        if limits.max_blocks < 2:
            raise StateAssemblyLimit("interpolation block cap reached; construction incomplete")
        nodes = nonzero_nodes(self.field, 2) if nodes is None else tuple(nodes)
        if len(nodes) > limits.max_blocks:
            raise StateAssemblyLimit("interpolation block cap reached; construction incomplete")
        if len(nodes) < 2:
            raise ValueError("determinant filtration requires two distinct nonzero nodes")
        source_shape = tuple(len(nodes)*size for size in self.shape)
        limits.tensor(source_shape)
        limits.tensor(self.shape)
        limits.maps(source_shape, self.shape)
        weights = constant_weights(self.field, nodes)
        self.require_valid()
        maps = []
        for axis, matrix in enumerate(self.polynomial_maps):
            rows = []
            for row in matrix:
                blocks = []
                for node, weight in zip(nodes, weights):
                    scale = self.field.div(weight, node) if axis == 2 else 1
                    blocks.extend(self.field.mul(scale, p.evaluate(node)) for p in row)
                rows.append(tuple(blocks))
            maps.append(LinearMap(self.field, source_shape[axis], rows))
        row = WitnessedOrder((source_key,)*len(nodes), (graded_key,), tuple(maps), "determinant-filtration")
        export = WitnessedOrderExport(row, ((source_key, self.source), (graded_key, self.graded)), nodes, 1)
        export.require_valid(limits)
        return export
