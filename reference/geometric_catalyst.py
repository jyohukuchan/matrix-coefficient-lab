"""An explicit geometric auxiliary tensor from a supplied matrix rank gap.

This is a finite algebraic bridge, not a method for finding a useful gap.
A supplied exact r-term M_(d^b) scheme with r<k^b gives an actual D=0
positive-gain catalyst. Intermediate matrix blocks match by coordinate
bijections; no tensor summands are cancelled.
"""

from dataclasses import dataclass

from .catalyst_compiler import DEFAULT_LIMITS, CompilerResourceLimit
from .catalysts import CatalyticCertificate
from .maps import LinearMap
from .schemes import BilinearScheme, naive_scheme
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor


@dataclass(frozen=True)
class GeometricCatalystConstruction:
    certificate: CatalyticCertificate
    auxiliary_scheme: TensorScheme
    exponent: int
    supplied_terms: int
    block_counts: tuple[int, ...]
    matrix_sizes: tuple[int, ...]
    block_offsets: tuple[int, ...]

    @property
    def source(self):
        return self.certificate.source

    @property
    def target(self):
        return self.certificate.target

    @property
    def auxiliary(self):
        return self.certificate.S

    @property
    def gain(self):
        return self.certificate.m


def build_geometric_catalyst(d, k, exponent, scheme, *, limits=DEFAULT_LIMITS):
    """Build S=sum_i k^(b-1-i) copies M_(d^i), m=k^b-r, and actual maps.

    In k*S, k^b scalar blocks supply m gain units and the supplied r-term
    top matrix decomposition. Each remaining middle M_(d^j) block is
    reindexed to M_d tensor M_(d^(j-1)), then all blocks are packed into
    M_d tensor S. The auxiliary decomposition returned alongside the
    certificate is an explicit conservative coordinate decomposition.
    """
    if type(d) is not int or d < 2 or type(k) is not int or k < 1:
        raise ValueError("require integer d>=2 and k>=1")
    if type(exponent) is not int or exponent < 1:
        raise ValueError("geometric exponent must be a positive integer")
    if not isinstance(scheme, BilinearScheme):
        raise ValueError("supply an exact square matrix scheme")
    # The target contains M_(d^b), hence at least 2^(6b) coefficients.
    # Reject obviously huge requests before powers or exponent-length loops.
    if 6*exponent >= limits.max_tensor_coefficients.bit_length():
        raise CompilerResourceLimit("geometric tensor limit reached; construction incomplete")
    top_size = d**exponent
    if not scheme.n == scheme.m == scheme.k == top_size:
        raise ValueError("supplied scheme must target precisely M_(d^exponent)")
    scalar_count = k**exponent
    if scheme.terms >= scalar_count:
        raise ValueError("supplied rank upper bound must be strictly smaller than k^exponent")
    gain = scalar_count-scheme.terms
    counts = tuple(k**(exponent-1-i) for i in range(exponent))
    sizes = tuple(d**i for i in range(exponent))
    offsets, dimension, auxiliary_terms = [], 0, 0
    for count, size in zip(counts, sizes):
        offsets.append(dimension)
        dimension += count*size*size
        auxiliary_terms += count*size**3
    source_shape = (k*dimension,)*3
    target_shape = (gain+d*d*dimension,)*3
    limits.maps(source_shape, target_shape)
    limits.scheme((dimension,)*3, auxiliary_terms)
    scheme.require_exact()
    f = scheme.field
    components = []
    for count, size in zip(counts, sizes):
        limits.scheme((size*size,)*3, size**3)
        component = naive_scheme(f, size, size, size).as_tensor_scheme()
        components.extend((component,)*count)
    auxiliary_scheme = TensorScheme.direct_sum(tuple(components))
    auxiliary_scheme.require_exact()
    S = auxiliary_scheme.target
    if S.shape != (dimension,)*3:
        raise ValueError("geometric auxiliary layout dimensions do not match")
    # Scalar block coordinates across the k independent source S copies.
    units = tuple(copy*dimension+coordinate for copy in range(k) for coordinate in range(counts[0]))
    if len(units) != scalar_count:
        raise ValueError("geometric scalar-copy count is inconsistent")
    families, maps = (scheme.a, scheme.b, scheme.c), []
    for family in families:
        rows = [[0]*source_shape[0] for _ in range(target_shape[0])]
        for scalar in range(gain):
            rows[scalar][units[scalar]] = 1
        for i, (count, size, offset) in enumerate(zip(counts, sizes, offsets)):
            next_size = d*size
            for copy in range(count):
                for matrix_digit in range(d*d):
                    row_d, column_d = divmod(matrix_digit, d)
                    for coordinate in range(size*size):
                        row_old, column_old = divmod(coordinate, size)
                        combined = (row_d*size+row_old)*next_size+column_d*size+column_old
                        output = gain+matrix_digit*dimension+offset+copy*size*size+coordinate
                        if i == exponent-1:
                            for term, vector in enumerate(family):
                                rows[output][units[gain+term]] = vector[combined]
                        else:
                            # Match the target's i block to the source's i+1
                            # matrix block, across all k independent S copies.
                            outer, inner = divmod(copy, counts[i+1])
                            old = outer*dimension+offsets[i+1]+inner*next_size*next_size+combined
                            rows[output][old] = 1
        maps.append(LinearMap(f, source_shape[0], rows))
    certificate = CatalyticCertificate(d, k, gain, FiniteTensor.zero(f, (0, 0, 0)), S, tuple(maps))
    certificate.require_valid()
    return GeometricCatalystConstruction(certificate, auxiliary_scheme, exponent, scheme.terms,
                                         counts, sizes, tuple(offsets))
