"""The concrete three-sector convolution degeneration from the Lean proof.

Direct counterparts of AuxiliarySeparation/Convolution/Basic.lean,
Sector/Weights.lean, and Sector/Degeneration.lean. These are finite coefficient
constructions, not the spectral argument or a generator for 9/4 matrix schemes.
"""

from dataclasses import dataclass
from functools import cached_property
from itertools import product

from .constructions import constant_weights, nonzero_nodes
from .fields import FiniteField
from .polynomials import Polynomial
from .tensor_degenerations import TensorDegeneration
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor


BRANCHES = ("left", "middle", "right")


def _natural(value):
    if type(value) is not int or value < 0:
        raise ValueError("dimensions and parameters must be nonnegative integers")
    return value


def convolution_tensor(field, a, b):
    """C(a,b)[i,j,k] = 1 when i+j=k, including empty input axes."""
    _natural(a)
    _natural(b)
    return FiniteTensor.from_function(field, (a, b, max(0, a+b-1)),
                                      lambda i, j, k: int(i+j == k))


@dataclass(frozen=True)
class ThreeSectorConstruction:
    """Weights and diagonal polynomial maps for C(a,3*h+a-1).

    Natural-number subtraction is truncated, as in the Lean definitions.
    All branch tensors retain the original first input axis of size a.
    Positive a,h give three nonempty branches; zero parameters are supported
    as boundary cases, without asserting three nonempty regions.
    """

    field: FiniteField
    a: int
    h: int

    def __post_init__(self):
        if not isinstance(self.field, FiniteField):
            raise ValueError("field must be an explicitly represented FiniteField")
        _natural(self.a)
        _natural(self.h)

    @property
    def source_width(self):
        return max(0, 3*self.h+self.a-1)

    @property
    def right_start(self):
        return max(0, 2*self.h+self.a-1)

    @property
    def middle_z_start(self):
        return max(0, self.h+self.a-1)

    @property
    def shape(self):
        return (self.a, self.source_width, max(0, self.a+self.source_width-1))

    @cached_property
    def source(self):
        return convolution_tensor(self.field, self.a, self.source_width)

    @cached_property
    def weights(self):
        """Unshifted signed integer weights on X, Y, Z, respectively."""
        return ((0,)*self.a,
                tuple(int(self.h <= j < self.right_start) for j in range(self.shape[1])),
                tuple(-int(self.middle_z_start <= k < self.right_start) for k in range(self.shape[2])))

    def total_weight(self, i, j, k):
        self.source.coefficient(i, j, k)  # validates the ambient coordinates
        return sum(weights[index] for weights, index in zip(self.weights, (i, j, k)))

    @cached_property
    def retained(self):
        return FiniteTensor.from_function(self.field, self.shape,
                                          lambda i, j, k: self.source.coefficient(i, j, k)
                                          if self.total_weight(i, j, k) == 0 else 0)

    @cached_property
    def erased(self):
        return FiniteTensor.from_function(self.field, self.shape,
                                          lambda i, j, k: self.source.coefficient(i, j, k)
                                          if self.total_weight(i, j, k) == 1 else 0)

    def branch_embedding(self, branch, u, r, s):
        """Embed coordinates of C(a,h), with its output coordinate s.

        The middle branch exchanges the last two legs and reverses the first.
        This method also accepts off-support coordinates s != u+r, so the
        entire local coefficient space, not just its nonzero support, is mapped.
        """
        if branch not in BRANCHES:
            raise ValueError("branch must be left, middle, or right")
        for index, size in ((u, self.a), (r, self.h), (s, max(0, self.a+self.h-1))):
            if type(index) is not int or not 0 <= index < size:
                raise ValueError("branch coordinate is out of range")
        if branch == "left":
            return u, r, s
        if branch == "middle":
            return self.a-1-u, self.h+s, self.middle_z_start+r
        return u, self.right_start+r, self.right_start+s

    @cached_property
    def _branch_supports(self):
        return tuple(frozenset(self.branch_embedding(branch, u, r, u+r)
                               for u, r in product(range(self.a), range(self.h)))
                     for branch in BRANCHES)

    def branch_support(self, branch):
        if branch not in BRANCHES:
            raise ValueError("branch must be left, middle, or right")
        return self._branch_supports[BRANCHES.index(branch)]

    @cached_property
    def branches(self):
        """Left, middle, right tensors in the SAME original ambient spaces."""
        return tuple(FiniteTensor.from_function(self.field, self.shape,
                                                lambda i, j, k, support=support: int((i, j, k) in support))
                     for support in self._branch_supports)

    @cached_property
    def diagonal_maps(self):
        """Diagonal entries of the local polynomial maps, after shifting Z by 1."""
        shifted = (self.weights[0], self.weights[1], tuple(w+1 for w in self.weights[2]))
        return tuple(tuple(Polynomial.monomial(self.field, w) for w in weights)
                     for weights in shifted)

    @cached_property
    def local_maps(self):
        """Full local map matrices L[output][input], with degree bounds (0,1,1)."""
        zero = Polynomial(self.field)
        return tuple(tuple(tuple(diagonal[i] if i == j else zero for j in range(len(diagonal)))
                           for i in range(len(diagonal))) for diagonal in self.diagonal_maps)

    def polynomial_coefficient(self, i, j, k):
        """Apply the diagonal local maps to the source, without using R or E."""
        coefficient = Polynomial.constant(self.field, self.source.coefficient(i, j, k))
        for diagonal, index in zip(self.diagonal_maps, (i, j, k)):
            coefficient = coefficient*diagonal[index]
        return coefficient

    @cached_property
    def _first_mismatch(self):
        for leg, (matrix, bound) in enumerate(zip(self.local_maps, (0, 1, 1))):
            if any(p.degree > bound for row in matrix for p in row):
                return f"local map on leg {leg} exceeds its degree bound"
        for i, j, k in product(*(range(d) for d in self.shape)):
            coordinate = i, j, k
            source = self.source.coefficient(*coordinate)
            retained = self.retained.coefficient(*coordinate)
            erased = self.erased.coefficient(*coordinate)
            if source and self.total_weight(*coordinate) not in (0, 1):
                return f"supported weight is not zero or one at {coordinate}"
            # Check disjointness as integers, even in characteristic two.
            branch_count = sum(int(coordinate in support) for support in self._branch_supports)
            if branch_count > 1 or retained != branch_count:
                return f"retained support is not the disjoint union of the branches at {coordinate}"
            if source != self.field.add(retained, erased) or (retained and erased):
                return f"source is not partitioned into retained and erased at {coordinate}"
            expected = Polynomial(self.field, (0, retained, erased))
            if self.polynomial_coefficient(*coordinate) != expected:
                return f"polynomial is not t*retained+t^2*erased at {coordinate}"
        return None

    def verify(self):
        """Check the region partition, support weights, and formal identity everywhere."""
        return self._first_mismatch is None

    def require_valid(self):
        if self._first_mismatch is not None:
            raise ValueError(f"invalid three-sector construction: {self._first_mismatch}")

    def generate_degeneration(self, source_scheme=None):
        """Apply the proof's diagonal maps to an exact decomposition of the source.

        Defaults to a trivial coordinate-pair decomposition. This preserves its
        number of terms; it does not find the best convolution rank. The result
        is a certified TensorDegeneration with target R and leading order 1.
        Use recover_power() to power this family before interpolated recovery.
        """
        self.require_valid()
        if source_scheme is None:
            source_scheme = TensorScheme.from_tensor(self.source)
        if not isinstance(source_scheme, TensorScheme) or source_scheme.target != self.source:
            raise ValueError("source scheme must target this construction's exact convolution tensor")
        source_scheme.require_exact()
        families = tuple(tuple(tuple(Polynomial.constant(self.field, x)*diagonal[i]
                                     for i, x in enumerate(row)) for row in family)
                         for family, diagonal in zip((source_scheme.a, source_scheme.b, source_scheme.c),
                                                     self.diagonal_maps))
        degeneration = TensorDegeneration(self.retained, 1, *families)
        degeneration.require_valid()
        return degeneration

    def recover_power(self, exponent, nodes=None, source_scheme=None):
        """Recover an exact decomposition of R^exponent over this fixed field.

        Power the polynomial family first, then interpolate its normalized
        leading coefficient. The proof maps give degree at most 2*exponent
        and leading order exponent, so exponent+1 nonzero nodes suffice.
        The actual source-family bound can be smaller in empty cases.

        Nodes are checked before allocating the tensor power. No field is
        automatically enlarged. To return to the prime field, descend the
        returned exact scheme once, after powering and recovery.
        """
        _natural(exponent)
        degeneration = self.generate_degeneration(source_scheme)
        count = max(0, degeneration.degree_bound-degeneration.leading)*exponent+1
        nodes = nonzero_nodes(self.field, count) if nodes is None else tuple(nodes)
        constant_weights(self.field, nodes)  # validates distinct nonzero field elements
        if len(nodes) < count:
            raise ValueError("not enough nodes for the normalized powered tensor degree bound")
        return degeneration.tensor_power(exponent).recover_scheme(nodes)

    @cached_property
    def side_labels(self):
        """Recover zero-based left/middle/right tags from ambient Y and Z."""
        return (tuple(0 if j < self.h else 1 if j < self.right_start else 2 for j in range(self.shape[1])),
                tuple(0 if k < self.middle_z_start else 1 if k < self.right_start else 2 for k in range(self.shape[2])))

    def tag_scheme(self, scheme=None):
        """Insert the two side tags of R, preserving its original first input."""
        from .tagging import tag_shared_scheme

        self.require_valid()
        if scheme is None:
            scheme = TensorScheme.from_tensor(self.retained)
        return tag_shared_scheme(scheme, self.branches, *self.side_labels)
