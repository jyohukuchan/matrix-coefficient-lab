"""Exact finite controls for three singleton orientations of length-t dots.

These constructors and coefficient witnesses are independent of the analytic
catalyst exclusions in research/variable-auxiliaries.md. Slice ranks below are
computed over the represented field; none is a generic-rank oracle or a
positive-gain catalytic restriction.
"""

from dataclasses import dataclass
from itertools import product

from .catalyst_search_experiment import field_matrix_rank
from .fields import FiniteField
from .subrank_constraints import export_selected_subrank_order
from .tensors import FiniteTensor, matrix_multiplication_tensor
from .witnessed_constraints import DEFAULT_LIMITS, StateAssemblyLimit, StateAssemblyLimits


def _parameters(field, t, limits):
    if not isinstance(field, FiniteField) or type(t) is not int or t < 1:
        raise ValueError("supply a represented field and positive integer dot length")
    if not isinstance(limits, StateAssemblyLimits):
        raise ValueError("supply explicit state assembly limits")


def _tensor_cap(shape, limits):
    limits.tensor(shape)
    if sum(shape) > limits.max_constraint_entries:
        raise StateAssemblyLimit("mixed-dot coordinate cap reached; controls incomplete")


def _names(first, second):
    if any(not isinstance(key, str) or not key for key in (first, second)) or first == second:
        raise ValueError("source and unit IDs must be distinct nonempty strings")


def mixed_dot_length_tensor(field, t, *, limits=DEFAULT_LIMITS):
    """S_t=DotX_t direct_sum DotY_t direct_sum DotZ_t, shape (1+2t)^3."""
    _parameters(field, t, limits)
    _tensor_cap((1+2*t,)*3, limits)
    if limits.max_blocks < 3:
        raise StateAssemblyLimit("mixed-dot orientation block cap reached; controls incomplete")
    dot = FiniteTensor.from_function(field, (1, t, t), lambda i, j, k: int(j == k))
    return FiniteTensor.direct_sum((dot, dot.permute_axes((1, 0, 2)), dot.permute_axes((1, 2, 0))))


def export_mixed_dot_length_subrank_order(field, t, source_key, *, unit_key="unit", limits=DEFAULT_LIMITS):
    """Three independent units <= S_t by coordinate selectors, all coefficients checked."""
    _parameters(field, t, limits)
    _names(source_key, unit_key)
    shape = (1+2*t,)*3
    _tensor_cap(shape, limits)
    limits.tensor((3,)*3)
    limits.maps(shape, (3,)*3)
    if sum(shape)+9 > limits.max_constraint_entries:
        raise StateAssemblyLimit("subrank selector bookkeeping cap reached; construction incomplete")
    if limits.max_blocks < 3:
        raise StateAssemblyLimit("subrank block cap reached; construction incomplete")
    return export_selected_subrank_order(mixed_dot_length_tensor(field, t, limits=limits),
                                         ((0, 1, t+1), (0, t, t+1), (0, t, 2*t)), source_key,
                                         unit_key=unit_key, limits=limits)


def export_mixed_dot_length_square_subrank_order(field, t, source_key, *, unit_key="unit", limits=DEFAULT_LIMITS):
    """(3+6t) units <= S_t^2 on disjoint orientation products; not optimality."""
    _parameters(field, t, limits)
    _names(source_key, unit_key)
    size, count = 1+2*t, 3+6*t
    shape = (size*size,)*3
    _tensor_cap(shape, limits)
    limits.tensor((count,)*3)
    limits.maps(shape, (count,)*3)
    if limits.max_blocks < count:
        raise StateAssemblyLimit("square subrank block cap reached; construction incomplete")
    if sum(shape)+3*count > limits.max_constraint_entries:
        raise StateAssemblyLimit("square selector bookkeeping cap reached; construction incomplete")
    offsets = ((0, 1, t+1), (0, t, t+1), (0, t, 2*t))
    selections = [[], [], []]
    for first, second in product(range(3), repeat=2):
        # Equal orientations have a singleton common leg: take one unit.
        # Unequal orientations permit all t diagonal pairs (r,r).
        for r in range(1 if first == second else t):
            for axis in range(3):
                a = offsets[axis][first]+(0 if axis == first else r)
                b = offsets[axis][second]+(0 if axis == second else r)
                selections[axis].append(a*size+b)
    source = mixed_dot_length_tensor(field, t, limits=limits).tensor_power(2)
    return export_selected_subrank_order(source, selections, source_key, unit_key=unit_key, limits=limits)


def _matrix(field, matrix, shape):
    matrix = tuple(tuple(row) for row in matrix)
    if len(matrix) != shape[0] or any(len(row) != shape[1] for row in matrix):
        raise ValueError(f"core matrix must have shape {shape}")
    for row in matrix:
        for value in row:
            field.check(value)
    return matrix


def mixed_dot_length_core_input(field, t, A, B, C, *, limits=DEFAULT_LIMITS):
    """Encode A(2x2), B((2t)x2), C(2x(2t)) in M2 tensor S_t first coordinates."""
    _parameters(field, t, limits)
    size = 1+2*t
    if 4*size > limits.max_constraint_entries:
        raise StateAssemblyLimit("core input bookkeeping cap reached; controls incomplete")
    A, B, C = (_matrix(field, matrix, shape) for matrix, shape in
               zip((A, B, C), ((2, 2), (2*t, 2), (2, 2*t))))
    vector = [0]*(4*size)
    for i, j in product(range(2), repeat=2):
        coordinate = (2*i+j)*size
        vector[coordinate] = A[i][j]
        for r in range(t):
            vector[coordinate+1+r] = B[t*i+r][j]
            vector[coordinate+1+t+r] = C[i][t*j+r]
    return tuple(vector)


def exact_length_slice_rank(tensor, vector, *, limits=DEFAULT_LIMITS):
    """Coefficient contraction followed by exact field matrix elimination."""
    if not isinstance(tensor, FiniteTensor) or not isinstance(limits, StateAssemblyLimits):
        raise ValueError("supply an exact tensor and state assembly limits")
    _tensor_cap(tensor.shape, limits)
    x, y, z = tensor.shape
    if y*z > limits.max_map_entries:
        raise StateAssemblyLimit("slice matrix cap reached; controls incomplete")
    vector = tuple(vector)
    if len(vector) != x:
        raise ValueError("slice input has the wrong coordinate count")
    for value in vector:
        tensor.field.check(value)
    f = tensor.field
    support = tuple((i, value) for i, value in enumerate(vector) if value)
    rows = tuple(tuple(f.sum(f.mul(value, coefficient) for i, value in support
                             if (coefficient := tensor.coefficient(i, j, k)))
                       for k in range(z)) for j in range(y))
    return field_matrix_rank(f, rows, z)


def _auxiliary(field, t, supplied, limits):
    expected = mixed_dot_length_tensor(field, t, limits=limits)
    if supplied is not None and supplied != expected:
        raise ValueError("supplied auxiliary does not match the actual three-oriented dots")
    return expected if supplied is None else supplied


def source_slice_control(field, t, vector, *, auxiliary=None, limits=DEFAULT_LIMITS):
    """Check rank(S_t(s,u,v))=t*[s!=0]+[u!=0]+[v!=0] on an actual input."""
    _parameters(field, t, limits)
    tensor = _auxiliary(field, t, auxiliary, limits)
    vector = tuple(vector)
    actual = exact_length_slice_rank(tensor, vector, limits=limits)
    expected = t*bool(vector[0])+bool(any(vector[1:1+t]))+bool(any(vector[1+t:]))
    if actual != expected:
        raise ValueError("actual source slice disagrees with the block rank formula")
    return actual


def core_slice_control(field, t, A, B, C, *, auxiliary=None, limits=DEFAULT_LIMITS):
    """Check rank=2t*rank(A)+2*rank(B)+2*rank(C) coefficientwise on this input."""
    _parameters(field, t, limits)
    size = 4*(1+2*t)
    _tensor_cap((size,)*3, limits)
    if size*size > limits.max_map_entries:
        raise StateAssemblyLimit("core slice matrix cap reached; controls incomplete")
    A, B, C = (_matrix(field, matrix, shape) for matrix, shape in
               zip((A, B, C), ((2, 2), (2*t, 2), (2, 2*t))))
    vector = mixed_dot_length_core_input(field, t, A, B, C, limits=limits)
    tensor = matrix_multiplication_tensor(field, 2, 2, 2).tensor_product(_auxiliary(field, t, auxiliary, limits))
    actual = exact_length_slice_rank(tensor, vector, limits=limits)
    expected = (2*t*field_matrix_rank(field, A, 2)+2*field_matrix_rank(field, B, 2)
                +2*field_matrix_rank(field, C, 2*t))
    if actual != expected:
        raise ValueError("actual core slice disagrees with the matrix rank formula")
    return actual


def mixed_dot_length_rank_profiles(field, t, *, limits=DEFAULT_LIMITS):
    """27 deterministic actual inputs realizing every rank triple (0,1,2)^3.

    This is a finite coefficient check, not a bound on arbitrary matrix spaces.
    Extension-field controls use a genuinely non-prime-field diagonal scalar.
    """
    _parameters(field, t, limits)
    if 27 > limits.max_constraint_entries:
        raise StateAssemblyLimit("rank-profile case cap reached; controls incomplete")
    size = 4*(1+2*t)
    _tensor_cap((size,)*3, limits)
    if size*size > limits.max_map_entries:
        raise StateAssemblyLimit("profile slice matrix cap reached; controls incomplete")
    tensor = matrix_multiplication_tensor(field, 2, 2, 2).tensor_product(mixed_dot_length_tensor(field, t, limits=limits))
    scalar = field.p if field.degree > 1 else 1
    profiles = []
    for a, b, c in product(range(3), repeat=3):
        matrix = lambda rows, columns, rank: tuple(tuple(scalar if i == j < rank else 0
                                                        for j in range(columns)) for i in range(rows))
        A, B, C = matrix(2, 2, a), matrix(2*t, 2, b), matrix(2, 2*t, c)
        vector = mixed_dot_length_core_input(field, t, A, B, C, limits=limits)
        actual = exact_length_slice_rank(tensor, vector, limits=limits)
        expected = 2*t*a+2*b+2*c
        if actual != expected:
            raise ValueError("actual rank-profile slice fails its exact coefficient formula")
        profiles.append(((a, b, c), actual))
    return tuple(profiles)


@dataclass(frozen=True)
class MixedDotLengthControl:
    length: int
    auxiliary: FiniteTensor
    catalyst: FiniteTensor
    source: FiniteTensor
    target: FiniteTensor
    catalyst_input: tuple[int, ...]
    catalyst_slice_rank: int
    source_vectors: tuple[tuple[int, ...], ...]
    source_slice_ranks: tuple[int, ...]
    target_vector: tuple[int, ...]
    target_slice_rank: int

    def require_valid(self, limits=DEFAULT_LIMITS):
        field, t = self.auxiliary.field, self.length
        _parameters(field, t, limits)
        size, d = 1+2*t, self.catalyst.shape
        _tensor_cap(tuple(x+5*size for x in d), limits)
        _tensor_cap(tuple(x+1+4*size for x in d), limits)
        if limits.max_blocks < 6:
            raise StateAssemblyLimit("control direct-sum block cap reached; controls incomplete")
        _auxiliary(field, t, self.auxiliary, limits)
        expected_source = FiniteTensor.direct_sum((self.catalyst,)+(self.auxiliary,)*5)
        expected_target = FiniteTensor.direct_sum((self.catalyst, FiniteTensor.unit(field),
            matrix_multiplication_tensor(field, 2, 2, 2).tensor_product(self.auxiliary)))
        if self.source != expected_source or self.target != expected_target:
            raise ValueError("control source/target coefficient tensors were corrupted")
        rho = exact_length_slice_rank(self.catalyst, self.catalyst_input, limits=limits)
        actual = tuple(exact_length_slice_rank(self.source, vector, limits=limits) for vector in self.source_vectors)
        target_rank = exact_length_slice_rank(self.target, self.target_vector, limits=limits)
        if (self.catalyst_slice_rank != rho or self.source_slice_ranks != actual
                or actual != (rho+5*t+10, rho+3*t+10, rho+10)
                or target_rank != self.target_slice_rank or target_rank != rho+4*t+9):
            raise ValueError("control ranks disagree with actual coefficient slices")


def mixed_dot_length_coefficient_control(field, t, *, catalyst=None, catalyst_input=None,
                                         auxiliary=None, limits=DEFAULT_LIMITS):
    """Full/two-killed/all-killed source slices and a full-rank target slice.

    rho denotes the supplied catalyst input's actual slice rank. No assertion
    that this is the catalyst's generic slice rank is made.
    """
    _parameters(field, t, limits)
    catalyst = FiniteTensor.zero(field, (0, 0, 0)) if catalyst is None else catalyst
    if not isinstance(catalyst, FiniteTensor) or catalyst.field != field:
        raise ValueError("catalyst must be an exact tensor over the supplied field")
    size = 1+2*t
    source_shape = tuple(x+5*size for x in catalyst.shape)
    target_shape = tuple(x+1+4*size for x in catalyst.shape)
    for shape in (source_shape, target_shape, (4*size,)*3):
        _tensor_cap(shape, limits)
        if shape[1]*shape[2] > limits.max_map_entries:
            raise StateAssemblyLimit("control slice matrix cap reached; controls incomplete")
    if limits.max_blocks < 6:
        raise StateAssemblyLimit("control direct-sum block cap reached; controls incomplete")
    S = _auxiliary(field, t, auxiliary, limits)
    d_vector = (1,)*catalyst.shape[0] if catalyst_input is None else tuple(catalyst_input)
    rho = exact_length_slice_rank(catalyst, d_vector, limits=limits)
    source = FiniteTensor.direct_sum((catalyst,)+(S,)*5)
    target = FiniteTensor.direct_sum((catalyst, FiniteTensor.unit(field),
                                    matrix_multiplication_tensor(field, 2, 2, 2).tensor_product(S)))
    vectors = tuple(d_vector+tuple(value for copy in range(5)
                                   for value in (int(copy >= killed),)+(1,)*(2*t)) for killed in (0, 2, 5))
    A = ((1, 0), (0, 1))
    B = tuple(tuple(int(i == j) for j in range(2)) for i in range(2*t))
    C = tuple(tuple(int(i == j) for j in range(2*t)) for i in range(2))
    target_vector = d_vector+(1,)+mixed_dot_length_core_input(field, t, A, B, C, limits=limits)
    ranks = tuple(exact_length_slice_rank(source, vector, limits=limits) for vector in vectors)
    result = MixedDotLengthControl(t, S, catalyst, source, target, d_vector, rho, vectors, ranks,
                                  target_vector, exact_length_slice_rank(target, target_vector, limits=limits))
    result.require_valid(limits)
    return result
