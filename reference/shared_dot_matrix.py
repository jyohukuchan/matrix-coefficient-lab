"""Exact matrix connections for shared dots; these are ordinary restrictions."""

from .constraint_generators import WitnessedOrderExport
from .fields import FiniteField
from .maps import LinearMap
from .sector_matrix import export_star_matrix_order
from .shared_dot_orders import shared_dot_tensor
from .tensors import matrix_multiplication_tensor
from .witnessed_constraints import DEFAULT_LIMITS, StateAssemblyLimit, StateAssemblyLimits, WitnessedOrder


def _preflight(field, q, source_key, target_key, source_shape, target_shape, limits, minimum):
    if not isinstance(field, FiniteField) or type(q) is not int or q < minimum:
        raise ValueError(f"supply a represented field and shared-dot length q>={minimum}")
    if any(not isinstance(key, str) or not key for key in (source_key, target_key)) or source_key == target_key:
        raise ValueError("tensor IDs must be distinct nonempty strings")
    if not isinstance(limits, StateAssemblyLimits):
        raise ValueError("supply explicit state assembly limits")
    limits.tensor(source_shape)
    limits.tensor(target_shape)
    limits.maps(source_shape, target_shape)
    if sum(source_shape)+sum(target_shape) > limits.max_constraint_entries:
        raise StateAssemblyLimit("shared-dot matrix bookkeeping cap reached")


def export_shared_dot_star_order(field, q, source_key, star_key, *, limits=DEFAULT_LIMITS):
    """The first coordinates 1,2 and other coordinates 0,1,2 give the star."""
    size = q+1
    _preflight(field, q, source_key, star_key, (size,)*3, (2, 3, 3), limits, 2)
    source = shared_dot_tensor(field, q, limits=limits)
    maps = tuple(LinearMap.selection(field, size, indices) for indices in ((1, 2), (0, 1, 2), (0, 1, 2)))
    # Reuse the independently checked star support through its square API.
    from .tensors import FiniteTensor
    support = {(0, 0, 1), (0, 1, 0), (1, 0, 2), (1, 2, 0)}
    star = FiniteTensor.from_function(field, (2, 3, 3), lambda i, j, k: int((i, j, k) in support))
    row = WitnessedOrder((source_key,), (star_key,), maps, "shared-dot-star")
    exported = WitnessedOrderExport(row, ((source_key, source), (star_key, star)))
    exported.require_valid(limits)
    return exported


def export_shared_dot_matrix_order(field, q, square_key, matrix_key, *, limits=DEFAULT_LIMITS):
    """C_q squared restricts to M2 via the signed star-square construction."""
    size = q+1
    _preflight(field, q, square_key, matrix_key, (size*size,)*3, (4,)*3, limits, 2)
    limits.maps((size*size,)*3, (4, 9, 9))
    if 3*size*size+22 > limits.max_constraint_entries:
        raise StateAssemblyLimit("shared-dot intermediate square bookkeeping cap reached")
    internal = "__shared_dot_star_square"
    while internal in (square_key, matrix_key):
        internal += "_"
    star = export_star_matrix_order(field, internal, matrix_key, limits=limits)
    selections = ((1, 2), (0, 1, 2), (0, 1, 2))
    projections = tuple(LinearMap.selection(field, size*size,
        tuple(i*size+j for i in indices for j in indices)) for indices in selections)
    source = shared_dot_tensor(field, q, limits=limits).tensor_power(2)
    if source.restrict(*projections) != star.registry[internal]:
        raise ValueError("shared-dot square does not project to the checked star square")
    maps = tuple(after.compose(before) for after, before in zip(star.row.maps, projections))
    target = matrix_multiplication_tensor(field, 2, 2, 2)
    row = WitnessedOrder((square_key,), (matrix_key,), maps, "shared-dot-square-matrix")
    exported = WitnessedOrderExport(row, ((square_key, source), (matrix_key, target)), exponent=2)
    exported.require_valid(limits)
    return exported


def export_shared_dot_two_copy_matrix_order(field, q, source_key, matrix_key, *, limits=DEFAULT_LIMITS):
    """Two C_q copies compute independent rows against one shared B matrix.

    In each copy X=(a0,a1,0,...), Y=(B11,B10,B00,B01-B10,0,...).
    Add output coordinates Z0+Z2 and Z1+Z3 to get the two row entries.
    Copy occurrences are integers; the minus sign is a field scalar.
    """
    size = q+1
    _preflight(field, q, source_key, matrix_key, (2*size,)*3, (4,)*3, limits, 3)
    if limits.max_blocks < 2:
        raise StateAssemblyLimit("shared-dot two-copy block cap reached")
    first = LinearMap.selection(field, 2*size, (0, 1, size, size+1))
    y_rows = ({2: 1}, {3: 1}, {1: 1, 3: -1}, {0: 1})
    second = LinearMap(field, 2*size, tuple(tuple(field.embed(row.get(i % size, 0))
        for i in range(2*size)) for row in y_rows))
    supports = ((0, 2), (1, 3), (size, size+2), (size+1, size+3))
    third = LinearMap(field, 2*size, tuple(tuple(int(i in support) for i in range(2*size))
        for support in supports))
    source = shared_dot_tensor(field, q, limits=limits)
    target = matrix_multiplication_tensor(field, 2, 2, 2)
    row = WitnessedOrder((source_key, source_key), (matrix_key,), (first, second, third),
                         "shared-dot-two-copy-matrix")
    exported = WitnessedOrderExport(row, ((source_key, source), (matrix_key, target)))
    exported.require_valid(limits)
    return exported
