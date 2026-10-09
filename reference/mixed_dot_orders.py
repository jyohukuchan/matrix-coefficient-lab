"""Coefficient-checked ordinary orders for three differently oriented dots.

These are actual subrank and matrix restrictions, distinct from the analytic
catalyst screen. They neither supply a rank-seven matrix scheme nor construct
a positive-gain catalytic certificate.
"""

from .constraint_generators import WitnessedOrderExport
from .fields import FiniteField
from .maps import LinearMap
from .mixed_dot_catalyst_screen import mixed_dot_tensor
from .sector_matrix import export_star_matrix_order
from .subrank_constraints import export_selected_subrank_order
from .tensors import FiniteTensor, matrix_multiplication_tensor
from .witnessed_constraints import (DEFAULT_LIMITS, StateAssemblyLimit,
                                    StateAssemblyLimits, WitnessedOrder)


def _names(*keys):
    if any(not isinstance(key, str) or not key for key in keys) or len(set(keys)) != len(keys):
        raise ValueError('tensor IDs must be distinct nonempty strings')


def _preflight(field, source_shape, target_shape, blocks, limits):
    if not isinstance(field, FiniteField):
        raise ValueError('supply an explicitly represented finite field')
    if not isinstance(limits, StateAssemblyLimits):
        raise ValueError('supply explicit state assembly limits')
    if limits.max_blocks < blocks:
        raise StateAssemblyLimit('mixed-dot block cap reached; construction incomplete')
    limits.tensor(source_shape)
    limits.tensor(target_shape)
    limits.maps(source_shape, target_shape)
    if sum(source_shape)+sum(target_shape) > limits.max_constraint_entries:
        raise StateAssemblyLimit('mixed-dot coordinate bookkeeping cap reached')


def export_mixed_dot_subrank_order(field, source_key, *, unit_key='unit', limits=DEFAULT_LIMITS):
    """Export three independent units <= S by actual coordinate selections."""
    _names(source_key, unit_key)
    _preflight(field, (5, 5, 5), (3, 3, 3), 3, limits)
    return export_selected_subrank_order(mixed_dot_tensor(field),
        ((0, 1, 3), (0, 2, 3), (0, 2, 4)), source_key,
        unit_key=unit_key, limits=limits)


def export_mixed_dot_square_subrank_order(field, square_key, *, unit_key='unit', limits=DEFAULT_LIMITS):
    """Export fifteen units <= S^2; no claim of optimal square subrank is made.

    The three equal-orientation products supply one selected unit each. The
    six unequal-orientation products supply two units each on their disjoint
    coordinate blocks. All off-diagonal coefficients are checked as well.
    """
    _names(square_key, unit_key)
    _preflight(field, (25, 25, 25), (15, 15, 15), 15, limits)
    selectors = (
        (0, 1, 2, 3, 4, 5, 10, 6, 8, 14, 15, 20, 16, 22, 18),
        (0, 2, 7, 3, 9, 10, 11, 12, 13, 14, 15, 21, 17, 22, 18),
        (0, 2, 8, 4, 9, 10, 16, 12, 14, 19, 20, 21, 22, 23, 24),
    )
    return export_selected_subrank_order(mixed_dot_tensor(field).tensor_power(2),
        selectors, square_key, unit_key=unit_key, limits=limits)


def _star_tensor(field):
    support = {(0, 0, 1), (0, 1, 0), (1, 0, 2), (1, 2, 0)}
    return FiniteTensor.from_function(field, (2, 3, 3),
        lambda i, j, k: int((i, j, k) in support))


def _star_projection_maps(field):
    first = LinearMap(field, 5, ((0, 1, 0, 1, 0), (0, 0, 1, 0, 1)))
    return (first, LinearMap.selection(field, 5, (2, 3, 4)),
            LinearMap.selection(field, 5, (4, 2, 3)))


def export_mixed_dot_star_sum_order(field, source_key, star_key, dot_key, *, limits=DEFAULT_LIMITS):
    """Export star + singleton-first dot_2 <= S without identifying their first legs."""
    _names(source_key, star_key, dot_key)
    _preflight(field, (5, 5, 5), (3, 5, 5), 2, limits)
    source = mixed_dot_tensor(field)
    star, dot = _star_tensor(field), FiniteTensor.from_function(field, (1, 2, 2),
        lambda i, j, k: int(j == k))
    first = LinearMap(field, 5,
        ((0, 1, 0, 1, 0), (0, 0, 1, 0, 1), (1, 0, 0, 0, 0)))
    maps = (first, LinearMap.selection(field, 5, (2, 3, 4, 0, 1)),
            LinearMap.selection(field, 5, (4, 2, 3, 0, 1)))
    row = WitnessedOrder((source_key,), (star_key, dot_key), maps, 'mixed-dots-star-plus-dot')
    exported = WitnessedOrderExport(row, ((source_key, source), (star_key, star), (dot_key, dot)))
    exported.require_valid(limits)
    return exported


def export_mixed_dot_two_copy_matrix_order(field, source_key, matrix_key, *, limits=DEFAULT_LIMITS):
    """Export M2 <= two copies of S by sharing the second matrix factor.

    Each copy computes one row of a 2-by-2 product. Their first matrix rows
    are independent and their second matrix coordinates are duplicated.
    """
    _names(source_key, matrix_key)
    _preflight(field, (10, 10, 10), (4, 4, 4), 2, limits)
    source = mixed_dot_tensor(field)
    target = matrix_multiplication_tensor(field, 2, 2, 2)

    def local_map(supports):
        return LinearMap(field, 10,
            tuple(tuple(int(index in support) for index in range(10)) for support in supports))

    maps = (local_map(((0, 3), (1, 4), (5, 8), (6, 9))),
            local_map(((3, 8), (0, 5), (4, 9), (2, 7))),
            local_map(((4,), (0, 2), (9,), (5, 7))))
    row = WitnessedOrder((source_key, source_key), (matrix_key,), maps,
        'mixed-dots-two-copy-matrix')
    exported = WitnessedOrderExport(row, ((source_key, source), (matrix_key, target)))
    exported.require_valid(limits)
    return exported


def _square_map(mapping):
    field, size = mapping.field, mapping.input_size
    return LinearMap(field, size*size,
        tuple(tuple(field.mul(a, b) for a in first for b in second)
              for first in mapping.rows for second in mapping.rows))


def export_mixed_dot_matrix_order(field, square_key, matrix_key, *, limits=DEFAULT_LIMITS):
    """Export M2 <= S^2 by composing S -> star with the signed star-square maps."""
    _names(square_key, matrix_key)
    _preflight(field, (25, 25, 25), (4, 4, 4), 1, limits)
    # The intermediate square projection allocates more entries than the
    # final maps, so bound it before constructing any tensor or local map.
    limits.maps((25, 25, 25), (4, 9, 9))
    if 75+22 > limits.max_constraint_entries:
        raise StateAssemblyLimit('mixed-dot square coordinate bookkeeping cap reached')
    source = mixed_dot_tensor(field).tensor_power(2)
    projections = tuple(_square_map(mapping) for mapping in _star_projection_maps(field))
    internal_key = '__mixed_dot_star_square'
    while internal_key in (square_key, matrix_key):
        internal_key += '_'
    star = export_star_matrix_order(field, internal_key, matrix_key, limits=limits)
    if source.restrict(*projections) != star.registry[internal_key]:
        raise ValueError('mixed-dot square does not project to the supplied star square')
    maps = tuple(after.compose(before) for after, before in zip(star.row.maps, projections))
    target = matrix_multiplication_tensor(field, 2, 2, 2)
    row = WitnessedOrder((square_key,), (matrix_key,), maps, 'mixed-dots-square-matrix')
    exported = WitnessedOrderExport(row, ((square_key, source), (matrix_key, target)), exponent=2)
    exported.require_valid(limits)
    return exported
