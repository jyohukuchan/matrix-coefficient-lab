"""Actual core projections for the analytic all-catalyst rank obstruction.

The substitution and replicated rectangular-rank arguments are in
research/mixed-dot-rank-catalysts.md. This module checks finite coefficients
and first concision; it does not compute tensor rank or formalize that proof.
"""

from dataclasses import dataclass

from .catalyst_search_experiment import field_matrix_rank
from .fields import FiniteField
from .maps import LinearMap
from .mixed_dot_lengths import mixed_dot_length_tensor
from .tensors import FiniteTensor, matrix_multiplication_tensor
from .witnessed_constraints import DEFAULT_LIMITS, StateAssemblyLimit, StateAssemblyLimits


def _preflight(field, length, repetitions, limits):
    if not isinstance(field, FiniteField) or any(type(v) is not int or v < 1
                                               for v in (length, repetitions)):
        raise ValueError("supply a represented field and positive length/copy count")
    if not isinstance(limits, StateAssemblyLimits):
        raise ValueError("supply explicit assembly limits")
    size = 4*(1+2*length)*repetitions
    target_shape = (4*repetitions, 4*length*repetitions, 4*length*repetitions)
    limits.tensor((size,)*3)
    limits.tensor(target_shape)
    limits.maps((size,)*3, target_shape)
    if size**3 > limits.max_map_entries or 3*size+sum(target_shape) > limits.max_constraint_entries:
        raise StateAssemblyLimit("core first-flattening or coordinate cap reached")
    if limits.max_blocks < max(3, repetitions):
        raise StateAssemblyLimit("mixed-dot rank-control block cap reached")


def _data(field, length, repetitions, limits):
    _preflight(field, length, repetitions, limits)
    size = 1+2*length
    auxiliary = mixed_dot_length_tensor(field, length, limits=limits)
    core = matrix_multiplication_tensor(field, 2, 2, 2).tensor_product(auxiliary)
    core = FiniteTensor.direct_sum((core,)*repetitions)
    rectangle = matrix_multiplication_tensor(field, 2, 2, 2*length)
    rectangle = FiniteTensor.direct_sum((rectangle,)*repetitions)
    selectors = [[], [], []]
    for copy in range(repetitions):
        offset = copy*4*size
        selectors[0].extend(offset+a*size for a in range(4))
        selectors[1].extend(offset+(2*j+k)*size+r
                            for j in range(2) for k in range(2) for r in range(length))
        selectors[2].extend(offset+(2*i+k)*size+r
                            for i in range(2) for k in range(2) for r in range(length))
    maps = tuple(LinearMap.selection(field, core.shape[axis], selector)
                 for axis, selector in enumerate(selectors))
    rows = tuple(tuple(core.coefficient(i, j, k)
                       for j in range(core.shape[1]) for k in range(core.shape[2]))
                 for i in range(core.shape[0]))
    rank = field_matrix_rank(field, rows, core.shape[1]*core.shape[2])
    return core, rectangle, maps, rank


@dataclass(frozen=True)
class MixedDotRankControl:
    core: FiniteTensor
    rectangle: FiniteTensor
    maps: tuple[LinearMap, ...]
    length: int
    repetitions: int
    first_flattening_rank: int

    @property
    def analytic_core_rank_lower(self):
        return 15*self.length*self.repetitions

    def require_valid(self, limits=DEFAULT_LIMITS):
        core, rectangle, maps, rank = _data(self.core.field, self.length, self.repetitions, limits)
        if (self.core != core or self.rectangle != rectangle or self.maps != maps
                or self.first_flattening_rank != rank or rank != core.shape[0]
                or core.restrict(*self.maps) != rectangle):
            raise ValueError("mixed-dot rank-control coefficients, maps or first concision disagree")


def mixed_dot_rank_control(field, length=2, *, repetitions=1, limits=DEFAULT_LIMITS):
    """Coefficient-checked projection to n independent M(2,2,2t) blocks.

    The analytic lower 15tn uses substitution and a replicated kernel-map
    theorem; finite first-flattening ranks alone do not establish it.
    """
    core, rectangle, maps, rank = _data(field, length, repetitions, limits)
    control = MixedDotRankControl(core, rectangle, maps, length, repetitions, rank)
    control.require_valid(limits)
    return control
