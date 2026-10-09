"""Bounded coordinate diagonal searches that return actual lower maps.

An induced matching in a tensor's nonzero coefficient support gives a
diagonal restriction after rescaling its output coordinates. Exhaustion
concerns this coordinate family only; general linear maps may do better.
The search is useful for cutting finite-state counterexamples without
mistaking a tensor-rank lower bound for an additive-state lower bound.
"""

from dataclasses import dataclass
from math import isfinite
from time import monotonic

from .maps import LinearMap
from .subrank_constraints import export_subrank_order
from .tensors import FiniteTensor
from .witnessed_constraints import DEFAULT_LIMITS, StateAssemblyLimit, StateAssemblyLimits


@dataclass(frozen=True)
class CoordinateSearchBudget:
    max_nodes: int = 100_000
    max_support: int = 2048
    max_diagonal: int = 64
    seconds: float = 1.0

    def __post_init__(self):
        for name in ('max_nodes', 'max_support', 'max_diagonal'):
            if type(getattr(self, name)) is not int or getattr(self, name) < 0:
                raise ValueError('coordinate search counts must be nonnegative integers')
        if (type(self.seconds) not in (int, float) or not isfinite(self.seconds)
                or self.seconds <= 0):
            raise ValueError('coordinate search time must be positive and finite')


@dataclass(frozen=True)
class CoordinateSubrankResult:
    status: str
    reason: str
    checked_nodes: int
    selections: tuple | None = None
    export: object | None = None


def find_coordinate_subrank(source, source_key, goal, *, unit_key='unit',
                            budget=CoordinateSearchBudget(), limits=DEFAULT_LIMITS):
    """Find goal independent units by coordinate selection and normalization.

    Every accepted diagonal is checked at all coefficients. The node/time
    budgets cover proposal search; exact witness verification follows it.
    finite_exhausted excludes only coordinate induced matchings of this size.
    A cap never establishes a lower or upper bound on general subrank.
    """
    if not isinstance(source, FiniteTensor):
        raise ValueError('supply an exact finite source tensor')
    if type(goal) is not int or goal < 1:
        raise ValueError('diagonal goal must be a positive integer')
    if not isinstance(budget, CoordinateSearchBudget):
        raise ValueError('supply explicit coordinate search budgets')
    if not isinstance(limits, StateAssemblyLimits):
        raise ValueError('supply explicit state assembly limits')
    if any(not isinstance(key, str) or not key for key in (source_key, unit_key)):
        raise ValueError('tensor IDs must be nonempty strings')
    if goal > min(source.shape):
        return CoordinateSubrankResult('finite_exhausted',
            'coordinate diagonal exceeds a source axis dimension', 0)
    if goal > budget.max_diagonal:
        return CoordinateSubrankResult('resource_cap', 'diagonal size cap reached', 0)
    try:
        limits.tensor(source.shape)
        limits.tensor((goal,)*3)
        limits.maps(source.shape, (goal,)*3)
        if goal > limits.max_blocks or sum(source.shape)+3*goal > limits.max_constraint_entries:
            raise StateAssemblyLimit('coordinate diagonal bookkeeping cap reached')
    except StateAssemblyLimit as error:
        return CoordinateSubrankResult('resource_cap', str(error), 0)
    support = []
    y, z = source.shape[1:]
    for offset, scalar in enumerate(source.coefficients):
        if scalar:
            if len(support) >= budget.max_support:
                return CoordinateSubrankResult('resource_cap', 'support bitset cap reached', 0)
            ij, k = divmod(offset, z)
            i, j = divmod(ij, y)
            support.append((i, j, k, scalar))
    if len(support) < goal:
        return CoordinateSubrankResult('finite_exhausted',
            'not enough nonzero coordinate terms for this diagonal', 0)
    covers = ({}, {}, {})
    for index, term in enumerate(support):
        for axis, coordinate in enumerate(term[:3]):
            covers[axis][coordinate] = covers[axis].get(coordinate, 0) | (1 << index)
    if any(len(cover) < goal for cover in covers):
        return CoordinateSubrankResult('finite_exhausted',
            'not enough supported coordinates on one axis', 0)
    deadline = monotonic()+budget.seconds
    checked, capped = 0, False

    def search(remaining, selected, chosen, axis_covers):
        nonlocal checked, capped
        if len(selected) == goal:
            return selected
        available = remaining & ~(axis_covers[0] | axis_covers[1] | axis_covers[2])
        if available.bit_count() < goal-len(selected):
            return None
        while available:
            if checked >= budget.max_nodes or monotonic() >= deadline:
                capped = True
                return None
            bit = available & -available
            index = bit.bit_length()-1
            available ^= bit
            checked += 1
            term = support[index]
            updated = tuple(axis_covers[axis] | covers[axis][term[axis]] for axis in range(3))
            wanted = chosen | bit
            # All cross coordinates are checked through support bitsets.
            # A selected box may contain exactly its diagonal nonzero entries.
            if updated[0] & updated[1] & updated[2] != wanted:
                continue
            answer = search(available, selected+(index,), wanted, updated)
            if answer is not None or capped:
                return answer
        return None

    answer = search((1 << len(support))-1, (), 0, (0, 0, 0))
    if answer is None:
        return CoordinateSubrankResult('resource_cap' if capped else 'finite_exhausted',
            'coordinate proposal budget reached' if capped else
            'no induced matching in the complete coordinate family', checked)
    terms = tuple(support[index] for index in answer)
    selections = tuple(tuple(term[axis] for term in terms) for axis in range(3))
    maps = [LinearMap.selection(source.field, size, indices)
            for size, indices in zip(source.shape, selections)]
    maps[2] = LinearMap(source.field, source.shape[2],
        tuple(tuple(source.field.inv(term[3]) if coordinate == term[2] else 0
                    for coordinate in range(source.shape[2])) for term in terms))
    exported = export_subrank_order(source, tuple(maps), source_key,
                                   unit_key=unit_key, limits=limits)
    return CoordinateSubrankResult('found',
        'selected diagonal passes full coefficient equality', checked, selections, exported)
