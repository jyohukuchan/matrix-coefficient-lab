"""Budgeted, seeded exact decomposition search by flips and plus moves.

Every move is a tensor identity. Reductions and final outputs are independently
coefficient checked. Search caps are never tensor-rank lower bounds. The seed
is supplied explicitly; this is not an unseeded coefficient solver.
"""

from collections import defaultdict, deque
from dataclasses import dataclass
from itertools import combinations
from random import Random
from time import monotonic

from .tensor_schemes import TensorScheme


@dataclass(frozen=True)
class FlipBudget:
    max_steps: int = 100_000
    timeout_seconds: float = 60.0
    probes_per_step: int = 12
    max_extra_terms: int = 1
    recent_states: int = 256
    max_tensor_entries: int = 100_000
    max_seed_terms: int = 128

    def __post_init__(self):
        for name in ('max_steps', 'max_extra_terms', 'recent_states', 'max_tensor_entries', 'max_seed_terms'):
            value = getattr(self, name)
            if type(value) is not int or value < 0:
                raise ValueError(f'{name} must be a nonnegative integer')
        if type(self.probes_per_step) is not int or self.probes_per_step < 1:
            raise ValueError('probe count must be a positive integer')
        if not isinstance(self.timeout_seconds, (int, float)) or isinstance(self.timeout_seconds, bool) or not 0 < self.timeout_seconds < float('inf'):
            raise ValueError('timeout must be a finite positive number')


class TermAlgebra:
    """Small finite-field tables; term normalization preserves the tensor."""

    def __init__(self, field):
        if field.order > 16:
            raise ValueError('flip experiment supports represented fields of order at most 16')
        self.field = field
        q = field.order
        self.add = tuple(tuple(field.add(a, b) for b in range(q)) for a in range(q))
        self.mul = tuple(tuple(field.mul(a, b) for b in range(q)) for a in range(q))
        self.neg = tuple(field.neg(a) for a in range(q))
        self.inv = (0,)+tuple(field.inv(a) for a in range(1, q))

    def scale(self, vector, scalar):
        return tuple(self.mul[scalar][v] for v in vector)

    def vector_add(self, a, b):
        return tuple(self.add[x][y] for x, y in zip(a, b))

    def normalized(self, vector):
        scalar = next((v for v in vector if v), 0)
        return (self.scale(vector, self.inv[scalar]), scalar) if scalar else (vector, 0)

    def canonical(self, terms):
        result = []
        for a, b, c in terms:
            a, sa = self.normalized(a)
            b, sb = self.normalized(b)
            if sa and sb and any(c):
                result.append((a, b, self.scale(c, self.mul[sa][sb])))
        return tuple(sorted(result))

    def reduce(self, terms):
        """Merge proportional terms agreeing on any two legs, until stable."""
        terms = self.canonical(terms)
        while True:
            count = len(terms)
            for first, second, remaining in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
                merged = {}
                for term in terms:
                    a, sa = self.normalized(term[first])
                    b, sb = self.normalized(term[second])
                    value = self.scale(term[remaining], self.mul[sa][sb])
                    key = a, b
                    merged[key] = self.vector_add(merged[key], value) if key in merged else value
                rebuilt = []
                for (a, b), value in merged.items():
                    if any(value):
                        term = [None]*3
                        term[first], term[second], term[remaining] = a, b, value
                        rebuilt.append(tuple(term))
                terms = self.canonical(rebuilt)
            if len(terms) == count:
                return terms

    def flip(self, terms, axis, first, second, scalar=1):
        """u*v + u'*v' = (u+t*u')*v + u'*(v'-t*v), common leg."""
        if (any(type(v) is not int for v in (axis, first, second, scalar))
                or first == second or not 0 <= axis < 3
                or not 0 <= first < len(terms) or not 0 <= second < len(terms)
                or not 1 <= scalar < self.field.order):
            raise ValueError('invalid flip indices or scalar')
        a, b = list(terms[first]), list(terms[second])
        common_a, sa = self.normalized(a[axis])
        common_b, sb = self.normalized(b[axis])
        if not sa or not sb or common_a != common_b:
            raise ValueError('flip requires proportional nonzero shared factors')
        u, v = (axis+1) % 3, (axis+2) % 3
        a[axis] = b[axis] = common_a
        a[u], b[u] = self.scale(a[u], sa), self.scale(b[u], sb)
        a[u] = self.vector_add(a[u], self.scale(b[u], scalar))
        b[v] = self.vector_add(b[v], self.scale(a[v], self.neg[scalar]))
        result = list(terms)
        result[first], result[second] = tuple(a), tuple(b)
        return self.reduce(result)

    def plus(self, terms, first, second, axis=0):
        """Two arbitrary terms -> three exact terms, with no shared-leg premise.

        abc + ABC = a(b-B)c + (a+A)BC + aB(c-C).
        """
        if (any(type(v) is not int for v in (axis, first, second))
                or first == second or not 0 <= axis < 3
                or not 0 <= first < len(terms) or not 0 <= second < len(terms)):
            raise ValueError('invalid plus indices')
        order = axis, (axis+1) % 3, (axis+2) % 3
        a, b, c = (terms[first][i] for i in order)
        A, B, C = (terms[second][i] for i in order)
        values = ((a, self.vector_add(b, self.scale(B, self.neg[1])), c),
                  (self.vector_add(a, A), B, C),
                  (a, B, self.vector_add(c, self.scale(C, self.neg[1]))))
        rebuilt = [term for i, term in enumerate(terms) if i not in (first, second)]
        for value in values:
            term = [None]*3
            for index, vector in zip(order, value):
                term[index] = vector
            rebuilt.append(tuple(term))
        return self.reduce(rebuilt)

    def shared_pairs(self, terms):
        pairs = []
        for axis in range(3):
            groups = defaultdict(list)
            for i, term in enumerate(terms):
                groups[self.normalized(term[axis])[0]].append(i)
            for group in groups.values():
                pairs.extend((axis, a, b) for a, b in combinations(group, 2))
        return pairs


def _scheme(target, terms):
    return TensorScheme(target, *(tuple(term[axis] for term in terms) for axis in range(3)))


@dataclass(frozen=True)
class FlipResult:
    status: str
    reason: str
    initial_terms: int
    goal_terms: int
    best_terms: int
    steps: int
    proposals: int
    plus_moves: int
    seed: int
    elapsed_seconds: float
    improvements: tuple[tuple[int, int], ...]
    scheme: TensorScheme | None


def search_flips(scheme, *, goal_terms, budget=FlipBudget(), seed=0):
    if not isinstance(scheme, TensorScheme) or not isinstance(budget, FlipBudget):
        raise ValueError('supply an exact tensor scheme and a flip budget')
    if type(goal_terms) is not int or goal_terms < 0 or type(seed) is not int:
        raise ValueError('goal terms must be nonnegative; seed must be an integer')
    start = monotonic()
    if (len(scheme.target.coefficients) > budget.max_tensor_entries or scheme.terms > budget.max_seed_terms):
        return FlipResult('resource_cap', 'input construction budget exceeded', scheme.terms, goal_terms,
                          scheme.terms, 0, 0, 0, seed, monotonic()-start, (), None)
    scheme.require_exact()
    algebra, rng = TermAlgebra(scheme.field), Random(seed)
    current = algebra.reduce(tuple(zip(scheme.a, scheme.b, scheme.c)))
    best = current
    _scheme(scheme.target, best).require_exact()
    improvements = [(0, len(best))] if len(best) < scheme.terms else []
    recent, seen = deque(), set()
    proposals, plus_moves, steps = 0, 0, 0
    reason = 'step cap reached; search incomplete'
    for steps in range(1, budget.max_steps+1):
        if len(best) <= goal_terms:
            steps -= 1
            reason = 'requested term bound reached with exact coefficients'
            break
        if monotonic()-start >= budget.timeout_seconds:
            steps -= 1
            reason = 'time cap reached; search incomplete'
            break
        if budget.recent_states:
            if len(recent) == budget.recent_states:
                seen.discard(recent.popleft())
            recent.append(current)
            seen.add(current)
        pairs = algebra.shared_pairs(current)
        use_plus = (len(current) < scheme.terms+budget.max_extra_terms and len(current) > 1
                    and (not pairs or rng.randrange(32) == 0))
        if use_plus:
            a, b = rng.sample(range(len(current)), 2)
            current = algebra.plus(current, a, b, rng.randrange(3))
            proposals += 1
            plus_moves += 1
        elif pairs:
            choices = []
            for _ in range(budget.probes_per_step):
                axis, a, b = rng.choice(pairs)
                if rng.randrange(2):
                    a, b = b, a
                proposed = algebra.flip(current, axis, a, b, rng.randrange(1, scheme.field.order))
                proposals += 1
                if proposed not in seen:
                    choices.append(proposed)
                if len(proposed) < len(best):
                    choices = [proposed]
                    break
            if choices:
                minimum = min(map(len, choices))
                current = rng.choice([choice for choice in choices if len(choice) == minimum])
            else:
                # Revisit a state rather than misreport a tabu dead end as UNSAT.
                seen.clear()
                recent.clear()
        else:
            reason = 'no allowed local moves; this says nothing about other decompositions'
            break
        if len(current) < len(best):
            _scheme(scheme.target, current).require_exact()
            best = current
            improvements.append((steps, len(best)))
    else:
        steps = budget.max_steps
    result = _scheme(scheme.target, best)
    result.require_exact()
    reached = len(best) <= goal_terms
    status = 'found' if reached else 'partial_improvement' if len(best) < scheme.terms else 'resource_cap'
    if reached:
        reason = 'requested term bound reached with exact coefficients'
    return FlipResult(status, reason, scheme.terms, goal_terms, len(best), steps, proposals,
                      plus_moves, seed, monotonic()-start, tuple(improvements), result)
