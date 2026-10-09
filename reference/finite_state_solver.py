"""Optional numeric LP proposals with independent exact rational certificates.

The solver proves either a negative-unit dual or feasibility of this FINITE
witnessed row family. Floating-point statuses alone prove neither. SciPy is an
optional experiment dependency, never imported by the ordinary reference API.
"""

from dataclasses import dataclass
from fractions import Fraction
from math import isfinite
from time import monotonic

from .witnessed_constraints import (DEFAULT_LIMITS, StateAssemblyLimit, WitnessedStateSystem,
                                    clear_rational_dual)


@dataclass(frozen=True)
class RationalStateCertificate:
    system: WitnessedStateSystem
    values: tuple[Fraction, ...]

    def __post_init__(self):
        if not isinstance(self.system, WitnessedStateSystem):
            raise ValueError('state certificate requires an actual finite witnessed system')
        values = tuple(self.values)
        if len(values) != len(self.system.keys) or any(not isinstance(v, Fraction) for v in values):
            raise ValueError('one exact Fraction is required per state coordinate')
        object.__setattr__(self, 'values', values)

    @property
    def assignment(self):
        return dict(zip(self.system.keys, self.values))

    def verify_inequalities(self):
        values = self.assignment
        return (values[self.system.unit_key] == 1 and
                all(sum(count*values[key] for key, count in row) >= 0
                    for row in self.system.sparse_vectors))

    def require_valid(self, limits=DEFAULT_LIMITS):
        self.system.require_valid(limits)
        if not self.verify_inequalities():
            raise ValueError('state fails exact normalization or a prescribed inequality')


@dataclass(frozen=True)
class StateSolverResult:
    status: str
    reason: str
    dual: object | None = None
    state: RationalStateCertificate | None = None
    numeric_statuses: tuple[int, ...] = ()


def _fractions(values, denominator_cap):
    if any(not isfinite(float(v)) for v in values):
        return None
    return tuple(Fraction(str(float(v))).limit_denominator(denominator_cap) for v in values)


def solve_finite_states(system, *, timeout_seconds=10.0, denominator_cap=1_000_000,
                        limits=DEFAULT_LIMITS):
    """Bound the combined numeric proposal phase; exact input checks precede it.

    The solver's native time limit is advisory, not a subprocess wall deadline.
    A success flag is accepted only after exact rational verification.
    """
    if not isinstance(system, WitnessedStateSystem):
        raise ValueError('supply a finite witnessed system')
    if type(denominator_cap) is not int or denominator_cap < 1:
        raise ValueError('denominator cap must be a positive integer')
    if (not isinstance(timeout_seconds, (int, float)) or isinstance(timeout_seconds, bool)
            or not isfinite(timeout_seconds) or timeout_seconds <= 0):
        raise ValueError('timeout must be finite and positive')
    try:
        system.require_valid(limits)
        vectors = system.dense_vectors(limits)
    except StateAssemblyLimit as error:
        return StateSolverResult('resource_cap', str(error))
    if any(abs(v) > 2**52 for row in vectors for v in row):
        return StateSolverResult('resource_cap', 'integer coefficient exceeds the proposal encoding budget')
    n, r = len(system.keys), len(vectors)
    unit = tuple(int(key == system.unit_key) for key in system.keys)
    # No rows: explicit normalized state suffices, without relying on an LP.
    if not r:
        state = RationalStateCertificate(system, tuple(Fraction(v) for v in unit))
        return StateSolverResult('finite_feasible', 'exact normalized state for the empty finite family', state=state)
    try:
        import numpy as np
        from scipy.optimize import linprog
    except ImportError:
        return StateSolverResult('unavailable', 'optional scipy/numpy experiment dependency is absent')
    deadline = monotonic()+timeout_seconds
    matrix = np.asarray(vectors, dtype=float).reshape(r, n)
    options = {'time_limit': timeout_seconds, 'primal_feasibility_tolerance': 1e-9,
               'dual_feasibility_tolerance': 1e-9}
    proposed_dual = linprog(np.zeros(r), A_eq=matrix.T, b_eq=-np.asarray(unit, dtype=float),
                           bounds=(0, None), method='highs', options=options)
    statuses = [int(proposed_dual.status)]
    if proposed_dual.success:
        values = _fractions(proposed_dual.x, denominator_cap)
        if values is not None and all(v >= 0 for v in values):
            try:
                dual = clear_rational_dual(system, values)
            except ValueError:
                pass
            else:
                return StateSolverResult('found', 'rationalized integer dual passes exact balance',
                                         dual=dual, numeric_statuses=tuple(statuses))
    remaining = deadline-monotonic()
    if remaining <= 0:
        return StateSolverResult('resource_cap', 'combined numeric proposal time cap reached',
                                 numeric_statuses=tuple(statuses))
    options['time_limit'] = remaining
    proposed_state = linprog(np.zeros(n), A_ub=-matrix, b_ub=np.zeros(r),
                            A_eq=np.asarray([unit], dtype=float), b_eq=[1.0],
                            bounds=(None, None), method='highs', options=options)
    statuses.append(int(proposed_state.status))
    if proposed_state.success:
        values = _fractions(proposed_state.x, denominator_cap)
        if values is not None:
            try:
                state = RationalStateCertificate(system, values)
            except ValueError:
                state = None
            if state is not None and state.verify_inequalities():
                return StateSolverResult('finite_feasible',
                                         'exact rational normalized state satisfies every selected row',
                                         state=state, numeric_statuses=tuple(statuses))
    return StateSolverResult('resource_cap',
                             'numeric proposals did not produce an exact certificate; finite problem unresolved',
                             numeric_statuses=tuple(statuses))
