"""Optional LP proposals must yield independently exact finite certificates."""

from contextlib import contextmanager
from dataclasses import replace
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from .fields import FiniteField
from .finite_state_solver import RationalStateCertificate, _fractions, solve_finite_states
from .maps import LinearMap
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, matrix_multiplication_tensor
from .witnessed_constraints import (DetectorConstraint, StateAssemblyLimits,
                                    WitnessedStateSystem, assemble_catalyst,
                                    rank_upper_constraint)


F2 = FiniteField(2)
HAS_SCIPY = importlib.util.find_spec("scipy") is not None and importlib.util.find_spec("numpy") is not None


def control_system(k=5):
    """Known exact rank7 INPUT certificate, not a new rank search."""
    artifact = Path(__file__).with_name("research") / "unseeded-f2-rank7.json"
    data = json.loads(artifact.read_text())
    matrix = matrix_multiplication_tensor(F2, 2, 2, 2)
    scheme = TensorScheme(matrix, data["a"], data["b"], data["c"])
    scheme.require_exact()
    return WitnessedStateSystem(2, k, (("unit", FiniteTensor.unit(F2)), ("matrix", matrix)),
                                (rank_upper_constraint("unit", "matrix", scheme),
                                 DetectorConstraint("unit", "matrix")))


class _Array:
    """Only structural NumPy operations used before mocked LP proposals."""
    def __init__(self, data):
        self.data = tuple(data)

    def reshape(self, rows, columns):
        flat = tuple(value for row in self.data for value in row)
        return _Array(tuple(flat[i*columns:(i+1)*columns] for i in range(rows)))

    @property
    def T(self):
        return _Array(tuple(zip(*self.data)))

    def __neg__(self):
        if self.data and isinstance(self.data[0], (tuple, list)):
            return _Array(tuple(tuple(-x for x in row) for row in self.data))
        return _Array(tuple(-x for x in self.data))


@contextmanager
def mocked_proposals(*proposals):
    """Exercise rational/status handling even when optional SciPy is absent."""
    numpy, scipy, optimize = ModuleType("numpy"), ModuleType("scipy"), ModuleType("scipy.optimize")
    numpy.asarray = lambda data, dtype=None: _Array(data)
    numpy.zeros = lambda size: _Array((0.0,)*size)
    optimize.linprog = Mock(side_effect=proposals)
    scipy.optimize, scipy.__path__ = optimize, []
    with patch.dict(sys.modules, {"numpy": numpy, "scipy": scipy, "scipy.optimize": optimize}):
        yield optimize.linprog


def proposal(success, values=(), status=0):
    return SimpleNamespace(success=success, x=values, status=status)


class ExactStateCertificateTests(unittest.TestCase):
    def test_exact_normalized_finite_state_acceptance_and_rejection(self):
        system = control_system(5)
        good = RationalStateCertificate(system, (Fraction(1), Fraction(6)))
        self.assertTrue(good.verify_inequalities())
        good.require_valid()
        self.assertEqual(good.assignment, {"unit": Fraction(1), "matrix": Fraction(6)})
        for values in ((Fraction(0), Fraction(6)), (Fraction(1), Fraction(4)),
                       (Fraction(1), Fraction(8)), (Fraction(1), Fraction(-1))):
            bad = RationalStateCertificate(system, values)
            self.assertFalse(bad.verify_inequalities())
            with self.assertRaises(ValueError):
                bad.require_valid()
        for values in ((1, 6), (Fraction(1),), (Fraction(1), 6.0)):
            with self.assertRaises(ValueError):
                RationalStateCertificate(system, values)

    def test_exact_assignment_cannot_override_corrupted_coefficient_witness(self):
        system = control_system(5)
        upper = system.constraints[0]
        old = upper.maps[0]
        rows = list(old.rows)
        rows[0] = (0,)*old.input_size
        bad_row = replace(upper, maps=(LinearMap(F2, old.input_size, rows),)+upper.maps[1:])
        bad_system = replace(system, constraints=(bad_row, system.constraints[1]))
        state = RationalStateCertificate(bad_system, (Fraction(1), Fraction(6)))
        self.assertTrue(state.verify_inequalities())  # scalar balance alone is insufficient
        with self.assertRaisesRegex(ValueError, "coefficient witness"):
            state.require_valid()
        with self.assertRaises(ValueError):
            solve_finite_states(bad_system)

    def test_nonfinite_and_negative_rational_approximations(self):
        self.assertEqual(_fractions((1.0, -0.5), 10), (Fraction(1), Fraction(-1, 2)))
        self.assertEqual(_fractions((.33333333,), 100), (Fraction(1, 3),))
        for value in (float("nan"), float("inf"), -float("inf")):
            self.assertIsNone(_fractions((value,), 100))


class NumericProposalTests(unittest.TestCase):
    def test_rounded_exact_dual_is_accepted_and_assembled(self):
        system = control_system(8)
        with mocked_proposals(proposal(True, (0.9999999, 1.0000001))) as solver:
            result = solve_finite_states(system, denominator_cap=100)
        self.assertEqual(solver.call_count, 1)
        self.assertEqual(result.status, "found")
        self.assertTrue(result.dual.verify_balance())
        result.dual.require_valid()
        catalyst = assemble_catalyst(result.dual)
        self.assertEqual((catalyst.d, catalyst.k, catalyst.m), (2, 8, 1))
        self.assertTrue(catalyst.verify())

    def test_incorrect_numeric_success_cannot_claim_dual_or_state(self):
        system = control_system(8)
        for dual_values in ((0.0, 0.0), (-1.0, 1.0), (1.0, -0.0001), (float("nan"), 1.0)):
            with mocked_proposals(proposal(True, dual_values), proposal(True, (1.0, 7.0))) as solver:
                result = solve_finite_states(system)
            self.assertEqual(solver.call_count, 2)
            self.assertEqual(result.status, "resource_cap")
            self.assertIsNone(result.dual)
            self.assertIsNone(result.state)
            self.assertIn("unresolved", result.reason)

    def test_invalid_dual_can_fall_back_to_an_independently_exact_finite_state(self):
        system = control_system(5)
        with mocked_proposals(proposal(True, (-1.0, 1.0)), proposal(True, (1.0, 6.0))):
            result = solve_finite_states(system)
        self.assertEqual(result.status, "finite_feasible")
        result.state.require_valid()
        self.assertIsNone(result.dual)

    def test_status_codes_alone_prove_neither_direction(self):
        system = control_system(5)
        for statuses in ((2, 2), (1, 1), (0, 0)):
            with mocked_proposals(proposal(False, (1.0, 1.0), statuses[0]),
                                  proposal(False, (1.0, 6.0), statuses[1])):
                result = solve_finite_states(system)
            self.assertEqual(result.status, "resource_cap")
            self.assertEqual(result.numeric_statuses, statuses)
            self.assertIsNone(result.state)
        # A contradictory reported status cannot invalidate a genuinely exact
        # certificate: only the rational witness supplies mathematical evidence.
        with mocked_proposals(proposal(True, (1.0, 1.0), 2)):
            result = solve_finite_states(control_system(8))
        self.assertEqual(result.status, "found")
        self.assertTrue(result.dual.verify_balance())

    def test_wrong_normalization_and_nonfinite_state_proposals_remain_unresolved(self):
        for values in ((0.0, 6.0), (1.0, 4.0), (1.0, float("inf"))):
            with mocked_proposals(proposal(False, status=2), proposal(True, values)):
                result = solve_finite_states(control_system(5))
            self.assertEqual(result.status, "resource_cap")
            self.assertIsNone(result.state)

    def test_malformed_proposal_lengths_cannot_claim_a_certificate(self):
        for values in ((), (1.0,), (1.0, 6.0, 0.0)):
            with mocked_proposals(proposal(True, values), proposal(True, values)):
                result = solve_finite_states(control_system(5))
            self.assertEqual(result.status, "resource_cap")
            self.assertIsNone(result.dual)
            self.assertIsNone(result.state)
            self.assertIn("unresolved", result.reason)

    def test_scipy_absent_and_resource_caps_have_honest_statuses(self):
        system = control_system(5)
        with patch.dict(sys.modules, {"scipy": None, "scipy.optimize": None}):
            missing = solve_finite_states(system)
        self.assertEqual(missing.status, "unavailable")
        self.assertIn("absent", missing.reason)
        for limits in (StateAssemblyLimits(max_tensor_entries=10),
                       StateAssemblyLimits(max_map_entries=1),
                       StateAssemblyLimits(max_constraint_entries=1)):
            with mocked_proposals() as solver:
                result = solve_finite_states(system, limits=limits)
            self.assertEqual(result.status, "resource_cap")
            solver.assert_not_called()
        with mocked_proposals() as solver:
            large_integer = solve_finite_states(control_system(2**53))
        self.assertEqual(large_integer.status, "resource_cap")
        self.assertIn("integer coefficient", large_integer.reason)
        solver.assert_not_called()

    def test_empty_finite_family_has_exact_state_without_optional_dependencies(self):
        system = replace(control_system(5), constraints=())
        with patch.dict(sys.modules, {"numpy": None, "scipy": None, "scipy.optimize": None}):
            result = solve_finite_states(system)
        self.assertEqual(result.status, "finite_feasible")
        self.assertEqual(result.numeric_statuses, ())
        self.assertEqual(result.state.values, (Fraction(1), Fraction(0)))
        result.state.require_valid()

    def test_combined_time_cap_stops_before_second_numeric_proposal(self):
        with mocked_proposals(proposal(False, status=1)) as solver:
            with patch("reference.finite_state_solver.monotonic", side_effect=(10.0, 11.0)):
                result = solve_finite_states(control_system(5), timeout_seconds=1)
        self.assertEqual(solver.call_count, 1)
        self.assertEqual(result.status, "resource_cap")
        self.assertIn("combined numeric proposal time cap", result.reason)
        self.assertEqual(result.numeric_statuses, (1,))
        self.assertIsNone(result.dual)
        self.assertIsNone(result.state)

    def test_second_proposal_receives_only_remaining_time_budget(self):
        with mocked_proposals(proposal(False, status=2), proposal(True, (1.0, 6.0))) as solver:
            with patch("reference.finite_state_solver.monotonic", side_effect=(10.0, 10.25)):
                result = solve_finite_states(control_system(5), timeout_seconds=1)
        self.assertEqual(solver.call_count, 2)
        self.assertEqual(solver.call_args.kwargs["options"]["time_limit"], .75)
        self.assertEqual(result.status, "finite_feasible")
        result.state.require_valid()

    def test_invalid_solver_options(self):
        system = control_system()
        for options in ({"denominator_cap": 0}, {"denominator_cap": True}, {"denominator_cap": 1.0},
                        {"timeout_seconds": 0}, {"timeout_seconds": True},
                        {"timeout_seconds": float("nan")}, {"timeout_seconds": float("inf")}):
            with self.assertRaises(ValueError):
                solve_finite_states(system, **options)
        with self.assertRaises(ValueError):
            solve_finite_states("not a system")


@unittest.skipUnless(HAS_SCIPY, "optional scipy/numpy unavailable")
class GenuineNumericControlTests(unittest.TestCase):
    def test_rank_seven_control_k5_is_exact_finite_feasible(self):
        result = solve_finite_states(control_system(5), timeout_seconds=5)
        self.assertEqual(result.status, "finite_feasible")
        result.state.require_valid()
        value = result.state.assignment["matrix"]
        self.assertGreaterEqual(value, Fraction(5))
        self.assertLessEqual(value, Fraction(7))

    def test_rank_seven_control_k8_has_exact_dual_and_actual_catalyst(self):
        result = solve_finite_states(control_system(8), timeout_seconds=5)
        self.assertEqual(result.status, "found")
        result.dual.require_valid()
        catalyst = assemble_catalyst(result.dual)
        self.assertTrue(catalyst.verify())
        self.assertEqual(catalyst.m, 1)


if __name__ == "__main__":
    unittest.main()
