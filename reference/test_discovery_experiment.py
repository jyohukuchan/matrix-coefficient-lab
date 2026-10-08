"""Replay checks exact maps, scalar normalization and integer balances."""

from copy import deepcopy
from fractions import Fraction
import unittest

from .constraint_generators import WitnessedOrderExport
from .discovery_experiment import replay_finite, run_experiment, serialize_finite
from .fields import FiniteField
from .finite_state_solver import RationalStateCertificate, StateSolverResult
from .proof_pipeline import ConstraintInventory
from .schemes import naive_scheme
from .witnessed_constraints import IntegerStateCertificate, rank_upper_constraint


class FiniteReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.book = ConstraintInventory(FiniteField(2))
        matrix = naive_scheme(FiniteField(2), 2, 2, 2).as_tensor_scheme()
        key = cls.book.register('M2', matrix.target)
        cls.book.add_export(WitnessedOrderExport(rank_upper_constraint('unit', key, matrix),
            (('unit', cls.book.inventory['unit']), (key, matrix.target))), 'supplied_upper')
        cls.book.add_detector('unit', 2)

    def test_state_replay_checks_all_inequalities_and_maps(self):
        system = self.book.system(2, 5)
        state = RationalStateCertificate(system, (Fraction(1), Fraction(5)))
        data = serialize_finite(self.book, system, StateSolverResult('finite_feasible', '', state=state))
        restored, checked = replay_finite(data)
        self.assertEqual(checked, ('finite_feasible',))
        self.assertEqual(restored.sparse_vectors, system.sparse_vectors)
        self.assertEqual(restored.registry, system.registry)
        for wrong in (['0', '5'], ['1', '4'], ['1', '9']):
            with self.assertRaises(ValueError):
                replay_finite({**data, 'state': wrong})
        corrupted = deepcopy(data)
        corrupted['rows'][0]['maps'] = [0, 0, 0]  # tensor leaves are not maps
        with self.assertRaises(ValueError):
            replay_finite(corrupted)

    def test_negative_unit_dual_replay_is_over_integers(self):
        system = self.book.system(2, 9)
        dual = IntegerStateCertificate(system, (1, 1), 1)
        data = serialize_finite(self.book, system, StateSolverResult('found', '', dual=dual))
        restored, checked = replay_finite(data)
        self.assertEqual(checked, ('negative_unit_dual',))
        self.assertEqual(restored.k, 9)
        for change in ({'weights': [1, 0], 'gain': 1}, {'weights': [1, 1], 'gain': 3}):
            with self.assertRaises(ValueError):
                replay_finite({**data, 'dual': change})

    def test_parameter_validation_precedes_expensive_construction(self):
        for kwargs in ({'k': True}, {'seeds': (0, 0)}, {'seeds': ()}, {'core_seconds': 0},
                       {'flip_seconds': float('nan')}, {'lp_seconds': True}, {'sat_seconds': -1}):
            with self.assertRaises(ValueError):
                run_experiment(**kwargs)

    def test_selected_and_reordered_rows_keep_their_actual_maps(self):
        system = self.book.system(2, 5, (1, 0))
        state = RationalStateCertificate(system, (Fraction(1), Fraction(5)))
        result = StateSolverResult('finite_feasible', '', state=state)
        restored, checked = replay_finite(serialize_finite(self.book, system, result))
        self.assertEqual(checked, ('finite_feasible',))
        self.assertEqual(restored.constraints, system.constraints)
        other = self.book.system(2, 5)
        with self.assertRaisesRegex(ValueError, 'different finite system'):
            serialize_finite(self.book, other, result)


if __name__ == '__main__':
    unittest.main()
