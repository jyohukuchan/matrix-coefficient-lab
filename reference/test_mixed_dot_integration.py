"""Saved coefficient families remain replayable after the mixed-dot pass."""

from contextlib import redirect_stdout
from fractions import Fraction
import importlib.util
import io
import json
from pathlib import Path
import runpy
import tempfile
import unittest

from .constraint_generators import WitnessedOrderExport
from .fields import FiniteField
from .finite_state_solver import RationalStateCertificate, StateSolverResult
from .proof_pipeline import ConstraintInventory
from .schemes import naive_scheme
from .witnessed_constraints import rank_upper_constraint


@unittest.skipUnless(importlib.util.find_spec('scipy'), 'optional SciPy LP proposals unavailable')
class MixedDotIntegrationTests(unittest.TestCase):
    def test_shared_dot_rank_schemes_and_matrix_contexts_replay(self):
        runner = runpy.run_path(str(Path(__file__).parent/'research/variable-proof-experiments.py'))
        field = FiniteField(2, (1, 1, 1))
        book = ConstraintInventory(field)
        matrix = naive_scheme(field, 2, 2, 2).as_tensor_scheme()
        book.register('M2', matrix.target)
        book.add_export(WitnessedOrderExport(rank_upper_constraint('unit', 'M2', matrix),
            (('unit', book.inventory['unit']), ('M2', matrix.target))), 'naive_matrix_upper')
        book.add_detector('unit', 2)
        system = book.system(2, 5)
        state = RationalStateCertificate(system, tuple(Fraction(1 if key == 'unit' else 5)
                                                      for key in system.keys))
        initial = runner['literal_export'](book, StateSolverResult('finite_feasible', 'small input', state=state))
        with tempfile.TemporaryDirectory() as folder, redirect_stdout(io.StringIO()):
            root = Path(folder)
            source = root/'input.json'
            source.write_text(json.dumps(initial)+'\n')
            runner['mixed_dots'](source, root/'first', 120, include_shared=True)
            data = json.loads((root/'first/finite-replay.json').read_text())
            report = json.loads((root/'first/report.json').read_text())
            self.assertTrue(report['shared_dot_connections'])
            self.assertFalse(report['supplied_rank_seven_matrix_row'])
            self.assertEqual(runner['replay'](data)['checked'], ['finite_feasible'])
            shared = next(row['positive'][0] for row in data['rows'] if row.get('label') == 'shared-dot-star')
            uppers = [row for row in data['rows'] if row.get('negative') == [shared]
                      and row.get('positive') == ['unit']*7]
            self.assertEqual(len(uppers), 1)
            self.assertTrue(any(row.get('family') == 'context:M2:shared_dot_matrix_rows' for row in data['rows']))
            corrupt = json.loads(json.dumps(data))
            corrupt['state'][shared] = '8'
            with self.assertRaisesRegex(ValueError, 'inequalities'):
                runner['replay'](corrupt)
            runner['mixed_dots'](root/'first/finite-replay.json', root/'second', 120, include_shared=True)
            repeated = json.loads((root/'second/finite-replay.json').read_text())
            self.assertEqual(len(data['rows']), len(repeated['rows']))
            self.assertEqual(set(data['inventory']), set(repeated['inventory']))

    def test_portable_extension_replay_rejects_corruption_and_is_idempotent(self):
        runner = runpy.run_path(str(Path(__file__).parent/'research/variable-proof-experiments.py'))
        field = FiniteField(2, (1, 1, 1))
        book = ConstraintInventory(field)
        matrix = naive_scheme(field, 2, 2, 2).as_tensor_scheme()
        key = book.register('M2', matrix.target)
        book.add_export(WitnessedOrderExport(rank_upper_constraint('unit', key, matrix),
            (('unit', book.inventory['unit']), (key, matrix.target))), 'naive_matrix_upper')
        book.add_detector('unit', 2)
        system = book.system(2, 5)
        state = RationalStateCertificate(system, tuple(Fraction(1 if key == 'unit' else 5)
                                                      for key in system.keys))
        result = StateSolverResult('finite_feasible', 'exact small input state', state=state)
        initial = runner['literal_export'](book, result)
        self.assertEqual(runner['replay'](initial)['checked'], ['finite_feasible'])
        with tempfile.TemporaryDirectory() as folder, redirect_stdout(io.StringIO()):
            root = Path(folder)
            source = root/'input.json'
            source.write_text(json.dumps(initial)+'\n')
            first, second = root/'first', root/'second'
            runner['mixed_dots'](source, first, 120)
            data = json.loads((first/'finite-replay.json').read_text())
            report = json.loads((first/'report.json').read_text())
            self.assertEqual(report['status'], 'finite_feasible')
            self.assertFalse(report['supplied_rank_seven_matrix_row'])
            self.assertEqual(runner['replay'](data)['checked'], ['finite_feasible'])
            square_key = next(row['positive'][0] for row in data['rows']
                if row.get('label') == 'mixed-dots-square-matrix')
            self.assertGreaterEqual(Fraction(data['state'][square_key]), 15)
            self.assertTrue(any(row.get('label') == 'mixed-dots-two-copy-matrix'
                                for row in data['rows']))
            self.assertTrue(any(row.get('family') == 'context:M2:mixed_dot_matrix_rows'
                                for row in data['rows']))
            swaps = [row for row in data['rows'] if row.get('family') == 'factor_swap:mixed_dots']
            self.assertEqual(len(swaps), 6)
            for row in swaps:
                self.assertEqual(data['state'][row['positive'][0]], data['state'][row['negative'][0]])
            self.assertTrue(any(item.get('source_entries') == 1_000_000
                                for item in report['skipped']))
            corrupt = json.loads(json.dumps(data))
            corrupt['state'][square_key] = '14'
            with self.assertRaisesRegex(ValueError, 'inequalities'):
                runner['replay'](corrupt)
            runner['mixed_dots'](first/'finite-replay.json', second, 120)
            repeated = json.loads((second/'finite-replay.json').read_text())
            self.assertEqual(len(data['rows']), len(repeated['rows']))
            self.assertEqual(set(data['inventory']), set(repeated['inventory']))
            self.assertEqual(runner['replay'](repeated)['checked'], ['finite_feasible'])


if __name__ == '__main__':
    unittest.main()
