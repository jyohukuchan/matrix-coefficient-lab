"""Independent actual-coefficient, all-characteristic and corruption controls."""
from copy import deepcopy
from itertools import combinations
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from .catalyst_search_experiment import field_matrix_rank
from .fields import FiniteField
from .star_degeneration import star_pencil_tensor
from .star_koszul import (StarKoszulLimit, StarKoszulLimits, star_integer_target_support,
                         star_integer_koszul, verify_star_integer_components)
from .tensors import FiniteTensor, matrix_multiplication_tensor


def certificate():
    return json.loads(Path(__file__).with_name('star_koszul.json').read_text())


class StarKoszulTests(unittest.TestCase):
    def test_integer_component_certificate(self):
        control = verify_star_integer_components()
        self.assertEqual((control.rank_lower, control.rank_one_factor, control.border_rank_lower), (1054, 70, 16))
        self.assertEqual((control.component_count, control.operation_count), (550, 1244))
        self.assertLessEqual(control.largest_integer_bits, 4)
        self.assertLess(Path(__file__).with_name('star_koszul.json').stat().st_size, 100_000)

    def test_complete_actual_tensor_and_independent_wedge_reconstruction(self):
        actual = star_integer_koszul()
        inputs = {subset: i for i, subset in enumerate(combinations(range(9), 4))}
        outputs = tuple(combinations(range(9), 5))
        for field in (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)), FiniteField(5)):
            target = FiniteTensor.direct_sum((FiniteTensor.unit(field),
                matrix_multiplication_tensor(field, 2, 2, 2).tensor_product(star_pencil_tensor(field))))
            support = set(star_integer_target_support())
            for x in range(9):
                for y in range(13):
                    for z in range(13):
                        self.assertEqual(target.coefficient(x, y, z), int((x, y, z) in support))
            expected = {}
            # Construct by deleting the added axis from the output wedge,
            # independently of the implementation's input-wedge insertion.
            for row, subset in enumerate(outputs):
                for position, x in enumerate(subset):
                    column = inputs[subset[:position]+subset[position+1:]]
                    sign = field.embed((-1)**position)
                    for y in range(13):
                        for z in range(13):
                            value = field.mul(sign, target.coefficient(x, y, z))
                            if value:
                                expected[(row*13+z, column*13+y)] = value
            self.assertEqual(expected, {key: field.embed(value) for key, value in actual.items()})

    def test_field_component_ranks_including_extension_scalar(self):
        sparse, data = star_integer_koszul(), certificate()
        for field in (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)), FiniteField(5)):
            scalar = field.p if field.degree > 1 else 1
            total = 0
            for component in data['components']:
                block = tuple(tuple(field.mul(scalar, field.embed(sparse.get((i, j), 0)))
                                    for j in component['columns']) for i in component['rows'])
                rank = field_matrix_rank(field, block, len(component['columns']))
                self.assertGreaterEqual(rank, component['unit_rank'])
                total += rank
            self.assertGreaterEqual(total, 1054)

    def test_corrupt_indices_operations_and_assertions_rejected(self):
        data = certificate()
        bad = deepcopy(data)
        bad['components'][1]['rows'][0] = bad['components'][0]['rows'][0]
        with self.assertRaisesRegex(ValueError, 'disjoint'):
            verify_star_integer_components(certificate=bad)
        bad = deepcopy(data)
        component = next(c for c in bad['components'] if c['operations'])
        component['operations'][0] = ['row_add', 0, 0, 1]
        with self.assertRaisesRegex(ValueError, 'distinct'):
            verify_star_integer_components(certificate=bad)
        bad = deepcopy(data)
        bad['components'][0]['determinant'] *= -1
        with self.assertRaisesRegex(ValueError, 'determinant'):
            verify_star_integer_components(certificate=bad)
        bad = deepcopy(data)
        bad['rank_one_factor'] = 69
        with self.assertRaises(ValueError):
            verify_star_integer_components(certificate=bad)

    def test_caps_precede_allocations_and_bit_growth(self):
        for limits in (StarKoszulLimits(max_nonzero_entries=0), StarKoszulLimits(max_coordinates=0),
                       StarKoszulLimits(max_integer_bits=0)):
            with patch('reference.star_koszul._certificate', side_effect=AssertionError('certificate allocation')):
                with self.assertRaises(StarKoszulLimit):
                    verify_star_integer_components(limits=limits)
        for limits in (StarKoszulLimits(max_component_entries=1), StarKoszulLimits(max_components=1),
                       StarKoszulLimits(max_operations=1243), StarKoszulLimits(max_certificate_bytes=1)):
            with patch('reference.star_koszul.star_integer_koszul', side_effect=AssertionError('matrix allocation')):
                with self.assertRaises(StarKoszulLimit):
                    verify_star_integer_components(limits=limits)
        data = certificate()
        component = next(c for c in data['components'] if c['operations'])
        component['operations'][0] = ['row_add', 0, 1, 1 << 100]
        with self.assertRaises(StarKoszulLimit):
            verify_star_integer_components(certificate=data, limits=StarKoszulLimits(max_integer_bits=8))
        with self.assertRaises(ValueError):
            verify_star_integer_components(limits=None)
