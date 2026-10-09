"""Independent coefficient, modular, unimodular and corruption controls."""

from itertools import combinations
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from .catalyst_search_experiment import field_matrix_rank
from .fields import FiniteField
from .maps import LinearMap
from .shared_dot_koszul import (SharedDotKoszulLimit, SharedDotKoszulLimits,
    shared_dot_three_field_koszul, shared_dot_three_integer_koszul, verify_shared_dot_three_integer_minor,
    shared_dot_two_integer_koszul, verify_shared_dot_two_integer_transform)
from .shared_dot_orders import shared_dot_tensor
from .tensors import FiniteTensor, matrix_multiplication_tensor


class SharedDotKoszulTests(unittest.TestCase):
    def test_two_dot_unimodular_operations_and_actual_field_matrices(self):
        control = verify_shared_dot_two_integer_transform()
        self.assertEqual((control.minor_size, control.minor_determinant, control.rank_one_factor), (124, -1, 6))
        self.assertEqual(control.border_rank_lower, 21)
        data = json.loads(Path(__file__).with_name('shared_dot_two_koszul.json').read_text())
        matrix = shared_dot_two_integer_koszul()
        for field in (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)), FiniteField(11)):
            target = FiniteTensor.direct_sum((FiniteTensor.unit(field),
                matrix_multiplication_tensor(field, 2, 2, 2).tensor_product(shared_dot_tensor(field, 2))))
            projected = target.restrict(LinearMap(field, 13, data['first_map']),
                                       LinearMap.identity(field, 13), LinearMap.identity(field, 13))
            domains, codomains = tuple(combinations(range(5), 2)), tuple(combinations(range(5), 3))
            for i, out in enumerate(codomains):
                for j, inp in enumerate(domains):
                    added = set(out)-set(inp)
                    for y in range(13):
                        for z in range(13):
                            expected = 0
                            if set(inp) <= set(out) and len(added) == 1:
                                axis, = added
                                sign = -1 if sum(v < axis for v in inp) % 2 else 1
                                expected = field.mul(field.embed(sign), projected.coefficient(axis, y, z))
                            self.assertEqual(field.embed(matrix[i*13+z][j*13+y]), expected)
            embedded = tuple(tuple(field.embed(value) for value in row) for row in matrix)
            self.assertGreaterEqual(field_matrix_rank(field, embedded, 130), 124)

    def test_two_dot_corrupt_operations_and_integer_caps(self):
        data = json.loads(Path(__file__).with_name('shared_dot_two_koszul.json').read_text())
        data['operations'][0] = ['row_add', 0, 0, 1]
        with self.assertRaises(ValueError):
            verify_shared_dot_two_integer_transform(certificate=data)
        with self.assertRaises(SharedDotKoszulLimit):
            verify_shared_dot_two_integer_transform(limits=SharedDotKoszulLimits(max_operations=902))
        with self.assertRaises(SharedDotKoszulLimit):
            verify_shared_dot_two_integer_transform(limits=SharedDotKoszulLimits(max_integer_bits=10))
        data['first_map'] = [[0]*13 for _ in range(5)]
        data['operations'] = []
        with self.assertRaisesRegex(ValueError, 'unit diagonal'):
            verify_shared_dot_two_integer_transform(certificate=data)

    def test_integer_minor_and_coefficient_reconstruction_over_control_fields(self):
        control = verify_shared_dot_three_integer_minor()
        self.assertEqual((control.minor_size, control.minor_determinant, control.rank_one_factor), (153, 1, 6))
        self.assertEqual(control.border_rank_lower, 26)
        data = json.loads(Path(__file__).with_name('shared_dot_three_koszul.json').read_text())
        for field in (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)), FiniteField(11)):
            target = FiniteTensor.direct_sum((FiniteTensor.unit(field),
                matrix_multiplication_tensor(field, 2, 2, 2).tensor_product(shared_dot_tensor(field, 3))))
            mapping = LinearMap(field, 17, data['first_map'])
            projected = target.restrict(mapping, LinearMap.identity(field, 17), LinearMap.identity(field, 17))
            matrix = shared_dot_three_field_koszul(field)
            domain = tuple(combinations(range(5), 2))
            codomain = tuple(combinations(range(5), 3))
            for i, out in enumerate(codomain):
                for j, inp in enumerate(domain):
                    added = set(out)-set(inp)
                    for y in range(17):
                        for z in range(17):
                            expected = 0
                            if set(inp) <= set(out) and len(added) == 1:
                                axis, = added
                                sign = -1 if sum(v < axis for v in inp) % 2 else 1
                                expected = field.mul(field.embed(sign), projected.coefficient(axis, y, z))
                            self.assertEqual(matrix[i*17+z][j*17+y], expected)
            minor = tuple(tuple(matrix[i][j] for j in data['minor_columns']) for i in data['minor_rows'])
            self.assertEqual(field_matrix_rank(field, minor, 153), 153)

    def test_corruption_and_caps_do_not_claim_exclusion(self):
        data = json.loads(Path(__file__).with_name('shared_dot_three_koszul.json').read_text())
        data['first_map'] = [[0]*17 for _ in range(5)]
        with self.assertRaisesRegex(ValueError, 'unit pivot'):
            verify_shared_dot_three_integer_minor(certificate=data)
        data['minor_rows'][1] = data['minor_rows'][0]
        with self.assertRaises(ValueError):
            verify_shared_dot_three_integer_minor(certificate=data)
        with patch('reference.shared_dot_koszul._certificate', side_effect=AssertionError('certificate read')):
            with self.assertRaises(SharedDotKoszulLimit):
                shared_dot_three_integer_koszul(limits=SharedDotKoszulLimits(max_matrix_entries=1))
        with self.assertRaises(ValueError):
            shared_dot_three_field_koszul(2)
