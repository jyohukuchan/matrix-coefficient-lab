"""Mixed-dot orders check complete coefficient arrays and preallocation caps."""

from dataclasses import replace
from unittest.mock import patch
import unittest

from .fields import FiniteField
from .maps import LinearMap
from .mixed_dot_orders import (export_mixed_dot_matrix_order,
    export_mixed_dot_square_subrank_order, export_mixed_dot_star_sum_order,
    export_mixed_dot_subrank_order, export_mixed_dot_two_copy_matrix_order, mixed_dot_tensor)
from .subrank_constraints import export_subrank_product_order
from .tensors import FiniteTensor, matrix_multiplication_tensor
from .witnessed_constraints import StateAssemblyLimit, StateAssemblyLimits


FIELDS = (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)), FiniteField(11))


class MixedDotOrderTests(unittest.TestCase):
    def test_square_supplies_fifteen_independent_units(self):
        for field in FIELDS:
            with self.subTest(field=field):
                exported = export_mixed_dot_square_subrank_order(field, 'square')
                self.assertEqual(exported.row.negative, ('unit',)*15)
                self.assertTrue(exported.verify())
                indices = [row.index(1) for row in exported.row.maps[0].rows]
                indices[1] = indices[0]
                wrong = replace(exported, row=replace(exported.row,
                    maps=(LinearMap.selection(field, 25, indices), *exported.row.maps[1:])))
                self.assertFalse(wrong.verify())

    def test_three_units_are_independent_in_every_field(self):
        for field in FIELDS:
            with self.subTest(field=field):
                exported = export_mixed_dot_subrank_order(field, 'mixed')
                self.assertEqual(exported.row.negative, ('unit',)*3)
                self.assertTrue(exported.verify())
                actual = exported.registry['mixed'].restrict(*exported.row.maps)
                expected = FiniteTensor.from_function(field, (3, 3, 3),
                    lambda i, j, k: int(i == j == k))
                self.assertEqual(actual, expected)
                square = export_subrank_product_order(exported, exported, 'square')
                self.assertEqual(square.row.negative, ('unit',)*9)
                self.assertTrue(square.verify())

    def test_star_plus_dot_keeps_the_dot_first_leg_independent(self):
        for field in FIELDS:
            with self.subTest(field=field):
                exported = export_mixed_dot_star_sum_order(field, 'mixed', 'star', 'dot')
                target = FiniteTensor.direct_sum((exported.registry['star'], exported.registry['dot']))
                self.assertEqual(target.shape, (3, 5, 5))
                self.assertEqual(exported.registry['mixed'].restrict(*exported.row.maps), target)
                first = LinearMap(field, 5,
                    (*exported.row.maps[0].rows[:2], (0, 0, 0, 0, 0)))
                wrong = replace(exported, row=replace(exported.row,
                    maps=(first, *exported.row.maps[1:])))
                self.assertFalse(wrong.verify())

    def test_square_matrix_full_coefficients_and_trace_correction(self):
        for field in FIELDS:
            with self.subTest(field=field):
                exported = export_mixed_dot_matrix_order(field, 'square', 'M2')
                source = mixed_dot_tensor(field).tensor_power(2)
                self.assertEqual(exported.registry['square'], source)
                self.assertEqual(source.restrict(*exported.row.maps),
                    matrix_multiplication_tensor(field, 2, 2, 2))
                rows = list(exported.row.maps[2].rows)
                rows[1] = (0,)*25
                wrong = replace(exported, row=replace(exported.row,
                    maps=(*exported.row.maps[:2], LinearMap(field, 25, rows))))
                self.assertFalse(wrong.verify())

    def test_two_copies_compute_independent_rows_with_shared_second_factor(self):
        for field in FIELDS:
            with self.subTest(field=field):
                exported = export_mixed_dot_two_copy_matrix_order(field, 'mixed', 'M2')
                source = FiniteTensor.direct_sum((mixed_dot_tensor(field),)*2)
                self.assertEqual(exported.row.positive, ('mixed', 'mixed'))
                self.assertEqual(source.restrict(*exported.row.maps),
                    matrix_multiplication_tensor(field, 2, 2, 2))
                rows = list(exported.row.maps[1].rows)
                rows[0] = (0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
                wrong = replace(exported, row=replace(exported.row,
                    maps=(exported.row.maps[0], LinearMap(field, 10, rows), exported.row.maps[2])))
                self.assertFalse(wrong.verify())

    def test_caps_precede_tensor_allocation_including_intermediate_maps(self):
        cases = (
            (export_mixed_dot_subrank_order, ('mixed',), StateAssemblyLimits(max_blocks=2)),
            (export_mixed_dot_star_sum_order, ('mixed', 'star', 'dot'), StateAssemblyLimits(max_blocks=1)),
            (export_mixed_dot_matrix_order, ('square', 'M2'), StateAssemblyLimits(max_tensor_entries=15000)),
            (export_mixed_dot_matrix_order, ('square', 'M2'), StateAssemblyLimits(max_map_entries=400)),
            (export_mixed_dot_matrix_order, ('square', 'M2'), StateAssemblyLimits(max_constraint_entries=90)),
            (export_mixed_dot_square_subrank_order, ('square',), StateAssemblyLimits(max_blocks=14)),
            (export_mixed_dot_square_subrank_order, ('square',), StateAssemblyLimits(max_map_entries=1000)),
            (export_mixed_dot_two_copy_matrix_order, ('mixed', 'M2'), StateAssemblyLimits(max_tensor_entries=999)),
            (export_mixed_dot_two_copy_matrix_order, ('mixed', 'M2'), StateAssemblyLimits(max_map_entries=119)),
            (export_mixed_dot_two_copy_matrix_order, ('mixed', 'M2'), StateAssemblyLimits(max_blocks=1)),
        )
        for function, keys, limits in cases:
            with self.subTest(function=function.__name__, limits=limits):
                with patch('reference.mixed_dot_orders.mixed_dot_tensor',
                           side_effect=AssertionError('allocation before preflight')):
                    with self.assertRaises(StateAssemblyLimit):
                        function(FiniteField(2), *keys, limits=limits)

    def test_invalid_fields_names_and_limits(self):
        for function, keys in ((export_mixed_dot_subrank_order, ('mixed',)),
                (export_mixed_dot_square_subrank_order, ('square',)),
                (export_mixed_dot_star_sum_order, ('mixed', 'star', 'dot')),
                (export_mixed_dot_two_copy_matrix_order, ('mixed', 'M2')),
                (export_mixed_dot_matrix_order, ('square', 'M2'))):
            with self.subTest(function=function.__name__):
                with self.assertRaises(ValueError):
                    function(None, *keys)
                with self.assertRaises(ValueError):
                    function(FiniteField(2), *keys, limits=None)
        for keys in (('same', 'same', 'dot'), ('mixed', '', 'dot'), ('mixed', 'star', None)):
            with self.assertRaises(ValueError):
                export_mixed_dot_star_sum_order(FiniteField(2), *keys)
        with self.assertRaises(ValueError):
            export_mixed_dot_subrank_order(FiniteField(2), 'unit')
        with self.assertRaises(ValueError):
            export_mixed_dot_matrix_order(FiniteField(2), 'same', 'same')


if __name__ == '__main__':
    unittest.main()
