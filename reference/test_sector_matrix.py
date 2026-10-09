"""Matrix extraction checks full identities, including off-support zeros."""

from dataclasses import replace
import unittest

from .constraint_generators import WitnessedOrderExport
from .extraction import simplify_scheme
from .fields import FiniteField
from .maps import LinearMap
from .sector_matrix import (export_matrix_column_sum_order, export_retained_matrix_order,
                            export_retained_two_copy_matrix_order, export_star_matrix_order)
from .sectors import ThreeSectorConstruction
from .tensor_schemes import TensorScheme
from .witnessed_constraints import StateAssemblyLimit, StateAssemblyLimits


class SectorMatrixTests(unittest.TestCase):
    def test_two_retained_copies_supply_independent_matrix_rows(self):
        for field in (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)), FiniteField(11)):
            for a, h in ((2, 1), (3, 1), (2, 2), (3, 2), (4, 3)):
                with self.subTest(field=field, a=a, h=h):
                    construction = ThreeSectorConstruction(field, a, h)
                    exported = export_retained_two_copy_matrix_order(construction, 'sector', 'M2')
                    self.assertEqual(exported.row.positive, ('sector', 'sector'))
                    self.assertTrue(exported.verify())
                    supplied = TensorScheme.from_tensor(construction.retained)
                    direct = TensorScheme.direct_sum((supplied, supplied))
                    recovered = simplify_scheme(direct.restrict(*exported.row.maps))
                    self.assertEqual(recovered.terms, 8)
                    recovered.require_exact()

    def test_two_copy_row_rejects_wrong_copy_placement_and_preflights(self):
        field = FiniteField(2)
        construction = ThreeSectorConstruction(field, 2, 1)
        exported = export_retained_two_copy_matrix_order(construction, 'sector', 'M2')
        first = LinearMap.selection(field, 4, (0, 1, 0, 1))
        self.assertFalse(replace(exported, row=replace(exported.row,
                                                     maps=(first, *exported.row.maps[1:]))).verify())
        for limits in (StateAssemblyLimits(max_blocks=1), StateAssemblyLimits(max_tensor_entries=100),
                       StateAssemblyLimits(max_map_entries=10)):
            with self.assertRaises(StateAssemblyLimit):
                export_retained_two_copy_matrix_order(construction, 'sector', 'M2', limits=limits)
        for construction in (None, ThreeSectorConstruction(field, 1, 1),
                             ThreeSectorConstruction(field, 2, 0)):
            with self.assertRaises(ValueError):
                export_retained_two_copy_matrix_order(construction, 'sector', 'M2')

    def test_column_sums_form_full_rectangular_matrix_tensors(self):
        for field in (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1))):
            for a, b, columns in ((2, 2, 2), (1, 3, 4), (2, 3, 3), (3, 2, 1)):
                with self.subTest(field=field, dimensions=(a, b, columns)):
                    exported = export_matrix_column_sum_order(field, a, b, columns, 'column', 'matrix')
                    self.assertTrue(exported.verify())
                    self.assertEqual(exported.row.positive, ('column',)*columns)
                    self.assertEqual(exported.registry['matrix'].shape, (a*b, b*columns, a*columns))

    def test_column_occurrences_do_not_collapse_in_characteristic_two(self):
        field = FiniteField(2)
        exported = export_matrix_column_sum_order(field, 2, 2, 2, 'column', 'M2')
        self.assertEqual(field.embed(len(exported.row.positive)), 0)
        self.assertEqual(exported.row.maps[0].rows[0], (1, 0, 0, 0, 1, 0, 0, 0))
        self.assertTrue(exported.verify())
        first = LinearMap(field, 8, tuple((0,)*8 for _ in range(4)))
        self.assertFalse(replace(exported, row=replace(exported.row,
                                                     maps=(first, *exported.row.maps[1:]))).verify())

    def test_column_sum_caps_preflight_large_dimensions(self):
        field = FiniteField(2)
        for limits in (StateAssemblyLimits(max_blocks=1), StateAssemblyLimits(max_tensor_entries=100),
                       StateAssemblyLimits(max_map_entries=10)):
            with self.assertRaises(StateAssemblyLimit):
                export_matrix_column_sum_order(field, 2, 2, 2, 'column', 'matrix', limits=limits)
        with self.assertRaises(StateAssemblyLimit):
            export_matrix_column_sum_order(field, 10**9, 10**9, 2, 'column', 'matrix')
        for value in (0, -1, True, 1.5):
            with self.assertRaises(ValueError):
                export_matrix_column_sum_order(field, 2, 2, value, 'column', 'matrix')

    def test_star_trace_correction_full_identity_and_signs(self):
        for field in (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)),
                      FiniteField(5), FiniteField(11)):
            with self.subTest(field=field):
                exported = export_star_matrix_order(field, 'square', 'M2')
                self.assertTrue(exported.verify())
                self.assertEqual(exported.row.maps[1].rows[3][3], field.embed(-1))
                self.assertEqual(exported.row.maps[2].rows[1][2], field.embed(-1))
                # Removing the trace correction loses some coefficient identities.
                third = list(exported.row.maps[2].rows)
                third[1] = tuple(field.embed(-1) if i == 2 else 0 for i in range(9))
                wrong = replace(exported.row, maps=(*exported.row.maps[:2], LinearMap(field, 9, third)))
                self.assertFalse(WitnessedOrderExport(wrong, exported.inventory).verify())

    def test_full_rectangular_identity_across_fields_and_parameters(self):
        for field in (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)), FiniteField(11)):
            for a, h in ((1, 1), (2, 1), (2, 2), (3, 1), (3, 2)):
                with self.subTest(field=field, a=a, h=h):
                    exported = export_retained_matrix_order(ThreeSectorConstruction(field, a, h),
                                                           'square', 'matrix')
                    self.assertTrue(exported.verify())
                    self.assertEqual(exported.registry['matrix'].shape, (a*a, 2*a, 2*a))

    def test_ordinary_row_extracts_a_matrix_scheme_from_sector_coefficients(self):
        field = FiniteField(2, (1, 1, 1))
        construction = ThreeSectorConstruction(field, 2, 1)
        exported = export_retained_matrix_order(construction, 'square', 'M2')
        sector = TensorScheme.from_tensor(construction.retained)
        matrix = simplify_scheme(sector.tensor_product(sector).restrict(*exported.row.maps))
        matrix.require_exact()
        self.assertEqual(matrix.terms, 8)
        self.assertEqual(matrix.target, exported.registry['M2'])
        self.assertEqual(matrix.to_matrix_scheme(2, 2, 2).multiply(((1, 2), (3, 1)),
                                                                  ((2, 1), (0, 3))),
                         [[2, 0], [1, 0]])

    def test_reversing_middle_input_is_necessary(self):
        field = FiniteField(2)
        exported = export_retained_matrix_order(ThreeSectorConstruction(field, 2, 1),
                                               'square', 'M2')
        wrong = replace(exported.row, maps=(LinearMap.identity(field, 4), *exported.row.maps[1:]))
        self.assertFalse(WitnessedOrderExport(wrong, exported.inventory).verify())

    def test_allocation_caps_and_invalid_parameters(self):
        field = FiniteField(2)
        for limits in (StateAssemblyLimits(max_tensor_entries=100),
                       StateAssemblyLimits(max_map_entries=100)):
            with self.assertRaises(StateAssemblyLimit):
                export_retained_matrix_order(ThreeSectorConstruction(field, 2, 1),
                                             'square', 'M2', limits=limits)
        with self.assertRaises(StateAssemblyLimit):
            export_retained_matrix_order(ThreeSectorConstruction(field, 10**9, 10**9),
                                         'square', 'matrix')
        for construction in (None, ThreeSectorConstruction(field, 0, 1),
                             ThreeSectorConstruction(field, 2, 0)):
            with self.assertRaises(ValueError):
                export_retained_matrix_order(construction, 'square', 'matrix')
        for keys in (('', 'matrix'), ('square', None), ('same', 'same')):
            with self.assertRaises(ValueError):
                export_retained_matrix_order(ThreeSectorConstruction(field, 2, 1), *keys)
            with self.assertRaises(ValueError):
                export_star_matrix_order(field, *keys)
        with self.assertRaises(ValueError):
            export_star_matrix_order(None, 'square', 'matrix')
        with self.assertRaises(StateAssemblyLimit):
            export_star_matrix_order(field, 'square', 'matrix',
                                     limits=StateAssemblyLimits(max_map_entries=10))


if __name__ == '__main__':
    unittest.main()
