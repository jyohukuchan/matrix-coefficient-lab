"""Subrank lower rows retain ordinary coefficient maps and scalar occurrences."""

from dataclasses import replace
from fractions import Fraction
import unittest
from unittest.mock import patch

from .constraint_generators import WitnessedOrderExport
from .fields import FiniteField
from .finite_state_solver import RationalStateCertificate
from .maps import LinearMap
from .proof_pipeline import ConstraintInventory
from .sectors import convolution_tensor
from .subrank_constraints import (export_convolution_subrank_order,
                                  export_matrix_subrank_order,
                                  export_projective_convolution_subrank_order,
                                  export_selected_subrank_order,
                                  export_subrank_order,
                                  export_subrank_product_order)
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, matrix_multiplication_tensor
from .witnessed_constraints import StateAssemblyLimit, StateAssemblyLimits


F2 = FiniteField(2)
F3 = FiniteField(3)
F4 = FiniteField(2, (1, 1, 1))


def diagonal(field, n):
    return FiniteTensor.from_function(field, (n,)*3, lambda i, j, k: int(i == j == k))


class ConvolutionSubrankTests(unittest.TestCase):
    def test_all_available_small_field_sizes_and_nonsquare_inputs(self):
        for field in (F2, F3, F4):
            for n in range(1, field.order+1):
                with self.subTest(field=field.order, n=n):
                    a, b = n+1, n+2
                    exported = export_convolution_subrank_order(field, a, b, n, "C")
                    self.assertIsInstance(exported, WitnessedOrderExport)
                    self.assertEqual(exported.row.positive, ("C",))
                    self.assertEqual(exported.row.negative, ("unit",)*n)
                    self.assertEqual(exported.nodes, tuple(range(n)))
                    self.assertEqual(exported.registry["C"], convolution_tensor(field, a, b))
                    # Checks the entire diagonal coefficient array, including
                    # all off-diagonal entries, rather than a sample of inputs.
                    self.assertEqual(exported.registry["C"].restrict(*exported.row.maps),
                                     diagonal(field, n))
                    for mapping in exported.row.maps[:2]:
                        self.assertTrue(all(not any(row[n:]) for row in mapping.rows))
                        for u in range(n):
                            for w, node in enumerate(exported.nodes):
                                evaluated = field.sum(field.mul(mapping.rows[u][i], field.pow(node, i))
                                                      for i in range(n))
                                self.assertEqual(evaluated, int(u == w))
                    self.assertTrue(exported.verify())

    def test_extension_nodes_use_genuine_encodings_and_reordering(self):
        exported = export_convolution_subrank_order(F4, 4, 5, 4, "C", nodes=(3, 1, 2, 0))
        self.assertEqual(exported.nodes, (3, 1, 2, 0))
        self.assertNotEqual(F4.embed(2), 2)
        self.assertTrue(any(value >= F4.p for mapping in exported.row.maps
                            for row in mapping.rows for value in row))
        self.assertEqual(exported.registry["C"].restrict(*exported.row.maps), diagonal(F4, 4))

    def test_inverse_vandermonde_transpose_is_essential(self):
        exported = export_convolution_subrank_order(F3, 3, 3, 3, "C")
        first, second, output = exported.row.maps
        self.assertNotEqual(first.rows, tuple(zip(*first.rows)))
        wrong = LinearMap(F3, 3, tuple(zip(*first.rows)))
        damaged = replace(exported, row=replace(exported.row, maps=(wrong, second, output)))
        self.assertFalse(damaged.verify())
        # Output uses evaluations, not a second interpolation matrix.
        wrong_output = LinearMap(F3, output.input_size,
                                 tuple(row+(0,)*(output.input_size-3) for row in first.rows))
        self.assertFalse(replace(exported, row=replace(exported.row,
                         maps=(first, second, wrong_output))).verify())

    def test_corrupted_sign_and_high_input_degree_are_rejected(self):
        exported = export_convolution_subrank_order(F3, 4, 5, 3, "C")
        first, second, output = exported.row.maps
        wrong_sign = LinearMap(F3, first.input_size,
                               (tuple(F3.neg(x) for x in first.rows[0]),)+first.rows[1:])
        self.assertFalse(replace(exported, row=replace(exported.row,
                         maps=(wrong_sign, second, output))).verify())
        # Padding by zero restricts to the degree < n subspace; allowing the
        # fourth coefficient changes actual formal tensor coefficients.
        rows = (first.rows[0][:-1]+(1,),)+first.rows[1:]
        wrong_degree = LinearMap(F3, 4, rows)
        self.assertFalse(replace(exported, row=replace(exported.row,
                         maps=(wrong_degree, second, output))).verify())

    def test_missing_nodes_and_invalid_parameters_are_rejected(self):
        for field, a, b, n in ((F2, 3, 3, 3), (F3, 2, 4, 3),
                                (F3, -1, 3, 1), (F3, 3, True, 1),
                                (F3, 3, 3, True), (3, 3, 3, 1)):
            with self.subTest(parameters=(field, a, b, n)), self.assertRaises(ValueError):
                export_convolution_subrank_order(field, a, b, n, "C")
        for nodes in ((0,), (0, 0), (0, 3), (0, True), (0, 1, 2)):
            with self.subTest(nodes=nodes), self.assertRaises(ValueError):
                export_convolution_subrank_order(F3, 2, 4, 2, "C", nodes=nodes)
        for key in ("", 1):
            with self.assertRaises(ValueError):
                export_convolution_subrank_order(F3, 2, 2, 2, key)

    def test_zero_diagonal_and_empty_inputs(self):
        for a, b in ((0, 3), (2, 0), (0, 0), (2, 3)):
            exported = export_convolution_subrank_order(F2, a, b, 0, "C")
            self.assertEqual(exported.row.negative, ())
            self.assertEqual(exported.nodes, ())
            self.assertEqual(exported.registry["C"].restrict(*exported.row.maps), diagonal(F2, 0))

    def test_caps_precede_tensor_polynomial_and_map_allocation(self):
        caps = (StateAssemblyLimits(max_tensor_entries=1),
                StateAssemblyLimits(max_map_entries=1),
                StateAssemblyLimits(max_blocks=1),
                StateAssemblyLimits(max_constraint_entries=1))
        for limits in caps:
            with self.subTest(limits=limits), \
                 patch("reference.subrank_constraints.convolution_tensor", side_effect=AssertionError("tensor allocation")), \
                 patch("reference.subrank_constraints.lagrange_basis", side_effect=AssertionError("polynomial allocation")), \
                 patch.object(LinearMap, "__post_init__", side_effect=AssertionError("map allocation")), \
                 self.assertRaises(StateAssemblyLimit):
                export_convolution_subrank_order(F3, 3, 4, 3, "C", limits=limits)
        # Empty axes do not bypass the cost of map column bookkeeping.
        with patch("reference.subrank_constraints.convolution_tensor", side_effect=AssertionError("tensor allocation")):
            with self.assertRaisesRegex(StateAssemblyLimit, "bookkeeping"):
                export_convolution_subrank_order(F2, 0, 10**9, 0, "C")


class GenericSubrankTests(unittest.TestCase):
    def test_matrix_diagonal_works_even_with_more_units_than_field_elements(self):
        exported = export_matrix_subrank_order(F2, 3, "M3")
        self.assertEqual(exported.row.negative, ("unit",)*3)
        expected = LinearMap.selection(F2, 9, (0, 4, 8))
        self.assertEqual(exported.row.maps, (expected,)*3)
        self.assertEqual(exported.registry["M3"].restrict(*exported.row.maps), diagonal(F2, 3))
        self.assertTrue(export_matrix_subrank_order(F2, 3, "M3", n=2).verify())

    def test_existing_scheme_and_generic_selectors(self):
        tensor = matrix_multiplication_tensor(F4, 2, 2, 2)
        scheme = TensorScheme.from_tensor(tensor)
        from_scheme = export_selected_subrank_order(scheme, ((0, 3),)*3, "M2")
        from_maps = export_subrank_order(tensor, from_scheme.row.maps, "M2")
        self.assertEqual(from_scheme, from_maps)
        self.assertEqual(tensor.restrict(*from_maps.row.maps), diagonal(F4, 2))
        bad = replace(scheme, a=((0,)*4,)+scheme.a[1:])
        with self.assertRaisesRegex(ValueError, "invalid tensor coefficient"):
            export_selected_subrank_order(bad, ((0, 3),)*3, "M2")

    def test_selector_that_creates_cross_terms_fails_full_coefficient_check(self):
        tensor = convolution_tensor(F2, 2, 2)
        with self.assertRaisesRegex(ValueError, "full coefficient equality"):
            export_selected_subrank_order(tensor, ((0, 1), (0, 1), (0, 1)), "C")

    def test_generic_validation_of_map_field_size_keys_and_selectors(self):
        tensor = FiniteTensor.unit(F3)
        valid = (LinearMap.identity(F3, 1),)*3
        for maps in (valid[:2], (LinearMap.identity(F2, 1),)*3,
                     (LinearMap.identity(F3, 2),)*3,
                     (valid[0], valid[1], LinearMap(F3, 1, ()))):
            with self.assertRaises(ValueError):
                export_subrank_order(tensor, maps, "source")
        for selections in (((1,),)*3, ((True,),)*3, ((0,),)*2, ((0,), (), (0,))):
            with self.assertRaises(ValueError):
                export_selected_subrank_order(tensor, selections, "source")
        self.assertTrue(export_subrank_order(tensor, valid, "unit").verify())
        with self.assertRaisesRegex(ValueError, "unit ID"):
            export_convolution_subrank_order(F3, 2, 2, 2, "unit")
        for d, n in ((0, None), (True, None), (2, 3), (2, True), (2, -1)):
            with self.assertRaises(ValueError):
                export_matrix_subrank_order(F2, d, "M", n=n)

    def test_matrix_and_selector_caps_precede_allocation(self):
        limits = StateAssemblyLimits(max_map_entries=1)
        with patch("reference.subrank_constraints.matrix_multiplication_tensor", side_effect=AssertionError("tensor allocation")):
            with self.assertRaises(StateAssemblyLimit):
                export_matrix_subrank_order(F2, 2, "M", limits=limits)
        with patch.object(LinearMap, "selection", side_effect=AssertionError("map allocation")):
            with self.assertRaises(StateAssemblyLimit):
                export_selected_subrank_order(FiniteTensor.unit(F2), ((0,),)*3,
                                               "U", limits=limits)


class ProjectiveConvolutionSubrankTests(unittest.TestCase):
    def test_all_available_projective_sizes_including_nonsquare_inputs(self):
        for field in (F2, F3, F4):
            for n in range(field.order+2):
                with self.subTest(field=field.order, n=n):
                    exported = export_projective_convolution_subrank_order(field, n+1, n+2, n, "C")
                    self.assertEqual(exported.row.negative, ("unit",)*n)
                    self.assertEqual(exported.nodes, tuple(range(max(0, n-1))))
                    self.assertEqual(exported.registry["C"].restrict(*exported.row.maps), diagonal(field, n))
                    self.assertTrue(exported.verify())
                    if n:
                        self.assertIn("projective", exported.row.label)
                        self.assertIn("infinity", exported.row.label)
                        output = exported.row.maps[2]
                        self.assertEqual(output.rows[-1], tuple(int(k == 2*(n-1))
                                                               for k in range(output.input_size)))
                        for mapping in exported.row.maps[:2]:
                            self.assertTrue(all(not any(row[n:]) for row in mapping.rows))
                            self.assertEqual(tuple(row[n-1] for row in mapping.rows), (0,)*(n-1)+(1,))

    def test_projective_node_supplies_f2_c33_and_f3_c44_extra_unit(self):
        for field, n in ((F2, 3), (F3, 4)):
            with self.subTest(field=field.order):
                exported = export_projective_convolution_subrank_order(field, n, n, n, "C")
                self.assertEqual(len(exported.row.negative), field.order+1)
                self.assertEqual(exported.registry["C"].restrict(*exported.row.maps), diagonal(field, n))
                with self.assertRaisesRegex(ValueError, "not enough field elements"):
                    export_convolution_subrank_order(field, n, n, n, "C")
        exported = export_projective_convolution_subrank_order(F2, 3, 3, 3, "C33")
        book = ConstraintInventory(F2)
        book.add_export(exported, "projective_subrank_lower")
        system = book.system(2, 5)
        self.assertEqual(system.sparse_vectors, ((("C33", 1), ("unit", -3)),))
        with self.assertRaises(ValueError):
            RationalStateCertificate(system, (Fraction(1), Fraction(2))).require_valid()
        RationalStateCertificate(system, (Fraction(1), Fraction(3))).require_valid()

    def test_extension_finite_nodes_and_coefficients_are_actual_field_encodings(self):
        exported = export_projective_convolution_subrank_order(F4, 5, 6, 5, "C", nodes=(3, 1, 2, 0))
        self.assertEqual(exported.nodes, (3, 1, 2, 0))
        self.assertTrue(all(type(node) is int for node in exported.nodes))
        self.assertNotEqual(F4.embed(2), 2)
        self.assertTrue(any(value >= F4.p for mapping in exported.row.maps
                            for row in mapping.rows for value in row))
        self.assertEqual(exported.registry["C"].restrict(*exported.row.maps), diagonal(F4, 5))

    def test_wrong_infinity_degree_sign_and_input_padding_are_rejected(self):
        exported = export_projective_convolution_subrank_order(F3, 4, 5, 3, "C", nodes=(1, 2))
        first, second, output = exported.row.maps
        # The old source top degree is 7, whereas the subspace top product
        # degree is 4. Infinity must extract degree 4 in this nonsquare case.
        self.assertEqual(output.input_size, 8)
        wrong_top = LinearMap(F3, output.input_size, output.rows[:-1]+(
            tuple(int(k == output.input_size-1) for k in range(output.input_size)),))
        wrong_adjacent = LinearMap(F3, output.input_size, output.rows[:-1]+(
            tuple(int(k == 2*(3-1)-1) for k in range(output.input_size)),))
        for damaged_output in (wrong_top, wrong_adjacent):
            self.assertFalse(replace(exported, row=replace(exported.row,
                             maps=(first, second, damaged_output))).verify())
        wrong_sign = LinearMap(F3, first.input_size,
                               first.rows[:-1]+(tuple(F3.neg(x) for x in first.rows[-1]),))
        self.assertFalse(replace(exported, row=replace(exported.row,
                         maps=(wrong_sign, second, output))).verify())
        bad_padding = LinearMap(F3, first.input_size,
                                 (first.rows[0][:-1]+(1,),)+first.rows[1:])
        self.assertFalse(replace(exported, row=replace(exported.row,
                         maps=(bad_padding, second, output))).verify())

    def test_parameter_and_finite_node_validation(self):
        for field, a, b, n in ((F2, 4, 4, 4), (F3, 3, 5, 4),
                                (F3, 4, -1, 1), (F3, 4, 4, True), (3, 4, 4, 1)):
            with self.subTest(parameters=(field, a, b, n)), self.assertRaises(ValueError):
                export_projective_convolution_subrank_order(field, a, b, n, "C")
        for nodes in ((0,), (0, 0), (0, 3), (0, None), (0, True), (0, 1, 2)):
            with self.subTest(nodes=nodes), self.assertRaises(ValueError):
                export_projective_convolution_subrank_order(F3, 3, 4, 3, "C", nodes=nodes)
        with self.assertRaises(ValueError):
            export_projective_convolution_subrank_order(F3, 3, 4, 3, "")

    def test_zero_delegates_affine_and_one_uses_only_constant_infinity(self):
        affine_zero = export_convolution_subrank_order(F2, 0, 3, 0, "C")
        self.assertEqual(export_projective_convolution_subrank_order(F2, 0, 3, 0, "C"), affine_zero)
        one = export_projective_convolution_subrank_order(F2, 2, 4, 1, "C")
        self.assertEqual(one.nodes, ())
        self.assertEqual(one.row.maps[0].rows, ((1, 0),))
        self.assertEqual(one.row.maps[1].rows, ((1, 0, 0, 0),))
        self.assertEqual(one.row.maps[2].rows, ((1, 0, 0, 0, 0),))
        self.assertEqual(one.registry["C"].restrict(*one.row.maps), FiniteTensor.unit(F2))

    def test_caps_precede_polynomial_tensor_and_map_allocation(self):
        for limits in (StateAssemblyLimits(max_blocks=2),
                       StateAssemblyLimits(max_tensor_entries=1),
                       StateAssemblyLimits(max_map_entries=1),
                       StateAssemblyLimits(max_constraint_entries=1)):
            with self.subTest(limits=limits), \
                 patch("reference.subrank_constraints.convolution_tensor", side_effect=AssertionError("tensor allocation")), \
                 patch("reference.subrank_constraints.lagrange_basis", side_effect=AssertionError("polynomial allocation")), \
                 patch.object(LinearMap, "__post_init__", side_effect=AssertionError("map allocation")), \
                 self.assertRaises(StateAssemblyLimit):
                export_projective_convolution_subrank_order(F2, 3, 4, 3, "C", limits=limits)


class ProductAndStateTests(unittest.TestCase):
    def test_product_distributes_four_binary_units_with_unequal_source_axes(self):
        first = export_convolution_subrank_order(F2, 2, 3, 2, "A")
        second = export_matrix_subrank_order(F2, 2, "B")
        exported = export_subrank_product_order(first, second, "A*B")
        self.assertEqual(exported.row.negative, ("unit",)*4)
        source = first.registry["A"].tensor_product(second.registry["B"])
        self.assertEqual(source.shape, (8, 12, 16))
        self.assertEqual(exported.registry["A*B"], source)
        self.assertEqual(source.restrict(*exported.row.maps), diagonal(F2, 4))
        self.assertEqual(sum(diagonal(F2, 4).coefficients), 4)
        for axis, mapping in enumerate(exported.row.maps):
            for u in range(2):
                for v in range(2):
                    expected = tuple(F2.mul(x, y) for x in first.row.maps[axis].rows[u]
                                     for y in second.row.maps[axis].rows[v])
                    self.assertEqual(mapping.rows[u*2+v], expected)

    def test_product_preserves_extension_coefficients_and_nonmatching_counts(self):
        first = export_convolution_subrank_order(F4, 2, 3, 2, "A", nodes=(2, 3))
        second = export_convolution_subrank_order(F4, 3, 3, 3, "B", nodes=(3, 0, 2))
        exported = export_subrank_product_order(first, second, "A*B")
        self.assertEqual(len(exported.row.negative), 6)
        self.assertTrue(any(value >= 2 for mapping in exported.row.maps
                            for row in mapping.rows for value in row))
        self.assertEqual(exported.registry["A*B"].restrict(*exported.row.maps), diagonal(F4, 6))

    def test_product_rejects_corrupted_witness_and_wrong_field_or_upper_row(self):
        first = export_convolution_subrank_order(F3, 2, 2, 2, "A")
        second = export_matrix_subrank_order(F3, 2, "B")
        damaged_map = LinearMap(F3, 2, ((0, 0), (0, 0)))
        damaged = replace(first, row=replace(first.row, maps=(damaged_map,)+first.row.maps[1:]))
        with self.assertRaisesRegex(ValueError, "full coefficient equality"):
            export_subrank_product_order(damaged, second, "A*B")
        with self.assertRaises(ValueError):
            export_subrank_product_order(first, export_matrix_subrank_order(F2, 2, "B"), "A*B")
        upper = replace(first, row=replace(first.row, positive=("unit", "unit"), negative=("A",)))
        with self.assertRaises(ValueError):
            export_subrank_product_order(upper, second, "A*B")

    def test_product_caps_precede_tensor_and_map_allocation(self):
        first = export_matrix_subrank_order(F2, 2, "A")
        second = export_matrix_subrank_order(F2, 2, "B")
        for limits in (StateAssemblyLimits(max_blocks=3),
                       StateAssemblyLimits(max_tensor_entries=100),
                       StateAssemblyLimits(max_map_entries=100)):
            with self.subTest(limits=limits), \
                 patch.object(FiniteTensor, "tensor_product", side_effect=AssertionError("tensor allocation")), \
                 patch.object(LinearMap, "__post_init__", side_effect=AssertionError("map allocation")), \
                 self.assertRaises(StateAssemblyLimit):
                export_subrank_product_order(first, second, "A*B", limits=limits)

    def test_zero_product_and_custom_unit_name(self):
        first = export_convolution_subrank_order(F3, 0, 2, 0, "A", unit_key="U")
        second = export_matrix_subrank_order(F3, 2, "B", unit_key="V")
        exported = export_subrank_product_order(first, second, "A*B", unit_key="scalar")
        self.assertEqual(exported.row.negative, ())
        self.assertEqual(exported.registry["scalar"], FiniteTensor.unit(F3))
        self.assertTrue(exported.verify())

    def test_book_integration_rejects_state_one_without_rank_additivity(self):
        book = ConstraintInventory(F2)
        exported = export_convolution_subrank_order(F2, 2, 3, 2, "C23")
        book.add_export(exported, "subrank_lower")
        system = book.system(2, 5)
        system.require_valid()
        self.assertEqual(system.sparse_vectors, ((("C23", 1), ("unit", -2)),))
        with self.assertRaises(ValueError):
            RationalStateCertificate(system, (Fraction(1), Fraction(1))).require_valid()
        RationalStateCertificate(system, (Fraction(1), Fraction(2))).require_valid()


if __name__ == "__main__":
    unittest.main()
