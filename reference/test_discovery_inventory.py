"""Finite discovery families retain coefficient witnesses and honest caps."""

from collections import Counter
from dataclasses import replace
from fractions import Fraction
import unittest

from .discovery_inventory import (DiscoveryBudget, commute_product_export,
                                  generate_discovery_inventory, scalar_lower_export,
                                  star_tensor)
from .fields import FiniteField
from .finite_state_solver import RationalStateCertificate
from .maps import LinearMap
from .tensors import FiniteTensor, matrix_multiplication_tensor
from .witnessed_constraints import (DetectorConstraint, StateAssemblyLimit,
                                    StateAssemblyLimits)


F2 = FiniteField(2)
F4 = FiniteField(2, (1, 1, 1))


class DiscoveryExportTests(unittest.TestCase):
    def test_scalar_lower_uses_actual_extension_inverse_and_selected_coordinate(self):
        tensor = FiniteTensor.from_function(F4, (2, 3, 2),
                                           lambda i, j, k: 2 if (i, j, k) == (1, 2, 1) else 0)
        exported = scalar_lower_export("extension", tensor)
        exported.require_valid()
        self.assertEqual(exported.row.positive, ("extension",))
        self.assertEqual(exported.row.negative, ("unit",))
        self.assertEqual(exported.row.maps[0].rows, ((0, 1),))
        self.assertEqual(exported.row.maps[1].rows, ((0, 0, 1),))
        self.assertEqual(exported.row.maps[2].rows, ((0, F4.inv(2)),))
        self.assertNotEqual(F4.inv(2), 1)
        self.assertEqual(tensor.restrict(*exported.row.maps), FiniteTensor.unit(F4))
        # Encoding 2 is an extension element, whereas scalar 2 vanishes.
        self.assertEqual(F4.embed(2), 0)
        broken = replace(exported.row, maps=exported.row.maps[:2]+(LinearMap(F4, 2, ((0, 1),)),))
        self.assertFalse(replace(exported, row=broken).verify())

    def test_scalar_lower_rejects_zero_and_wrong_unit_alias(self):
        for shape in ((2, 3, 2), (0, 2, 1)):
            with self.assertRaisesRegex(ValueError, "zero tensor"):
                scalar_lower_export("zero", FiniteTensor.zero(F4, shape))
        with self.assertRaisesRegex(ValueError, "unit key"):
            scalar_lower_export("unit", FiniteTensor(F4, (1, 1, 1), (2,)))
        actual_unit = scalar_lower_export("unit", FiniteTensor.unit(F4))
        actual_unit.require_valid()
        self.assertEqual(len(actual_unit.inventory), 1)

    def test_commuting_unequal_axis_factors_has_exact_maps_in_both_directions(self):
        first = FiniteTensor.from_function(F4, (2, 1, 3),
                                          lambda i, j, k: (i+2*k+1) % 4)
        second = FiniteTensor.from_function(F4, (1, 3, 2),
                                           lambda i, j, k: (2*j+k+1) % 4)
        forward = commute_product_export("first", first, "second", second)
        backward = commute_product_export("first", first, "second", second, reverse=True)
        forward.require_valid()
        backward.require_valid()
        self.assertNotEqual(forward.registry["first@second"], forward.registry["second@first"])
        self.assertEqual(forward.row.positive, backward.row.negative)
        self.assertEqual(forward.row.negative, backward.row.positive)
        for axis, (a, b) in enumerate(zip(first.shape, second.shape)):
            self.assertEqual(backward.row.maps[axis].compose(forward.row.maps[axis]),
                             LinearMap.identity(F4, a*b))
            self.assertEqual(forward.row.maps[axis].compose(backward.row.maps[axis]),
                             LinearMap.identity(F4, a*b))
        identity_maps = tuple(LinearMap.identity(F4, a*b) for a, b in zip(first.shape, second.shape))
        self.assertFalse(replace(forward, row=replace(forward.row, maps=identity_maps)).verify())

    def test_matrix_context_commutation_is_an_actual_permutation(self):
        matrix = matrix_multiplication_tensor(F2, 2, 2, 2)
        test = star_tensor(F2)
        forward = commute_product_export("M2", matrix, "star", test)
        backward = commute_product_export("M2", matrix, "star", test, reverse=True)
        forward.require_valid()
        backward.require_valid()
        self.assertEqual(forward.registry["M2@star"].shape, (8, 12, 12))
        self.assertNotEqual(forward.registry["M2@star"], forward.registry["star@M2"])
        with self.assertRaises(ValueError):
            commute_product_export("M2", matrix, "different-field", FiniteTensor.unit(F4))


class SmallDiscoveryInventoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.book, cls.skipped = generate_discovery_inventory(F2, contexts=())

    def test_small_f2_inventory_and_selected_finite_families_are_exact(self):
        self.book.system(2, 5).require_valid()
        families = Counter(self.book.families)
        self.assertGreater(families["scalar_lower"], 0)
        self.assertGreater(families["rank_upper"], 0)
        self.assertGreater(families["detector"], 0)
        self.assertNotIn("known_control", families)
        for family in ("scalar_lower", "rank_upper", "detector"):
            indices = tuple(i for i, value in enumerate(self.book.families) if value == family)
            selected = self.book.system(2, 5, indices)
            self.assertEqual(len(selected.constraints), families[family])
            selected.require_valid()
        skipped = {entry["operation"]: entry for entry in self.skipped}
        self.assertIn("sector", skipped)
        self.assertIn("determinant", skipped)
        self.assertIn("no field extension is implicit", skipped["exact_type_fourier"]["reason"])
        self.assertEqual(skipped["projective-C33"]["fallback"], "checked coordinate-pair scheme")

    def test_interned_aliases_link_upper_rows_and_exact_detector_products(self):
        matrix = matrix_multiplication_tensor(F2, 2, 2, 2)
        matrix_key = self.book._canonical[matrix]
        self.assertEqual(matrix_key, "M2*unit")
        self.assertNotIn("M2", self.book.inventory)
        for row in self.book.rows:
            if isinstance(row, DetectorConstraint):
                expected = matrix.tensor_product(self.book.inventory[row.test])
                self.assertEqual(self.book.inventory[row.product], expected)
                self.assertEqual(row.product, self.book._canonical[expected])
        upper = next(row for row, family in zip(self.book.rows, self.book.families)
                     if family == "rank_upper" and row.negative == (matrix_key,))
        self.assertEqual(upper.positive, ("unit",)*8)
        unit_detector = next(row for row in self.book.rows
                             if isinstance(row, DetectorConstraint) and row.test == "unit")
        self.assertEqual(unit_detector.product, upper.negative[0])

    def test_diagonal_lowers_reject_the_previous_unit_valued_convolution_state(self):
        indices = tuple(i for i, family in enumerate(self.book.families)
                        if family in ('subrank_lower', 'subrank_product'))
        system = self.book.system(2, 5, indices)
        system.require_valid()
        old_weak_values = RationalStateCertificate(
            system, (Fraction(1),)*len(system.keys))
        self.assertFalse(old_weak_values.verify_inequalities())
        lower = next(row for row, family in zip(self.book.rows, self.book.families)
                     if family == 'subrank_lower' and row.positive == ('C33',))
        self.assertEqual(lower.negative, ('unit',)*3)
        # Projective infinity supplies the third unit over the two-element field.
        self.assertEqual(len(lower.negative), F2.order+1)

    def test_extension_nodes_strengthen_the_c33_diagonal_lower(self):
        book, _ = generate_discovery_inventory(
            F4, contexts=(), budget=DiscoveryBudget(max_upper_terms=1))
        lower = next(row for row, family in zip(book.rows, book.families)
                     if family == 'subrank_lower' and row.positive == ('C33',))
        self.assertEqual(lower.negative, ('unit',)*3)
        row_index = book.rows.index(lower)
        book.system(2, 5, (row_index,)).require_valid()

    def test_row_and_decomposition_caps_return_explicit_incomplete_reports(self):
        capped, skipped = generate_discovery_inventory(F2, contexts=(), budget=DiscoveryBudget(max_rows=1))
        self.assertEqual(len(capped.rows), 1)
        self.assertTrue(any(entry["reason"] == "row cap reached" for entry in skipped))
        capped.system(2, 5).require_valid()
        sparse, skipped = generate_discovery_inventory(F2, contexts=(), budget=DiscoveryBudget(max_upper_terms=1))
        self.assertNotIn("rank_upper", sparse.families)
        self.assertTrue(any(entry["reason"] == "supplied decomposition term cap reached" for entry in skipped))
        sparse.system(2, 5).require_valid()

    def test_context_lift_cap_does_not_claim_a_missing_lift_or_impossibility(self):
        book, skipped = generate_discovery_inventory(
            F2, include_known_control=True, contexts=("M2",),
            budget=DiscoveryBudget(max_lift_entries=1, max_upper_terms=1))
        self.assertNotIn("context:M2:known_control", book.families)
        self.assertIn("known_control", book.families)
        self.assertIn("factor_swap", book.families)
        self.assertTrue(any(entry["operation"].startswith("lift-") and
                            entry["reason"] == "lift construction cap reached" for entry in skipped))
        # Selected upper/swap rows still carry full coefficient evidence.
        indices = tuple(i for i, family in enumerate(book.families)
                        if family in ("known_control", "factor_swap"))
        book.system(2, 5, indices).require_valid()

    def test_atom_budget_and_dense_budget_fail_honestly(self):
        for maximum in (1, 2, 5):
            with self.assertRaises(StateAssemblyLimit):
                generate_discovery_inventory(F2, contexts=(), budget=DiscoveryBudget(max_atoms=maximum))
        with self.assertRaises(StateAssemblyLimit):
            generate_discovery_inventory(F2, contexts=(), limits=StateAssemblyLimits(max_tensor_entries=1))

    def test_invalid_discovery_arguments(self):
        for kwargs in ({"include_known_control": 1}, {"contexts": ("M2", "M2")},
                       {"contexts": ("unsupported",)}):
            with self.assertRaises(ValueError):
                generate_discovery_inventory(F2, **kwargs)
        with self.assertRaises(ValueError):
            generate_discovery_inventory("F2", contexts=())
        for kwargs in ({"max_rows": 0}, {"max_atoms": True}, {"max_upper_terms": 1.0},
                       {"max_lift_entries": -1}):
            with self.assertRaises(ValueError):
                DiscoveryBudget(**kwargs)


if __name__ == "__main__":
    unittest.main()
