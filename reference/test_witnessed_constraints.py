"""Independent full-coefficient checks for finite dual -> catalytic maps."""

from dataclasses import replace
from fractions import Fraction
import unittest
from unittest.mock import patch

from .fields import FiniteField
from .maps import LinearMap
from .schemes import naive_scheme, strassen_scheme
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor
from .witnessed_constraints import (DetectorConstraint, IntegerStateCertificate,
                                    StateAssemblyLimit, StateAssemblyLimits,
                                    WitnessedOrder, WitnessedStateSystem,
                                    assemble_catalyst, clear_rational_dual,
                                    find_nonnegative_dual, rank_upper_constraint)


F2 = FiniteField(2)


def toy_system(k=9, scheme=None):
    scheme = naive_scheme(F2, 2, 2, 2).as_tensor_scheme() if scheme is None else scheme
    unit = FiniteTensor.unit(scheme.field)
    return WitnessedStateSystem(2, k, (("unit", unit), ("matrix", scheme.target)),
                               (rank_upper_constraint("unit", "matrix", scheme),
                                DetectorConstraint("unit", "matrix")))


class WitnessedStateTests(unittest.TestCase):
    def test_exact_dual_assembles_real_positive_gain_maps(self):
        system = toy_system()
        result = find_nonnegative_dual(system)
        self.assertEqual((result.status, result.checked_subsets), ("found", 3))
        self.assertEqual(result.coefficients, (Fraction(1), Fraction(1)))
        integer = clear_rational_dual(system, result.coefficients)
        catalyst = assemble_catalyst(integer)
        self.assertTrue(catalyst.verify())
        self.assertEqual((integer.m, integer.weights), (1, (1, 1)))
        self.assertEqual((catalyst.D.shape, catalyst.S.shape), ((8,)*3, (1,)*3))
        self.assertEqual((catalyst.source.shape, catalyst.target.shape), ((17,)*3, (13,)*3))
        self.assertEqual(catalyst.source.restrict(*catalyst.maps), catalyst.target)
        self.assertEqual(catalyst.source.restrict(*catalyst.strip_gain()), catalyst.core_target)

    def test_denominator_clearing_keeps_gain_copies_in_characteristic_two(self):
        system = toy_system(10)
        integer = clear_rational_dual(system, (Fraction(1, 2), Fraction(1, 2)))
        self.assertEqual((integer.m, integer.weights), (2, (1, 1)))
        self.assertEqual(F2.embed(integer.m), 0)
        catalyst = assemble_catalyst(integer)
        self.assertTrue(catalyst.verify())
        for i in (8, 9):
            self.assertEqual(catalyst.target.coefficient(i, i, i), 1)
        # A discrepancy by two units vanishes modulo 2, but is NOT a dual.
        wrong_gain = replace(integer, m=4)
        self.assertFalse(wrong_gain.verify_balance())
        with self.assertRaisesRegex(ValueError, "over the integers"):
            assemble_catalyst(wrong_gain)

    def test_ragged_multiple_detector_blocks_and_product_distribution(self):
        unit = FiniteTensor.unit(F2)
        two_units = FiniteTensor.direct_sum((unit, unit))
        matrix = naive_scheme(F2, 2, 2, 2).as_tensor_scheme()
        product = matrix.tensor_product(TensorScheme.from_tensor(two_units))
        # 8 copies of X -> M2 tensor X, ordered by matrix term then X coords.
        maps = []
        for family, size in zip((matrix.a, matrix.b, matrix.c), matrix.shape):
            rows = []
            for out in range(size):
                for x in range(2):
                    rows.append(tuple(vector[out] if coord == x else 0
                                      for vector in family for coord in range(2)))
            maps.append(LinearMap(F2, 16, tuple(rows)))
        upper_product = WitnessedOrder(("two",)*8, ("product",), tuple(maps), "tensor_upper")
        selection = (LinearMap.selection(F2, 2, (0,)),)*3
        rows = (rank_upper_constraint("unit", "matrix", matrix),
                DetectorConstraint("unit", "matrix"),
                upper_product, DetectorConstraint("two", "product"),
                WitnessedOrder(("two",), ("unit",), selection))
        system = WitnessedStateSystem(2, 9, (("unit", unit), ("two", two_units),
                                            ("matrix", matrix.target), ("product", product.target)), rows)
        integer = IntegerStateCertificate(system, (1, 1, 1, 1, 1), 2)
        catalyst = assemble_catalyst(integer)
        self.assertTrue(catalyst.verify())
        self.assertEqual((catalyst.D.shape, catalyst.S.shape), ((26,)*3, (3,)*3))
        self.assertEqual(catalyst.target.shape, (40,)*3)
        self.assertEqual(catalyst.matrix_factor, matrix.target.tensor_product(catalyst.S))
        # Swapping two product rows on ONLY one axis breaks the exact tensor.
        old = catalyst.maps[0]
        rows = list(old.rows)
        rows[29], rows[31] = rows[31], rows[29]
        corrupted = replace(catalyst, maps=(LinearMap(F2, old.input_size, rows),)+catalyst.maps[1:])
        self.assertFalse(corrupted.verify())

    def test_exact_map_rejection_despite_correct_integer_balance(self):
        system = toy_system()
        original = system.constraints[0]
        old = original.maps[0]
        rows = list(old.rows)
        rows[0] = (0,)*old.input_size
        broken = replace(original, maps=(LinearMap(F2, old.input_size, rows),)+original.maps[1:])
        system = replace(system, constraints=(broken, system.constraints[1]))
        certificate = IntegerStateCertificate(system, (1, 1), 1)
        self.assertTrue(certificate.verify_balance())
        with self.assertRaisesRegex(ValueError, "exact coefficient witness"):
            assemble_catalyst(certificate)

    def test_detector_does_not_accept_matching_shapes_only(self):
        system = toy_system()
        tensor = system.registry["matrix"]
        false_product = replace(tensor, coefficients=(0,)*len(tensor.coefficients))
        with self.assertRaisesRegex(ValueError, "exact registered"):
            replace(system, inventory=(("unit", system.registry["unit"]),
                                       ("matrix", false_product)),
                    constraints=(system.constraints[1],)).require_valid()

    def test_finite_exhaustion_caps_and_redundant_rows(self):
        system = toy_system(5, strassen_scheme(F2).as_tensor_scheme())
        exhausted = find_nonnegative_dual(system)
        self.assertEqual((exhausted.status, exhausted.checked_subsets), ("finite_exhausted", 3))
        self.assertIn("no global impossibility", exhausted.reason)
        self.assertEqual(find_nonnegative_dual(system, max_subsets=1).status, "resource_cap")
        self.assertEqual(find_nonnegative_dual(system, max_support=1).status, "resource_cap")
        system = toy_system()
        unit = system.registry["unit"]
        zero_row = WitnessedOrder(("unit",), ("unit",), (LinearMap.identity(F2, 1),)*3)
        system = replace(system, constraints=(zero_row,)+system.constraints)
        result = find_nonnegative_dual(system)
        self.assertEqual(result.status, "found")
        self.assertEqual(result.coefficients, (Fraction(0), Fraction(1), Fraction(1)))
        self.assertTrue(assemble_catalyst(clear_rational_dual(system, result.coefficients)).verify())

    def test_lower_and_additive_directions_have_real_map_witnesses(self):
        system = toy_system()
        unit = FiniteTensor.unit(F2)
        two = FiniteTensor.direct_sum((unit, unit))
        inventory = system.inventory+(("two", two),)
        add_positive = WitnessedOrder(("two",), ("unit", "unit"), (LinearMap.identity(F2, 2),)*3, "addPos")
        add_negative = WitnessedOrder(("unit", "unit"), ("two",), (LinearMap.identity(F2, 2),)*3, "addNeg")
        lower = WitnessedOrder(("two",), (), (LinearMap(F2, 2, ()),)*3, "lower")
        system = replace(system, inventory=inventory, constraints=system.constraints+(add_positive, add_negative, lower))
        integer = IntegerStateCertificate(system, (1, 1, 1, 1, 0), 1)
        self.assertTrue(integer.verify_balance())
        catalyst = assemble_catalyst(integer)
        self.assertEqual(catalyst.D.shape, (12,)*3)
        self.assertTrue(catalyst.verify())

    def test_extension_coefficients_preserved_in_witnessed_maps(self):
        extension = FiniteField(2, (1, 1, 1))
        matrix = strassen_scheme(extension).rescale_terms(2, 1).as_tensor_scheme()
        system = toy_system(8, matrix)
        certificate = assemble_catalyst(clear_rational_dual(system, (1, 1)))
        self.assertEqual(certificate.field, extension)
        self.assertTrue(certificate.verify())
        self.assertTrue(any(value >= 2 for mapping in certificate.maps for row in mapping.rows for value in row))

    def test_resource_caps_apply_before_expanding_integer_weights(self):
        system = toy_system()
        integer = IntegerStateCertificate(system, (10**9, 10**9), 10**9)
        self.assertTrue(integer.verify_balance())
        with self.assertRaisesRegex(StateAssemblyLimit, "incomplete"):
            assemble_catalyst(integer)
        with self.assertRaises(StateAssemblyLimit):
            assemble_catalyst(IntegerStateCertificate(system, (1, 1), 1),
                              limits=StateAssemblyLimits(max_tensor_entries=100))
        with self.assertRaises(StateAssemblyLimit):
            assemble_catalyst(IntegerStateCertificate(system, (1, 1), 1),
                              limits=StateAssemblyLimits(max_map_entries=10))

    def test_bad_registry_weights_and_rational_dual_types(self):
        system = toy_system()
        for change in ({"d": True}, {"k": 0}, {"unit_key": "matrix"},
                       {"inventory": system.inventory+system.inventory[:1]},
                       {"constraints": (DetectorConstraint("missing", "matrix"),)}):
            with self.assertRaises(ValueError):
                replace(system, **change)
        for weights, m in (((True, 1), 1), ((1,), 1), ((-1, 1), 1), ((1, 1), 0)):
            with self.assertRaises(ValueError):
                IntegerStateCertificate(system, weights, m)
        for coefficients in ((1.0, 1), (-1, 1), (Fraction(1, 2), 1), ()):
            with self.assertRaises(ValueError):
                clear_rational_dual(system, coefficients)

    def test_sparse_balance_and_dense_constraint_preflight(self):
        system = toy_system()
        extras = tuple((f"unused{i}", FiniteTensor.unit(F2)) for i in range(1000))
        system = replace(system, inventory=system.inventory+extras)
        integer = IntegerStateCertificate(system, (1, 1), 1)
        with patch.object(WitnessedStateSystem, "dense_vectors", side_effect=AssertionError("dense allocation")):
            self.assertTrue(integer.verify_balance())
            self.assertEqual(clear_rational_dual(system, (1, 1)).weights, (1, 1))
            self.assertEqual(find_nonnegative_dual(system, max_subsets=0).status, "resource_cap")
        # Unused IDs do not become LP coordinates, but active rows still
        # count against the declared dense scalar cap before matrix creation.
        self.assertEqual(find_nonnegative_dual(system, limits=StateAssemblyLimits(max_constraint_entries=4)).status,
                         "found")
        with self.assertRaises(StateAssemblyLimit):
            find_nonnegative_dual(system, limits=StateAssemblyLimits(max_constraint_entries=3))


if __name__ == "__main__":
    unittest.main()
