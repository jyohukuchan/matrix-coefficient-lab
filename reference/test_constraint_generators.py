"""Actual proof-sector restrictions and context-distributed ordinary rows."""

from dataclasses import replace
import unittest
from unittest.mock import patch

from .constraint_generators import export_sector_order, lift_witnessed_order
from .fields import FiniteField
from .maps import LinearMap
from .schemes import naive_scheme
from .sectors import ThreeSectorConstruction
from .tensors import FiniteTensor
from .witnessed_constraints import (StateAssemblyLimit, StateAssemblyLimits,
                                    WitnessedStateSystem, find_nonnegative_dual,
                                    rank_upper_constraint)


F2 = FiniteField(2)
F4 = FiniteField(2, (1, 1, 1))
F16 = FiniteField(2, (1, 1, 0, 0, 1))


def altered_third(export, scalars):
    """Perturb individual node blocks, maintaining valid matrix dimensions."""
    mapping = export.row.maps[2]
    size = mapping.output_size
    rows = tuple(tuple(mapping.field.mul(value, scalars[column//size])
                       for column, value in enumerate(row)) for row in mapping.rows)
    return replace(export, row=replace(export.row, maps=export.row.maps[:2]+(LinearMap(mapping.field, mapping.input_size, rows),)))


class ConstraintGeneratorTests(unittest.TestCase):
    def test_powers_zero_one_two_actual_full_coefficients_and_atom_counts(self):
        for field in (F4, F16):
            c = ThreeSectorConstruction(field, 2, 1)
            for exponent in (0, 1, 2):
                export = export_sector_order(c, "source", "retained", exponent=exponent)
                self.assertEqual(export.nodes, tuple(range(1, exponent+2)))
                self.assertEqual(export.row.positive, ("source",)*(exponent+1))
                self.assertEqual(export.row.negative, ("retained",))
                self.assertEqual(export.registry["source"], c.source.tensor_power(exponent))
                self.assertEqual(export.registry["retained"], c.retained.tensor_power(exponent))
                source = FiniteTensor.direct_sum((export.registry["source"],)*(exponent+1))
                self.assertEqual(source.restrict(*export.row.maps), export.registry["retained"])
                self.assertTrue(export.verify())
        # Canonical node 2 is alpha in F4, despite scalar cast 2=0.
        self.assertEqual(F4.embed(2), 0)
        self.assertTrue(export_sector_order(ThreeSectorConstruction(F4, 2, 1), "s", "r").verify())

    def test_more_distinct_nodes_and_strict_degree_contract_even_without_erased_part(self):
        c = ThreeSectorConstruction(F16, 2, 1)
        export = export_sector_order(c, "s", "r", nodes=(2, 5, 9, 13))
        self.assertEqual(len(export.row.positive), 4)
        self.assertTrue(export.verify())
        no_erased = ThreeSectorConstruction(F2, 1, 1)
        self.assertFalse(any(no_erased.erased.coefficients))
        with self.assertRaises(ValueError):
            export_sector_order(no_erased, "s", "r", exponent=1)
        self.assertTrue(export_sector_order(no_erased, "s", "r", exponent=0).verify())

    def test_missing_weights_normalization_and_shift_are_detected(self):
        c = ThreeSectorConstruction(F4, 2, 1)
        export = export_sector_order(c, "s", "r", nodes=(1, 2))
        from .constructions import constant_weights
        weights = constant_weights(F4, export.nodes)
        missing_weight = altered_third(export, tuple(F4.inv(w) for w in weights))
        missing_normalization = altered_third(export, export.nodes)
        reversed_shift = altered_third(export, tuple(F4.inv(node) for node in export.nodes))
        for bad in (missing_weight, missing_normalization, reversed_shift):
            self.assertFalse(bad.verify())
            with self.assertRaisesRegex(ValueError, "coefficient equality"):
                bad.require_valid()

    def test_invalid_nodes_parameters_and_formal_construction(self):
        c = ThreeSectorConstruction(F4, 2, 1)
        for nodes in ((1,), (0, 1), (1, 1), (1, 4)):
            with patch.object(FiniteTensor, "tensor_power") as powers:
                with self.assertRaises(ValueError):
                    export_sector_order(c, "s", "r", nodes=nodes)
                powers.assert_not_called()
        for exponent in (-1, True, 1.0):
            with self.assertRaises(ValueError):
                export_sector_order(c, "s", "r", exponent=exponent)
        for a, h in ((0, 1), (2, 0)):
            with self.assertRaisesRegex(ValueError, "positive a,h"):
                export_sector_order(ThreeSectorConstruction(F4, a, h), "s", "r")
        with self.assertRaises(ValueError):
            export_sector_order(c, "same", "same")
        class WrongShift(ThreeSectorConstruction):
            @property
            def diagonal_maps(self):
                from .polynomials import Polynomial
                first, second, third = super().diagonal_maps
                return first, second, tuple(Polynomial.constant(self.field, 1) for _ in third)
        with self.assertRaisesRegex(ValueError, "three-sector"):
            export_sector_order(WrongShift(F4, 2, 1), "s", "r")

    def test_sector_resource_preflight_before_powers_and_seed_map_clones(self):
        c = ThreeSectorConstruction(F16, 2, 1)
        for arguments in ({"exponent": 10**9},
                          {"limits": StateAssemblyLimits(max_map_entries=1)},
                          {"limits": StateAssemblyLimits(max_blocks=1)},
                          {"limits": StateAssemblyLimits(max_tensor_entries=100)}):
            with patch.object(FiniteTensor, "tensor_power") as powers:
                with self.assertRaisesRegex(StateAssemblyLimit, "incomplete"):
                    export_sector_order(c, "s", "r", **arguments)
                powers.assert_not_called()

    def test_context_lift_preserves_unequal_axes_occurrences_and_metadata(self):
        original = export_sector_order(ThreeSectorConstruction(F4, 2, 1), "s", "r")
        context = FiniteTensor(F4, (1, 2, 3), (2, 0, 3, 0, 1, 0))
        lifted = lift_witnessed_order(original.row, original.inventory, context, {"s": "sx", "r": "rx"})
        self.assertEqual(lifted.row.positive, ("sx", "sx"))
        self.assertEqual(lifted.row.negative, ("rx",))
        self.assertEqual(lifted.key_pairs, (("s", "sx"), ("r", "rx")))
        self.assertEqual(lifted.registry["sx"], original.registry["s"].tensor_product(context))
        self.assertEqual(lifted.registry["rx"], original.registry["r"].tensor_product(context))
        self.assertTrue(lifted.verify())
        source = FiniteTensor.direct_sum((lifted.registry["sx"],)*2)
        self.assertEqual(source.restrict(*lifted.row.maps), lifted.registry["rx"])

    def test_context_lift_nondiagonal_upper_relation_and_empty_context(self):
        scheme = naive_scheme(F4, 2, 2, 2).as_tensor_scheme()
        inventory = (("unit", FiniteTensor.unit(F4)), ("matrix", scheme.target))
        row = rank_upper_constraint("unit", "matrix", scheme)
        context = FiniteTensor(F4, (2, 1, 2), (2, 0, 1, 3))
        lifted = lift_witnessed_order(row, inventory, context, {"unit": "x", "matrix": "mx"})
        self.assertEqual(lifted.row.positive, ("x",)*8)
        self.assertTrue(lifted.verify())
        mapping = lifted.row.maps[0]
        # Swap context-coordinate output digits; dimensions remain valid.
        rows = tuple(mapping.rows[q^1] for q in range(mapping.output_size))
        bad = replace(lifted, row=replace(lifted.row, maps=(LinearMap(F4, mapping.input_size, rows),)+lifted.row.maps[1:]))
        self.assertFalse(bad.verify())
        empty = lift_witnessed_order(row, inventory, FiniteTensor.zero(F4, (2, 0, 3)),
                                    {"unit": "z", "matrix": "mz"})
        self.assertTrue(empty.verify())

    def test_context_input_witness_validation_and_preallocation_caps(self):
        export = export_sector_order(ThreeSectorConstruction(F4, 2, 1), "s", "r")
        context = FiniteTensor.unit(F4)
        for key_map in ({"s": "sx"}, {"s": "s", "r": "rx"}, {"s": "x", "r": "x"}):
            with self.assertRaises(ValueError):
                lift_witnessed_order(export.row, export.inventory, context, key_map)
        with self.assertRaises(ValueError):
            lift_witnessed_order(export.row, export.inventory, FiniteTensor.unit(F2), {"s": "sx", "r": "rx"})
        bad = altered_third(export, (0, 0))
        with self.assertRaises(ValueError):
            lift_witnessed_order(bad.row, bad.inventory, context, {"s": "sx", "r": "rx"})
        with patch.object(FiniteTensor, "tensor_product") as products:
            with self.assertRaises(StateAssemblyLimit):
                lift_witnessed_order(export.row, export.inventory, context, {"s": "sx", "r": "rx"},
                                      limits=StateAssemblyLimits(max_map_entries=1))
            products.assert_not_called()

    def test_generated_sector_row_enters_finite_state_system_without_spurious_dual(self):
        export = export_sector_order(ThreeSectorConstruction(F4, 2, 1), "s", "r")
        inventory = (("unit", FiniteTensor.unit(F4)),)+export.inventory
        system = WitnessedStateSystem(2, 5, inventory, (export.row,))
        system.require_valid()
        result = find_nonnegative_dual(system)
        self.assertEqual(result.status, "finite_exhausted")
        self.assertIn("no global impossibility", result.reason)


if __name__ == "__main__":
    unittest.main()
