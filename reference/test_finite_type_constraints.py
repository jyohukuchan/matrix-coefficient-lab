"""Actual exact-type pulls and finite Fourier rows with all overhead retained."""

from dataclasses import replace
import unittest
from unittest.mock import patch

from .finite_type_constraints import export_exact_type_order, export_finite_separation_order
from .fields import FiniteField
from .separation import FiniteSeparation
from .tensors import FiniteTensor, shared_first_tensor
from .witnessed_constraints import StateAssemblyLimit, StateAssemblyLimits


class FiniteTypeConstraintsTests(unittest.TestCase):
    def test_exact_type_keeps_ordered_different_branch_products(self):
        f = FiniteField(3)
        branches = (FiniteTensor(f, (2, 1, 2), (1, 0, 0, 1)),
                    FiniteTensor(f, (2, 1, 2), (0, 1, 2, 0)))
        exported = export_exact_type_order(branches, (1, 1), "power", "type")
        self.assertEqual(exported.words, ((0, 1), (1, 0)))
        self.assertEqual(exported.branches, (branches[0].tensor_product(branches[1]),
                                            branches[1].tensor_product(branches[0])))
        self.assertNotEqual(*exported.branches)
        self.assertTrue(exported.verify())
        self.assertEqual(exported.restriction.registry["type"], shared_first_tensor(exported.branches))

    def test_counts_with_zero_and_empty_word(self):
        f = FiniteField(2)
        branches = (FiniteTensor.unit(f),)*3
        for counts, words in (((0, 2, 0), ((1, 1),)), ((0, 0, 0), ((),))):
            exported = export_exact_type_order(branches, counts, "source", "target")
            self.assertEqual(exported.words, words)
            self.assertTrue(exported.verify())

    def test_actual_binary_words_do_not_cancel_into_a_direct_sum(self):
        f = FiniteField(2)
        exported = export_exact_type_order((FiniteTensor.unit(f),)*2, (1, 1), "s", "t")
        self.assertEqual(exported.restriction.registry["t"].shape, (1, 2, 2))
        self.assertEqual(exported.restriction.registry["t"].coefficients, (1, 0, 0, 1))
        # Erasing the word label on just one split axis creates forbidden cross terms.
        row = exported.restriction.row
        from .maps import LinearMap
        bad = LinearMap.selection(f, row.maps[1].input_size, (1, 1))
        self.assertFalse(replace(exported.restriction, row=replace(row, maps=(row.maps[0], bad, row.maps[2]))).verify())

    def test_finite_fourier_row_retains_period_times_interpolation_nodes(self):
        f = FiniteField(11)
        branches = (FiniteTensor(f, (1, 1, 1), (1,)), FiniteTensor(f, (1, 1, 1), (3,)))
        separation = FiniteSeparation(branches)
        exported = export_finite_separation_order(separation, "shared", "separated")
        self.assertEqual((separation.period, len(exported.nodes), len(exported.row.positive)), (10, 2, 20))
        self.assertTrue(exported.verify())
        self.assertEqual(exported.registry["separated"], separation.direct_sum_target)

    def test_extension_fourier_row_preserves_a_nonbase_coefficient(self):
        f = FiniteField(2, (1, 1, 0, 0, 1))
        separation = FiniteSeparation((FiniteTensor(f, (1, 1, 1), (2,)),))
        exported = export_finite_separation_order(separation, "s", "t")
        self.assertTrue(exported.verify())
        self.assertEqual(exported.registry["t"].coefficients, (2,))

    def test_caps_are_checked_before_word_power_allocation(self):
        branches = (FiniteTensor.unit(FiniteField(2)),)*2
        with patch.object(FiniteTensor, "tensor_power", side_effect=AssertionError("must not allocate")):
            with self.assertRaises(StateAssemblyLimit):
                export_exact_type_order(branches, (10**9, 1), "s", "t")
            with self.assertRaises(StateAssemblyLimit):
                export_exact_type_order(branches, (4, 4), "s", "t", limits=StateAssemblyLimits(max_blocks=8))

    def test_missing_nodes_and_invalid_counts_are_rejected(self):
        f = FiniteField(11)
        branches = (FiniteTensor.unit(f),)*2
        with self.assertRaises(ValueError):
            export_exact_type_order(branches, (1, True), "s", "t")
        with self.assertRaises(ValueError):
            export_finite_separation_order(FiniteSeparation(branches), "s", "t", nodes=(1,))
        with self.assertRaises(ValueError):
            export_finite_separation_order(FiniteSeparation(branches), "s", "t", nodes=(1, 1))

    def test_empty_axes_do_not_bypass_the_word_position_cap(self):
        branches = (FiniteTensor.zero(FiniteField(2), (1, 0, 0)),)*2
        with patch.object(FiniteTensor, "tensor_power", side_effect=AssertionError("must not allocate")):
            with self.assertRaisesRegex(StateAssemblyLimit, "position cap"):
                export_exact_type_order(branches, (4095, 1), "s", "t")


if __name__ == "__main__":
    unittest.main()
