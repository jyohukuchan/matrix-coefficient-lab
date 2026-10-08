"""Actual geometric positive-gain catalysts from supplied coefficient gaps."""

from dataclasses import replace
import unittest
from unittest.mock import patch

from .catalyst_compiler import CompilerLimits, CompilerResourceLimit
from .fields import FiniteField
from .geometric_catalyst import build_geometric_catalyst
from .maps import LinearMap
from .schemes import naive_scheme, strassen_scheme
from .tensors import FiniteTensor


F2 = FiniteField(2)
F4 = FiniteField(2, (1, 1, 1))


class GeometricCatalystTests(unittest.TestCase):
    def test_supplied_seven_term_b1_bridge_with_actual_gain(self):
        result = build_geometric_catalyst(2, 8, 1, strassen_scheme(F2))
        self.assertEqual(result.gain, 1)
        self.assertEqual(result.auxiliary, FiniteTensor.unit(F2))
        self.assertEqual(result.auxiliary_scheme.terms, 1)
        self.assertEqual(result.source.shape, (8, 8, 8))
        self.assertEqual(result.target.shape, (5, 5, 5))
        self.assertTrue(result.certificate.verify())
        self.assertEqual(result.source.restrict(*result.certificate.maps), result.target)

    def test_b2_matches_middle_blocks_and_reindexes_actual_matrix_products(self):
        scheme = naive_scheme(F2, 4, 4, 4)
        result = build_geometric_catalyst(2, 9, 2, scheme)
        self.assertEqual(result.gain, 17)
        self.assertEqual((result.block_counts, result.matrix_sizes, result.block_offsets), ((9, 1), (1, 2), (0, 9)))
        self.assertEqual(result.auxiliary.shape, (13, 13, 13))
        self.assertEqual(result.auxiliary_scheme.terms, 17)
        self.assertEqual(result.source.shape, (117, 117, 117))
        self.assertEqual(result.target.shape, (69, 69, 69))
        self.assertTrue(result.certificate.verify())
        self.assertEqual(result.source.restrict(*result.certificate.maps), result.target)
        # The top block contains Md tensor M2 in generic axis-pair order.
        # For q=(1,0) and s=(0,1), the row-major M4 coordinate is (2,1).
        top_generic = 17+2*13+9+1
        vector = result.certificate.maps[0].rows[top_generic]
        expected = tuple(scheme.a[t][2*4+1] for t in range(64))
        unit_positions = tuple(copy*13+j for copy in range(9) for j in range(9))
        self.assertEqual(tuple(vector[pos] for pos in unit_positions[17:]), expected)
        self.assertEqual(sum(result.target.coefficients), 17+8*(9+8))

    def test_extension_coefficients_and_even_gain_remain_actual_copies(self):
        scheme = strassen_scheme(F4).rescale_terms(2, 3)
        result = build_geometric_catalyst(2, 9, 1, scheme)
        self.assertEqual(result.gain, 2)
        self.assertEqual(F4.embed(result.gain), 0)
        self.assertEqual(result.target.coefficient(0, 0, 0), 1)
        self.assertEqual(result.target.coefficient(1, 1, 1), 1)
        self.assertTrue(result.certificate.verify())
        self.assertTrue(any(value >= 2 for mapping in result.certificate.maps for row in mapping.rows for value in row))

    def test_wrong_product_index_bijection_fails_full_coefficient_check(self):
        result = build_geometric_catalyst(2, 9, 2, naive_scheme(F2, 4, 4, 4))
        mapping = result.certificate.maps[0]
        rows = list(mapping.rows)
        # Swap two top generic product coordinates on only one leg.
        first, second = 17+9+1, 17+13+9
        rows[first], rows[second] = rows[second], rows[first]
        bad = replace(result.certificate, maps=(LinearMap(F2, mapping.input_size, rows),)+result.certificate.maps[1:])
        self.assertFalse(bad.verify())
        with self.assertRaises(ValueError):
            bad.require_valid()

    def test_missing_gap_bad_sizes_nonexact_scheme_and_invalid_parameters(self):
        for d, k, exponent, scheme in ((2, 7, 1, strassen_scheme(F2)),
                                       (2, 8, 1, naive_scheme(F2, 2, 2, 2)),
                                       (2, 9, 2, strassen_scheme(F2)),
                                       (2, 9, 1, "not a scheme")):
            with self.assertRaises(ValueError):
                build_geometric_catalyst(d, k, exponent, scheme)
        source = strassen_scheme(F2)
        bad = replace(source, c=((0, 0, 0, 0),)*source.terms)
        with self.assertRaisesRegex(ValueError, "tensor coefficient"):
            build_geometric_catalyst(2, 8, 1, bad)
        for d, k, exponent in ((1, 8, 1), (True, 8, 1), (2, 0, 1), (2, 8, 0), (2, 8, True)):
            with self.assertRaises(ValueError):
                build_geometric_catalyst(d, k, exponent, source)

    def test_huge_exponent_and_dense_resources_rejected_before_allocations(self):
        source = strassen_scheme(F2)
        for exponent, limits in ((10**9, CompilerLimits()),
                                 (1, CompilerLimits(max_map_entries=1)),
                                 (1, CompilerLimits(max_tensor_coefficients=10)),
                                 (1, CompilerLimits(max_scheme_scalars=1))):
            with patch.object(FiniteTensor, "direct_sum") as allocations:
                with self.assertRaisesRegex(CompilerResourceLimit, "incomplete"):
                    build_geometric_catalyst(2, 8, exponent, source, limits=limits)
                allocations.assert_not_called()


if __name__ == "__main__":
    unittest.main()
