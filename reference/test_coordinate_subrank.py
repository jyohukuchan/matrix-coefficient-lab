"""Independent finite matching checks and actual nonmultiplicative subrank cuts."""

from itertools import permutations, product
import unittest

from .coordinate_subrank import CoordinateSearchBudget, find_coordinate_subrank
from .fields import FiniteField
from .sectors import ThreeSectorConstruction, convolution_tensor
from .subrank_constraints import export_projective_convolution_subrank_order
from .tensors import FiniteTensor
from .witnessed_constraints import StateAssemblyLimits


F2 = FiniteField(2)
F4 = FiniteField(2, (1, 1, 1))


def brute_coordinate_diagonal(tensor, n):
    # Enumerate local coordinate substitutions independently, not support subsets.
    for xs, ys, zs in product(*(permutations(range(d), n) for d in tensor.shape)):
        if all(bool(tensor.coefficient(xs[i], ys[j], zs[k])) == (i == j == k)
               for i, j, k in product(range(n), repeat=3)):
            return True
    return False


class CoordinateSubrankTests(unittest.TestCase):
    def test_complete_binary_cube_family_matches_independent_map_enumeration(self):
        for mask in range(256):
            tensor = FiniteTensor(F2, (2, 2, 2), tuple((mask >> i) & 1 for i in range(8)))
            expected = brute_coordinate_diagonal(tensor, 2)
            result = find_coordinate_subrank(tensor, 'cube', 2)
            self.assertEqual(result.status == 'found', expected, mask)
            self.assertIn(result.status, ('found', 'finite_exhausted'))
            if result.export is not None:
                result.export.require_valid()

    def test_rectangular_supports_and_extension_weights_have_actual_normalization(self):
        for seed in range(32):
            coefficients = tuple((i*13+seed*7+i*i+seed*seed) % 4
                                 if (i+seed) % 3 == 0 else 0 for i in range(18))
            tensor = FiniteTensor(F4, (2, 3, 3), coefficients)
            result = find_coordinate_subrank(tensor, 'rectangular', 2)
            self.assertEqual(result.status == 'found', brute_coordinate_diagonal(tensor, 2))
            if result.export is not None:
                result.export.require_valid()
        tensor = FiniteTensor.from_function(F4, (2, 2, 2),
            lambda i, j, k: (2 if i == 0 else 3) if i == j == k else 0)
        result = find_coordinate_subrank(tensor, 'weighted', 2)
        result.export.require_valid()
        self.assertEqual(result.export.row.maps[2].rows,
                         ((F4.inv(2), 0), (0, F4.inv(3))))

    def test_retained_three_region_gap_has_three_checked_units(self):
        tensor = ThreeSectorConstruction(F2, 3, 1).retained
        result = find_coordinate_subrank(tensor, 'retained', 3)
        self.assertEqual(result.status, 'found')
        self.assertEqual(result.export.row.negative, ('unit',)*3)
        result.export.require_valid()

    def test_product_can_exceed_product_of_individual_subrank_bounds(self):
        left = convolution_tensor(F2, 2, 1)
        right = convolution_tensor(F2, 1, 2)
        for tensor in (left, right):
            self.assertEqual(find_coordinate_subrank(tensor, 'factor', 2).status,
                             'finite_exhausted')
        combined = left.tensor_product(right)
        self.assertEqual(combined.shape, (2, 2, 4))
        result = find_coordinate_subrank(combined, 'Fourier-component', 2)
        self.assertEqual(result.status, 'found')
        result.export.require_valid()

    def test_coordinate_exhaustion_is_weaker_than_general_linear_subrank(self):
        tensor = convolution_tensor(F2, 3, 3)
        self.assertFalse(brute_coordinate_diagonal(tensor, 3))
        self.assertEqual(find_coordinate_subrank(tensor, 'C33', 3).status, 'finite_exhausted')
        export_projective_convolution_subrank_order(F2, 3, 3, 3, 'C33').require_valid()

    def test_caps_never_claim_complete_coordinate_exhaustion(self):
        tensor = convolution_tensor(F2, 2, 2)
        for budget in (CoordinateSearchBudget(max_nodes=0),
                       CoordinateSearchBudget(max_support=1),
                       CoordinateSearchBudget(max_diagonal=1)):
            self.assertEqual(find_coordinate_subrank(tensor, 'C22', 2, budget=budget).status,
                             'resource_cap')
        self.assertEqual(find_coordinate_subrank(tensor, 'C22', 2,
            limits=StateAssemblyLimits(max_map_entries=1)).status, 'resource_cap')
        empty = FiniteTensor.zero(F2, (2, 2, 2))
        self.assertEqual(find_coordinate_subrank(empty, 'zero', 1).status, 'finite_exhausted')

    def test_invalid_arguments(self):
        tensor = FiniteTensor.unit(F2)
        for goal in (0, -1, True, 1.0):
            with self.assertRaises(ValueError):
                find_coordinate_subrank(tensor, 'unit', goal)
        for value in (float('nan'), float('inf'), 0, True):
            with self.assertRaises(ValueError):
                CoordinateSearchBudget(seconds=value)
        with self.assertRaises(ValueError):
            CoordinateSearchBudget(max_nodes=True)
        with self.assertRaises(ValueError):
            find_coordinate_subrank(tensor, '', 1)


if __name__ == '__main__':
    unittest.main()
