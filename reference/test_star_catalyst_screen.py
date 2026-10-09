"""Screening uses exact necessary conditions, never finite generic-rank guesses."""

import unittest

from .fields import FiniteField
from .star_catalyst_screen import StarScreenBudget, screen_star_catalyst
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor


F2 = FiniteField(2)


def random_boundary_tensor():
    masks = (306381, 47558)
    return FiniteTensor.from_function(F2, (4, 4, 5),
        lambda i, j, k: int(k == j+i) if i < 2 else (masks[i-2] >> (5*j+k)) & 1)


def larger_unexcluded_tensor():
    # A screen control only: no catalytic maps or generic-rank claim.
    return FiniteTensor.from_function(F2, (6, 8, 8),
        lambda i, j, k: int(k == (j+i) % 8))


class StarScreenTests(unittest.TestCase):
    def test_previous_random_boundary_is_excluded_by_new_shape_refinement(self):
        result = screen_star_catalyst(random_boundary_tensor())
        self.assertEqual(result['status'], 'excluded')
        self.assertEqual(result['flattening_ranks'], (4, 4, 5))
        self.assertIn('three-plane refinement requires other flattening ranks at least six', result['reasons'])
        self.assertNotIn('base_field_slices', result)
        self.assertEqual(result['required_same_field_tensor_rank'], 7)

    def test_first_dimension_five_and_low_core_concision_are_excluded(self):
        diagonal = FiniteTensor.from_function(F2, (3, 6, 6),
            lambda i, j, k: int(j == k and ((1, 5, 2, 6, 3, 7)[j] >> i) & 1))
        result = screen_star_catalyst(diagonal, budget=StarScreenBudget(max_slices=0))
        self.assertEqual(result['status'], 'excluded')
        self.assertEqual(result['flattening_ranks'], (3, 6, 6))
        self.assertIn('slice-space refinement requires first flattening rank at least six', result['reasons'])
        small = FiniteTensor.from_function(F2, (2, 3, 3), lambda i, j, k: int(j == k))
        self.assertEqual(screen_star_catalyst(small)['status'], 'excluded')
        dimension_five = FiniteTensor.from_function(F2, (5, 10, 10),
            lambda i, j, k: int(k == (j+i) % 10))
        self.assertEqual(screen_star_catalyst(dimension_five)['flattening_ranks'], (5, 10, 10))
        self.assertEqual(screen_star_catalyst(dimension_five)['status'], 'excluded')
        six_diagonal = FiniteTensor.from_function(F2, (6, 6, 6), lambda i, j, k: int(i == j == k))
        result = screen_star_catalyst(six_diagonal)
        self.assertIn('generic slice rank at most seven requires another flattening rank at least ten', result['reasons'])

    def test_nonconcise_first_leg_uses_actual_slice_basis(self):
        # Seventh slice repeats the first. A flattening-kernel vector must not
        # be mistaken for a nonzero slice in the concise first space.
        tensor = FiniteTensor.from_function(F2, (7, 8, 8),
            lambda i, j, k: int(k == (j+(i if i < 6 else 0)) % 8))
        result = screen_star_catalyst(tensor)
        self.assertEqual(result['flattening_ranks'], (6, 8, 8))
        self.assertEqual(result['status'], 'not_excluded')
        self.assertEqual(result['base_field_slices']['total_nonzero'], 63)
        self.assertEqual(result['base_field_slices']['checked'], 63)
        self.assertNotIn(0, result['base_field_slices']['histogram'])
        self.assertTrue(result['base_field_slices']['complete'])
        self.assertEqual(result['base_field_slices']['generic_rank'], 'not_computed')

    def test_slice_cap_does_not_create_an_exclusion(self):
        result = screen_star_catalyst(larger_unexcluded_tensor(), budget=StarScreenBudget(max_slices=0))
        self.assertEqual(result['status'], 'not_excluded')
        self.assertEqual(result['base_field_slices']['checked'], 0)
        self.assertFalse(result['base_field_slices']['complete'])
        self.assertEqual(result['base_field_slices']['maximum_observed'], 0)

    def test_rank_upper_is_checked_and_resource_caps_do_not_assert_impossibility(self):
        target = larger_unexcluded_tensor()
        supplied = TensorScheme.from_tensor(target)
        result = screen_star_catalyst(target, supplied_scheme=supplied)
        self.assertEqual(result['supplied_upper_terms'], supplied.terms)
        with self.assertRaises(ValueError):
            screen_star_catalyst(target, supplied_scheme=TensorScheme.from_tensor(FiniteTensor.unit(F2)))
        broken = TensorScheme(target, (), (), ())
        with self.assertRaises(ValueError):
            screen_star_catalyst(target, supplied_scheme=broken)
        self.assertEqual(screen_star_catalyst(target, budget=StarScreenBudget(max_tensor_entries=1))['status'], 'resource_cap')
        huge_empty = FiniteTensor.zero(F2, (1_000_000, 0, 0))
        self.assertEqual(screen_star_catalyst(huge_empty)['status'], 'resource_cap')

    def test_invalid_budgets_and_empty_tensor(self):
        for cap in (True, -1, 1.5):
            with self.assertRaises(ValueError):
                StarScreenBudget(max_slices=cap)
        self.assertEqual(screen_star_catalyst(FiniteTensor.zero(F2, (0, 0, 0)))['status'], 'excluded')


if __name__ == '__main__':
    unittest.main()
