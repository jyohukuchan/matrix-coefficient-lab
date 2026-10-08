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


class StarScreenTests(unittest.TestCase):
    def test_unexcluded_random_boundary_is_not_a_catalyst_and_no_generic_rank_inference(self):
        result = screen_star_catalyst(random_boundary_tensor())
        self.assertEqual(result['status'], 'not_excluded')
        self.assertEqual(result['flattening_ranks'], (4, 4, 5))
        self.assertEqual(result['base_field_slices']['histogram'], {4: 10, 3: 4, 2: 1})
        self.assertTrue(result['base_field_slices']['complete'])
        self.assertEqual(result['base_field_slices']['generic_rank'], 'not_computed')
        self.assertEqual(result['required_same_field_tensor_rank'], 6)

    def test_small_shape_and_three_dimensional_diagonal_are_excluded(self):
        diagonal = FiniteTensor.from_function(F2, (3, 6, 6),
            lambda i, j, k: int(j == k and ((1, 5, 2, 6, 3, 7)[j] >> i) & 1))
        result = screen_star_catalyst(diagonal, budget=StarScreenBudget(max_slices=0))
        self.assertEqual(result['status'], 'excluded')
        self.assertEqual(result['flattening_ranks'], (3, 6, 6))
        self.assertIn('first-concise dimension three with diagonal slice support', result['reasons'])
        small = FiniteTensor.from_function(F2, (2, 3, 3), lambda i, j, k: int(j == k))
        self.assertEqual(screen_star_catalyst(small)['status'], 'excluded')

    def test_nonconcise_first_leg_uses_basis_and_reports_a_real_low_rank_slice(self):
        # Fourth slice repeats the first. A flattening-kernel vector must not
        # be mistaken for a nonzero slice in the concise first space.
        tensor = FiniteTensor.from_function(F2, (4, 5, 5),
            lambda i, j, k: int((j == k and (i if i < 3 else 0) == j//2)
                                or ((i if i < 3 else 0) == 1 and (j, k) == (0, 1))))
        result = screen_star_catalyst(tensor)
        self.assertEqual(result['flattening_ranks'], (3, 5, 5))
        self.assertEqual(result['status'], 'excluded')
        vector = result['counterexample_input']
        self.assertEqual(vector, (0, 0, 1, 0))
        matrix = tuple(tuple(F2.sum(F2.mul(vector[i], tensor.coefficient(i, j, k)) for i in range(4))
                             for k in range(5)) for j in range(5))
        self.assertEqual(matrix[-1][-1], 1)
        self.assertEqual(sum(any(row) for row in matrix), 1)
        self.assertEqual(result['base_field_slices']['histogram'], {1: 1})
        self.assertFalse(result['base_field_slices']['complete'])

    def test_slice_cap_does_not_create_an_exclusion(self):
        result = screen_star_catalyst(random_boundary_tensor(), budget=StarScreenBudget(max_slices=0))
        self.assertEqual(result['status'], 'not_excluded')
        self.assertEqual(result['base_field_slices']['checked'], 0)
        self.assertFalse(result['base_field_slices']['complete'])
        self.assertEqual(result['base_field_slices']['maximum_observed'], 0)

    def test_rank_upper_is_checked_and_resource_caps_do_not_assert_impossibility(self):
        target = random_boundary_tensor()
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
