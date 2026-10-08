"""Coefficient identities and bounded behavior of the seeded flip experiment."""

from dataclasses import replace
from itertools import product
import unittest
from unittest.mock import patch

from .fields import FiniteField
from .flip_search import FlipBudget, TermAlgebra, search_flips
from .schemes import naive_multiply, naive_scheme
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor


FIELDS = (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)))
SHAPE = (2, 3, 2)
TERMS = (((1, 0), (1, 0, 1), (0, 1)),
         ((1, 1), (0, 1, 1), (1, 1)),
         ((0, 1), (1, 1, 0), (1, 0)))


def coefficient_tensor(field, terms, shape=SHAPE):
    """Independent direct sum of triple products at every coordinate."""
    return FiniteTensor.from_function(field, shape, lambda i, j, k:
        field.sum(field.mul(field.mul(a[i], b[j]), c[k]) for a, b, c in terms))


def term_scheme(field, terms, shape=SHAPE):
    target = coefficient_tensor(field, terms, shape)
    return TensorScheme(target, *(tuple(term[axis] for term in terms) for axis in range(3)))


class FlipIdentityTests(unittest.TestCase):
    def test_plus_identity_all_axes_fields_and_both_term_orders(self):
        for field, axis, first, second in product(FIELDS, range(3), (0, 1), (0, 1)):
            if first == second:
                continue
            algebra = TermAlgebra(field)
            original = coefficient_tensor(field, TERMS)
            result = algebra.plus(TERMS, first, second, axis)
            self.assertEqual(coefficient_tensor(field, result), original)
            self.assertLessEqual(len(result), len(TERMS)+1)
            self.assertTrue(term_scheme(field, result).verify())

    def test_flip_identity_all_axes_scalars_and_proportional_common_factors(self):
        for field, axis in product(FIELDS, range(3)):
            algebra, scale = TermAlgebra(field), min(2, field.order-1)
            first, second = list(TERMS[0]), list(TERMS[1])
            common = tuple(1 if i == 0 else scale if i == SHAPE[axis]-1 else 0 for i in range(SHAPE[axis]))
            first[axis], second[axis] = tuple(field.mul(scale, x) for x in common), common
            original_terms = (tuple(first), tuple(second), TERMS[2])
            expected = coefficient_tensor(field, original_terms)
            for scalar, left, right in product(range(1, field.order), (0, 1), (0, 1)):
                if left == right:
                    continue
                result = algebra.flip(original_terms, axis, left, right, scalar)
                self.assertEqual(coefficient_tensor(field, result), expected)
                self.assertLessEqual(len(result), len(original_terms))

    def test_canonical_absorbs_genuine_extension_scalars_into_third_leg(self):
        field = FIELDS[-1]
        algebra = TermAlgebra(field)
        a, b, c = TERMS[0]
        scaled = (tuple(field.mul(2, x) for x in a), tuple(field.mul(3, x) for x in b), c)
        canonical = algebra.canonical((scaled,))
        self.assertEqual(canonical[0][:2], (a, b))
        self.assertEqual(canonical[0][2], tuple(field.mul(field.mul(2, 3), x) for x in c))
        self.assertEqual(coefficient_tensor(field, canonical), coefficient_tensor(field, (scaled,)))
        self.assertEqual(algebra.canonical(canonical), canonical)

    def test_reduce_zero_terms_and_cancellation_over_every_field(self):
        for field in FIELDS:
            algebra = TermAlgebra(field)
            a, b, c = TERMS[0]
            cancelled = ((a, b, c), (a, b, tuple(field.neg(x) for x in c)))
            zeros = (((0, 0), b, c), (a, (0, 0, 0), c), (a, b, (0, 0)))
            self.assertEqual(algebra.reduce(cancelled+zeros), ())
            original = TERMS+cancelled+zeros
            reduced = algebra.reduce(original)
            self.assertEqual(coefficient_tensor(field, reduced), coefficient_tensor(field, TERMS))
            self.assertEqual(algebra.reduce(reduced), reduced)

    def test_reduce_merges_proportional_pairs_on_each_pair_of_legs(self):
        for field, remaining in product(FIELDS, range(3)):
            algebra, scalar = TermAlgebra(field), min(2, field.order-1)
            a, b, c = TERMS[0]
            first, second = list((a, b, c)), list((a, b, c))
            scaled_axis = next(axis for axis in range(3) if axis != remaining)
            second[scaled_axis] = tuple(field.mul(scalar, x) for x in second[scaled_axis])
            second[remaining] = tuple(int(i == 0) for i in range(SHAPE[remaining]))
            original = (tuple(first), tuple(second))
            reduced = algebra.reduce(original)
            self.assertLessEqual(len(reduced), 1)
            self.assertEqual(coefficient_tensor(field, reduced), coefficient_tensor(field, original))

    def test_invalid_move_premises_and_ranges(self):
        algebra = TermAlgebra(FIELDS[0])
        for call in (lambda: algebra.flip(TERMS, -1, 0, 1),
                     lambda: algebra.flip(TERMS, 3, 0, 1),
                     lambda: algebra.flip(TERMS, 0, 0, 0),
                     lambda: algebra.flip(TERMS, 0, 0, 1, 0),
                     lambda: algebra.flip(TERMS, 0, 0, 1, 2),
                     lambda: algebra.flip(TERMS, 0, 0, 1),  # no common proportional factor
                     lambda: algebra.plus(TERMS, 0, 0),
                     lambda: algebra.plus(TERMS, 0, 1, 3)):
            with self.assertRaises(ValueError):
                call()
        with self.assertRaises(ValueError):
            TermAlgebra(FiniteField(17))

    def test_negative_out_of_range_and_noninteger_move_indices_are_rejected(self):
        algebra = TermAlgebra(FIELDS[0])
        for first, second in ((-1, 0), (0, -1), (3, 0), (0, 3), (True, 0), (0, 1.0)):
            for call in (lambda: algebra.flip(TERMS, 0, first, second),
                         lambda: algebra.plus(TERMS, first, second, 0)):
                with self.assertRaises(ValueError):
                    call()
        for axis in (True, 1.0, "0"):
            with self.assertRaises(ValueError):
                algebra.flip(TERMS, axis, 0, 1)
            with self.assertRaises(ValueError):
                algebra.plus(TERMS, 0, 1, axis)
        for scalar in (True, 1.0, -1):
            with self.assertRaises(ValueError):
                algebra.flip(TERMS, 0, 0, 1, scalar)


class FlipSearchTests(unittest.TestCase):
    def test_seeded_naive_matrix_eight_to_seven_with_small_budget(self):
        field = FIELDS[0]
        source = naive_scheme(field, 2, 2, 2).as_tensor_scheme()
        budget = FlipBudget(max_steps=100, timeout_seconds=10, probes_per_step=12)
        result = search_flips(source, goal_terms=7, budget=budget, seed=0)
        self.assertEqual(result.status, "found")
        self.assertEqual((result.initial_terms, result.best_terms), (8, 7))
        self.assertLessEqual(result.steps, 100)
        self.assertGreater(result.plus_moves, 0)
        self.assertTrue(result.scheme.verify())
        self.assertEqual(result.scheme.target, source.target)
        matrix = result.scheme.to_matrix_scheme(2, 2, 2)
        for seed in range(5):
            left = [[(i+j+seed) % 2 for j in range(2)] for i in range(2)]
            right = [[(i+2*j+seed+1) % 2 for j in range(2)] for i in range(2)]
            self.assertEqual(matrix.multiply(left, right), naive_multiply(field, left, right))
        replay = search_flips(source, goal_terms=7, budget=budget, seed=0)
        self.assertEqual(replay.scheme, result.scheme)
        self.assertEqual(replay.improvements, result.improvements)

    def test_step_and_time_caps_return_exact_best_without_impossibility_claim(self):
        source = naive_scheme(FIELDS[0], 2, 2, 2).as_tensor_scheme()
        result = search_flips(source, goal_terms=6, budget=FlipBudget(max_steps=0), seed=0)
        self.assertEqual(result.status, "resource_cap")
        self.assertEqual((result.steps, result.proposals), (0, 0))
        self.assertIn("search incomplete", result.reason)
        self.assertTrue(result.scheme.verify())
        with patch("reference.flip_search.monotonic", side_effect=(0.0, 1.0, 1.1)):
            timed = search_flips(source, goal_terms=6, budget=FlipBudget(timeout_seconds=.1), seed=0)
        self.assertEqual(timed.status, "resource_cap")
        self.assertIn("time cap", timed.reason)
        self.assertIn("search incomplete", timed.reason)
        self.assertTrue(timed.scheme.verify())

    def test_input_construction_caps_and_goal_already_met(self):
        source = naive_scheme(FIELDS[0], 2, 2, 2).as_tensor_scheme()
        for budget in (FlipBudget(max_seed_terms=7), FlipBudget(max_tensor_entries=63)):
            result = search_flips(source, goal_terms=7, budget=budget)
            self.assertEqual(result.status, "resource_cap")
            self.assertIn("budget exceeded", result.reason)
            self.assertIsNone(result.scheme)
        already = search_flips(source, goal_terms=8, budget=FlipBudget(max_steps=0))
        self.assertEqual(already.status, "found")
        self.assertEqual(already.steps, 0)
        self.assertTrue(already.scheme.verify())

    def test_initial_exact_reduction_counts_as_partial_improvement(self):
        unit = FiniteTensor.unit(FIELDS[0])
        # Three unit terms sum to one over F2; exact reduction needs no moves.
        source = TensorScheme(unit, ((1,),)*3, ((1,),)*3, ((1,),)*3)
        result = search_flips(source, goal_terms=0, budget=FlipBudget(max_steps=0))
        self.assertEqual(result.status, "partial_improvement")
        self.assertEqual((result.initial_terms, result.best_terms), (3, 1))
        self.assertEqual(result.improvements, ((0, 1),))
        self.assertIn("search incomplete", result.reason)
        self.assertTrue(result.scheme.verify())

    def test_invalid_budgets_search_inputs_and_corrupt_certificates(self):
        for change in ({"max_steps": -1}, {"max_extra_terms": True}, {"recent_states": 1.0},
                       {"max_tensor_entries": -1}, {"max_seed_terms": -1},
                       {"probes_per_step": 0}, {"probes_per_step": True},
                       {"timeout_seconds": 0}, {"timeout_seconds": True},
                       {"timeout_seconds": float("nan")}, {"timeout_seconds": float("inf")}):
            with self.assertRaises(ValueError):
                FlipBudget(**change)
        source = naive_scheme(FIELDS[0], 2, 2, 2).as_tensor_scheme()
        bad = replace(source, c=((0,)*4,)*8)
        for options in ({"goal_terms": -1}, {"goal_terms": True}, {"goal_terms": 7, "seed": False},
                        {"goal_terms": 7, "budget": "invalid"}):
            with self.assertRaises(ValueError):
                search_flips(source, **options)
        with self.assertRaisesRegex(ValueError, "tensor coefficient"):
            search_flips(bad, goal_terms=7)
        with self.assertRaises(ValueError):
            search_flips("not a scheme", goal_terms=7)


if __name__ == "__main__":
    unittest.main()
