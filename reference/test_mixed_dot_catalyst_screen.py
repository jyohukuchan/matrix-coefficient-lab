"""Durable coefficient/profile controls; finite samples are not generic proofs."""

from dataclasses import replace
from itertools import product
from random import Random
import unittest
from unittest.mock import patch

from .catalyst_search_experiment import field_matrix_rank
from .fields import FiniteField
from .mixed_dot_catalyst_screen import (MixedDotScreenBudget, MixedDotScreenLimit,
                                       exact_slice_rank, mixed_dot_coefficient_control,
                                       mixed_dot_core_input, mixed_dot_tensor,
                                       screen_mixed_dot_catalyst,
                                       weighted_diagonal_dimension_maximum,
                                       weighted_linear_dimension_maximum)
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, matrix_multiplication_tensor


FIELDS = (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)))


def matrix_of_rank(rows, columns, rank):
    return tuple(tuple(int(i == j and i < rank) for j in range(columns)) for i in range(rows))


class MixedDotCoefficientTests(unittest.TestCase):
    def test_actual_auxiliary_and_all_small_field_slice_inputs(self):
        for field in FIELDS:
            with self.subTest(field=field.order):
                S = mixed_dot_tensor(field)
                self.assertEqual(S.shape, (5, 5, 5))
                scheme = TensorScheme.from_tensor(S)
                scheme.require_exact()
                self.assertEqual(scheme.terms, 6)
                for vector in product(range(field.order), repeat=5):
                    expected = 2*bool(vector[0])+bool(any(vector[1:3]))+bool(any(vector[3:5]))
                    self.assertEqual(exact_slice_rank(S, vector), expected)

    def test_core_matrix_formula_for_every_rank_profile(self):
        for field in FIELDS:
            S = mixed_dot_tensor(field)
            core = matrix_multiplication_tensor(field, 2, 2, 2).tensor_product(S)
            for a, b, c in product(range(3), repeat=3):
                with self.subTest(field=field.order, ranks=(a, b, c)):
                    vector = mixed_dot_core_input(field, matrix_of_rank(2, 2, a),
                                                   matrix_of_rank(4, 2, b), matrix_of_rank(2, 4, c))
                    self.assertEqual(exact_slice_rank(core, vector), 4*a+2*b+2*c)

    def test_genuine_extension_encoded_inputs(self):
        field, random = FIELDS[2], Random(20261009)
        core = matrix_multiplication_tensor(field, 2, 2, 2).tensor_product(mixed_dot_tensor(field))
        nonbase_seen = False
        for trial in range(12):
            matrices = tuple(tuple(tuple(random.randrange(field.order) for j in range(columns))
                                   for i in range(rows)) for rows, columns in ((2, 2), (4, 2), (2, 4)))
            vector = mixed_dot_core_input(field, *matrices)
            nonbase_seen |= any(value >= field.p for value in vector)
            ranks = tuple(field_matrix_rank(field, matrix, len(matrix[0])) for matrix in matrices)
            self.assertEqual(exact_slice_rank(core, vector), 4*ranks[0]+2*ranks[1]+2*ranks[2])
        self.assertTrue(nonbase_seen)
        self.assertEqual(field.embed(2), 0)

    def test_source_target_and_known_upper_controls_in_every_field(self):
        for field, kind in product(FIELDS, ("zero", "unit", "matrix")):
            with self.subTest(field=field.order, kind=kind):
                control = mixed_dot_coefficient_control(field, kind)
                extra = {"zero": 0, "unit": 1, "matrix": 4}[kind]
                self.assertEqual(control.source.shape, (25+extra,)*3)
                self.assertEqual(control.target.shape, (21+extra,)*3)
                self.assertEqual(control.source_slice_ranks, (20+extra, 16+extra, 10+extra))
                self.assertEqual(control.target_slice_rank, 17+extra)
                self.assertEqual(control.auxiliary_scheme.terms, 6)
                self.assertEqual(control.target_scheme.terms, {"zero": 43, "unit": 44, "matrix": 50}[kind])
                self.assertTrue(control.auxiliary_scheme.verify())
                self.assertTrue(control.target_scheme.verify())
                self.assertLessEqual(control.target_slice_rank, control.target_scheme.terms)
                corrupted = replace(control.auxiliary_scheme,
                                    a=((0,)*5,)+control.auxiliary_scheme.a[1:])
                self.assertFalse(corrupted.verify())

    def test_unit_exception_has_dimension_eighteen_and_actual_rank_twelve(self):
        for field in FIELDS:
            control = mixed_dot_coefficient_control(field, "unit")
            selected = control.unit_exception_coordinates
            self.assertEqual(len(set(selected)), 18)
            rows = tuple(tuple(int(i == j) for j in range(22)) for i in selected)
            self.assertEqual(field_matrix_rank(field, rows, 22), 18)
            self.assertEqual(control.unit_exception_vector[:2], (0, 0))
            self.assertTrue(all(not control.unit_exception_vector[i] for i in range(22) if i not in selected))
            self.assertEqual(control.unit_exception_rank, 12)
            self.assertGreater(control.unit_exception_rank, 11)


class MixedDotDimensionAndScreenTests(unittest.TestCase):
    def test_weighted_dimension_maxima_and_unique_maximizing_profiles(self):
        expected = {"zero": (10, 17, (0, 2, 2, 1), 54),
                    "unit": (11, 18, (0, 2, 2, 1, 1), 108),
                    "matrix": (14, 21, (0, 2, 2, 1, 2), 162)}
        for kind, (bound, dimension, profile, cases) in expected.items():
            result = weighted_linear_dimension_maximum(kind)
            self.assertEqual(result["rank_bound"], bound)
            self.assertEqual(result["maximum_dimension"], dimension)
            self.assertEqual(result["maximizing_profiles"], (profile,))
            self.assertEqual(result["checked_profiles"], cases)
        # Recompute independently, using explicit projection blocks and their
        # dimension/rank weights instead of relying on saved maxima.
        for kind in expected:
            blocks = (((0, 0), (2, 4), (4, 8)),
                      ((0, 0), (4, 2), (8, 4)),
                      ((0, 0), (4, 2), (8, 4)), ((0, 0), (1, 1)))
            if kind == "unit":
                blocks += (((0, 0), (1, 1)),)
            elif kind == "matrix":
                blocks += (((0, 0), (2, 2), (4, 4)),)
            bound = expected[kind][0]
            maximum = max(sum(d for d, r in profile) for profile in product(*blocks)
                          if sum(r for d, r in profile) <= bound)
            self.assertEqual(maximum, expected[kind][1])

    def test_supported_exact_cases_are_analytically_excluded(self):
        for field in FIELDS:
            for tensor, kind in ((FiniteTensor.zero(field, (0, 0, 0)), "zero"),
                                 (FiniteTensor.unit(field), "unit"),
                                 (matrix_multiplication_tensor(field, 2, 2, 2), "matrix")):
                report = screen_mixed_dot_catalyst(tensor, m=3)
                self.assertEqual(report["status"], "analytically_excluded")
                self.assertEqual(report["catalyst_kind"], kind)
                self.assertEqual(report["gain"], 3)
                self.assertFalse(report["finite_sampling_proves_all_fields"])
                self.assertIn("analytic argument", report["verification"])
                if kind == "unit":
                    self.assertEqual(report["unit_exception"], {"dimension": 18,
                                     "witnessed_slice_rank": 12, "source_slice_bound": 11})
                else:
                    self.assertGreater(report["common_kernel_dimension_lower"],
                                       report["dimension_control"]["maximum_dimension"])

    def test_arbitrary_catalysts_are_excluded_but_unsupported_auxiliaries_are_not(self):
        field = FIELDS[0]
        diagonal = FiniteTensor.from_function(field, (2,)*3, lambda i, j, k: int(i == j == k))
        corrupt_diagonal = FiniteTensor(field, diagonal.shape, (0,)+diagonal.coefficients[1:])
        report = screen_mixed_dot_catalyst(corrupt_diagonal)
        self.assertEqual(report["status"], "analytically_excluded")
        self.assertEqual(report["catalyst_kind"], "arbitrary_finite")
        S = mixed_dot_tensor(field)
        damaged = FiniteTensor(field, S.shape, (0,)+S.coefficients[1:])
        self.assertEqual(screen_mixed_dot_catalyst(FiniteTensor.unit(field), auxiliary=damaged)["status"], "unclassified")
        # A same-shaped corrupt matrix is not treated as M2 merely by shape.
        matrix = matrix_multiplication_tensor(field, 2, 2, 2)
        corrupt_matrix = FiniteTensor(field, matrix.shape, (0,)+matrix.coefficients[1:])
        report = screen_mixed_dot_catalyst(corrupt_matrix, m=3)
        self.assertEqual(report["status"], "analytically_excluded")
        self.assertEqual(report["catalyst_kind"], "arbitrary_finite")
        control = report['all_catalyst_rank_control']
        self.assertFalse(control['iteration_materialized'])
        self.assertGreater(control['target_rank_lower'], control['source_rank_upper'])
        self.assertTrue(report['arbitrary_finite_catalyst'])

    def test_same_orientation_dot_families_and_scope_controls(self):
        for field, n, axes in product(FIELDS, (1, 2, 3, 5),
                                      ((0, 1, 2), (1, 0, 2), (1, 2, 0))):
            with self.subTest(field=field.order, n=n, axes=axes):
                dot = FiniteTensor.from_function(field, (1, 2, 2), lambda i, j, k: int(j == k))
                tensor = FiniteTensor.direct_sum((dot,)*n).permute_axes(axes)
                report = screen_mixed_dot_catalyst(tensor)
                self.assertEqual(report["status"], "analytically_excluded")
                self.assertEqual(report["dot_copies"], n)
                self.assertEqual(report["singleton_leg"], axes.index(0))
                self.assertEqual(report["coefficient_controls"]["source_slice_ranks_full_two_zero_all_zero"],
                                 (2*n+20, 2*n+16, 2*n+10))
                self.assertEqual(report["coefficient_controls"]["target_full_rank_slice"], 2*n+17)
                corrupted = FiniteTensor(field, tensor.shape, (0,)+tensor.coefficients[1:])
                result = screen_mixed_dot_catalyst(corrupted)
                self.assertEqual(result["status"], "analytically_excluded")
                self.assertEqual(result["catalyst_kind"], "arbitrary_finite")
        mixed = FiniteTensor.direct_sum((dot, dot.permute_axes((1, 0, 2))))
        self.assertEqual(screen_mixed_dot_catalyst(mixed)["status"], "analytically_excluded")
        self.assertEqual(screen_mixed_dot_catalyst(dot,
            budget=MixedDotScreenBudget(max_slice_cases=0))["status"], "resource_cap")

    def test_mixed_orientation_all_count_exclusion_and_coefficient_recognition(self):
        for field, counts in product(FIELDS, ((1, 1, 0), (2, 3, 4), (3, 0, 4), (1, 5, 3), (3, 3, 3))):
            with self.subTest(field=field.order, counts=counts):
                dot = FiniteTensor.from_function(field, (1, 2, 2), lambda i, j, k: int(j == k))
                oriented = (dot, dot.permute_axes((1, 0, 2)), dot.permute_axes((1, 2, 0)))
                D = FiniteTensor.direct_sum(tuple(tensor for tensor, n in zip(oriented, counts) for _ in range(n)))
                report = screen_mixed_dot_catalyst(D)
                self.assertEqual(report['dot_counts'], counts)
                self.assertEqual(report['status'], 'analytically_excluded')
                axis = counts.index(min(counts))
                rho = sum(counts)+counts[axis]
                self.assertEqual(report['coefficient_controls']['source_slice_ranks_full_two_zero_all_zero'],
                                 (rho+20, rho+16, rho+10))
                self.assertEqual(report['coefficient_controls']['target_full_rank_slice'], rho+17)
                self.assertEqual(report['rank_obstruction']['source_rank'], 2*sum(counts)+30)
                self.assertEqual(report['rank_obstruction']['target_rank_lower'], 2*sum(counts)+31)
                corrupt = FiniteTensor(field, D.shape, (0,)+D.coefficients[1:])
                result = screen_mixed_dot_catalyst(corrupt)
                self.assertEqual(result['status'], 'analytically_excluded')
                self.assertEqual(result['catalyst_kind'], 'arbitrary_finite')

    def test_exact_diagonal_catalysts_are_excluded_with_actual_n2_n3_controls(self):
        for field, n in product(FIELDS, range(4)):
            with self.subTest(field=field.order, n=n):
                diagonal = FiniteTensor.from_function(field, (n,)*3, lambda i, j, k: int(i == j == k))
                result = screen_mixed_dot_catalyst(diagonal)
                self.assertEqual(result["status"], "analytically_excluded")
                if n >= 2:
                    self.assertEqual(result["diagonal_size"], n)
                    self.assertIn("analytic all-n", result["diagonal_argument"]["scope"])
                    self.assertEqual(result["coefficient_controls"]["source_slice_ranks_full_two_zero_all_zero"],
                                     (n+20, n+16, n+10))
                    self.assertEqual(result["coefficient_controls"]["target_full_rank_slice"], n+17)
                    self.assertEqual(result["coefficient_controls"]["target_upper_terms"], n+43)
                    control = mixed_dot_coefficient_control(field, "diagonal", diagonal_size=n)
                    self.assertEqual(control.catalyst, diagonal)
                    self.assertTrue(control.target_scheme.verify())
        n = 8
        diagonal = FiniteTensor.from_function(FIELDS[0], (n,)*3, lambda i, j, k: int(i == j == k))
        self.assertEqual(screen_mixed_dot_catalyst(diagonal)["status"], "analytically_excluded")

    def test_diagonal_profile_maximum_is_n_plus_seventeen_without_exponential_enumeration(self):
        for n in (*range(8), 10**9):
            result = weighted_diagonal_dimension_maximum(n)
            self.assertEqual(result["rank_bound"], n+10)
            self.assertEqual(result["maximum_dimension"], n+17)
            self.assertEqual(result["maximizing_profiles"], ((0, 2, 2, n+1),))
            self.assertEqual(result["checked_profiles"], 27)
        with self.assertRaises(MixedDotScreenLimit):
            weighted_diagonal_dimension_maximum(3, budget=MixedDotScreenBudget(max_profile_cases=26))

    def test_diagonal_allocation_caps_and_arbitrary_coefficient_scope(self):
        with patch("reference.mixed_dot_catalyst_screen.mixed_dot_tensor", side_effect=AssertionError("tensor allocation")):
            with self.assertRaises(MixedDotScreenLimit):
                mixed_dot_coefficient_control(FIELDS[0], "diagonal", diagonal_size=10**9)
        # Existing large input is checked before any canonical tensor is built.
        diagonal = FiniteTensor.from_function(FIELDS[0], (50,)*3, lambda i, j, k: int(i == j == k))
        with patch("reference.mixed_dot_catalyst_screen.mixed_dot_tensor", side_effect=AssertionError("tensor allocation")):
            self.assertEqual(screen_mixed_dot_catalyst(diagonal)["status"], "resource_cap")
        for field in FIELDS:
            diagonal = FiniteTensor.from_function(field, (3,)*3, lambda i, j, k: int(i == j == k))
            corrupted = FiniteTensor(field, diagonal.shape, (0,)+diagonal.coefficients[1:])
            self.assertEqual(screen_mixed_dot_catalyst(corrupted)["status"], "analytically_excluded")

    def test_caps_are_honest_and_precede_dense_control_allocations(self):
        caps = (MixedDotScreenBudget(max_tensor_entries=1),
                MixedDotScreenBudget(max_matrix_entries=1),
                MixedDotScreenBudget(max_slice_cases=0))
        for budget in caps:
            with self.subTest(budget=budget), \
                 patch("reference.mixed_dot_catalyst_screen.mixed_dot_tensor", side_effect=AssertionError("tensor allocation")), \
                 patch("reference.mixed_dot_catalyst_screen.strassen_scheme", side_effect=AssertionError("scheme allocation")), \
                 self.assertRaises(MixedDotScreenLimit):
                mixed_dot_coefficient_control(FIELDS[0], "unit", budget=budget)
        for budget in (*caps, MixedDotScreenBudget(max_profile_cases=0)):
            report = screen_mixed_dot_catalyst(FiniteTensor.unit(FIELDS[0]), budget=budget)
            self.assertEqual(report["status"], "resource_cap")
        with self.assertRaises(MixedDotScreenLimit):
            weighted_linear_dimension_maximum("zero", budget=MixedDotScreenBudget(max_profile_cases=53))
        huge_empty = FiniteTensor.zero(FIELDS[0], (10**9, 0, 0))
        self.assertEqual(screen_mixed_dot_catalyst(huge_empty)["status"], "resource_cap")

    def test_parameter_validation(self):
        for value in (True, -1, 1.5):
            with self.assertRaises(ValueError):
                MixedDotScreenBudget(max_profile_cases=value)
        with self.assertRaises(ValueError):
            screen_mixed_dot_catalyst(FiniteTensor.unit(FIELDS[0]), m=0)
        with self.assertRaises(ValueError):
            screen_mixed_dot_catalyst(FiniteTensor.unit(FIELDS[0]), auxiliary=FiniteTensor.unit(FIELDS[1]))
        with self.assertRaises(ValueError):
            mixed_dot_core_input(FIELDS[0], ((1,),), ((0, 0),)*4, ((0,)*4,)*2)
        with self.assertRaises(ValueError):
            mixed_dot_core_input(FIELDS[0], ((True, 0), (0, 1)), ((0, 0),)*4, ((0,)*4,)*2)
        with self.assertRaises(ValueError):
            weighted_linear_dimension_maximum("other")


if __name__ == "__main__":
    unittest.main()
