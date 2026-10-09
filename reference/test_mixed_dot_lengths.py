"""Exact, capped length controls; finite witnesses do not prove exclusions."""

from dataclasses import replace
from itertools import product
import unittest
from unittest.mock import patch

from .fields import FiniteField
from .maps import LinearMap
from .mixed_dot_catalyst_screen import mixed_dot_core_input, mixed_dot_tensor
from .mixed_dot_lengths import (core_slice_control, exact_length_slice_rank,
                                export_mixed_dot_length_square_subrank_order,
                                export_mixed_dot_length_subrank_order,
                                mixed_dot_length_coefficient_control,
                                mixed_dot_length_core_input, mixed_dot_length_rank_profiles,
                                mixed_dot_length_tensor, source_slice_control)
from .tensors import FiniteTensor, matrix_multiplication_tensor
from .witnessed_constraints import StateAssemblyLimit, StateAssemblyLimits


FIELDS = (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)))


class MixedDotLengthTests(unittest.TestCase):
    def test_lengths_and_three_unit_selectors_all_coefficients(self):
        for field in FIELDS:
            for t in (1, 2, 3, 4):
                with self.subTest(field=field, length=t):
                    tensor = mixed_dot_length_tensor(field, t)
                    self.assertEqual(tensor.shape, (1+2*t,)*3)
                    self.assertEqual(sum(bool(x) for x in tensor.coefficients), 3*t)
                    exported = export_mixed_dot_length_subrank_order(field, t, "mixed")
                    exported.require_valid()
                    self.assertEqual(exported.row.negative, ("unit",)*3)
                    selected = tensor.restrict(*exported.row.maps)
                    self.assertEqual(selected, FiniteTensor.from_function(field, (3,)*3,
                                                                         lambda i, j, k: int(i == j == k)))
                    self.assertEqual(sum(bool(x) for x in selected.coefficients), 3)

    def test_length_two_matches_existing_constructor_and_encoder(self):
        for field in FIELDS:
            self.assertEqual(mixed_dot_length_tensor(field, 2), mixed_dot_tensor(field))
            q = field.order
            A = ((1, 0), (q-1, 1))
            B = tuple(tuple((i+j+1) % q for j in range(2)) for i in range(4))
            C = tuple(tuple((2*i+j+1) % q for j in range(4)) for i in range(2))
            self.assertEqual(mixed_dot_length_core_input(field, 2, A, B, C),
                             mixed_dot_core_input(field, A, B, C))

    def test_source_slice_zero_nonzero_blocks_include_extension_coefficients(self):
        for field in FIELDS:
            scalar = field.p if field.degree > 1 else field.order-1
            for t in (1, 2, 3, 4):
                for s, u, v in product(range(2), repeat=3):
                    vector = (scalar*s,)+(scalar*u,)*t+(scalar*v,)*t
                    self.assertEqual(source_slice_control(field, t, vector), t*s+u+v)

    def test_all_27_actual_core_rank_profiles_per_length_and_field(self):
        for field in FIELDS:
            for t in (1, 2, 3, 4):
                profiles = mixed_dot_length_rank_profiles(field, t)
                self.assertEqual(len(profiles), 27)
                self.assertEqual(tuple(ranks for ranks, _ in profiles), tuple(product(range(3), repeat=3)))
                for (a, b, c), actual in profiles:
                    self.assertEqual(actual, 2*t*a+2*b+2*c)

    def test_general_encoder_offsets_and_unbalanced_extension_matrices(self):
        field, t = FIELDS[-1], 3
        A = ((2, 3), (0, 1))
        B = ((0, 0), (2, 1), (3, 0), (0, 0), (1, 2), (0, 3))
        C = ((3, 0, 1, 0, 0, 2), (0, 2, 0, 1, 3, 0))
        vector = mixed_dot_length_core_input(field, t, A, B, C)
        for i, j in product(range(2), repeat=2):
            base = (2*i+j)*(1+2*t)
            self.assertEqual(vector[base], A[i][j])
            self.assertEqual(vector[base+1:base+1+t], tuple(B[t*i+r][j] for r in range(t)))
            self.assertEqual(vector[base+1+t:base+1+2*t], tuple(C[i][t*j+r] for r in range(t)))
        self.assertIn(2, vector)
        self.assertIn(3, vector)
        self.assertEqual(core_slice_control(field, t, A, B, C), 4*t+8)
        # This public control accepts one-pass coefficient iterables as well.
        self.assertEqual(core_slice_control(field, t, iter(A), iter(B), iter(C)), 4*t+8)

    def test_repeated_source_and_target_rank_controls_with_actual_catalyst_slices(self):
        for field in FIELDS:
            for t in (1, 2, 3, 4):
                scalar = field.p if field.degree > 1 else 1
                D = FiniteTensor.from_function(field, (2, 2, 2), lambda i, j, k: int(i == j == k))
                control = mixed_dot_length_coefficient_control(field, t, catalyst=D,
                                                               catalyst_input=(scalar, scalar))
                self.assertEqual(control.catalyst_slice_rank, 2)
                self.assertEqual(control.source_slice_ranks, (2+5*t+10, 2+3*t+10, 2+10))
                self.assertEqual(control.target_slice_rank, 2+4*t+9)
                self.assertEqual(control.source.shape, (2+5*(1+2*t),)*3)
                self.assertEqual(control.target.shape, (3+4*(1+2*t),)*3)
        control = mixed_dot_length_coefficient_control(FIELDS[-1], 4,
            catalyst=matrix_multiplication_tensor(FIELDS[-1], 2, 2, 2), catalyst_input=(2, 0, 0, 2))
        self.assertEqual(control.catalyst_slice_rank, 4)
        self.assertEqual(control.source_slice_ranks, (34, 26, 14))
        self.assertEqual(control.target_slice_rank, 29)

    def test_square_subrank_witnesses_are_actual_diagonals_not_optimality_claims(self):
        for field in FIELDS:
            for t in (1, 2, 3, 4):
                exported = export_mixed_dot_length_square_subrank_order(field, t, "square")
                self.assertEqual(len(exported.row.negative), 3+6*t)
                source = exported.registry["square"]
                self.assertEqual(source, mixed_dot_length_tensor(field, t).tensor_power(2))
                count = 3+6*t
                diagonal = source.restrict(*exported.row.maps)
                self.assertEqual(diagonal, FiniteTensor.from_function(field, (count,)*3,
                                                                      lambda i, j, k: int(i == j == k)))

    def test_corrupted_coefficients_selectors_and_recorded_ranks_are_rejected(self):
        field, t = FIELDS[-1], 2
        tensor = mixed_dot_length_tensor(field, t)
        coefficients = list(tensor.coefficients)
        coefficients[0] = 0
        damaged = replace(tensor, coefficients=tuple(coefficients))
        with self.assertRaises(ValueError):
            source_slice_control(field, t, (1,)*5, auxiliary=damaged)
        A, B, C = ((1, 0), (0, 1)), ((1, 0), (0, 1), (0, 0), (0, 0)), ((1, 0, 0, 0), (0, 1, 0, 0))
        with self.assertRaises(ValueError):
            core_slice_control(field, t, A, B, C, auxiliary=damaged)
        exported = export_mixed_dot_length_subrank_order(field, t, "mixed")
        bad_first = LinearMap.selection(field, 5, (0, 1, 4))
        self.assertFalse(replace(exported, row=replace(exported.row,
            maps=(bad_first,)+exported.row.maps[1:])).verify())
        control = mixed_dot_length_coefficient_control(field, t)
        with self.assertRaises(ValueError):
            replace(control, source_slice_ranks=(20, 16, 0)).require_valid()
        with self.assertRaises(ValueError):
            replace(control, target=FiniteTensor.zero(field, control.target.shape)).require_valid()

    def test_resource_caps_precede_tensor_and_map_allocations(self):
        field = FIELDS[0]
        with patch("reference.mixed_dot_lengths.FiniteTensor.from_function", side_effect=AssertionError("allocation")):
            with self.assertRaises(StateAssemblyLimit):
                mixed_dot_length_tensor(field, 10**9)
            with self.assertRaises(StateAssemblyLimit):
                mixed_dot_length_tensor(field, 2, limits=StateAssemblyLimits(max_blocks=2))
        operations = (
            (lambda: export_mixed_dot_length_subrank_order(field, 2, "s", limits=StateAssemblyLimits(max_map_entries=1))),
            (lambda: export_mixed_dot_length_square_subrank_order(field, 4, "s", limits=StateAssemblyLimits(max_tensor_entries=100))),
            (lambda: mixed_dot_length_rank_profiles(field, 2, limits=StateAssemblyLimits(max_constraint_entries=26))),
            (lambda: mixed_dot_length_rank_profiles(field, 2, limits=StateAssemblyLimits(max_map_entries=1))),
            (lambda: mixed_dot_length_coefficient_control(field, 4, limits=StateAssemblyLimits(max_map_entries=1))),
        )
        with patch("reference.mixed_dot_lengths.mixed_dot_length_tensor", side_effect=AssertionError("allocation")):
            for operation in operations:
                with self.assertRaisesRegex(StateAssemblyLimit, "incomplete"):
                    operation()
        with self.assertRaises(StateAssemblyLimit):
            mixed_dot_length_core_input(field, 10**9, None, None, None)
        huge_empty = FiniteTensor.zero(field, (10**9, 0, 0))
        with self.assertRaises(StateAssemblyLimit):
            exact_length_slice_rank(huge_empty, ())

    def test_invalid_parameters_and_coefficients(self):
        for t in (0, -1, True, 1.0):
            with self.assertRaises(ValueError):
                mixed_dot_length_tensor(FIELDS[0], t)
        with self.assertRaises(ValueError):
            mixed_dot_length_tensor("F2", 2)
        with self.assertRaises(ValueError):
            mixed_dot_length_tensor(FIELDS[0], 2, limits=None)
        for names in (("", "unit"), ("unit", "unit")):
            with self.assertRaises(ValueError):
                export_mixed_dot_length_subrank_order(FIELDS[0], 2, names[0], unit_key=names[1])
        with self.assertRaises(ValueError):
            source_slice_control(FIELDS[0], 2, (1,))
        with self.assertRaises(ValueError):
            source_slice_control(FIELDS[0], 2, (True, 0, 0, 0, 0))
        with self.assertRaises(ValueError):
            mixed_dot_length_core_input(FIELDS[0], 2, ((1,),), ((0, 0),)*4, ((0,)*4,)*2)
        with self.assertRaises(ValueError):
            mixed_dot_length_coefficient_control(FIELDS[0], 2, catalyst=FiniteTensor.unit(FIELDS[1]))


if __name__ == "__main__":
    unittest.main()
