"""Exact complete-convolution controls for the audited rank obstruction."""

from itertools import product
import unittest
from unittest.mock import patch

from .fields import FiniteField
from .maps import LinearMap
from .projective_convolution import projective_convolution_scheme
from .research.convolution_commutator_experiment import (convolution_commutator_control,
                                                       mixed_output_quotient,
                                                       output_flattening_rank,
                                                       slice_map)
from .tensors import FiniteTensor
from .schemes import strassen_scheme
from .tensor_schemes import TensorScheme
from .witnessed_constraints import StateAssemblyLimit, StateAssemblyLimits


FIELDS = (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)), FiniteField(5))


class ConvolutionCommutatorTests(unittest.TestCase):
    def test_all_fields_dimensions_repetitions_and_both_input_orientations(self):
        cases = 0
        for field, a, b, q, swap in product(FIELDS, range(2, 5), range(2, 5), (1, 2), (False, True)):
            with self.subTest(field=field.order, a=a, b=b, q=q, swap=swap):
                control = convolution_commutator_control(field, a, b, q=q, swapped_inputs=swap)
                degree = a if swap else b
                L = a+b-1
                n, z = q*(1+4*degree), q*(1+4*L)
                self.assertEqual(control.output_rank, z)
                self.assertEqual(control.identity_rank, n)
                self.assertEqual(control.commutator_rank, 4*q*degree)
                self.assertEqual(control.normalized_slices[0], LinearMap.identity(field, n))
                self.assertEqual(control.lower_bound, q*(1+4*L+2*degree))
                self.assertEqual(control.supplied_scheme.terms, q*(1+8*a*b))
                self.assertLessEqual(control.lower_bound, control.supplied_scheme.terms)
                self.assertTrue(control.supplied_scheme.verify())
                outside = set(range(z))-set(control.image_coordinates)
                self.assertTrue(all(not any(mapping.rows[k])
                                    for mapping in control.slices for k in outside))
                if degree == max(a, b):
                    self.assertEqual(control.comparison_gap, q*(1+2*max(a, b)-L))
                    self.assertGreater(control.comparison_gap, 0)
                cases += 1
        self.assertEqual(cases, 144)

    def test_characteristic_two_commutator_does_not_vanish(self):
        for field in (FIELDS[0], FIELDS[2]):
            control = convolution_commutator_control(field, 2, 4, m=2, q=2)
            n = control.identity_rank
            expected_diagonal = []
            for copy in range(2):
                expected_diagonal.extend((0, 0))
                expected_diagonal.extend((1,)*(4*4))
            self.assertEqual(control.commutator.rows,
                             tuple(tuple(expected_diagonal[i] if i == j else 0 for j in range(n))
                                   for i in range(n)))
            self.assertEqual(control.commutator_rank, 32)

    def test_mixed_complement_preserves_slices_and_full_tensor_restriction(self):
        for field in FIELDS:
            with self.subTest(field=field.order):
                control = convolution_commutator_control(field, 2, 4, q=2, mixed_complement=True)
                n, z = control.identity_rank, control.output_rank
                self.assertEqual(len(control.complement), z-n)
                self.assertTrue(all(not any(control.quotient.apply(vector)) for vector in control.complement))
                selected = set(control.image_coordinates)
                discarded = tuple(k for k in range(z) if k not in selected)
                self.assertEqual(tuple(tuple(column[k] for k in discarded) for column in control.complement),
                                 tuple(tuple(int(i == j) for j in range(z-n)) for i in range(z-n)))
                # The complement is actually mixed: the quotient has nonzero
                # coefficients on coordinates outside V, not just projection.
                self.assertTrue(any(control.quotient.rows[i][k]
                                    for i in range(n) for k in discarded))
                x, y, _ = control.tensor.shape
                compressed = control.tensor.restrict(LinearMap.identity(field, x),
                                                       LinearMap.identity(field, y), control.quotient)
                for vector, expected in zip(control.slice_vectors, control.normalized_slices):
                    self.assertEqual(slice_map(compressed, vector), expected)
                if field.degree > 1:
                    self.assertTrue(any(value >= field.p for row in control.quotient.rows for value in row))

    def test_input_swap_reaches_the_stronger_bound_for_nonsquare_convolution(self):
        left = convolution_commutator_control(FIELDS[1], 4, 2)
        right = convolution_commutator_control(FIELDS[1], 4, 2, swapped_inputs=True)
        self.assertEqual(left.tensor.permute_axes((1, 0, 2)), right.tensor)
        self.assertEqual(left.lower_bound, 25)
        self.assertEqual(right.lower_bound, 29)
        self.assertEqual(max(left.lower_bound, right.lower_bound), 1+4*5+2*4)

    def test_no_scalar_gain_and_singleton_convolutions(self):
        for field, a, b, m in product(FIELDS, (1, 2), (1, 2), (0, 2)):
            control = convolution_commutator_control(field, a, b, m=m,
                                                       swapped_inputs=a > b)
            self.assertEqual(control.lower_bound, m+4*(a+b-1)+2*max(a, b))
            self.assertGreater(control.comparison_gap, 0)
            self.assertLessEqual(control.lower_bound, control.supplied_scheme.terms)

    def test_rank_computation_uses_actual_tensor_coefficients(self):
        control = convolution_commutator_control(FIELDS[2], 2, 3)
        self.assertEqual(output_flattening_rank(control.tensor), control.output_rank)
        empty = FiniteTensor.zero(FIELDS[2], control.tensor.shape)
        self.assertEqual(output_flattening_rank(empty), 0)
        self.assertEqual(slice_map(empty, control.slice_vectors[0]).rows,
                         ((0,)*empty.shape[1],)*empty.shape[2])

    def test_lower_bounds_remain_below_strassen_interpolation_upper_controls(self):
        for field in FIELDS:
            with self.subTest(field=field.order):
                control = convolution_commutator_control(field, 2, 2, m=2, q=2)
                # Projective interpolation supplies rank 3 over F2 as well;
                # no extra finite nodes or field extension are inserted.
                convolution = projective_convolution_scheme(field, 2, 2)
                matrix = strassen_scheme(field).as_tensor_scheme()
                scalar = TensorScheme.from_tensor(FiniteTensor.unit(field))
                block = TensorScheme.direct_sum((scalar, scalar, matrix.tensor_product(convolution)))
                supplied = TensorScheme.direct_sum((block, block))
                supplied.require_exact()
                self.assertEqual(supplied.target, control.tensor)
                self.assertEqual(supplied.terms, 2*(2+7*3))
                self.assertLessEqual(control.lower_bound, supplied.terms)

    def test_parameter_validation(self):
        for arguments in ((2, 2, 2), (FIELDS[0], 0, 2), (FIELDS[0], 2, True)):
            with self.assertRaises(ValueError):
                convolution_commutator_control(*arguments)
        for options in ({"m": -1}, {"m": True}, {"q": 0}, {"q": True},
                        {"swapped_inputs": 1}, {"mixed_complement": 1}, {"limits": True}):
            with self.assertRaises(ValueError):
                convolution_commutator_control(FIELDS[0], 2, 2, **options)
        for selected in ((0, 0), (3,), (True,)):
            with self.assertRaises(ValueError):
                mixed_output_quotient(FIELDS[0], selected, 3)
        with self.assertRaises(ValueError):
            mixed_output_quotient(FIELDS[0], (0,), True)

    def test_resource_caps_precede_all_tensor_and_map_allocations(self):
        for options in ({"q": 10**9}, {"m": 10**9},
                        {"limits": StateAssemblyLimits(max_tensor_entries=1)},
                        {"limits": StateAssemblyLimits(max_map_entries=1)},
                        {"limits": StateAssemblyLimits(max_constraint_entries=1)}):
            with self.subTest(options=options), \
                 patch("reference.research.convolution_commutator_experiment.convolution_tensor",
                       side_effect=AssertionError("tensor allocation")), \
                 patch.object(FiniteTensor, "unit", side_effect=AssertionError("unit allocation")), \
                 patch.object(LinearMap, "__post_init__", side_effect=AssertionError("map allocation")), \
                 self.assertRaises(StateAssemblyLimit):
                convolution_commutator_control(FIELDS[0], 2, 2, **options)
        with patch.object(LinearMap, "__post_init__", side_effect=AssertionError("map allocation")):
            with self.assertRaises(StateAssemblyLimit):
                mixed_output_quotient(FIELDS[0], (), 10**9)
        tensor = FiniteTensor.unit(FIELDS[0])
        with self.assertRaises(StateAssemblyLimit):
            output_flattening_rank(tensor, limits=StateAssemblyLimits(max_tensor_entries=0))
        with self.assertRaises(StateAssemblyLimit):
            slice_map(tensor, (1,), limits=StateAssemblyLimits(max_map_entries=0))


if __name__ == "__main__":
    unittest.main()
