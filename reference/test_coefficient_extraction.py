"""Formal, node-free extraction with independently checked finite targets."""

from dataclasses import replace
from math import comb
import unittest
from unittest.mock import patch

from .catalyst_compiler import CompilerLimits, CompilerResourceLimit
from .coefficient_extraction import extract_leading_coefficient
from .fields import FiniteField
from .polynomials import Polynomial
from .shared_dot_orders import shared_dot_degeneration, shared_dot_tensor
from .tensor_degenerations import TensorDegeneration
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor


FIELDS = (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)))


class CoefficientExtractionTests(unittest.TestCase):
    def test_shared_cq3_powers_over_all_fields_without_nodes(self):
        for field in FIELDS:
            for exponent in (1, 2):
                with self.subTest(field=field, exponent=exponent):
                    degeneration = shared_dot_degeneration(field, 3).tensor_power(exponent)
                    with self.assertRaises(ValueError):
                        degeneration.recover_scheme()  # insufficient nonzero field elements
                    with patch("reference.tensor_degenerations.constant_weights", side_effect=AssertionError("interpolation")):
                        scheme = extract_leading_coefficient(degeneration)
                    scheme.require_exact()
                    self.assertEqual(scheme.field, field)
                    self.assertEqual(scheme.target, shared_dot_tensor(field, 3).tensor_power(exponent))
                    self.assertLessEqual(scheme.terms, comb(degeneration.leading+2, 2)*degeneration.terms)
                    self.assertEqual(scheme.terms, {1: 9, 2: 135}[exponent])
                    left = tuple((i+1) % field.order for i in range(scheme.shape[0]))
                    right = tuple((2*i+1) % field.order for i in range(scheme.shape[1]))
                    self.assertEqual(scheme.apply(left, right), scheme.target.contract(left, right))

    def test_rectangular_tensor_with_noise_on_multiple_legs_and_extension_coefficients(self):
        for field in FIELDS:
            target = FiniteTensor.from_function(field, (2, 3, 4),
                                               lambda i, j, k: (i+2*j+k+1) % field.order)
            exact = TensorScheme.from_tensor(target)
            a = tuple(tuple(Polynomial(field, (0, value, 0, 0, row[(i+1) % len(row)]))
                            for i, value in enumerate(row)) for row in exact.a)
            b = tuple(tuple(Polynomial(field, (0, value, 0, row[(i+1) % len(row)]))
                            for i, value in enumerate(row)) for row in exact.b)
            c = tuple(tuple(Polynomial.constant(field, value) for value in row) for row in exact.c)
            degeneration = TensorDegeneration(target, 2, a, b, c)
            scheme = extract_leading_coefficient(degeneration)
            self.assertEqual(scheme.shape, (2, 3, 4))
            self.assertEqual(scheme.terms, exact.terms)
            scheme.require_exact()
            if field.degree > 1:
                self.assertIn(2, target.coefficients)
                self.assertIn(3, target.coefficients)

    def test_leading_zero_cancellation_zero_vectors_and_higher_noise(self):
        for field in FIELDS:
            unit = FiniteTensor.unit(field)
            one, zero = Polynomial.constant(field, 1), Polynomial(field)
            degeneration = TensorDegeneration(unit, 2,
                ((one,), (-one,), (Polynomial.monomial(field, 2),), (one,), (Polynomial.monomial(field, 1),)),
                ((one,), (one,), (one,), (Polynomial.monomial(field, 3),), (zero,)),
                ((one,),)*5)
            scheme = extract_leading_coefficient(degeneration)
            self.assertEqual(scheme.terms, 1)
            self.assertEqual(scheme.target, unit)
            scheme.require_exact()
            constant = TensorDegeneration.from_scheme(scheme)
            self.assertEqual(extract_leading_coefficient(constant), scheme)

    def test_zero_targets_and_empty_axes_have_exact_empty_output(self):
        for shape in ((2, 3, 4), (0, 2, 3)):
            tensor = FiniteTensor.zero(FIELDS[0], shape)
            degeneration = TensorDegeneration(tensor, 0, (), (), ())
            scheme = extract_leading_coefficient(degeneration)
            self.assertEqual(scheme.terms, 0)
            self.assertEqual(scheme.target, tensor)
            scheme.require_exact()

    def test_lower_degree_noise_rejected_even_when_extracted_coefficient_is_correct(self):
        field = FIELDS[0]
        one = Polynomial.constant(field, 1)
        bad = TensorDegeneration(FiniteTensor.unit(field), 1,
                                ((Polynomial(field, (1, 1)),),), ((one,),), ((one,),))
        self.assertEqual(bad.tensor_polynomial(0, 0, 0).coefficient(1), 1)
        with self.assertRaisesRegex(ValueError, "invalid degeneration"):
            extract_leading_coefficient(bad)

    def test_sampling_indistinguishable_forgery_and_target_corruption_are_rejected(self):
        field = FIELDS[0]
        degeneration = shared_dot_degeneration(field, 3)
        a = list(degeneration.a)
        a[-1] = (a[-1][0]+Polynomial(field, (0, 1, 1)),)+a[-1][1:]
        forged = replace(degeneration, a=a)
        for node in range(field.order):
            self.assertEqual(forged.evaluate_at(node), degeneration.evaluate_at(node))
        # The leading coefficient is still C3; lower formal noise is the fault.
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    self.assertEqual(forged.tensor_polynomial(i, j, k).coefficient(3),
                                     degeneration.target.coefficient(i, j, k))
        with self.assertRaises(ValueError):
            extract_leading_coefficient(forged)
        coefficients = list(degeneration.target.coefficients)
        coefficients[0] = 1
        wrong_target = replace(degeneration, target=replace(degeneration.target, coefficients=tuple(coefficients)))
        with self.assertRaises(ValueError):
            extract_leading_coefficient(wrong_target)

    def test_actual_support_count_caps_precede_polynomial_validation_and_output(self):
        degeneration = shared_dot_degeneration(FIELDS[0], 3).tensor_power(2)
        with patch("reference.coefficient_extraction.TensorDegeneration.require_valid", side_effect=AssertionError("validation")):
            with self.assertRaises(CompilerResourceLimit):
                extract_leading_coefficient(degeneration, limits=CompilerLimits(max_terms=134))
        # The combinatorial bound is 700; the actual nonzero-support count is
        # 135, so this smaller resource cap still admits the real extraction.
        scheme = extract_leading_coefficient(degeneration, limits=CompilerLimits(max_terms=135))
        self.assertEqual(scheme.terms, 135)
        scheme.require_exact()

    def test_huge_leading_degrees_coordinates_and_formal_work_are_capped_early(self):
        field = FIELDS[0]
        huge_leading = TensorDegeneration(FiniteTensor.zero(field, (1, 1, 1)), 10**9, (), (), ())
        huge_axis = TensorDegeneration(FiniteTensor.zero(field, (10**9, 0, 0)), 0, (), (), ())
        # Small leading degree can still have costly dense high-degree noise.
        dense = Polynomial(field, (1,)*201)
        expensive = TensorDegeneration(FiniteTensor.unit(field), 0,
                                      ((dense,),), ((dense,),), ((Polynomial.constant(field, 1),),))
        holes = TensorDegeneration(FiniteTensor.unit(field), 0,
            ((Polynomial(field, (1,)*101),),), ((Polynomial(field, (1,)+(0,)*99+(1,)),),),
            ((Polynomial.constant(field, 1),),))
        with patch("reference.coefficient_extraction.TensorDegeneration.require_valid", side_effect=AssertionError("validation")):
            for degeneration in (huge_leading, huge_axis):
                with self.assertRaisesRegex(CompilerResourceLimit, "incomplete"):
                    extract_leading_coefficient(degeneration)
            with self.assertRaisesRegex(CompilerResourceLimit, "formal coefficient operation"):
                extract_leading_coefficient(expensive, limits=CompilerLimits(max_scheme_scalars=1000))
            with self.assertRaisesRegex(CompilerResourceLimit, "formal coefficient operation"):
                extract_leading_coefficient(holes, limits=CompilerLimits(max_scheme_scalars=1000))

    def test_dense_input_formal_storage_and_target_caps_precede_validation(self):
        field = FIELDS[0]
        noise = Polynomial(field, (1,)+(0,)*99+(1,))
        one = Polynomial.constant(field, 1)
        degeneration = TensorDegeneration(FiniteTensor.unit(field), 0, ((noise,),), ((one,),), ((one,),))
        shared = shared_dot_degeneration(field, 3)
        with patch("reference.coefficient_extraction.TensorDegeneration.require_valid", side_effect=AssertionError("validation")):
            with self.assertRaisesRegex(CompilerResourceLimit, "input polynomial"):
                extract_leading_coefficient(degeneration, limits=CompilerLimits(max_scheme_scalars=100))
            with self.assertRaisesRegex(CompilerResourceLimit, "formal polynomial validation"):
                extract_leading_coefficient(degeneration, limits=CompilerLimits(max_scheme_scalars=200))
            with self.assertRaises(CompilerResourceLimit):
                extract_leading_coefficient(shared,
                                           limits=CompilerLimits(max_tensor_coefficients=63))
        scheme = extract_leading_coefficient(degeneration, limits=CompilerLimits(max_scheme_scalars=500))
        self.assertEqual(scheme.terms, 1)

    def test_invalid_arguments(self):
        for value in (None, FiniteTensor.unit(FIELDS[0]), TensorScheme.from_tensor(FiniteTensor.unit(FIELDS[0]))):
            with self.assertRaises(ValueError):
                extract_leading_coefficient(value)
        with self.assertRaises(ValueError):
            extract_leading_coefficient(shared_dot_degeneration(FIELDS[0], 3), limits=None)


if __name__ == "__main__":
    unittest.main()
