"""Formal star-pencil identities and base-field coefficient extraction."""

from dataclasses import replace
from itertools import product
import unittest
from unittest.mock import patch

from .catalyst_compiler import CompilerLimits, CompilerResourceLimit
from .discovery_inventory import star_tensor
from .fields import FiniteField
from .polynomials import Polynomial
from .star_degeneration import (extract_star_power, require_star_identity,
                                star_error_tensor,
                                star_polynomial_degeneration)
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor


FIELDS = (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)),
          FiniteField(3, (1, 0, 1)))


class StarDegenerationTests(unittest.TestCase):
    def test_three_term_formal_identity_all_coordinates_and_characteristics(self):
        for field in FIELDS:
            with self.subTest(field=field):
                degeneration = star_polynomial_degeneration(field)
                self.assertEqual(degeneration.target, star_tensor(field))
                self.assertEqual(degeneration.terms, 3)
                self.assertEqual(degeneration.degree_bound, 2)
                self.assertTrue(degeneration.verify())
                require_star_identity(degeneration)
                self.assertEqual(degeneration.a[2], (Polynomial.constant(field, field.neg(1)),)*2)
                self.assertTrue(all(p.degree <= 0 for row in degeneration.a for p in row))
                error = star_error_tensor(field)
                for coordinate in product(range(2), range(3), range(3)):
                    self.assertEqual(degeneration.tensor_polynomial(*coordinate),
                                     Polynomial(field, (0, degeneration.target.coefficient(*coordinate),
                                                        error.coefficient(*coordinate))))

    def test_wrong_constant_cancellation_and_higher_noise_are_detected(self):
        for field in FIELDS:
            degeneration = star_polynomial_degeneration(field)
            a = degeneration.a[:2]+((Polynomial(field),)*2,)
            without_cancellation = replace(degeneration, a=a)
            self.assertFalse(without_cancellation.verify())
            with self.assertRaises(ValueError):
                require_star_identity(without_cancellation)
            extra = Polynomial.monomial(field, 3)
            c = degeneration.c[:2]+((degeneration.c[2][0]+extra,)+degeneration.c[2][1:],)
            wrong_noise = replace(degeneration, c=c)
            self.assertTrue(wrong_noise.verify())  # generic leading equality still holds
            with self.assertRaisesRegex(ValueError, "formal identity"):
                require_star_identity(wrong_noise)

    def test_finite_field_sampling_cannot_substitute_for_formal_verification(self):
        field = FiniteField(2)
        degeneration = star_polynomial_degeneration(field)
        vanishing_function = Polynomial(field, (0, 1, 1))  # t^2-t on F2
        c = degeneration.c[:2]+((degeneration.c[2][0]+vanishing_function,)+degeneration.c[2][1:],)
        forged = replace(degeneration, c=c)
        for node in range(field.order):
            self.assertEqual(forged.evaluate_at(node), degeneration.evaluate_at(node))
        self.assertFalse(forged.verify())
        with self.assertRaises(ValueError):
            require_star_identity(forged)

    def test_power_one_and_two_extract_exact_schemes_without_field_extension(self):
        for field in FIELDS:
            for exponent in (1, 2):
                with self.subTest(field=field, exponent=exponent):
                    result = extract_star_power(field, exponent)
                    self.assertEqual(result.powered.terms, 3**exponent)
                    self.assertEqual(result.powered.leading, exponent)
                    self.assertEqual(result.scheme.terms, (exponent+1)*3**exponent)
                    self.assertEqual(result.term_bound, result.scheme.terms)
                    self.assertEqual(result.scheme.field, field)
                    self.assertEqual(result.scheme.target, star_tensor(field).tensor_power(exponent))
                    result.powered.require_valid()
                    result.scheme.require_exact()
                    left = tuple(i % field.order for i in range(result.scheme.shape[0]))
                    right = tuple((2*i+1) % field.order for i in range(result.scheme.shape[1]))
                    self.assertEqual(result.scheme.apply(left, right), result.scheme.target.contract(left, right))
        # F2 has only one nonzero interpolation node. Direct expansion still
        # produces both six- and twenty-seven-term exact certificates there.
        self.assertEqual(extract_star_power(FiniteField(2), 2).scheme.terms, 27)

    def test_full_power_two_polynomial_including_noise_and_lower_coefficients(self):
        for field in FIELDS:
            result = extract_star_power(field, 2)
            star, error = star_tensor(field), star_error_tensor(field)
            target = star.tensor_product(star)
            cross1, cross2 = star.tensor_product(error), error.tensor_product(star)
            top = error.tensor_product(error)
            for coordinate in product(*(range(size) for size in target.shape)):
                expected = Polynomial(field, (0, 0, target.coefficient(*coordinate),
                                               field.add(cross1.coefficient(*coordinate),
                                                         cross2.coefficient(*coordinate)),
                                               top.coefficient(*coordinate)))
                self.assertEqual(result.powered.tensor_polynomial(*coordinate), expected)

    def test_omitted_or_incorrect_degree_splits_fail_exact_coefficients(self):
        result = extract_star_power(FiniteField(3), 2)
        a, b, c = result.scheme.a, result.scheme.b, result.scheme.c
        kept = tuple(i for i in range(result.scheme.terms) if i % 3 != 1)
        omitted = TensorScheme(result.scheme.target, tuple(a[i] for i in kept),
                               tuple(b[i] for i in kept), tuple(c[i] for i in kept))
        self.assertFalse(omitted.verify())
        wrong_c = tuple(c[q+2-i] for q in range(0, len(c), 3) for i in range(3))
        incorrect_split = TensorScheme(result.scheme.target, a, b, wrong_c)
        self.assertFalse(incorrect_split.verify())

    def test_zero_power_and_preallocation_resource_caps(self):
        field = FiniteField(2)
        with patch("reference.star_degeneration.star_polynomial_degeneration",
                   side_effect=AssertionError("zero power needs no star allocation")):
            unit = extract_star_power(field, 0, limits=CompilerLimits(max_tensor_coefficients=1,
                                                                    max_scheme_scalars=3, max_terms=1))
        self.assertEqual(unit.scheme.target, FiniteTensor.unit(field))
        self.assertEqual(unit.scheme.terms, 1)
        self.assertEqual(unit.term_bound, 1)
        for limits in (CompilerLimits(max_tensor_coefficients=17), CompilerLimits(max_terms=5),
                       CompilerLimits(max_scheme_scalars=47)):
            with patch("reference.star_degeneration.star_polynomial_degeneration") as constructor:
                with self.assertRaises(CompilerResourceLimit):
                    extract_star_power(field, 1, limits=limits)
                constructor.assert_not_called()
        with patch("reference.star_degeneration.star_polynomial_degeneration") as constructor:
            with self.assertRaisesRegex(CompilerResourceLimit, "incomplete"):
                extract_star_power(field, 10**9)
            constructor.assert_not_called()

    def test_invalid_inputs(self):
        for exponent in (-1, True, 1.0):
            with self.assertRaises(ValueError):
                extract_star_power(FiniteField(2), exponent)
        for function in (extract_star_power, star_polynomial_degeneration):
            with self.assertRaises(ValueError):
                function("F2")
        with self.assertRaises(ValueError):
            extract_star_power(FiniteField(2), limits=None)
        with self.assertRaises(ValueError):
            require_star_identity(None)
        wrong_target = replace(star_polynomial_degeneration(FiniteField(2)),
                               target=FiniteTensor.zero(FiniteField(2), (2, 3, 3)))
        with self.assertRaises(ValueError):
            require_star_identity(wrong_target)


if __name__ == "__main__":
    unittest.main()
