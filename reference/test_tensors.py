"""Generic tensors, independent bilinear maps, and matrix API compatibility."""

from dataclasses import replace
from itertools import product
from random import Random
import unittest

from . import (FiniteField, FiniteTensor, Polynomial, TensorScheme, TensorDegeneration,
               descend, matrix_multiplication_tensor, naive_scheme, strassen_scheme)
from .degenerations import left_perturbation_fixture


F2 = FiniteField(2)
F3 = FiniteField(3)
F4 = FiniteField(2, (1, 1, 1))
F5 = FiniteField(5)
F9 = FiniteField(3, (1, 0, 1))


def convolution_tensor(field):
    # A non-matrix target: (dim X, dim Y, dim Z)=(2,3,4).
    return FiniteTensor.from_function(field, (2, 3, 4), lambda i, j, k: int(i+j == k))


def perturbed_tensor(scheme):
    f = scheme.field
    a = [[Polynomial(f, (0, x, row[(i+1) % len(row)])) for i, x in enumerate(row)]
         for row in scheme.a]
    constants = lambda family: [[Polynomial.constant(f, x) for x in row] for row in family]
    return TensorDegeneration(scheme.target, 1, a, constants(scheme.b), constants(scheme.c))


class TensorTests(unittest.TestCase):
    def test_arbitrary_coefficients_and_lexicographic_order(self):
        tensor = FiniteTensor(F5, (2, 2, 2), (0, 1, 2, 3, 4, 0, 1, 2))
        self.assertEqual(tensor.coefficient(0, 1, 0), 2)
        self.assertEqual(tensor.coefficient(1, 0, 0), 4)
        self.assertEqual(tensor.coefficient(1, 1, 1), 2)
        from_function = FiniteTensor.from_function(F5, (2, 2, 2),
                                                   lambda i, j, k: (4*i+2*j+k) % 5)
        self.assertEqual(tensor, from_function)

    def test_direct_convolution_is_polynomial_multiplication(self):
        target = convolution_tensor(F5)
        self.assertEqual(target.contract((1, 2), (3, 4, 1)), (3, 0, 4, 2))
        scheme = TensorScheme.from_tensor(target)
        self.assertEqual(scheme.terms, 6)
        self.assertTrue(scheme.verify())
        self.assertEqual(scheme.apply((1, 2), (3, 4, 1)), (3, 0, 4, 2))

    def test_exhaustive_non_matrix_f2_tensors_and_inputs(self):
        for coefficients in product(range(2), repeat=4):
            target = FiniteTensor(F2, (1, 2, 2), coefficients)
            scheme = TensorScheme.from_tensor(target)
            self.assertTrue(scheme.verify())
            for a, b, c in product(range(2), repeat=3):
                expected = ((coefficients[0]*a*b + coefficients[2]*a*c) % 2,
                            (coefficients[1]*a*b + coefficients[3]*a*c) % 2)
                self.assertEqual(target.contract((a,), (b, c)), expected)
                self.assertEqual(scheme.apply((a,), (b, c)), expected)

    def test_random_arbitrary_tensors_over_extensions(self):
        random = Random(401)
        for field in (F3, F4, F9):
            for shape in ((2, 3, 4), (3, 1, 2), (1, 2, 5)):
                coefficients = tuple(random.randrange(field.order)
                                     for _ in range(shape[0]*shape[1]*shape[2]))
                target = FiniteTensor(field, shape, coefficients)
                scheme = TensorScheme.from_tensor(target)
                self.assertTrue(scheme.verify())
                for _ in range(5):
                    left = tuple(random.randrange(field.order) for _ in range(shape[0]))
                    right = tuple(random.randrange(field.order) for _ in range(shape[1]))
                    self.assertEqual(scheme.apply(left, right), target.contract(left, right))

    def test_axis_pairing_for_tensor_product(self):
        first = FiniteTensor.from_function(F5, (2, 1, 3), lambda i, j, k: (i+k+1) % 5)
        second = FiniteTensor.from_function(F5, (1, 2, 2), lambda i, j, k: (2*j+k+2) % 5)
        combined = first.tensor_product(second)
        self.assertEqual(combined.shape, (2, 2, 6))
        for i, j, k in product(range(2), range(2), range(6)):
            k1, k2 = divmod(k, 2)
            expected = (first.coefficient(i, 0, k1)*second.coefficient(0, j, k2)) % 5
            self.assertEqual(combined.coefficient(i, j, k), expected)
        scheme = TensorScheme.from_tensor(first).tensor_product(TensorScheme.from_tensor(second))
        self.assertTrue(scheme.verify())
        self.assertEqual(scheme.target, combined)
        self.assertEqual(scheme.apply((1, 3), (4, 2)), combined.contract((1, 3), (4, 2)))

    def test_tensor_powers_use_independent_axis_pairs(self):
        target = convolution_tensor(F3)
        scheme = TensorScheme.from_tensor(target)
        self.assertEqual(target.tensor_power(0), FiniteTensor.unit(F3))
        self.assertEqual(scheme.tensor_power(0).apply((2,), (2,)), (1,))
        powered = scheme.tensor_power(2)
        self.assertEqual(powered.shape, (4, 9, 16))
        self.assertEqual(powered.terms, 36)
        for i, j, k in product(range(4), range(9), range(16)):
            i1, i2 = divmod(i, 2)
            j1, j2 = divmod(j, 3)
            k1, k2 = divmod(k, 4)
            self.assertEqual(powered.target.coefficient(i, j, k),
                             int(i1+j1 == k1 and i2+j2 == k2))
        self.assertTrue(powered.verify())

    def test_zero_tensor_empty_axes_and_zero_term_decompositions(self):
        for shape in ((2, 3, 4), (0, 2, 3), (2, 0, 3), (2, 3, 0), (0, 0, 0)):
            target = FiniteTensor.zero(F4, shape)
            scheme = TensorScheme.from_tensor(target)
            self.assertEqual(scheme.terms, 0)
            self.assertTrue(scheme.verify())
            left, right = (1,)*shape[0], (1,)*shape[1]
            self.assertEqual(scheme.apply(left, right), (0,)*shape[2])
            self.assertEqual(target.contract(left, right), (0,)*shape[2])
            self.assertEqual(target.tensor_product(FiniteTensor.unit(F4)), target)
            degeneration = TensorDegeneration.from_scheme(scheme, leading=3)
            self.assertTrue(degeneration.verify())
            recovered = degeneration.recover_scheme()
            self.assertTrue(recovered.verify())
            self.assertEqual(recovered.terms, 0)
            self.assertEqual(descend(scheme).target, target.project_to_prime_field())

    def test_bad_decomposition_is_rejected(self):
        scheme = TensorScheme.from_tensor(convolution_tensor(F5))
        c = [list(row) for row in scheme.c]
        c[0][0] = 2
        bad = replace(scheme, c=c)
        self.assertFalse(bad.verify())
        for action in (lambda: bad.apply((0, 0), (0, 0, 0)),
                       lambda: bad.tensor_power(0), lambda: descend(bad)):
            with self.assertRaisesRegex(ValueError, "invalid tensor coefficient"):
                action()

    def test_validation_and_copying_mutable_inputs(self):
        coefficients = [1, 2]
        target = FiniteTensor(F3, [1, 1, 2], coefficients)
        coefficients[0] = 0
        self.assertEqual(target.coefficients, (1, 2))
        for shape in ((1, 2), (-1, 1, 1), (1.0, 1, 1)):
            with self.assertRaises(ValueError):
                FiniteTensor.zero(F3, shape)
        with self.assertRaises(ValueError):
            FiniteTensor(F3, (1, 1, 2), (1,))
        with self.assertRaises(ValueError):
            FiniteTensor(F3, (1, 1, 1), (3,))
        for index in (-1, 1, 0.0):
            with self.assertRaises(ValueError):
                target.coefficient(index, 0, 0)
        with self.assertRaises(ValueError):
            target.contract((), (1,))
        with self.assertRaises(ValueError):
            target.tensor_product(FiniteTensor.unit(F4))
        with self.assertRaises(ValueError):
            target.tensor_power(-1)
        scheme = TensorScheme.from_tensor(target)
        a = [list(row) for row in scheme.a]
        copied = replace(scheme, a=a)
        self.assertTrue(copied.verify())
        a[0][0] = 0
        self.assertEqual(copied.a, scheme.a)
        with self.assertRaises(ValueError):
            replace(scheme, b=())
        with self.assertRaises(ValueError):
            replace(scheme, a=((3,),))


class GenericDegenerationTests(unittest.TestCase):
    def test_non_matrix_target_recovery_and_powers(self):
        random = Random(402)
        target = convolution_tensor(F5)
        d = perturbed_tensor(TensorScheme.from_tensor(target))
        self.assertTrue(d.verify())
        for exponent in (0, 1, 2):
            powered = d.tensor_power(exponent)
            self.assertTrue(powered.verify())
            self.assertEqual(powered.leading, exponent)
            self.assertEqual(powered.recovery_node_count, exponent+1)
            recovered = powered.recover_scheme()
            self.assertEqual(recovered.target, target.tensor_power(exponent))
            self.assertEqual(recovered.terms, (exponent+1)*6**exponent)
            self.assertTrue(recovered.verify())
            left = tuple(random.randrange(5) for _ in range(recovered.shape[0]))
            right = tuple(random.randrange(5) for _ in range(recovered.shape[1]))
            self.assertEqual(recovered.apply(left, right), recovered.target.contract(left, right))
        self.assertFalse(d.evaluate_at(0).verify())
        with self.assertRaises(ValueError):
            d.evaluate_at(0).apply((0, 0), (0, 0, 0))

    def test_arbitrary_extension_coefficients_are_preserved_by_recovery(self):
        target = FiniteTensor(F4, (1, 2, 2), (2, 3, 1, 2))
        d = TensorDegeneration.from_scheme(TensorScheme.from_tensor(target), leading=2)
        self.assertTrue(d.verify())
        recovered = d.recover_scheme()
        self.assertEqual(recovered.target, target)
        self.assertEqual(recovered.apply((3,), (1, 2)), target.contract((3,), (1, 2)))

    def test_target_not_just_shape_is_checked(self):
        target = convolution_tensor(F3)
        d = TensorDegeneration.from_scheme(TensorScheme.from_tensor(target), leading=1)
        bad = replace(d, target=FiniteTensor.zero(F3, target.shape))
        self.assertFalse(bad.verify())
        with self.assertRaises(ValueError):
            bad.recover_scheme()

    def test_formal_validation_nodes_and_fields(self):
        target = FiniteTensor.unit(F2)
        one = Polynomial.constant(F2, 1)
        bad = TensorDegeneration(target, 1, ((Polynomial(F2, (0, 0, 1)),),), ((one,),), ((one,),))
        self.assertTrue(bad.evaluate_at(1).verify())
        self.assertFalse(bad.verify())
        with self.assertRaises(ValueError):
            bad.recover_scheme()
        d = perturbed_tensor(TensorScheme.from_tensor(convolution_tensor(F3)))
        for nodes in ((), (0, 1), (1, 1), (1,)):
            with self.assertRaises(ValueError):
                d.recover_scheme(nodes)
        with self.assertRaises(ValueError):
            replace(d, leading=-1)
        with self.assertRaises(ValueError):
            replace(d, a=((one,),)*d.terms)
        with self.assertRaises(ValueError):
            d.tensor_power(-1)
        with self.assertRaises(ValueError):
            d.tensor_product(TensorDegeneration.from_scheme(TensorScheme.from_tensor(FiniteTensor.unit(F4))))


class ProjectionAndCompatibilityTests(unittest.TestCase):
    def test_generic_descent_with_one_fixed_extension(self):
        for field in (F4, F9):
            target = convolution_tensor(field)
            seed = TensorScheme.from_tensor(target).rescale_terms(field.p, 1)
            self.assertTrue(seed.verify())
            for exponent in (0, 1, 2):
                projected = descend(seed.tensor_power(exponent))
                self.assertEqual(projected.terms, 6**exponent * field.degree**2)
                self.assertTrue(projected.verify())
                self.assertEqual(projected.target, target.tensor_power(exponent).project_to_prime_field())
                left, right = (1,)*projected.shape[0], (1,)*projected.shape[1]
                self.assertEqual(projected.apply(left, right), projected.target.contract(left, right))

    def test_linear_projection_of_non_base_tensor_is_not_descent(self):
        target = FiniteTensor(F4, (1, 1, 2), (2, 3))
        scheme = TensorScheme.from_tensor(target)
        with self.assertRaisesRegex(ValueError, "non-base-field"):
            descend(scheme)
        projected = scheme.project_to_prime_field()
        self.assertEqual(projected.target.coefficients, (0, 1))
        self.assertTrue(projected.verify())
        self.assertEqual(projected.apply((1,), (1,)), (0, 1))
        # pi(alpha^2)=1 while pi(alpha)^2=0: pi is not a field homomorphism.
        alpha = TensorScheme.from_tensor(FiniteTensor(F4, (1, 1, 1), (2,)))
        self.assertEqual(alpha.tensor_power(2).project_to_prime_field().target.coefficients, (1,))
        self.assertEqual(alpha.project_to_prime_field().tensor_power(2).target.coefficients, (0,))

    def test_matrix_adapter_round_trip_and_recovery(self):
        for field in (F3, F4):
            for matrix_scheme in (strassen_scheme(field), naive_scheme(field, 2, 3, 1)):
                generic = matrix_scheme.as_tensor_scheme()
                self.assertTrue(generic.verify())
                self.assertEqual(generic.target,
                                 matrix_multiplication_tensor(field, matrix_scheme.n, matrix_scheme.m, matrix_scheme.k))
                self.assertEqual(generic.to_matrix_scheme(matrix_scheme.n, matrix_scheme.m, matrix_scheme.k), matrix_scheme)
                d = left_perturbation_fixture(matrix_scheme)
                recovered = d.as_tensor_degeneration().recover_scheme()
                self.assertEqual(recovered.to_matrix_scheme(matrix_scheme.n, matrix_scheme.m, matrix_scheme.k),
                                 d.recover_scheme())

    def test_matrix_conversion_rejects_same_shape_different_target(self):
        target = FiniteTensor.zero(F3, (4, 4, 4))
        scheme = TensorScheme.from_tensor(target)
        self.assertTrue(scheme.verify())
        with self.assertRaisesRegex(ValueError, "not the requested matrix"):
            scheme.to_matrix_scheme(2, 2, 2)
        self.assertTrue(strassen_scheme(F3).as_tensor_scheme().to_matrix_scheme(2, 2, 2).verify())


if __name__ == "__main__":
    unittest.main()
