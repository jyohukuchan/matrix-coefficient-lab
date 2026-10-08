"""Correctness checks, including exhaustive small cases and bad certificates."""

from dataclasses import replace
from itertools import product
from random import Random
import unittest

from .constructions import (constant_weights, fourier_filter, is_primitive_root,
                            nonzero_nodes, primitive_root, recover_leading_coefficient,
                            separation_period)
from .fields import FiniteField
from .schemes import descend, naive_multiply, naive_scheme, strassen_scheme


F2 = FiniteField(2)
F3 = FiniteField(3)
F4 = FiniteField(2, (1, 1, 1))
F9 = FiniteField(3, (1, 0, 1))
F16 = FiniteField(2, (1, 1, 0, 0, 1))
F25 = FiniteField(5, (2, 0, 1))


def random_matrix(field, n, m, random):
    return [[random.randrange(field.order) for _ in range(m)] for _ in range(n)]


class FieldTests(unittest.TestCase):
    def test_known_extension_products_and_scalar_embedding(self):
        # alpha^2 = alpha+1 in F_4; beta^2 = -1 in F_9.
        self.assertEqual(F4.mul(2, 2), 3)
        self.assertEqual(F9.mul(3, 3), 2)
        self.assertEqual(F25.mul(5, 5), 3)
        self.assertEqual(F4.embed(2), 0)
        self.assertEqual(F4.coordinates(2), (0, 1))
        self.assertEqual(F4.from_coordinates((0, 1)), 2)

    def test_small_field_laws_exhaustively(self):
        for f in (F2, F3, F4, F9):
            with self.subTest(order=f.order):
                for a in range(f.order):
                    self.assertEqual(f.add(a, 0), a)
                    self.assertEqual(f.mul(a, 1), a)
                    self.assertEqual(f.add(a, f.neg(a)), 0)
                    if a:
                        self.assertEqual(f.mul(a, f.inv(a)), 1)
                    for b, c in product(range(f.order), repeat=2):
                        self.assertEqual(f.add(f.add(a, b), c), f.add(a, f.add(b, c)))
                        self.assertEqual(f.mul(f.mul(a, b), c), f.mul(a, f.mul(b, c)))
                        self.assertEqual(f.mul(a, f.add(b, c)), f.add(f.mul(a, b), f.mul(a, c)))
                        self.assertEqual(f.mul(a, b), f.mul(b, a))

    def test_larger_field_inverses(self):
        for f in (F16, F25):
            for a in range(1, f.order):
                self.assertEqual(f.mul(a, f.inv(a)), 1)
                self.assertEqual(f.pow(a, -1), f.inv(a))

    def test_invalid_field_and_element_rejected(self):
        for p in (0, 1, 4, 9, 2.0):
            with self.assertRaises(ValueError):
                FiniteField(p)
        for p, modulus in ((2, (1, 0, 1)), (3, (2, 0, 1)),
                           (2, (1, 0, 0, 0, 1)), (2, (1,)), (3, (1, 2))):
            with self.assertRaises(ValueError):
                FiniteField(p, modulus)
        with self.assertRaises(ZeroDivisionError):
            F4.inv(0)
        for a in (-1, 4, 1.0):
            with self.assertRaises(ValueError):
                F4.mul(a, 1)
        with self.assertRaises(ValueError):
            F4.from_coordinates((1,))

    def test_irreducibility_against_quadratic_root_oracle(self):
        # A monic quadratic is irreducible iff it has no root in F_p.
        for p in (2, 3, 5, 7):
            for a, b in product(range(p), repeat=2):
                has_root = any((x*x + b*x + a) % p == 0 for x in range(p))
                if has_root:
                    with self.assertRaises(ValueError):
                        FiniteField(p, (a, b, 1))
                else:
                    self.assertEqual(FiniteField(p, (a, b, 1)).order, p*p)


class SchemeTests(unittest.TestCase):
    def test_seed_exact_coefficients_in_several_characteristics(self):
        for f in (F2, F3, FiniteField(5), F4, F9):
            self.assertTrue(strassen_scheme(f).verify())
            self.assertEqual(strassen_scheme(f).terms, 7)
            self.assertTrue(naive_scheme(f, 2, 3, 2).verify())

    def test_all_two_by_two_inputs_over_f2(self):
        scheme = strassen_scheme(F2)
        matrices = [[[a, b], [c, d]] for a, b, c, d in product(range(2), repeat=4)]
        for left, right in product(matrices, repeat=2):
            expected = naive_multiply(F2, left, right)
            self.assertEqual(scheme.multiply(left, right), expected)
            self.assertEqual(scheme.multiply_recursive(left, right), expected)

    def test_random_matrices_and_recursive_blocks(self):
        random = Random(20261008)
        for f in (F2, F3, F4, F9, F25):
            seed = strassen_scheme(f)
            powered = seed.tensor_power(2)
            self.assertTrue(powered.verify())
            for n in (1, 2, 4, 8):
                for _ in range(3):
                    left, right = (random_matrix(f, n, n, random) for _ in range(2))
                    expected = naive_multiply(f, left, right)
                    self.assertEqual(seed.multiply_recursive(left, right), expected)
                    if n == 4:
                        self.assertEqual(powered.multiply(left, right), expected)

    def test_tensor_product_rectangular_coordinate_order(self):
        scheme = naive_scheme(F3, 2, 3, 1).tensor_product(naive_scheme(F3, 1, 2, 3))
        self.assertEqual((scheme.n, scheme.m, scheme.k, scheme.terms), (2, 6, 3, 36))
        self.assertTrue(scheme.verify())
        random = Random(31)
        left = random_matrix(F3, 2, 6, random)
        right = random_matrix(F3, 6, 3, random)
        self.assertEqual(scheme.multiply(left, right), naive_multiply(F3, left, right))

    def test_fixed_extension_descent_before_and_after_powers(self):
        random = Random(51)
        for extension in (F4, F9, F16):
            seed = strassen_scheme(extension).rescale_terms(extension.p, 1)
            self.assertTrue(seed.verify())
            for exponent in (0, 1, 2):
                scheme = descend(seed.tensor_power(exponent))
                self.assertEqual(scheme.terms, 7**exponent * extension.degree**2)
                self.assertTrue(scheme.verify())
                left, right = (random_matrix(scheme.field, scheme.n, scheme.n, random)
                               for _ in range(2))
                self.assertEqual(scheme.multiply(left, right), naive_multiply(scheme.field, left, right))
            repeated = descend(seed).tensor_power(2)
            self.assertEqual(repeated.terms, (7 * extension.degree**2)**2)
            self.assertTrue(repeated.verify())

    def test_prime_field_descent_and_rectangular_descent(self):
        seed = naive_scheme(F3, 2, 3, 1)
        self.assertEqual(descend(seed), seed)
        extension_seed = naive_scheme(F9, 2, 3, 1).rescale_terms(3, 4)
        projected = descend(extension_seed)
        self.assertTrue(projected.verify())
        left, right = [[1, 2, 0], [0, 1, 1]], [[2], [1], [2]]
        self.assertEqual(projected.multiply(left, right), naive_multiply(F3, left, right))

    def test_corrupted_certificate_is_rejected(self):
        seed = strassen_scheme(F3)
        changed = [list(row) for row in seed.c]
        changed[0][1] = F3.add(changed[0][1], 1)
        bad = replace(seed, c=changed)
        self.assertFalse(bad.verify())
        with self.assertRaisesRegex(ValueError, "invalid tensor coefficient"):
            bad.multiply([[0, 0], [0, 0]], [[0, 0], [0, 0]])
        with self.assertRaises(ValueError):
            descend(bad)
        with self.assertRaises(ValueError):
            bad.tensor_power(0)

    def test_shape_and_parameter_checks(self):
        seed = strassen_scheme(F2)
        with self.assertRaises(ValueError):
            replace(seed, a=((1,),) * 7)
        with self.assertRaises(ValueError):
            seed.multiply([[1]], [[1]])
        with self.assertRaises(ValueError):
            seed.multiply_recursive([[0]*3 for _ in range(3)], [[0]*3 for _ in range(3)])
        with self.assertRaises(ValueError):
            seed.tensor_power(-1)
        with self.assertRaises(ValueError):
            seed.tensor_product(strassen_scheme(F3))
        with self.assertRaises(ZeroDivisionError):
            seed.rescale_terms(0, 1)


class ConstructionTests(unittest.TestCase):
    def test_period_avoids_characteristic_and_has_bounds(self):
        for f in (F2, F3, FiniteField(5), F4, F9):
            self.assertEqual(separation_period(f, 0), 1)
            for blocks in range(1, 30):
                period = separation_period(f, blocks)
                self.assertIn(period, (5*blocks, 5*blocks+1))
                self.assertNotEqual(f.embed(period), 0)
                self.assertTrue(5*blocks <= period <= 6*blocks)
        self.assertEqual(separation_period(FiniteField(5), 1), 6)

    def test_fourier_filter_exact_including_aliases_and_negatives(self):
        for f in (F16, F25, FiniteField(11)):
            period = separation_period(f, 1)
            root = primitive_root(f, period)
            self.assertTrue(is_primitive_root(f, root, period))
            for exponent in range(-3*period, 3*period+1):
                self.assertEqual(fourier_filter(f, period, root, exponent), int(exponent % period == 0))
        self.assertEqual(fourier_filter(F2, 1, 1, -42), 1)

    def test_fourier_invalid_field_or_root_rejected(self):
        with self.assertRaises(ValueError):
            primitive_root(F4, 5)
        with self.assertRaises(ValueError):
            fourier_filter(F16, 5, 1, 0)
        with self.assertRaises(ValueError):
            fourier_filter(F4, 2, 1, 0)
        with self.assertRaises(ValueError):
            primitive_root(F4, 0)
        with self.assertRaises(ValueError):
            separation_period(F4, -1)

    def test_interpolation_random_polynomials_with_vanishing_prefix(self):
        random = Random(92)
        for f in (F2, F3, F4, F9, F16, F25):
            for count in range(1, min(f.order, 6)):
                nodes = nonzero_nodes(f, count)
                self.assertEqual(f.sum(constant_weights(f, nodes)), 1)
                for leading in range(4):
                    for _ in range(5):
                        normalized = [random.randrange(f.order) for _ in range(count)]
                        coefficients = [0]*leading + normalized + [0, 0]
                        recovered = recover_leading_coefficient(f, coefficients, leading, nodes)
                        self.assertEqual(recovered, normalized[0])

    def test_interpolation_hypotheses_enforced(self):
        for nodes in ((), (0,), (1, 1)):
            with self.assertRaises(ValueError):
                constant_weights(F4, nodes)
        with self.assertRaises(ValueError):
            nonzero_nodes(F4, 4)
        with self.assertRaises(ValueError):
            recover_leading_coefficient(F4, [1, 1], 1, (1, 2))
        with self.assertRaises(ValueError):
            recover_leading_coefficient(F4, [0, 1, 1, 1], 1, (1, 2))
        self.assertEqual(recover_leading_coefficient(F4, [0, 0], 10, (1,)), 0)


if __name__ == "__main__":
    unittest.main()
