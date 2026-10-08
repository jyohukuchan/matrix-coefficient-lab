"""Formal polynomial certificates, exact recovery, and composition controls."""

from dataclasses import replace
from itertools import product
from random import Random
import unittest

from .constructions import recover_leading_coefficient
from .degenerations import PolynomialDegeneration, left_perturbation_fixture
from .fields import FiniteField
from .polynomials import Polynomial
from .schemes import descend, naive_multiply, naive_scheme, strassen_scheme


F2 = FiniteField(2)
F3 = FiniteField(3)
F4 = FiniteField(2, (1, 1, 1))
F9 = FiniteField(3, (1, 0, 1))


def scalar_degeneration(field, leading, a, b=None, c=None):
    one = Polynomial.constant(field, 1)
    return PolynomialDegeneration(field, 1, 1, 1, leading,
                                  [[p] for p in a], [[p] for p in (b or [one]*len(a))],
                                  [[p] for p in (c or [one]*len(a))])


def random_matrix(field, n, m, random):
    return [[random.randrange(field.order) for _ in range(m)] for _ in range(n)]


class PolynomialTests(unittest.TestCase):
    def test_formal_equality_differs_from_finite_field_function_equality(self):
        t = Polynomial(F2, (0, 1))
        self.assertNotEqual(t, t**2)
        for node in range(F2.order):
            self.assertEqual(t.evaluate(node), (t**2).evaluate(node))

    def test_known_coefficients_and_canonical_zero(self):
        # (1+t)^2 has no middle term in characteristic two.
        one_plus_t = Polynomial(F4, (1, 1, 0, 0))
        self.assertEqual(one_plus_t.coefficients, (1, 1))
        self.assertEqual((one_plus_t**2).coefficients, (1, 0, 1))
        self.assertEqual((one_plus_t-one_plus_t).degree, -1)
        self.assertEqual((one_plus_t**0).coefficients, (1,))
        self.assertEqual(Polynomial.monomial(F4, 4, 0), Polynomial(F4))
        self.assertEqual(Polynomial(F4).coefficient(19), 0)
        # alpha*(alpha*t) = (alpha+1)*t.
        self.assertEqual((Polynomial.constant(F4, 2)*Polynomial(F4, (0, 2))).coefficients, (0, 3))

    def test_evaluation_homomorphism(self):
        random = Random(301)
        for field in (F2, F3, F4, F9):
            for _ in range(15):
                a, b = [Polynomial(field, tuple(random.randrange(field.order) for _ in range(5)))
                        for _ in range(2)]
                for node in range(field.order):
                    self.assertEqual((a+b).evaluate(node), field.add(a.evaluate(node), b.evaluate(node)))
                    self.assertEqual((a*b).evaluate(node), field.mul(a.evaluate(node), b.evaluate(node)))
                    self.assertEqual((a**3).evaluate(node), field.pow(a.evaluate(node), 3))

    def test_polynomial_validation(self):
        with self.assertRaises(ValueError):
            Polynomial(F2, (-1,))
        with self.assertRaises(ValueError):
            Polynomial(F2) + Polynomial(F3)
        with self.assertRaises(TypeError):
            Polynomial(F2) * 1
        for degree in (-1, 1.0):
            with self.assertRaises(ValueError):
                Polynomial.monomial(F2, degree)
            with self.assertRaises(ValueError):
                Polynomial(F2).coefficient(degree)
            with self.assertRaises(ValueError):
                Polynomial(F2)**degree


class DegenerationTests(unittest.TestCase):
    def test_low_coefficients_cancel_across_terms(self):
        # (1+t) + (-1) = t: neither term alone vanishes at zero.
        for field in (F2, F3, F4, F9):
            d = scalar_degeneration(field, 1,
                                    [Polynomial(field, (1, 1)), Polynomial.constant(field, field.neg(1))])
            self.assertTrue(d.verify())
            self.assertEqual(d.tensor_polynomial(0, 0, 0), Polynomial(field, (0, 1)))
            self.assertEqual(d.recovery_node_count, 1)
            scheme = d.recover_scheme()
            self.assertEqual(scheme.terms, 2)
            for a, b in product(range(field.order), repeat=2):
                self.assertEqual(scheme.multiply([[a]], [[b]]), [[field.mul(a, b)]])

    def test_noise_is_not_proportional_to_target_and_evaluation_is_not_exact(self):
        d = left_perturbation_fixture(strassen_scheme(F4))
        self.assertTrue(d.verify())
        off_target = d.tensor_polynomial(0, 2, 0)
        self.assertEqual(off_target.coefficient(1), 0)
        self.assertEqual(off_target.coefficient(2), 1)
        for node in range(F4.order):
            self.assertFalse(d.evaluate_at(node).verify())
        with self.assertRaises(ValueError):
            d.evaluate_at(2).multiply([[1, 0], [0, 1]], [[1, 0], [0, 1]])
        self.assertTrue(d.recover_scheme().verify())

    def test_recovered_schemes_and_matrix_products(self):
        random = Random(302)
        for field in (F3, F4, F9):
            d = left_perturbation_fixture(strassen_scheme(field))
            scheme = d.recover_scheme()
            self.assertEqual((d.degree_bound, d.leading, d.recovery_node_count), (2, 1, 2))
            self.assertEqual(scheme.terms, 14)
            self.assertTrue(scheme.verify())
            for _ in range(10):
                a, b = [random_matrix(field, 2, 2, random) for _ in range(2)]
                self.assertEqual(scheme.multiply(a, b), naive_multiply(field, a, b))
            # Cross-check each tensor polynomial against scalar coefficient recovery.
            for i, j, h in product(range(4), repeat=3):
                recovered = recover_leading_coefficient(field, d.tensor_polynomial(i, j, h).coefficients,
                                                       d.leading, (1, 2))
                self.assertEqual(recovered, d.tensor_polynomial(i, j, h).coefficient(d.leading))

    def test_all_three_coefficient_families_can_be_polynomial(self):
        d = scalar_degeneration(F9, 1, [Polynomial(F9, (0, 1, 1))],
                                [Polynomial(F9, (1, 3))], [Polynomial(F9, (1, 4))])
        self.assertTrue(d.verify())
        self.assertEqual((d.degree_bound, d.recovery_node_count), (4, 4))
        scheme = d.recover_scheme()
        self.assertEqual(scheme.terms, 4)
        for a, b in product(range(F9.order), repeat=2):
            self.assertEqual(scheme.multiply([[a]], [[b]]), [[F9.mul(a, b)]])

    def test_power_then_recovery_has_linear_degree_overhead(self):
        random = Random(303)
        for field in (F4, F9):
            d = left_perturbation_fixture(strassen_scheme(field))
            for exponent in (0, 1, 2):
                powered = d.tensor_power(exponent)
                self.assertTrue(powered.verify())
                self.assertEqual(powered.leading, exponent)
                self.assertEqual(powered.degree_bound, 2*exponent)
                self.assertEqual(powered.recovery_node_count, exponent+1)
                recovered = powered.recover_scheme()
                self.assertEqual(recovered.terms, (exponent+1)*7**exponent)
                self.assertTrue(recovered.verify())
                a, b = [random_matrix(field, recovered.n, recovered.n, random) for _ in range(2)]
                self.assertEqual(recovered.multiply(a, b), naive_multiply(field, a, b))
            repeated = d.recover_scheme().tensor_power(2)
            self.assertEqual(repeated.terms, 196)
            self.assertTrue(repeated.verify())
            a, b = [random_matrix(field, 4, 4, random) for _ in range(2)]
            self.assertEqual(repeated.multiply(a, b), d.tensor_power(2).recover_scheme().multiply(a, b))
            self.assertEqual(d.recover_scheme().multiply_recursive(a, b), naive_multiply(field, a, b))

    def test_recovery_then_fixed_extension_descent(self):
        random = Random(304)
        for field in (F4, F9):
            seed = strassen_scheme(field).rescale_terms(field.p, 1)
            d = left_perturbation_fixture(seed)
            exact = d.tensor_power(2).recover_scheme()
            base = descend(exact)
            self.assertEqual(base.terms, 147*field.degree**2)
            self.assertTrue(base.verify())
            for _ in range(5):
                a, b = [random_matrix(base.field, 4, 4, random) for _ in range(2)]
                self.assertEqual(base.multiply(a, b), naive_multiply(base.field, a, b))

    def test_rectangular_polynomial_tensor_product(self):
        first = left_perturbation_fixture(naive_scheme(F4, 2, 3, 1))
        second = left_perturbation_fixture(naive_scheme(F4, 1, 2, 3))
        powered = first.tensor_product(second)
        self.assertEqual((powered.n, powered.m, powered.k, powered.leading), (2, 6, 3, 2))
        recovered = powered.recover_scheme()
        self.assertEqual(recovered.terms, 108)
        self.assertTrue(recovered.verify())
        random = Random(305)
        a, b = random_matrix(F4, 2, 6, random), random_matrix(F4, 6, 3, random)
        self.assertEqual(recovered.multiply(a, b), naive_multiply(F4, a, b))

    def test_monomial_lift_and_exhaustive_f2_matrix_inputs(self):
        d = PolynomialDegeneration.from_scheme(strassen_scheme(F2), leading=5)
        self.assertTrue(d.verify())
        self.assertEqual((d.degree_bound, d.recovery_node_count), (5, 1))
        scheme = d.recover_scheme()
        self.assertEqual(scheme.terms, 7)
        matrices = [[[a, b], [c, e]] for a, b, c, e in product(range(2), repeat=4)]
        for a, b in product(matrices, repeat=2):
            self.assertEqual(scheme.multiply(a, b), naive_multiply(F2, a, b))

    def test_custom_nodes_extra_nodes_and_permutation(self):
        d = left_perturbation_fixture(strassen_scheme(F4))
        for nodes in ((2, 1), (1, 2, 3), (3, 1, 2), iter((2, 3))):
            scheme = d.recover_scheme(nodes)
            self.assertTrue(scheme.verify())
            self.assertIn(scheme.terms, (14, 21))
        # Encoded node 2 is alpha, not the integer scalar 2=0 in this field.
        self.assertEqual(F4.embed(2), 0)
        self.assertTrue(d.recover_scheme((1, 2)).verify())

    def test_bad_formal_certificates_rejected_despite_matching_evaluations(self):
        # t^2 and t induce the same function on F_2, but [t]t^2=0.
        bad_leading = scalar_degeneration(F2, 1, [Polynomial(F2, (0, 0, 1))])
        # 1+t+t^2 evaluates to 1 at the only nonzero F_2 node, yet [1]P != 0.
        bad_lower = scalar_degeneration(F2, 1, [Polynomial(F2, (1, 1, 1))])
        for bad in (bad_leading, bad_lower):
            self.assertTrue(bad.evaluate_at(1).verify())
            self.assertFalse(bad.verify())
            with self.assertRaisesRegex(ValueError, "invalid degeneration"):
                bad.recover_scheme()
            with self.assertRaises(ValueError):
                bad.tensor_power(0)
        self.assertTrue(PolynomialDegeneration.from_scheme(naive_scheme(F2, 1, 1, 1)).verify())

    def test_corrupted_matrix_certificate_rejected(self):
        d = left_perturbation_fixture(strassen_scheme(F4))
        changed = [list(row) for row in d.c]
        changed[0][1] = changed[0][1] + Polynomial.constant(F4, 1)
        bad = replace(d, c=changed)
        self.assertFalse(bad.verify())
        with self.assertRaises(ValueError):
            bad.recover_scheme()

    def test_not_enough_nodes_and_invalid_nodes_rejected(self):
        d = left_perturbation_fixture(strassen_scheme(F4))
        for nodes in ((), (0, 1), (1, 1), (1,), (1, 4)):
            with self.assertRaises(ValueError):
                d.recover_scheme(nodes)
        with self.assertRaises(ValueError):
            left_perturbation_fixture(strassen_scheme(F2)).recover_scheme()
        # A fixed finite field cannot supply arbitrarily many interpolation nodes.
        scalar = scalar_degeneration(F4, 1, [Polynomial(F4, (0, 1, 1))])
        with self.assertRaises(ValueError):
            scalar.tensor_power(3).recover_scheme()

    def test_shapes_fields_indices_and_exponents_checked(self):
        d = left_perturbation_fixture(strassen_scheme(F4))
        with self.assertRaises(ValueError):
            replace(d, leading=-1)
        with self.assertRaises(ValueError):
            replace(d, a=((Polynomial(F4),),)*7)
        with self.assertRaises(ValueError):
            replace(d, b=())
        changed = [list(row) for row in d.c]
        changed[0][0] = Polynomial(F3)
        with self.assertRaises(ValueError):
            replace(d, c=changed)
        changed[0][0] = (1,)
        with self.assertRaises(ValueError):
            replace(d, c=changed)
        with self.assertRaises(ValueError):
            d.tensor_product(left_perturbation_fixture(strassen_scheme(F3)))
        with self.assertRaises(ValueError):
            d.tensor_polynomial(4, 0, 0)
        for exponent in (-1, 1.0):
            with self.assertRaises(ValueError):
                d.tensor_power(exponent)
        with self.assertRaises(ValueError):
            PolynomialDegeneration.from_scheme(strassen_scheme(F4), leading=-1)

    def test_nested_input_lists_do_not_mutate_certificates(self):
        source = left_perturbation_fixture(strassen_scheme(F4))
        a = [list(row) for row in source.a]
        d = replace(source, a=a)
        self.assertTrue(d.verify())
        a[0][0] = Polynomial(F4)
        self.assertEqual(d.a, source.a)
        self.assertTrue(d.recover_scheme().verify())


if __name__ == "__main__":
    unittest.main()
