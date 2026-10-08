"""Actual independent Fourier maps, square noise, and separated direct sums."""

from dataclasses import replace
from itertools import product
import unittest
from unittest.mock import patch

from . import (FiniteField, FiniteTensor, FiniteSeparation, LinearMap, Polynomial,
               TensorDegeneration, TensorScheme, ThreeSectorConstruction,
               convolution_scheme, descend)
from .separation import phase, square_weights


F11 = FiniteField(11)
F16 = FiniteField(2, (1, 1, 0, 0, 1))
F81 = FiniteField(3, (2, 1, 0, 0, 1))


def scalar_branches(field, values):
    return tuple(FiniteTensor(field, (1, 1, 1), (value,)) for value in values)


class SeparationTests(unittest.TestCase):
    def test_signed_weights_and_no_wrap_are_integer_identities(self):
        for m in range(1, 6):
            for g, h, u, v in product(range(1, m+1), repeat=4):
                e = phase(g, h, u, v)
                weight = sum(square_weights(g, h, u, v))
                self.assertLess(abs(e), 5*m)
                self.assertEqual(weight, (g-h)**2+h*e)
                if e == 0:
                    self.assertEqual(weight, (g-h)**2)
                    self.assertEqual(weight == 0, g == h and u == v)
                    self.assertLessEqual(weight, (m-1)**2)

    def test_local_fourier_maps_contract_real_copies_and_preserve_coefficients(self):
        branches = tuple(FiniteTensor.from_function(F11, (2, 2, 1),
                                                    lambda i, j, k, h=h: (i+3*j+h+1) % 11)
                         for h in range(2))
        s = FiniteSeparation(branches)
        copies = FiniteTensor.direct_sum((s.source,)*s.period)
        mapped = copies.restrict(*s.projection_maps)
        self.assertEqual(mapped, s.projected)
        scheme = s.projection_scheme()
        self.assertEqual(scheme.target, mapped)
        self.assertTrue(scheme.verify())
        for i, j, k in product(*(range(d) for d in s.shape)):
            x, g = divmod(i, 2)
            hy, u = divmod(j, 2)
            hz, v = divmod(k, 2)
            h, y = divmod(hy, 2)
            h2, z = divmod(hz, 1)
            e = (u+1)-(v+1)+2*((g+1)-(h+1))
            expected = branches[h].coefficient(x, y, z) if h == h2 and e == 0 else 0
            self.assertEqual(mapped.coefficient(i, j, k), expected)

    def test_real_square_noise_survives_above_the_leading_coefficient(self):
        s = FiniteSeparation(scalar_branches(F16, (1, 2, 7)))
        d = s.generate_degeneration()
        self.assertTrue(s.verify())
        self.assertEqual(s.shifts, (0, 6, 9))
        self.assertEqual((s.leading, d.degree_bound, s.normalized_degree_bound), (15, 25, 4))
        # g=1,h=2,u=3,v=1: phase zero, square one, extension coefficient 2.
        self.assertEqual(d.tensor_polynomial(0, 5, 3), Polynomial.monomial(F16, 16, 2))
        self.assertEqual(s.target.coefficient(0, 5, 3), 0)
        self.assertEqual(d.tensor_polynomial(1, 3, 3), Polynomial.monomial(F16, 15, 2))
        for coordinate in product(range(3), range(9), range(9)):
            self.assertEqual(d.tensor_polynomial(*coordinate), s.polynomial_coefficient(*coordinate))
        with self.assertRaisesRegex(ValueError, "formal coefficient degree"):
            d.recover_scheme((1,), degree_bound=s.leading)
        for bad in (-1, True, 19.0):
            with self.assertRaises(ValueError):
                d.recover_scheme(degree_bound=bad)
        # The supplied square bound is validated formally and saves nodes.
        self.assertEqual(d.recovery_node_count, 11)
        recovered = d.recover_scheme(degree_bound=s.leading+4)
        self.assertEqual(recovered.terms, 5*s.period*3)
        self.assertEqual(recovered.target, s.target)
        self.assertTrue(recovered.verify())

    def test_nonzero_unshifted_weighted_maps_match_polynomial_evaluation(self):
        s = FiniteSeparation(scalar_branches(F16, (1, 2, 7)))
        copies = FiniteTensor.direct_sum((s.source,)*s.period)
        for node in (1, 2, 5):
            mapped = copies.restrict(*s.weighted_maps(node))
            normalized = FiniteTensor.from_function(F16, s.shape, lambda i, j, k:
                         F16.div(s.polynomial_coefficient(i, j, k).evaluate(node), F16.pow(node, s.leading)))
            self.assertEqual(mapped, normalized)
        with self.assertRaises(ValueError):
            s.weighted_maps(0)

    def test_recovery_powers_and_full_branch_dot_product_direct_sum(self):
        s = FiniteSeparation(scalar_branches(F11, (1, 3)))
        exact = s.recover_scheme()
        direct = s.as_direct_sum(exact)
        self.assertTrue(direct.verify())
        self.assertEqual(exact.terms, 40)
        self.assertEqual(direct.target, s.direct_sum_target)
        for i, j, k in product(range(2), range(4), range(4)):
            hy, u = divmod(j, 2)
            hz, v = divmod(k, 2)
            expected = (1, 3)[i] if i == hy == hz and u == v else 0
            self.assertEqual(direct.target.coefficient(i, j, k), expected)
        self.assertEqual(exact.apply((1, 2), (3, 4, 5, 6)), s.target.contract((1, 2), (3, 4, 5, 6)))
        for exponent in (0, 2):
            recovered = s.recover_power(exponent)
            self.assertEqual(recovered.target, s.target.tensor_power(exponent))
            self.assertEqual(recovered.terms, (exponent+1)*20**exponent)
            self.assertTrue(recovered.verify())

    def test_adjusted_periods_in_characteristics_five_and_three(self):
        for field, values, period in ((FiniteField(5, (2, 0, 1)), (2,), 6),
                                     (F81, (1, 3, 7), 16)):
            s = FiniteSeparation(scalar_branches(field, values))
            self.assertEqual(s.period, period)
            self.assertEqual(field.embed(5*len(values)), 0)
            self.assertNotEqual(field.embed(period), 0)
            recovered = s.recover_scheme()
            self.assertTrue(recovered.verify())
            self.assertEqual(recovered.target, s.target)

    def test_optimized_convolution_regions_tags_separation_and_descent(self):
        c = ThreeSectorConstruction(F16, 2, 1)
        source = convolution_scheme(F16, 2, c.source_width)
        retained = c.recover_power(1, source_scheme=source)
        tagged = c.tag_scheme(retained)
        s = FiniteSeparation(c.branches)
        self.assertEqual(tagged.target, s.source)
        exact = s.recover_scheme(source_scheme=tagged)
        self.assertEqual((source.terms, retained.terms, tagged.terms, exact.terms), (5, 10, 10, 750))
        self.assertEqual(exact.shape, (6, 36, 45))
        self.assertTrue(exact.verify())
        direct = s.as_direct_sum(exact)
        self.assertEqual(direct.target, s.direct_sum_target)
        base = descend(direct)
        self.assertEqual(base.terms, 12000)
        self.assertTrue(base.verify())
        left = (1, 0, 0, 1, 1, 1)
        right = tuple(j % 2 for j in range(36))
        self.assertEqual(base.apply(left, right), base.target.contract(left, right))

    def test_invalid_periods_roots_sources_and_nodes(self):
        branches = scalar_branches(F11, (1, 3))
        for kwargs in ({"period": 5}, {"period": True}, {"period": 11}, {"root": 1}, {"root": 11}):
            with self.assertRaises(ValueError):
                FiniteSeparation(branches, **kwargs)
        with self.assertRaises(ValueError):
            FiniteSeparation(scalar_branches(FiniteField(2), (1, 1, 1)))
        for family in ((), (branches[0], FiniteTensor.unit(FiniteField(2))),
                       (branches[0], FiniteTensor.zero(F11, (2, 1, 1)))):
            with self.assertRaises(ValueError):
                FiniteSeparation(family)
        s = FiniteSeparation(branches)
        source = TensorScheme.from_tensor(s.source)
        bad = replace(source, c=((0, 0),)*source.terms)
        for scheme in (bad, TensorScheme.from_tensor(s.target), "not a scheme"):
            with self.assertRaises(ValueError):
                s.recover_power(0, source_scheme=scheme)
        for exponent, nodes in ((20, None), (1, (1,)), (1, (0, 1)), (1, (1, 1)), (1, (1, 11))):
            with patch.object(TensorScheme, "direct_sum") as copying:
                with self.assertRaises(ValueError):
                    s.recover_power(exponent, nodes)
                copying.assert_not_called()
        with self.assertRaises(ValueError):
            s.as_direct_sum(TensorScheme.from_tensor(s.source))
        for exponent in (-1, True, 1.0):
            with self.assertRaises(ValueError):
                s.recover_power(exponent)

    def test_missing_normalization_and_wrong_weights_are_detected(self):
        class MissingNormalization(FiniteSeparation):
            @property
            def projection_maps(self):
                first, second, third = super().projection_maps
                return first, second, LinearMap(self.field, third.input_size,
                                                tuple(tuple(self.field.mul(self.field.embed(self.period), x) for x in row)
                                                      for row in third.rows))

        class WrongWeight(FiniteSeparation):
            @property
            def weights(self):
                first, second, third = super().weights
                return tuple(w-1 for w in first), second, third

        branches = scalar_branches(F11, (1, 3))
        with self.assertRaisesRegex(ValueError, "local Fourier maps"):
            MissingNormalization(branches).generate_degeneration()
        self.assertFalse(WrongWeight(branches).verify())
        with self.assertRaises(ValueError):
            WrongWeight(branches).generate_degeneration()

    def test_empty_axes_and_zero_branches(self):
        for shape in ((0, 1, 1), (1, 0, 1), (1, 1, 0), (1, 1, 1)):
            s = FiniteSeparation((FiniteTensor.zero(F11, shape),)*2)
            recovered = s.recover_scheme()
            self.assertTrue(recovered.verify())
            self.assertFalse(any(recovered.target.coefficients))
            self.assertEqual(s.as_direct_sum(recovered).target, s.direct_sum_target)


if __name__ == "__main__":
    unittest.main()
