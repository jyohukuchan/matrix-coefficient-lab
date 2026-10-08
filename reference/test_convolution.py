"""Full coefficient interpolation and its connection to the proof's regions."""

from itertools import product
from random import Random
import unittest

from .convolution import convolution_scheme, lagrange_basis
from .fields import FiniteField
from .polynomials import Polynomial
from .schemes import descend
from .sectors import ThreeSectorConstruction


F8 = FiniteField(2, (1, 1, 0, 1))


class ConvolutionTests(unittest.TestCase):
    def test_basis_reconstructs_every_monomial_formally(self):
        for field, nodes in ((FiniteField(5), (0, 3, 1, 4)), (F8, (7, 0, 2, 4, 1))):
            basis = lagrange_basis(field, nodes)
            for i, j in product(range(len(nodes)), repeat=2):
                self.assertEqual(basis[i].evaluate(nodes[j]), int(i == j))
            for degree in range(len(nodes)):
                reconstructed = Polynomial(field)
                for node, p in zip(nodes, basis):
                    reconstructed = reconstructed + Polynomial.constant(field, field.pow(node, degree))*p
                self.assertEqual(reconstructed, Polynomial.monomial(field, degree))

    def test_convolution_in_several_characteristics_and_empty_axes(self):
        random = Random(1701)
        for field in (FiniteField(2), FiniteField(3), FiniteField(5), F8):
            for a, b in product(range(4), repeat=2):
                if a and b and a+b-1 > field.order:
                    continue
                scheme = convolution_scheme(field, a, b)
                self.assertEqual(scheme.terms, a+b-1 if a and b else 0)
                self.assertTrue(scheme.verify())
                for _ in range(3):
                    left = tuple(random.randrange(field.order) for _ in range(a))
                    right = tuple(random.randrange(field.order) for _ in range(b))
                    expected = [0]*max(0, a+b-1)
                    for i, j in product(range(a), range(b)):
                        expected[i+j] = field.add(expected[i+j], field.mul(left[i], right[j]))
                    self.assertEqual(scheme.apply(left, right), tuple(expected))

    def test_custom_extra_nodes_and_extension_encodings(self):
        scheme = convolution_scheme(F8, 2, 3, iter((0, 7, 2, 4, 1)))
        self.assertEqual(scheme.terms, 5)
        self.assertTrue(scheme.verify())
        self.assertEqual(scheme.a[2], (1, 2))
        self.assertEqual(F8.embed(2), 0)
        self.assertEqual(scheme.apply((2, 5), (3, 7, 6)), scheme.target.contract((2, 5), (3, 7, 6)))

    def test_interpolated_source_runs_through_regions_power_and_descent(self):
        c = ThreeSectorConstruction(F8, 2, 2)
        source = convolution_scheme(F8, 2, c.source_width)
        self.assertEqual(source.target, c.source)
        self.assertEqual(source.terms, 8)
        degeneration = c.generate_degeneration(source)
        self.assertEqual(degeneration.terms, 8)
        for coordinate in product(*(range(d) for d in c.shape)):
            self.assertEqual(degeneration.tensor_polynomial(*coordinate), c.polynomial_coefficient(*coordinate))
        exact = c.recover_power(2, source_scheme=source)
        self.assertEqual(exact.terms, 192)
        base = descend(exact)
        self.assertEqual(base.terms, 1728)
        self.assertTrue(base.verify())
        self.assertEqual(base.target, ThreeSectorConstruction(FiniteField(2), 2, 2).retained.tensor_power(2))
        left, right = (1, 0, 1, 1), tuple(i % 2 for i in range(49))
        self.assertEqual(base.apply(left, right), base.target.contract(left, right))

    def test_bad_nodes_dimensions_and_insufficient_fields(self):
        for nodes in ((0, 0, 1), (0, 1), (0, 1, 8), (0, 1, True), ()):
            with self.assertRaises(ValueError):
                convolution_scheme(F8, 2, 2, nodes)
        with self.assertRaises(ValueError):
            convolution_scheme(FiniteField(2), 2, 2)
        for bad in (-1, True, 1.0):
            with self.assertRaises(ValueError):
                convolution_scheme(F8, bad, 2)
        self.assertEqual(lagrange_basis(F8, ()), ())
        self.assertEqual(convolution_scheme(FiniteField(2), 0, 9).terms, 0)


if __name__ == "__main__":
    unittest.main()
