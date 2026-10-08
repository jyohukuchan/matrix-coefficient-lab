"""Formal projective interpolation and saved unseeded coefficient certificates."""

from itertools import product
import json
from pathlib import Path
from random import Random
import unittest

from .fields import FiniteField
from .polynomials import Polynomial
from .projective_convolution import projective_convolution_scheme
from .sectors import ThreeSectorConstruction, convolution_tensor
from .tensor_schemes import TensorScheme
from .tensors import matrix_multiplication_tensor


F4 = FiniteField(2, (1, 1, 1))
F16 = FiniteField(2, (1, 1, 0, 0, 1))


def direct_convolution(field, left, right):
    result = [0]*max(0, len(left)+len(right)-1)
    for i, j in product(range(len(left)), range(len(right))):
        result[i+j] = field.add(result[i+j], field.mul(left[i], right[j]))
    return tuple(result)


class ProjectiveConvolutionTests(unittest.TestCase):
    def test_every_monomial_reconstructed_as_formal_polynomial(self):
        for field, a, b, nodes in ((FiniteField(2), 2, 2, (0, 1)),
                                   (FiniteField(3), 2, 3, (2, 0, 1)),
                                   (F4, 3, 3, (2, 0, 3, 1)),
                                   (F16, 3, 4, (0, 2, 5, 8, 15))):
            scheme = projective_convolution_scheme(field, a, b, nodes)
            degree = a+b-2
            output = tuple(Polynomial(field, coefficients) for coefficients in scheme.c)
            self.assertEqual(output[-1].degree, degree)
            self.assertEqual(output[-1].coefficient(degree), 1)
            for monomial in range(degree+1):
                recovered = Polynomial(field)
                for node, polynomial in zip(nodes, output[:-1]):
                    recovered = recovered+Polynomial.constant(field, field.pow(node, monomial))*polynomial
                if monomial == degree:
                    recovered = recovered+output[-1]
                self.assertEqual(recovered, Polynomial.monomial(field, monomial))
            self.assertTrue(scheme.verify())

    def test_small_schemes_all_coefficients_and_all_inputs(self):
        for field, a, b in ((FiniteField(2), 2, 2), (FiniteField(3), 2, 2),
                            (FiniteField(3), 2, 3), (F4, 2, 3), (F4, 3, 3)):
            scheme = projective_convolution_scheme(field, a, b)
            self.assertEqual(scheme.terms, a+b-1)
            self.assertTrue(scheme.verify())
            for left, right in product(product(range(field.order), repeat=a),
                                       product(range(field.order), repeat=b)):
                self.assertEqual(scheme.apply(left, right), direct_convolution(field, left, right))

    def test_infinity_output_has_required_lower_coefficients(self):
        for field, a, b, expected in ((FiniteField(2), 2, 2, (0, 1, 1)),
                                      (FiniteField(3), 2, 3, (0, 2, 0, 1)),
                                      (F4, 3, 3, (0, 1, 0, 0, 1))):
            scheme = projective_convolution_scheme(field, a, b)
            self.assertEqual(scheme.c[-1], expected)
            self.assertNotEqual(scheme.c[-1], (0,)*(a+b-2)+(1,))
            self.assertEqual(scheme.a[-1], (0,)*(a-1)+(1,))
            self.assertEqual(scheme.b[-1], (0,)*(b-1)+(1,))
            # The nonzero formal polynomial vanishes on all finite nodes.
            polynomial = Polynomial(field, expected)
            for node in range(a+b-2):
                self.assertEqual(polynomial.evaluate(node), 0)

    def test_extension_node_encodings_and_actual_values(self):
        nodes = (15, 2, 0, 9, 5)
        scheme = projective_convolution_scheme(F16, 3, 4, iter(nodes))
        self.assertEqual(scheme.a[1], (1, 2, F16.pow(2, 2)))
        self.assertEqual(F16.embed(2), 0)
        random = Random(2233)
        for _ in range(24):
            left = tuple(random.randrange(F16.order) for _ in range(3))
            right = tuple(random.randrange(F16.order) for _ in range(4))
            self.assertEqual(scheme.apply(left, right), direct_convolution(F16, left, right))

    def test_scalar_polynomial_and_empty_boundaries(self):
        for a, b in ((1, 1), (1, 2), (1, 4), (3, 1)):
            field = FiniteField(3)
            scheme = projective_convolution_scheme(field, a, b)
            self.assertEqual(scheme.terms, a+b-1)
            self.assertTrue(scheme.verify())
            for left, right in product(product(range(3), repeat=a), product(range(3), repeat=b)):
                self.assertEqual(scheme.apply(left, right), direct_convolution(field, left, right))
        scalar = projective_convolution_scheme(FiniteField(2), 1, 1, ())
        self.assertEqual((scalar.a, scalar.b, scalar.c), (((1,),), ((1,),), ((1,),)))
        for a, b in ((0, 0), (0, 8), (9, 0)):
            scheme = projective_convolution_scheme(FiniteField(2), a, b)
            self.assertEqual(scheme.terms, 0)
            self.assertTrue(scheme.verify())
            self.assertEqual(scheme.apply((0,)*a, (0,)*b), (0,)*max(0, a+b-1))

    def test_bad_nodes_dimensions_and_insufficient_fields(self):
        for nodes in ((0,), (0, 1, 2), (0, 0), (0, True), (0, 4), (0, 1.0)):
            with self.assertRaises(ValueError):
                projective_convolution_scheme(F4, 2, 2, nodes)
        with self.assertRaises(ValueError):
            projective_convolution_scheme(FiniteField(2), 2, 3)
        with self.assertRaises(ValueError):
            projective_convolution_scheme(FiniteField(3), 3, 3)
        with self.assertRaises(ValueError):
            projective_convolution_scheme(F4, 0, 3, (0,))
        with self.assertRaises(ValueError):
            projective_convolution_scheme(F4, 1, 1, (0,))
        for bad in (-1, True, 1.0):
            with self.assertRaises(ValueError):
                projective_convolution_scheme(F4, bad, 2)
            with self.assertRaises(ValueError):
                projective_convolution_scheme(F4, 2, bad)
        with self.assertRaises(ValueError):
            projective_convolution_scheme(2, 2, 2)

    def test_projective_source_connects_to_sector_power_and_fixed_descent(self):
        construction = ThreeSectorConstruction(F4, 2, 1)
        source = projective_convolution_scheme(F4, 2, construction.source_width)
        self.assertEqual(source.terms, 5)
        self.assertEqual(source.target, construction.source)
        recovered = construction.recover_power(2, source_scheme=source)
        self.assertEqual(recovered.terms, 75)
        self.assertTrue(recovered.verify())
        base = recovered.descend()
        self.assertEqual(base.terms, 300)
        self.assertEqual(base.target, ThreeSectorConstruction(FiniteField(2), 2, 1).retained.tensor_power(2))
        base.require_exact()
        left, right = (1, 0, 1, 1), tuple(i % 2 for i in range(16))
        self.assertEqual(base.apply(left, right), base.target.contract(left, right))

    def test_saved_unseeded_matrix_certificate_without_solver(self):
        artifact = Path(__file__).parent/"research"/"unseeded-f2-rank7.json"
        data = json.loads(artifact.read_text(encoding="utf-8"))
        self.assertEqual(data["target"]["matrix_dimensions"], [2, 2, 2])
        field = FiniteField(data["target"]["field_characteristic"])
        target = matrix_multiplication_tensor(field, 2, 2, 2)
        scheme = TensorScheme(target, data["a"], data["b"], data["c"])
        self.assertEqual(scheme.terms, 7)
        scheme.require_exact()
        for left, right in product(product(range(2), repeat=4), repeat=2):
            self.assertEqual(scheme.apply(left, right), target.contract(left, right))
        # The same formal coefficient identity also holds after extension.
        extended = TensorScheme(matrix_multiplication_tensor(F4, 2, 2, 2), data["a"], data["b"], data["c"])
        extended.require_exact()

    def test_saved_unseeded_extension_convolution_without_solver(self):
        artifact = Path(__file__).parent/"research"/"unseeded-f4-c23.json"
        data = json.loads(artifact.read_text(encoding="utf-8"))
        field_data = data["target"]
        field = FiniteField(field_data["field_characteristic"], tuple(field_data["field_modulus"]))
        a, b = field_data["convolution_dimensions"]
        scheme = TensorScheme(convolution_tensor(field, a, b), data["a"], data["b"], data["c"])
        self.assertEqual(scheme.terms, 4)
        self.assertTrue(any(value >= field.p for family in (scheme.a, scheme.b, scheme.c)
                            for row in family for value in row))
        scheme.require_exact()
        for left, right in product(product(range(field.order), repeat=a), product(range(field.order), repeat=b)):
            self.assertEqual(scheme.apply(left, right), direct_convolution(field, left, right))


if __name__ == "__main__":
    unittest.main()
