"""Proof-derived regions, actual local restrictions, and coordinate controls."""

from dataclasses import replace
from itertools import product
from random import Random
import unittest

from . import FiniteField, FiniteTensor, Polynomial, TensorScheme, ThreeSectorConstruction, convolution_tensor
from .sectors import BRANCHES


F2 = FiniteField(2)
F3 = FiniteField(3)
F4 = FiniteField(2, (1, 1, 1))
F5 = FiniteField(5)
F9 = FiniteField(3, (1, 0, 1))


class SectorTests(unittest.TestCase):
    def test_explicit_small_example_weights_and_support(self):
        c = ThreeSectorConstruction(F2, 2, 2)
        self.assertEqual(c.shape, (2, 7, 8))
        self.assertEqual(c.weights, ((0, 0), (0, 0, 1, 1, 1, 0, 0),
                                    (0, 0, 0, -1, -1, 0, 0, 0)))
        self.assertEqual(c.branch_support("left"),
                         frozenset(((0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 2))))
        self.assertEqual(c.branch_support("middle"),
                         frozenset(((1, 2, 3), (1, 3, 4), (0, 3, 3), (0, 4, 4))))
        self.assertEqual(c.branch_support("right"),
                         frozenset(((0, 5, 5), (0, 6, 6), (1, 5, 6), (1, 6, 7))))
        erased = frozenset((i, j, k) for i, j, k in product(range(2), range(7), range(8))
                           if c.erased.coefficient(i, j, k))
        self.assertEqual(erased, frozenset(((0, 2, 2), (1, 4, 5))))
        self.assertEqual(c.polynomial_coefficient(0, 3, 3), Polynomial(F2, (0, 1)))
        self.assertEqual(c.polynomial_coefficient(0, 2, 2), Polynomial(F2, (0, 0, 1)))
        # A negative unshifted weight is possible OUTSIDE convolution support.
        self.assertEqual(c.total_weight(0, 0, 3), -1)
        self.assertEqual(c.polynomial_coefficient(0, 0, 3), Polynomial(F2))

    def test_partition_and_formal_identity_over_parameter_grid(self):
        for field in (F2, F3, F4, F5, F9):
            for a, h in product(range(5), range(4)):
                with self.subTest(field=field.order, a=a, h=h):
                    c = ThreeSectorConstruction(field, a, h)
                    self.assertTrue(c.verify())
                    source_count = sum(int(x != 0) for x in c.source.coefficients)
                    self.assertEqual(source_count, a*c.source_width)
                    self.assertEqual(sum(c.retained.coefficients), 3*a*h)
                    self.assertEqual(sum(c.erased.coefficients), a*max(0, a-1))
                    for branch in BRANCHES:
                        self.assertEqual(len(c.branch_support(branch)), a*h)
                    for i, j, k in product(*(range(d) for d in c.shape)):
                        polynomial = c.polynomial_coefficient(i, j, k)
                        self.assertEqual(polynomial.coefficient(0), 0)
                        self.assertEqual(polynomial.coefficient(1), c.retained.coefficient(i, j, k))
                        self.assertEqual(polynomial.coefficient(2), c.erased.coefficient(i, j, k))
                        self.assertLessEqual(polynomial.degree, 2)
                        self.assertEqual(c.source.coefficient(i, j, k),
                                         sum(branch.coefficient(i, j, k) for branch in c.branches)
                                         + c.erased.coefficient(i, j, k))

    def test_branch_embeddings_preserve_all_local_coefficients(self):
        for a, h in ((1, 1), (2, 2), (3, 1), (3, 2)):
            c = ThreeSectorConstruction(F3, a, h)
            for name, branch in zip(BRANCHES, c.branches):
                for u, r, s in product(range(a), range(h), range(a+h-1)):
                    coordinate = c.branch_embedding(name, u, r, s)
                    self.assertEqual(branch.coefficient(*coordinate), int(u+r == s))
            self.assertEqual(c.branch_embedding("middle", 0, 0, 0),
                             (a-1, h, h+a-1))

    def test_first_input_is_shared_by_all_three_branches(self):
        c = ThreeSectorConstruction(F5, 3, 2)
        self.assertEqual(c.shape[0], 3)
        for name, branch in zip(BRANCHES, c.branches):
            self.assertEqual(branch.shape, c.source.shape)
            self.assertEqual({i for i, _, _ in c.branch_support(name)}, set(range(3)))
        first_map = c.local_maps[0]
        self.assertEqual(len(first_map), 3)
        for i, j in product(range(3), repeat=2):
            self.assertEqual(first_map[i][j], Polynomial.constant(F5, int(i == j)))
        # The same impulse input produces output in all three ambient regions.
        x, y = (1, 0, 0), (1,)*c.shape[1]
        for branch in c.branches:
            self.assertTrue(any(branch.contract(x, y)))
        with self.assertRaises(ValueError):
            c.retained.contract(x*3, y)

    def test_branches_compute_convolutions_and_reversed_middle_correlation(self):
        random = Random(501)
        for field in (F2, F3, F4, F5):
            for a, h in ((1, 2), (2, 2), (3, 1)):
                c = ThreeSectorConstruction(field, a, h)
                for _ in range(5):
                    x = tuple(random.randrange(field.order) for _ in range(a))
                    y = tuple(random.randrange(field.order) for _ in range(c.shape[1]))
                    left, middle, right = [[0]*c.shape[2] for _ in range(3)]
                    for i, r in product(range(a), range(h)):
                        left[i+r] = field.add(left[i+r], field.mul(x[i], y[r]))
                        output = c.right_start+i+r
                        right[output] = field.add(right[output], field.mul(x[i], y[c.right_start+r]))
                    for r in range(h):
                        middle[c.middle_z_start+r] = field.sum(field.mul(x[i], y[h+a-1-i+r])
                                                              for i in range(a))
                    expected = tuple(tuple(row) for row in (left, middle, right))
                    actual = tuple(branch.contract(x, y) for branch in c.branches)
                    self.assertEqual(actual, expected)
                    self.assertEqual(c.retained.contract(x, y),
                                     tuple(field.sum(row[k] for row in expected) for k in range(c.shape[2])))

    def test_full_local_restriction_matches_diagonal_implementation(self):
        # Independently contract the full L[output][input] matrices with C.
        for field in (F2, F3, F4):
            for a, h in ((2, 1), (3, 1)):
                c = ThreeSectorConstruction(field, a, h)
                lx, ly, lz = c.local_maps
                for i, j, k in product(*(range(d) for d in c.shape)):
                    value = Polynomial(field)
                    for old_i, old_j in product(range(a), range(c.source_width)):
                        old_k = old_i+old_j
                        value = value + lx[i][old_i]*ly[j][old_j]*lz[k][old_k]
                    self.assertEqual(value, c.polynomial_coefficient(i, j, k))
                for matrix, bound in zip(c.local_maps, (0, 1, 1)):
                    self.assertTrue(all(p.degree <= bound for row in matrix for p in row))

    def test_generated_polynomial_families_match_the_actual_restriction(self):
        for field in (F2, F3, F4, F9):
            for a, h in ((1, 1), (2, 2), (3, 1), (2, 0), (0, 2)):
                c = ThreeSectorConstruction(field, a, h)
                source = TensorScheme.from_tensor(c.source)
                if field.degree > 1:
                    source = source.rescale_terms(field.p, 1)
                degeneration = c.generate_degeneration(source)
                self.assertTrue(degeneration.verify())
                self.assertEqual(degeneration.target, c.retained)
                self.assertEqual(degeneration.leading, 1)
                self.assertEqual(degeneration.terms, source.terms)
                self.assertLessEqual(degeneration.degree_bound, 2)
                for i, j, k in product(*(range(d) for d in c.shape)):
                    self.assertEqual(degeneration.tensor_polynomial(i, j, k), c.polynomial_coefficient(i, j, k))
                # First-leg map is the identity, including extension coefficients.
                self.assertEqual(degeneration.a,
                                 tuple(tuple(Polynomial.constant(field, x) for x in row) for row in source.a))

    def test_missing_shift_and_missing_reversal_are_detected(self):
        class MissingShift(ThreeSectorConstruction):
            @property
            def diagonal_maps(self):
                first, second, third = super().diagonal_maps
                return first, second, tuple(Polynomial.constant(self.field, 1) for _ in third)

        class MissingReversal(ThreeSectorConstruction):
            def branch_embedding(self, name, u, r, s):
                coordinate = super().branch_embedding(name, u, r, s)
                return (u, coordinate[1], coordinate[2]) if name == "middle" else coordinate

        for cls in (MissingShift, MissingReversal):
            bad = cls(F3, 2, 2)
            self.assertFalse(bad.verify())
            with self.assertRaisesRegex(ValueError, "invalid three-sector"):
                bad.generate_degeneration()

    def test_bad_source_schemes_are_rejected(self):
        c = ThreeSectorConstruction(F3, 2, 2)
        source = TensorScheme.from_tensor(c.source)
        changed = [list(row) for row in source.c]
        changed[0][0] = 2
        bad = replace(source, c=changed)
        for scheme in (bad, TensorScheme.from_tensor(c.retained), "not a scheme",
                       TensorScheme.from_tensor(ThreeSectorConstruction(F2, 2, 2).source)):
            with self.assertRaises(ValueError):
                c.generate_degeneration(scheme)

    def test_empty_axes_natural_subtraction_and_input_validation(self):
        for a, b, output in ((0, 0, 0), (0, 3, 2), (3, 0, 2)):
            tensor = convolution_tensor(F2, a, b)
            self.assertEqual(tensor.shape, (a, b, output))
            self.assertFalse(any(tensor.coefficients))
        for a, h in ((0, 0), (0, 2), (1, 0), (3, 0)):
            c = ThreeSectorConstruction(F2, a, h)
            self.assertTrue(c.verify())
            self.assertFalse(any(c.retained.coefficients))
            self.assertTrue(c.generate_degeneration().verify())
        for bad in (-1, 1.0, True):
            with self.assertRaises(ValueError):
                ThreeSectorConstruction(F2, bad, 1)
            with self.assertRaises(ValueError):
                convolution_tensor(F2, 1, bad)
        c = ThreeSectorConstruction(F2, 2, 2)
        for coordinate in ((-1, 0, 0), (2, 0, 0), (0, 7, 0), (0, 0, 8)):
            with self.assertRaises(ValueError):
                c.total_weight(*coordinate)
        with self.assertRaises(ValueError):
            c.branch_support("unknown")
        with self.assertRaises(ValueError):
            c.branch_embedding("middle", 0, 2, 0)
        with self.assertRaises(ValueError):
            c.branch_embedding("left", 0, 0, 3)
        with self.assertRaises(ValueError):
            ThreeSectorConstruction("not a field", 1, 1)


if __name__ == "__main__":
    unittest.main()
