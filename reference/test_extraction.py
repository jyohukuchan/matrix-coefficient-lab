"""Exact proof-branch extraction, boundary identification, and matrix execution."""

from dataclasses import replace
from itertools import product
import unittest

from .convolution import convolution_scheme
from .extraction import (branch_maps, extract_branch, factor_singleton_leg,
                         proof_matrix_pipeline, simplify_scheme, square_from_boundary)
from .fields import FiniteField
from .maps import LinearMap
from .schemes import naive_multiply
from .sectors import BRANCHES, ThreeSectorConstruction, convolution_tensor
from .separation import FiniteSeparation
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, matrix_multiplication_tensor


F16 = FiniteField(2, (1, 1, 0, 0, 1))


class ExtractionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pipeline = proof_matrix_pipeline()

    def test_real_proof_path_records_counts_and_extracts_all_branches(self):
        p = self.pipeline
        self.assertEqual(dict(p.stage_terms), {"convolution": 5, "retained": 10, "tagged": 10,
                                               "projected": 150, "separated": 750, "matrix": 8})
        self.assertEqual(p.period, 15)
        for raw, merged, factored in zip(p.branch_schemes, p.merged_branches, p.factored_branches):
            self.assertEqual(raw.target, convolution_tensor(F16, 2, 1))
            self.assertEqual(raw.terms, 750)
            self.assertLess(merged.terms, raw.terms)
            self.assertEqual(factored.terms, 2)
            for s in (raw, merged, factored):
                self.assertTrue(s.verify())
                self.assertEqual(s.apply((7, 13), (9,)), (F16.mul(7, 9), F16.mul(13, 9)))

    def test_auxiliary_dot_selection_and_middle_reverse_swap(self):
        c = ThreeSectorConstruction(F16, 2, 2)
        separation = FiniteSeparation(c.branches)
        for branch, dot in product(BRANCHES, range(3)):
            maps = branch_maps(c, separation, branch, dot)
            selected = separation.target.restrict(*maps)
            if branch == "middle":
                self.assertEqual(selected.shape, (2, 3, 2))
                selected = selected.permute_axes((0, 2, 1))
            self.assertEqual(selected, convolution_tensor(F16, 2, 2))
        correct = list(branch_maps(c, separation, "middle"))
        # Replacing reversed X with unreversed X changes actual coefficients.
        correct[0] = LinearMap.selection(F16, separation.shape[0], (1, 4))
        wrong = separation.target.restrict(*correct).permute_axes((0, 2, 1))
        self.assertNotEqual(wrong, convolution_tensor(F16, 2, 2))
        # Selecting different Y,Z dot indices annihilates the block.
        mismatch = list(branch_maps(c, separation, "left", 0))
        mismatch[2] = branch_maps(c, separation, "left", 1)[2]
        self.assertFalse(any(separation.target.restrict(*mismatch).coefficients))

    def test_proportional_merging_preserves_extension_scalars_and_cancellation(self):
        f = F16
        a, b = (1, 7), (1, 13)
        c, d = (9, 4), (3, 11)
        target = FiniteTensor.from_function(f, (2, 2, 2), lambda i, j, k:
                                            f.mul(f.mul(a[i], b[j]), f.add(c[k], d[k])))
        scaled_a = tuple(f.mul(2, x) for x in a)
        scaled_b = tuple(f.mul(3, x) for x in b)
        scaled_d = tuple(f.div(x, f.mul(2, 3)) for x in d)
        scheme = TensorScheme(target, (a, scaled_a, (0, 0)), (b, scaled_b, b), (c, scaled_d, c))
        simplified = simplify_scheme(scheme)
        self.assertEqual(simplified.terms, 1)
        self.assertTrue(simplified.verify())
        cancelled = TensorScheme(FiniteTensor.zero(f, (2, 2, 2)), (a, a), (b, b), (c, c))
        self.assertEqual(simplify_scheme(cancelled).terms, 0)

    def test_singleton_matrix_factorization_with_nontrivial_pivots(self):
        f = FiniteField(5)
        matrix = ((0, 0, 2), (3, 1, 4), (3, 1, 1), (0, 0, 0))
        target = FiniteTensor.from_function(f, (4, 1, 3), lambda i, j, k: matrix[i][k])
        # Generated families are nontrivially scaled, not target-coordinate seeds.
        source = TensorScheme.from_tensor(target).rescale_terms(2, 3)
        factored = factor_singleton_leg(source)
        self.assertEqual(factored.terms, 2)
        self.assertTrue(factored.verify())
        self.assertEqual(factored.apply((1, 2, 3, 4), (2,)), target.contract((1, 2, 3, 4), (2,)))
        zero = TensorScheme.from_tensor(FiniteTensor.zero(f, (3, 1, 2)))
        self.assertEqual(factor_singleton_leg(zero).terms, 0)

    def test_cyclic_tensor_product_and_nonidentity_side_reindexing(self):
        for a in (1, 2, 3):
            boundary = convolution_scheme(F16, a, 1)
            factors = tuple(boundary.permute_axes(order) for order in ((0, 1, 2), (1, 2, 0), (2, 0, 1)))
            unindexed = factors[0].target.tensor_product(factors[1].target).tensor_product(factors[2].target)
            for i, j, k in product(range(a), repeat=3):
                self.assertEqual(unindexed.coefficient(i*a+j, k*a+j, i*a+k), 1)
            if a > 1:
                self.assertNotEqual(unindexed, matrix_multiplication_tensor(F16, a, a, a))
            scheme = square_from_boundary((boundary,)*3)
            self.assertEqual(scheme.terms, a**3)
            self.assertEqual(scheme.as_tensor_scheme().target, matrix_multiplication_tensor(F16, a, a, a))
            self.assertTrue(scheme.verify())

    def test_generated_matrix_execution_and_fixed_field_descent_after_power(self):
        p = self.pipeline
        s = p.matrix_scheme
        for seed in range(8):
            left = [[(seed+3*i+5*j) % 16 for j in range(2)] for i in range(2)]
            right = [[(7*seed+2*i+9*j) % 16 for j in range(2)] for i in range(2)]
            self.assertEqual(s.multiply(left, right), naive_multiply(F16, left, right))
        for exponent in (0, 1, 2):
            base = p.descend_power(exponent)
            self.assertEqual(base.field, FiniteField(2))
            self.assertEqual(base.n, 2**exponent)
            self.assertEqual(base.terms, 16*8**exponent)
            self.assertTrue(base.verify())
            n = base.n
            left = [[(i+j) % 2 for j in range(n)] for i in range(n)]
            right = [[(i+2*j+1) % 2 for j in range(n)] for i in range(n)]
            self.assertEqual(base.multiply(left, right), naive_multiply(base.field, left, right))

    def test_invalid_targets_maps_and_corrupt_certificates_are_rejected(self):
        c = ThreeSectorConstruction(F16, 2, 1)
        sep = FiniteSeparation(c.branches)
        for branch, dot in (("other", 0), ("left", -1), ("left", 3), ("left", True)):
            with self.assertRaises(ValueError):
                branch_maps(c, sep, branch, dot)
        with self.assertRaises(ValueError):
            branch_maps(c, FiniteSeparation((c.branches[0],)*3), "left")
        with self.assertRaises(ValueError):
            extract_branch(c, sep, self.pipeline.branch_schemes[0], "left")
        s = self.pipeline.factored_branches[0]
        corrupted = replace(s, c=((0, 0),)*s.terms)
        for call in (lambda: simplify_scheme(corrupted), lambda: factor_singleton_leg(corrupted),
                     lambda: square_from_boundary((s, s, corrupted))):
            with self.assertRaises(ValueError):
                call()
        with self.assertRaises(ValueError):
            factor_singleton_leg(convolution_scheme(F16, 2, 2))
        for family in ((), (s, s), (s, s, convolution_scheme(F16, 3, 1))):
            with self.assertRaises(ValueError):
                square_from_boundary(family)
        for a in (-1, 0, True, 2.0):
            with self.assertRaises(ValueError):
                proof_matrix_pipeline(a)


if __name__ == "__main__":
    unittest.main()
