"""Connect the proof's actual regions to powers, interpolation, and descent."""

from dataclasses import replace
from itertools import product
from random import Random
import unittest
from unittest.mock import patch

from . import (FiniteField, FiniteTensor, Polynomial, TensorDegeneration,
               TensorScheme, ThreeSectorConstruction, descend)


F2 = FiniteField(2)
F3 = FiniteField(3)
F4 = FiniteField(2, (1, 1, 1))
F9 = FiniteField(3, (1, 0, 1))

# Independently enumerated a=2,h=1 regions in shape (2,4,5).
RETAINED = frozenset(((0, 0, 0), (1, 0, 1), (1, 1, 2),
                      (0, 2, 2), (0, 3, 3), (1, 3, 4)))
ERASED = frozenset(((0, 1, 1), (1, 2, 3)))


class SectorPipelineTests(unittest.TestCase):
    def assert_computes_target(self, scheme):
        self.assertTrue(scheme.verify())
        random = Random(719)
        for _ in range(3):
            left = tuple(random.randrange(scheme.field.order) for _ in range(scheme.shape[0]))
            right = tuple(random.randrange(scheme.field.order) for _ in range(scheme.shape[1]))
            self.assertEqual(scheme.apply(left, right), scheme.target.contract(left, right))

    def test_formal_power_uses_one_parameter_and_preserves_noise(self):
        for field in (F2, F3, F4):
            c = ThreeSectorConstruction(field, 2, 1)
            powered = c.generate_degeneration().tensor_power(2)
            self.assertEqual((powered.leading, powered.degree_bound, powered.recovery_node_count), (2, 4, 3))
            for i, j, k in product(range(4), range(16), range(25)):
                first = i//2, j//4, k//5
                second = i % 2, j % 4, k % 5
                r1, r2 = int(first in RETAINED), int(second in RETAINED)
                e1, e2 = int(first in ERASED), int(second in ERASED)
                expected = Polynomial(field, (0, 0, r1*r2, field.embed(r1*e2+e1*r2), e1*e2))
                self.assertEqual(powered.tensor_polynomial(i, j, k), expected)
            self.assertTrue(powered.verify())

    def test_recovered_powers_match_independent_region_product(self):
        for field in (F4, F9):
            c = ThreeSectorConstruction(field, 2, 1)
            for exponent in (0, 1, 2):
                with self.subTest(field=field.order, exponent=exponent):
                    scheme = c.recover_power(exponent)
                    self.assertEqual(scheme.shape, (2**exponent, 4**exponent, 5**exponent))
                    self.assertEqual(scheme.terms, (exponent+1)*8**exponent)
                    self.assertEqual(scheme.target, c.retained.tensor_power(exponent))
                    self.assert_computes_target(scheme)
                    base = descend(scheme)
                    self.assertEqual(base.terms, field.degree**2*scheme.terms)
                    self.assertEqual(base.target, ThreeSectorConstruction(FiniteField(field.p), 2, 1)
                                     .retained.tensor_power(exponent))
                    self.assert_computes_target(base)
            expected = FiniteTensor.from_function(field, (4, 16, 25), lambda i, j, k:
                                                 int((i//2, j//4, k//5) in RETAINED
                                                     and (i % 2, j % 4, k % 5) in RETAINED))
            self.assertEqual(scheme.target, expected)
            self.assertEqual(sum(scheme.target.coefficients), 36)
            # All nine ordered branch pairs still share a^2 input coordinates.
            self.assertEqual(scheme.shape[0], 4)

    def test_documented_a2_h2_example(self):
        c = ThreeSectorConstruction(F4, 2, 2)
        scheme = c.recover_power(2)
        self.assertEqual(scheme.shape, (4, 49, 64))
        self.assertEqual(scheme.terms, 588)
        self.assertEqual(sum(scheme.target.coefficients), 144)
        self.assert_computes_target(scheme)
        base = descend(scheme)
        self.assertEqual(base.terms, 2352)
        self.assertEqual(base.target, ThreeSectorConstruction(F2, 2, 2).retained.tensor_power(2))
        self.assert_computes_target(base)

    def test_recovery_and_descent_orders_have_fixed_overhead(self):
        c = ThreeSectorConstruction(F4, 2, 1)
        recovered = c.recover_power(2)
        repeated_recovery = c.recover_power(1).tensor_power(2)
        self.assertEqual((recovered.terms, repeated_recovery.terms), (192, 256))
        self.assertEqual(recovered.target, repeated_recovery.target)
        self.assert_computes_target(repeated_recovery)
        base = descend(recovered)
        repeated_descent = descend(c.recover_power(1)).tensor_power(2)
        self.assertEqual((base.terms, repeated_descent.terms), (768, 4096))
        self.assertEqual(base.target, repeated_descent.target)
        self.assert_computes_target(repeated_descent)

    def test_custom_extension_source_and_node_order(self):
        c = ThreeSectorConstruction(F9, 2, 1)
        source = TensorScheme.from_tensor(c.source).rescale_terms(3, 4)
        self.assertTrue(any(x >= F9.p for row in source.a for x in row))
        scheme = c.recover_power(2, nodes=iter((8, 3, 1, 5)), source_scheme=source)
        self.assertEqual(scheme.terms, 4*source.terms**2)
        self.assertEqual(scheme.target, c.retained.tensor_power(2))
        self.assert_computes_target(scheme)
        self.assert_computes_target(descend(scheme))

    def test_unavailable_or_bad_nodes_fail_before_power_allocation(self):
        cases = ((F2, 1, None), (F3, 2, None), (F4, 3, None), (F4, 10**6, None),
                 (F4, 2, (1, 2)), (F4, 2, (1, 1, 2)), (F4, 2, (0, 1, 2)),
                 (F4, 2, (1, 2, 4)), (F4, 2, ()), (F4, 2, (1, 2, True)))
        for field, exponent, nodes in cases:
            with self.subTest(field=field.order, exponent=exponent, nodes=nodes):
                with patch.object(TensorDegeneration, "tensor_power") as powering:
                    with self.assertRaises(ValueError):
                        ThreeSectorConstruction(field, 2, 1).recover_power(exponent, nodes)
                    powering.assert_not_called()

    def test_zeroth_power_empty_axes_and_zero_target(self):
        unit = ThreeSectorConstruction(F2, 2, 1).recover_power(0)
        self.assertEqual(unit.target, FiniteTensor.unit(F2))
        self.assertEqual(unit.terms, 1)
        self.assertEqual(unit.apply((1,), (1,)), (1,))
        for a, h in ((0, 0), (0, 2), (1, 0), (2, 0)):
            c = ThreeSectorConstruction(F4, a, h)
            for exponent in (0, 1, 2):
                scheme = c.recover_power(exponent)
                self.assertEqual(scheme.target, c.retained.tensor_power(exponent))
                self.assert_computes_target(scheme)

    def test_source_and_exponent_are_validated_even_for_zero_power(self):
        c = ThreeSectorConstruction(F4, 2, 1)
        source = TensorScheme.from_tensor(c.source)
        corrupted = [list(row) for row in source.c]
        corrupted[0][0] = 0
        for bad in (replace(source, c=corrupted), TensorScheme.from_tensor(c.retained),
                    TensorScheme.from_tensor(ThreeSectorConstruction(F2, 2, 1).source)):
            with self.assertRaises(ValueError):
                c.recover_power(0, source_scheme=bad)
        for exponent in (-1, True, 1.0):
            with self.assertRaises(ValueError):
                c.recover_power(exponent)

    def test_evaluations_and_corrupted_recovery_cannot_be_applied(self):
        c = ThreeSectorConstruction(F4, 2, 1)
        evaluation = c.generate_degeneration().tensor_power(2).evaluate_at(1)
        self.assertFalse(evaluation.verify())
        with self.assertRaises(ValueError):
            evaluation.apply((0,)*4, (0,)*16)
        recovered = c.recover_power(2)
        corrupted = [list(row) for row in recovered.c]
        corrupted[0][0] = F4.add(corrupted[0][0], 1)
        bad = replace(recovered, c=corrupted)
        self.assertFalse(bad.verify())
        with self.assertRaises(ValueError):
            bad.apply((0,)*4, (0,)*16)


if __name__ == "__main__":
    unittest.main()
