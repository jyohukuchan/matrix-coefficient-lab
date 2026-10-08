"""Exact scalar elimination and repeated rectangular matrix sharing."""

from dataclasses import replace
import unittest
from unittest.mock import patch

from .catalyst_compiler import CompilerLimits, CompilerResourceLimit
from .fields import FiniteField
from .maps import LinearMap
from .rectangular_restrictions import share_first_axis, share_matrix_left_operand
from .scalar_summands import eliminate_scalar_summands
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, matrix_multiplication_tensor, shared_first_tensor


class ResearchRestrictionTests(unittest.TestCase):
    def test_scalar_elimination_with_cross_mixed_decomposition(self):
        for field in (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1))):
            core = FiniteTensor(field, (1, 2, 2), (1, 0, 0, 1))
            target = FiniteTensor.direct_sum((core, FiniteTensor.unit(field)))
            # Mix the final scalar output with two old output coordinates.
            # This makes the eliminated scalar pivot's old C entries nonzero,
            # requiring the actual substitution rather than just dropping it.
            source = TensorScheme.from_tensor(target)
            a, b = (1, 1), (1, 0, 1)
            c = (1, 1, 1)
            # Insert a term and its additive inverse, keeping the exact target.
            scheme = TensorScheme(target, (a,)+source.a+(a,), (b,)+source.b+(b,),
                                  (c,)+source.c+(tuple(field.neg(value) for value in c),))
            scheme.require_exact()
            result = eliminate_scalar_summands(scheme, core)
            self.assertEqual(result.terms, scheme.terms-1)
            self.assertEqual(result.target, core)
            self.assertTrue(result.verify())
            self.assertEqual(result.apply((1,), (1, 1)), core.contract((1,), (1, 1)))

    def test_multiple_scalar_gains_and_empty_core(self):
        field = FiniteField(2)
        core = FiniteTensor.zero(field, (0, 0, 0))
        target = FiniteTensor.direct_sum((core,)+(FiniteTensor.unit(field),)*3)
        result = eliminate_scalar_summands(TensorScheme.from_tensor(target), core, 3)
        self.assertEqual((result.shape, result.terms), ((0, 0, 0), 0))
        self.assertTrue(result.verify())
        self.assertEqual(eliminate_scalar_summands(result, core, 0), result)

    def test_reject_general_summand_cancellation_and_wrong_target(self):
        field = FiniteField(2)
        core = FiniteTensor.unit(field)
        nonscalar = FiniteTensor(field, (2, 2, 2), (1, 0, 0, 0, 0, 0, 0, 0))
        scheme = TensorScheme.from_tensor(FiniteTensor.direct_sum((core, nonscalar)))
        with self.assertRaisesRegex(ValueError, "exactly the core"):
            eliminate_scalar_summands(scheme, core, 2)
        for copies in (-1, True, 1000000000):
            with self.assertRaises(ValueError):
                eliminate_scalar_summands(scheme, core, copies)

    def test_sharing_rectangular_operand_keeps_scalar_gains(self):
        for field in (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1))):
            restriction = share_matrix_left_operand(field, 4, 2, 4, 2, scalar_gain=1)
            self.assertEqual(restriction.source.shape, (18, 18, 34))
            rectangle = matrix_multiplication_tensor(field, 4, 2, 8)
            expected = FiniteTensor.direct_sum((FiniteTensor.unit(field),)*2+(rectangle,))
            self.assertEqual(restriction.target, expected)
            self.assertTrue(restriction.verify())
            self.assertEqual(restriction.target.shape, (10, 18, 34))
            # X must sum both original matrix blocks, rather than pick one.
            old = restriction.maps[0]
            rows = [list(row) for row in old.rows]
            rows[2][10] = 0
            bad = replace(restriction, maps=(LinearMap(field, old.input_size, rows),)+restriction.maps[1:])
            self.assertFalse(bad.verify())

    def test_general_rectangles_single_block_and_resource_caps(self):
        field = FiniteField(2)
        for n, inner, m, count in ((1, 3, 2, 1), (2, 3, 1, 3), (3, 1, 2, 2)):
            restriction = share_matrix_left_operand(field, n, inner, m, count)
            self.assertEqual(restriction.target, matrix_multiplication_tensor(field, n, inner, m*count))
            self.assertTrue(restriction.verify())
        with patch("reference.rectangular_restrictions.matrix_multiplication_tensor") as allocation:
            with self.assertRaises(CompilerResourceLimit):
                share_matrix_left_operand(field, 4, 2, 4, 10**9)
            allocation.assert_not_called()
        with self.assertRaises(CompilerResourceLimit):
            share_matrix_left_operand(field, 4, 2, 4, 2, limits=CompilerLimits(max_map_entries=1))

    def test_generic_sharing_retains_copy_labels_and_gains(self):
        field = FiniteField(2, (1, 1, 1))
        tensor = FiniteTensor(field, (2, 1, 3), (1, 2, 0, 0, 1, 3))
        restriction = share_first_axis(tensor, 3, scalar_gain=2)
        expected = FiniteTensor.direct_sum((FiniteTensor.unit(field),)*6+(shared_first_tensor((tensor,)*3),))
        self.assertEqual(restriction.target, expected)
        self.assertTrue(restriction.verify())
        self.assertEqual((restriction.source.shape, restriction.target.shape), ((12, 9, 15), (8, 9, 15)))
        # Two independent tensor labels still cannot contribute together.
        self.assertEqual(restriction.target.coefficient(6, 6, 10), 0)
        self.assertEqual(restriction.target.coefficient(6, 6, 7), 2)


if __name__ == "__main__":
    unittest.main()
