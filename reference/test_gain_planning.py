"""Exact recurrence choices, bounded generation, and one-time descent."""

from fractions import Fraction
import unittest
from unittest.mock import patch

from .catalyst_compiler import CompilerLimits, coordinate_catalyst
from .fields import FiniteField
from .gain_planning import compile_gain_target, plan_gain_catalytic_parameters
from .schemes import naive_multiply
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor


class GainPlanningTests(unittest.TestCase):
    def inputs(self, field=None, m=1):
        field = FiniteField(2) if field is None else field
        certificate = coordinate_catalyst(S=FiniteTensor.unit(field), m=m)
        return certificate, TensorScheme.from_tensor(certificate.S), TensorScheme.from_tensor(certificate.D)

    def test_gain_accepts_a_block_rejected_by_the_old_k_only_gap(self):
        inputs = self.inputs()
        tau = Fraction(31, 10)
        self.assertGreaterEqual(inputs[0].k**tau.denominator, inputs[0].d**tau.numerator)
        plan = plan_gain_catalytic_parameters(*inputs, tau, max_b=1)
        self.assertEqual((plan.status, plan.exponent, plan.iterations, plan.multiplier),
                         ("planned", 1, 1, 8))
        self.assertTrue(plan.numeric_inequality)

    def test_actual_generated_gap_and_independent_matrix_product(self):
        result = compile_gain_target(*self.inputs(), Fraction(31, 10))
        self.assertEqual(result.status, "verified")
        self.assertEqual((result.prime_scheme.n, result.prime_scheme.terms), (2, 8))
        left, right = [[1, 1], [0, 1]], [[0, 1], [1, 1]]
        self.assertEqual(result.prime_scheme.multiply(left, right),
                         naive_multiply(FiniteField(2), left, right))

    def test_strict_boundary_is_not_accepted(self):
        plan = plan_gain_catalytic_parameters(*self.inputs(), Fraction(3), max_b=3)
        self.assertEqual(plan.status, "resource_cap")
        self.assertFalse(plan.numeric_inequality)
        self.assertIsNone(plan.exponent)

    def test_caps_remain_incomplete(self):
        for caps in ({"max_b": 0}, {"max_iterations": 0}):
            result = compile_gain_target(*self.inputs(), Fraction(4), **caps)
            self.assertEqual(result.status, "resource_cap")
            self.assertIsNone(result.prime_scheme)

    def test_final_dense_cap_precedes_powered_compilation(self):
        with patch("reference.gain_planning.compile_gain_powered_catalyst",
                   side_effect=AssertionError("must not allocate a prefix")):
            result = compile_gain_target(*self.inputs(), Fraction(31, 10),
                                         limits=CompilerLimits(max_tensor_coefficients=1))
        self.assertEqual(result.plan.status, "planned")
        self.assertEqual(result.status, "resource_cap")
        self.assertIsNone(result.extension_scheme)

    def test_extension_and_characteristic_dividing_gain_descend_once(self):
        # k=10,m=2 gives Q=8 in characteristic two; two units do not disappear.
        inputs = self.inputs(FiniteField(2, (1, 1, 1)), m=2)
        result = compile_gain_target(*inputs, Fraction(5), max_iterations=4)
        self.assertEqual(result.status, "verified")
        self.assertEqual(result.plan.iterations, 2)
        self.assertEqual((result.extension_scheme.n, result.extension_scheme.terms), (4, 64))
        self.assertEqual(result.prime_scheme.terms, 4*64)
        self.assertTrue(result.prime_scheme.verify())

    def test_input_validation_is_not_a_synthetic_rank_oracle(self):
        certificate, S, D = self.inputs()
        with self.assertRaises(ValueError):
            plan_gain_catalytic_parameters(certificate, D, S, Fraction(4))
        with self.assertRaises(ValueError):
            plan_gain_catalytic_parameters(certificate, S, D, 4.0)
        with self.assertRaises(ValueError):
            plan_gain_catalytic_parameters(certificate, S, D, Fraction(4), max_b=True)


if __name__ == "__main__":
    unittest.main()
