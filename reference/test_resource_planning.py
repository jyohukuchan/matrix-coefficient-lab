"""Cost forecasts must distinguish hypotheses, actual inputs and outputs."""

from fractions import Fraction
import unittest

from .catalyst_compiler import coordinate_catalyst
from .fields import FiniteField
from .gain_planning import plan_gain_catalytic_parameters
from .resource_planning import (CandidateCosts, ExpansionBudget, forecast_candidate,
                                forecast_verified_inputs)
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor


class ResourcePlanningTests(unittest.TestCase):
    def test_hypothetical_first_sufficient_block_is_25_and_exceeds_dense_budget(self):
        costs = CandidateCosts(2, 5, 1, 4, 0, (2, 3, 3), (0, 0, 0))
        forecast = forecast_candidate(costs, Fraction(12, 5), max_b=25, max_iterations=1)
        sufficient = [b for b in forecast.blocks if b.strict_gap]
        self.assertEqual(len(sufficient), 1)
        block = sufficient[0]
        self.assertEqual((block.block_exponent, block.matrix_size, block.multiplier),
                         (25, 33_554_432, 1_092_751_820_893_346_269))
        self.assertEqual(block.tensor_entries, block.matrix_size**6)
        self.assertEqual(block.scheme_scalars, 3*block.prime_terms*block.matrix_size**2)
        self.assertIn("scheme_scalars", block.budget_failures)
        self.assertFalse(block.within_expansion_budget)
        self.assertEqual((forecast.input_status, forecast.output_status),
                         ("unverified_assumptions", "numerical_only"))

    def test_verified_forecast_agrees_with_actual_planner(self):
        certificate = coordinate_catalyst(S=FiniteTensor.unit(FiniteField(2)))
        S, D = TensorScheme.from_tensor(certificate.S), TensorScheme.from_tensor(certificate.D)
        tau = Fraction(31, 10)
        forecast = forecast_verified_inputs(certificate, S, D, tau, max_b=1)
        plan = plan_gain_catalytic_parameters(certificate, S, D, tau, max_b=1)
        block = forecast.blocks[0]
        self.assertEqual((block.matrix_size, block.prime_terms, block.multiplier),
                         (plan.predicted_size, plan.predicted_terms, plan.multiplier))
        self.assertEqual(forecast.input_status, "verified_inputs")
        self.assertEqual(forecast.output_status, "numerical_only")
        self.assertTrue(block.within_expansion_budget)

    def test_one_time_extension_overhead_and_integer_gain(self):
        field = FiniteField(2, (1, 1, 1))
        certificate = coordinate_catalyst(S=FiniteTensor.unit(field), m=2)
        S, D = TensorScheme.from_tensor(certificate.S), TensorScheme.from_tensor(certificate.D)
        block = forecast_verified_inputs(certificate, S, D, Fraction(5), max_b=1).blocks[0]
        self.assertEqual((block.powered_gain, block.iterations, block.matrix_size), (2, 2, 4))
        self.assertEqual((block.extension_terms, block.prime_terms), (64, 256))
        self.assertEqual(block.packed_coefficient_bytes, (block.scheme_scalars+7)//8)

    def test_strict_boundary_and_iteration_cap_are_not_success(self):
        costs = CandidateCosts(2, 9, 1, 1, 0, (1, 1, 1), (0, 0, 0))
        for tau, cap in ((Fraction(3), 4), (Fraction(4), 0)):
            forecast = forecast_candidate(costs, tau, max_b=1, max_iterations=cap)
            self.assertFalse(forecast.blocks[0].strict_gap)
        empty = forecast_candidate(costs, Fraction(4), max_b=0)
        self.assertEqual(empty.blocks, ())

    def test_budget_counts_are_not_an_output_certificate(self):
        costs = CandidateCosts(2, 9, 1, 1, 0, (1, 1, 1), (0, 0, 0))
        forecast = forecast_candidate(costs, Fraction(4), max_b=1,
                                      budget=ExpansionBudget(max_terms=0))
        self.assertEqual(forecast.blocks[0].budget_failures, ("terms",))
        with self.assertRaises(ValueError):
            ExpansionBudget(max_terms=True)
        with self.assertRaises(ValueError):
            forecast_candidate(costs, 2.4)
        with self.assertRaises(ValueError):
            forecast_candidate(costs, Fraction(4), max_b=True)
        certificate = coordinate_catalyst(S=FiniteTensor.unit(FiniteField(2)))
        S, D = TensorScheme.from_tensor(certificate.S), TensorScheme.from_tensor(certificate.D)
        with self.assertRaises(ValueError):
            forecast_verified_inputs(certificate, D, S, Fraction(4))


if __name__ == "__main__":
    unittest.main()
