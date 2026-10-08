"""Actual finite catalytic witnesses and explicitly conditional cost planning."""

from dataclasses import replace
from fractions import Fraction
import unittest

from .catalysts import (CatalyticCertificate, plan_catalytic_parameters,
                        plan_conditional_costs)
from .fields import FiniteField
from .maps import LinearMap
from .schemes import naive_scheme
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor


F2 = FiniteField(2)


def unit_catalyst_certificate(m=1):
    """An actual coefficient witness, intentionally unable to give tau=9/4.

    m scalar terms plus the eight ordinary 2x2 product terms give a
    restriction from (m+8) diagonal units. D has zero axes and S is unit.
    """
    matrix = naive_scheme(F2, 2, 2, 2).as_tensor_scheme()
    terms = TensorScheme.direct_sum((TensorScheme.from_tensor(FiniteTensor.unit(F2)),)*m+(matrix,))
    # Column q is the q-th vector of the decomposition; restricting a
    # diagonal unit sum by these matrices realizes its rank-one sum.
    maps = tuple(LinearMap(F2, terms.terms, tuple(tuple(row[out] for row in family)
                                                 for out in range(size)))
                 for family, size in zip((terms.a, terms.b, terms.c), terms.shape))
    return CatalyticCertificate(2, m+8, m, FiniteTensor.zero(F2, (0, 0, 0)), FiniteTensor.unit(F2), maps)


class CatalyticTests(unittest.TestCase):
    def test_actual_f2_nine_unit_witness_full_coefficients_and_gain_removal(self):
        certificate = unit_catalyst_certificate()
        self.assertEqual(certificate.source.shape, (9, 9, 9))
        self.assertEqual(certificate.target.shape, (5, 5, 5))
        self.assertTrue(certificate.verify())
        self.assertEqual(certificate.restricted, certificate.target)
        maps = certificate.strip_gain()
        self.assertEqual(tuple(mapping.output_size for mapping in maps), (4, 4, 4))
        self.assertEqual(certificate.source.restrict(*maps), certificate.core_target)
        self.assertEqual(sum(certificate.target.coefficients), 9)
        self.assertEqual(sum(certificate.core_target.coefficients), 8)

    def test_scalar_gain_is_an_integer_copy_count_in_characteristic_two(self):
        certificate = unit_catalyst_certificate(2)
        self.assertEqual(F2.embed(certificate.m), 0)
        self.assertEqual(certificate.source.shape, (10, 10, 10))
        self.assertEqual(certificate.target.shape, (6, 6, 6))
        self.assertTrue(certificate.verify())
        self.assertEqual(certificate.target.coefficient(0, 0, 0), 1)
        self.assertEqual(certificate.target.coefficient(1, 1, 1), 1)
        self.assertEqual(sum(certificate.target.coefficients), 10)
        self.assertEqual(certificate.source.restrict(*certificate.strip_gain()), certificate.core_target)

    def test_nonempty_D_is_preserved_when_stripping_gain(self):
        base = unit_catalyst_certificate()
        D = FiniteTensor(F2, (1, 2, 1), (1, 1))
        identity_maps = tuple(LinearMap.identity(F2, size) for size in D.shape)
        maps = []
        for identity, original in zip(identity_maps, base.maps):
            rows = tuple(row+(0,)*original.input_size for row in identity.rows)
            rows += tuple((0,)*identity.input_size+row for row in original.rows)
            maps.append(LinearMap(F2, identity.input_size+original.input_size, rows))
        certificate = replace(base, D=D, maps=tuple(maps))
        self.assertTrue(certificate.verify())
        self.assertEqual(certificate.source.restrict(*certificate.strip_gain()), certificate.core_target)
        self.assertEqual(certificate.core_target.coefficient(0, 0, 0), 1)
        self.assertEqual(certificate.core_target.coefficient(0, 1, 0), 1)

    def test_corrupt_local_map_fails_even_when_scalar_gain_is_correct(self):
        base = unit_catalyst_certificate()
        old = base.maps[0]
        rows = list(old.rows)
        rows[1] = (0,)*old.input_size
        corrupted = replace(base, maps=(LinearMap(F2, old.input_size, rows),)+base.maps[1:])
        self.assertEqual(corrupted.restricted.coefficient(0, 0, 0), 1)
        self.assertFalse(corrupted.verify())
        with self.assertRaisesRegex(ValueError, "invalid catalytic coefficient"):
            corrupted.require_valid()
        with self.assertRaises(ValueError):
            corrupted.strip_gain()

    def test_invalid_certificate_dimensions_fields_and_zero_auxiliary(self):
        base = unit_catalyst_certificate()
        for change in ({"d": 1}, {"d": True}, {"k": 0}, {"m": 0}, {"m": 1.0},
                       {"S": FiniteTensor.zero(F2, (1, 1, 1))},
                       {"D": FiniteTensor.zero(FiniteField(3), (0, 0, 0))},
                       {"maps": base.maps[:2]},
                       {"maps": (LinearMap.identity(F2, 9),)+base.maps[1:]}):
            with self.assertRaises(ValueError):
                replace(base, **change)

    def test_actual_certificate_cannot_certify_nine_quarters_gap(self):
        base = unit_catalyst_certificate()
        S = TensorScheme.from_tensor(base.S)
        D = TensorScheme.from_tensor(base.D)
        plan = plan_catalytic_parameters(base, S, D, Fraction(9, 4))
        self.assertEqual(plan.status, "gap_not_certified")
        self.assertIn("no conclusion", plan.reason)
        self.assertIsNone(plan.b)
        self.assertFalse(plan.numeric_inequality)
        self.assertGreaterEqual(9**4, 2**9)
        corrupted = replace(S, c=((0,),))
        with self.assertRaises(ValueError):
            plan_catalytic_parameters(base, corrupted, D, Fraction(9, 4))
        with self.assertRaises(ValueError):
            plan_catalytic_parameters(base, D, D, Fraction(9, 4))

    def test_certified_scheme_upper_bounds_are_used_without_rank_minimization(self):
        base = unit_catalyst_certificate()
        # Three unit terms sum to unit in characteristic two. This is an
        # honest supplied upper bound of three, rather than minimum rank one.
        S = TensorScheme(base.S, ((1,),)*3, ((1,),)*3, ((1,),)*3)
        D = TensorScheme.from_tensor(base.D)
        plan = plan_catalytic_parameters(base, S, D, Fraction(4))
        self.assertEqual((plan.status, plan.b, plan.iterations), ("conditional", 2, 1))
        self.assertEqual((plan.predicted_size, plan.predicted_rank), (4, 243))
        self.assertTrue(plan.numeric_inequality)

    def test_synthetic_cost_recurrence_exact_strict_inequality(self):
        # These costs have NO catalytic witness. Success is numerical only.
        plan = plan_conditional_costs(2, 4, 1, 0, 4, Fraction(9, 4))
        self.assertEqual((plan.status, plan.b, plan.iterations), ("conditional", 1, 17))
        self.assertEqual(plan.catalyst_cost, 0)
        self.assertEqual((plan.predicted_size, plan.predicted_rank), (2**17, 4**17))
        self.assertTrue(plan.numeric_inequality)
        self.assertIn("decomposition not constructed", plan.reason)
        # At iteration 16 there is equality, so strict inequality must fail.
        incomplete = plan_conditional_costs(2, 4, 1, 0, 4, Fraction(9, 4), max_iterations=16)
        self.assertEqual(incomplete.status, "resource_cap")
        self.assertFalse(incomplete.numeric_inequality)
        self.assertIn("incomplete", incomplete.reason)

    def test_synthetic_block_choice_and_additive_catalyst_cost(self):
        plan = plan_conditional_costs(2, 4, 16, 1, 1, Fraction(9, 4))
        self.assertEqual(plan.b, 17)
        self.assertEqual(plan.catalyst_cost, sum(8**i*4**(16-i) for i in range(17)))
        r, u = 1, 1
        for _ in range(plan.iterations):
            r = 16*(4**17*r+plan.catalyst_cost)
            u *= 2**17
        self.assertEqual((plan.predicted_rank, plan.predicted_size), (r, u))
        self.assertTrue(plan.numeric_inequality)
        for options in ({"max_b": 16}, {"max_b": 0}, {"max_iterations": 0}):
            incomplete = plan_conditional_costs(2, 4, 16, 1, 1, Fraction(9, 4), **options)
            self.assertEqual(incomplete.status, "resource_cap")
            self.assertIn("incomplete", incomplete.reason)

    def test_invalid_conditional_planner_inputs(self):
        good = dict(d=2, k=4, s_rank=1, catalyst_rank=0, extension_degree=1, tau=Fraction(9, 4))
        for change in ({"d": 1}, {"k": True}, {"s_rank": 0}, {"catalyst_rank": -1},
                       {"extension_degree": 0}, {"tau": 2.25}, {"tau": Fraction(0)},
                       {"tau": Fraction(-1)}, {"max_b": -1}, {"max_iterations": True}):
            with self.assertRaises(ValueError):
                plan_conditional_costs(**(good | change))


if __name__ == "__main__":
    unittest.main()
