"""Actual slice/commutator formulas and honest scope of shared-dot screens."""

from dataclasses import replace
from itertools import product
import unittest
from unittest.mock import patch

from .fields import FiniteField
from .shared_dot_catalyst_screen import screen_shared_dot_catalyst, shared_dot_commutator_control
from .shared_dot_orders import shared_dot_tensor
from .tensors import FiniteTensor
from .witnessed_constraints import StateAssemblyLimit, StateAssemblyLimits


FIELDS = (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)), FiniteField(11))


class SharedDotScreenTests(unittest.TestCase):
    def test_actual_invertible_slices_commutators_and_repetition(self):
        for field, q, gain in product(FIELDS, (2, 3, 4), (1, 2)):
            control = shared_dot_commutator_control(field, q, gain=gain)
            control.require_valid()
            self.assertEqual(control.identity_rank, gain+4*(q+1))
            self.assertEqual(control.commutator_rank, 4*(q+1))
            self.assertEqual(control.rank_lower_control, gain+6*(q+1))
            self.assertEqual(control.normalized[0].rows,
                tuple(tuple(int(i == j) for j in range(control.identity_rank)) for i in range(control.identity_rank)))
        repeated = shared_dot_commutator_control(FIELDS[2], 3, repetitions=2)
        self.assertEqual(repeated.rank_lower_control, 50)
        self.assertEqual(repeated.commutator_rank, 32)

    def test_all_field_exclusion_and_remaining_frontier(self):
        for field in FIELDS:
            for q, gain in ((1, 1), (2, 1), (4, 1), (3, 2), (2, 3)):
                auxiliary = shared_dot_tensor(field, q)
                for catalyst in (FiniteTensor.zero(field, (0,)*3), FiniteTensor.unit(field), auxiliary):
                    report = screen_shared_dot_catalyst(catalyst, auxiliary=auxiliary, m=gain)
                    self.assertEqual(report['status'], 'analytically_excluded')
                    self.assertEqual(report['exclusion'], 'all_field_integral_koszul' if q in (2, 3) else 'all_field_border_commutator')
                    self.assertTrue(report['arbitrary_finite_catalyst'])
                    self.assertFalse(report['finite_sampling_proves_all_fields'])
            shared = shared_dot_tensor(field, 3)
            zero = screen_shared_dot_catalyst(FiniteTensor.zero(field, (0,)*3), auxiliary=shared)
            self.assertEqual(zero['status'], 'analytically_excluded')
            self.assertEqual(zero['exclusion'], 'all_field_integral_koszul')
            self.assertEqual(zero['coefficient_controls']['rank_lower_control'], 25)
            self.assertEqual(screen_shared_dot_catalyst(FiniteTensor.unit(field), auxiliary=shared)['status'], 'analytically_excluded')
            c2 = shared_dot_tensor(field, 2)
            self.assertEqual(screen_shared_dot_catalyst(FiniteTensor.unit(field), auxiliary=c2)['status'], 'analytically_excluded')

    def test_corrupted_controls_and_auxiliaries_are_rejected(self):
        control = shared_dot_commutator_control(FIELDS[0], 3)
        for bad in (replace(control, commutator_rank=0), replace(control, identity_rank=16),
                    replace(control, vectors=((0,)*17,)+control.vectors[1:])):
            with self.assertRaises(ValueError):
                bad.require_valid()
        tensor = shared_dot_tensor(FIELDS[0], 3)
        coefficients = list(tensor.coefficients)
        coefficients[5] = 0
        damaged = FiniteTensor(FIELDS[0], tensor.shape, coefficients)
        self.assertEqual(screen_shared_dot_catalyst(FiniteTensor.unit(FIELDS[0]), auxiliary=damaged)['status'], 'unclassified')

    def test_caps_precede_allocations_and_invalid_parameters(self):
        with patch('reference.shared_dot_catalyst_screen.shared_dot_tensor', side_effect=AssertionError('allocation')):
            with self.assertRaises(StateAssemblyLimit):
                shared_dot_commutator_control(FIELDS[0], 10**9)
            with self.assertRaises(StateAssemblyLimit):
                shared_dot_commutator_control(FIELDS[0], 3, limits=StateAssemblyLimits(max_map_entries=1))
        shared = shared_dot_tensor(FIELDS[0], 3)
        self.assertEqual(screen_shared_dot_catalyst(FiniteTensor.unit(FIELDS[0]), auxiliary=shared,
            limits=StateAssemblyLimits(max_tensor_entries=1))['status'], 'resource_cap')
        self.assertEqual(screen_shared_dot_catalyst(FiniteTensor.unit(FIELDS[0]), auxiliary=shared,
            limits=StateAssemblyLimits(max_map_entries=1000))['status'], 'resource_cap')
        with self.assertRaises(ValueError):
            shared_dot_commutator_control(FIELDS[0], 3, gain=0)
