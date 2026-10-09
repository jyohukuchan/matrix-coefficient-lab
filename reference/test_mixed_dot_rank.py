"""Controls for the actual replicated rectangular component, not rank samples."""

from dataclasses import replace
import unittest
from unittest.mock import patch

from .fields import FiniteField
from .maps import LinearMap
from .mixed_dot_rank import mixed_dot_rank_control
from .tensors import FiniteTensor
from .witnessed_constraints import StateAssemblyLimit, StateAssemblyLimits


class MixedDotRankTests(unittest.TestCase):
    def test_replicated_component_and_first_concision(self):
        for field in (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1))):
            for length, copies in ((1, 1), (2, 1), (2, 2), (3, 1)):
                with self.subTest(field=field.order, length=length, copies=copies):
                    control = mixed_dot_rank_control(field, length, repetitions=copies)
                    self.assertEqual(control.first_flattening_rank, 4*(1+2*length)*copies)
                    self.assertEqual(control.rectangle.shape, (4*copies, 4*length*copies, 4*length*copies))
                    self.assertEqual(control.core.restrict(*control.maps), control.rectangle)
                    self.assertEqual(control.analytic_core_rank_lower, 15*length*copies)

    def test_coefficient_map_and_rank_corruption_fail(self):
        control = mixed_dot_rank_control(FiniteField(2))
        bad = FiniteTensor(control.rectangle.field, control.rectangle.shape,
                           (0,)+control.rectangle.coefficients[1:])
        for forged in (replace(control, rectangle=bad),
                       replace(control, first_flattening_rank=19),
                       replace(control, maps=(LinearMap.selection(control.core.field, 20, (5, 0, 10, 15)),)+control.maps[1:])):
            with self.assertRaises(ValueError):
                forged.require_valid()

    def test_resource_caps_precede_tensor_allocations(self):
        for limits in (StateAssemblyLimits(max_tensor_entries=1),
                       StateAssemblyLimits(max_map_entries=1),
                       StateAssemblyLimits(max_blocks=1),
                       StateAssemblyLimits(max_constraint_entries=1)):
            with patch('reference.mixed_dot_rank.mixed_dot_length_tensor', side_effect=AssertionError('allocation')):
                with self.assertRaises(StateAssemblyLimit):
                    mixed_dot_rank_control(FiniteField(2), limits=limits)
        with self.assertRaises(StateAssemblyLimit):
            mixed_dot_rank_control(FiniteField(2), 10**9)

    def test_parameter_validation(self):
        for kwargs in ({'length': True}, {'length': 0}, {'repetitions': -1}, {'limits': None}):
            with self.assertRaises(ValueError):
                mixed_dot_rank_control(FiniteField(2), **kwargs)
