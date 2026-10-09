"""Full coefficient and genuine-extension checks for shared-dot matrix maps."""

from dataclasses import replace
from itertools import product
from random import Random
import unittest
from unittest.mock import patch

from .fields import FiniteField
from .maps import LinearMap
from .shared_dot_matrix import (export_shared_dot_matrix_order, export_shared_dot_star_order,
                                export_shared_dot_two_copy_matrix_order)
from .tensors import FiniteTensor
from .witnessed_constraints import StateAssemblyLimit, StateAssemblyLimits


FIELDS = (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)), FiniteField(11))


class SharedDotMatrixTests(unittest.TestCase):
    def test_actual_star_square_and_two_copy_restrictions(self):
        for field, q in product(FIELDS, (3, 4)):
            for function in (export_shared_dot_star_order, export_shared_dot_matrix_order,
                             export_shared_dot_two_copy_matrix_order):
                exported = function(field, q, "shared", "target")
                exported.require_valid()
                source = FiniteTensor.direct_sum(tuple(exported.registry[key] for key in exported.row.positive))
                self.assertEqual(source.restrict(*exported.row.maps), exported.registry["target"])
        export_shared_dot_star_order(FIELDS[0], 2, "c2", "star").require_valid()
        export_shared_dot_matrix_order(FIELDS[0], 2, "c2square", "m2").require_valid()

    def test_two_copy_bilinear_evaluation_and_signs(self):
        field, random = FIELDS[2], Random(20261009)
        exported = export_shared_dot_two_copy_matrix_order(field, 3, "shared", "m2")
        source = FiniteTensor.direct_sum((exported.registry["shared"],)*2)
        target = exported.registry["m2"]
        xmap, ymap, zmap = exported.row.maps
        for trial in range(10):
            x = tuple(random.randrange(field.order) for i in range(4))
            y = tuple(random.randrange(field.order) for i in range(4))
            pulled_x = tuple(field.sum(field.mul(x[r], xmap.rows[r][i]) for r in range(4)) for i in range(8))
            pulled_y = tuple(field.sum(field.mul(y[r], ymap.rows[r][i]) for r in range(4)) for i in range(8))
            self.assertEqual(zmap.apply(source.contract(pulled_x, pulled_y)), target.contract(x, y))
        odd = export_shared_dot_two_copy_matrix_order(FIELDS[3], 3, "shared", "m2")
        rows = list(odd.row.maps[1].rows)
        rows[2] = tuple(0 if i in (3, 7) else value for i, value in enumerate(rows[2]))
        bad = LinearMap(FIELDS[3], 8, rows)
        self.assertFalse(replace(odd, row=replace(odd.row, maps=(odd.row.maps[0], bad, odd.row.maps[2]))).verify())

    def test_caps_and_invalid_inputs_precede_allocations(self):
        with patch("reference.shared_dot_matrix.shared_dot_tensor", side_effect=AssertionError("allocation")):
            for function in (export_shared_dot_star_order, export_shared_dot_matrix_order,
                             export_shared_dot_two_copy_matrix_order):
                with self.assertRaises(StateAssemblyLimit):
                    function(FIELDS[0], 10**9, "source", "target")
                with self.assertRaises(StateAssemblyLimit):
                    function(FIELDS[0], 3, "source", "target", limits=StateAssemblyLimits(max_map_entries=1))
                with self.assertRaises(ValueError):
                    function(FIELDS[0], 3, "same", "same")
            with self.assertRaises(StateAssemblyLimit):
                export_shared_dot_two_copy_matrix_order(FIELDS[0], 3, "source", "target",
                    limits=StateAssemblyLimits(max_blocks=1))
        with self.assertRaises(ValueError):
            export_shared_dot_two_copy_matrix_order(FIELDS[0], 2, "source", "target")
