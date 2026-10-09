"""Adjugate witness and actual rank-five controls, including extension fields."""

from random import Random
import unittest

from .fields import FiniteField
from .mixed_dot_border import adjugate_slice_relation, mixed_dot_border_control
from .mixed_dot_catalyst_screen import MixedDotScreenBudget, MixedDotScreenLimit
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor


FIELDS = (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)), FiniteField(11))


class MixedDotBorderTests(unittest.TestCase):
    def test_unit_cofactor_and_nonzero_obstruction_in_every_control_field(self):
        for field in FIELDS:
            report = mixed_dot_border_control(field)
            self.assertEqual(report['border_rank'], 6)
            self.assertEqual(report['obstruction_coefficient'], field.neg(1))
            middle = report['slices'][1].rows
            minor = tuple(tuple(middle[r][c] for c in range(5) if c != 2) for r in range(5) if r != 3)
            self.assertEqual(minor, tuple(tuple(int(i == j) for j in range(4)) for i in range(4)))

    def test_actual_rank_five_decompositions_obey_the_adjugate_identity(self):
        random = Random(20261009)
        for field in FIELDS:
            for trial in range(3):
                families = tuple(tuple(tuple(random.randrange(field.order) for i in range(5)) for r in range(5))
                                 for axis in range(3))
                target = FiniteTensor.from_function(field, (5,)*3, lambda i, j, k:
                    field.sum(field.mul(field.mul(a[i], b[j]), c[k]) for a, b, c in zip(*families)))
                TensorScheme(target, *families).require_exact()
                vectors = tuple(tuple(random.randrange(field.order) for i in range(5)) for r in range(3))
                _, _, bracket = adjugate_slice_relation(target, vectors)
                self.assertFalse(any(value for row in bracket.rows for value in row))

    def test_caps_and_invalid_parameters(self):
        for budget in (MixedDotScreenBudget(max_matrix_entries=24), MixedDotScreenBudget(max_slice_cases=2)):
            with self.assertRaises(MixedDotScreenLimit):
                mixed_dot_border_control(FIELDS[0], budget=budget)
        with self.assertRaises(ValueError):
            adjugate_slice_relation(FiniteTensor.unit(FIELDS[0]), ((1,),)*3)
