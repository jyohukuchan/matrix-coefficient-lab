"""Actual powered-catalyst coefficient maps and matrix absorption arrays."""

from dataclasses import replace
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from .catalyst_compiler import (CompilerLimits, CompilerResourceLimit,
                               CompiledRestriction, compile_absorption, compile_gain_absorption,
                               compile_gain_powered_catalyst,
                               compile_powered_catalyst, coordinate_catalyst,
                               unit_from_nonzero)
from .fields import FiniteField
from .maps import LinearMap
from .schemes import naive_multiply, naive_scheme
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, matrix_multiplication_tensor
from .witnessed_constraints import (DetectorConstraint, WitnessedStateSystem,
                                    assemble_catalyst, clear_rational_dual,
                                    find_nonnegative_dual, rank_upper_constraint)


F2 = FiniteField(2)
F4 = FiniteField(2, (1, 1, 1))


def cross_mixed_catalyst():
    """Swap one D input with an S input, so D is actually used in M2 output."""
    certificate = coordinate_catalyst(D=FiniteTensor.unit(F2))
    maps = []
    for mapping in certificate.maps:
        rows = tuple((row[-1],)+row[1:-1]+(row[0],) for row in mapping.rows)
        maps.append(LinearMap(F2, mapping.input_size, rows))
    result = replace(certificate, maps=tuple(maps))
    result.require_valid()
    return result


class CatalystCompilerTests(unittest.TestCase):
    def test_unit_selector_normalizes_an_extension_coefficient(self):
        S = FiniteTensor(F4, (2, 1, 2), (0, 0, 0, 2))
        maps = unit_from_nonzero(S)
        self.assertEqual(S.restrict(*maps), FiniteTensor.unit(F4))
        with self.assertRaises(ValueError):
            unit_from_nonzero(FiniteTensor.zero(F4, (1, 1, 1)))

    def test_actual_coordinate_certificate_and_b1_repeated_absorption(self):
        certificate = coordinate_catalyst()
        self.assertTrue(certificate.verify())
        self.assertEqual((certificate.d, certificate.k, certificate.m), (2, 9, 1))
        D = TensorScheme.from_tensor(certificate.D)
        S = TensorScheme.from_tensor(certificate.S)
        powered = compile_powered_catalyst(certificate, D)
        self.assertEqual((powered.matrix_size, powered.copies), (2, 9))
        self.assertTrue(powered.restriction.verify())
        seed = naive_scheme(F2, 1, 1, 1)
        for size, terms in ((2, 9), (4, 81)):
            step = compile_absorption(powered, seed, S)
            self.assertEqual((step.matrix_scheme.n, step.matrix_scheme.terms), (size, terms))
            self.assertEqual(step.auxiliary_copies, terms)
            self.assertEqual(step.catalyst_terms, 0)
            self.assertTrue(step.restriction.verify())
            self.assertEqual(step.restriction.target, matrix_multiplication_tensor(F2, size, size, size))
            left = [[(i+j) % 2 for j in range(size)] for i in range(size)]
            right = [[(i+1) % 2 for j in range(size)] for i in range(size)]
            self.assertEqual(step.matrix_scheme.multiply(left, right), naive_multiply(F2, left, right))
            seed = step.matrix_scheme

    def test_actual_b2_powered_maps_and_absorption(self):
        certificate = coordinate_catalyst()
        powered = compile_powered_catalyst(certificate, TensorScheme.from_tensor(certificate.D), 2)
        self.assertEqual((powered.exponent, powered.matrix_size, powered.copies), (2, 4, 81))
        self.assertEqual(powered.restriction.source.shape, (81, 81, 81))
        self.assertEqual(powered.restriction.target.shape, (16, 16, 16))
        self.assertTrue(powered.restriction.verify())
        step = compile_absorption(powered, naive_scheme(F2, 1, 1, 1), TensorScheme.from_tensor(certificate.S))
        self.assertEqual((step.matrix_scheme.n, step.matrix_scheme.terms), (4, 81))
        self.assertTrue(step.matrix_scheme.verify())

    def test_nonempty_catalyst_nontrivial_auxiliary_and_extension_descent(self):
        D = FiniteTensor(F4, (1, 2, 1), (2, 3))
        S = FiniteTensor(F4, (2, 1, 2), (0, 0, 0, 2))
        certificate = coordinate_catalyst(D=D, S=S, m=2)
        self.assertTrue(certificate.verify())
        self.assertEqual(certificate.k, 10)
        Ds = TensorScheme.from_tensor(D).rescale_terms(2, 3)
        Ss = TensorScheme.from_tensor(S).rescale_terms(3, 2)
        powered = compile_powered_catalyst(certificate, Ds)
        seed = naive_scheme(F4, 1, 1, 1)
        step = compile_absorption(powered, seed, Ss)
        self.assertEqual(step.catalyst_terms, Ds.terms)
        self.assertEqual(step.matrix_scheme.terms, Ss.terms*(certificate.k+Ds.terms))
        self.assertTrue(step.restriction.verify())
        left, right = [[1, 2], [3, 1]], [[2, 1], [1, 3]]
        self.assertEqual(step.matrix_scheme.multiply(left, right), naive_multiply(F4, left, right))
        base = step.descend()
        self.assertEqual(base.terms, 4*step.matrix_scheme.terms)
        self.assertTrue(base.verify())
        self.assertEqual(base.field, F2)

    def test_nonempty_powered_catalyst_formula_and_conservative_terms(self):
        D = FiniteTensor.unit(F2)
        certificate = cross_mixed_catalyst()
        powered = compile_powered_catalyst(certificate, TensorScheme.from_tensor(D), 2)
        # D_2 = (D tensor M_2) + nine copies D, with one same S.
        self.assertEqual(powered.catalyst.shape, (13, 13, 13))
        self.assertEqual(powered.catalyst_scheme.terms, 8+9)
        self.assertEqual(powered.restriction.source.shape, (94, 94, 94))
        self.assertEqual(powered.restriction.target.shape, (29, 29, 29))
        self.assertTrue(powered.restriction.verify())
        step = compile_absorption(powered, naive_scheme(F2, 1, 1, 1), TensorScheme.from_tensor(certificate.S))
        self.assertEqual((step.matrix_scheme.n, step.matrix_scheme.terms), (4, 98))
        self.assertTrue(step.matrix_scheme.verify())

    def test_cross_mixed_catalyst_reused_across_eight_seed_terms(self):
        certificate = cross_mixed_catalyst()
        powered = compile_powered_catalyst(certificate, TensorScheme.from_tensor(certificate.D))
        # The core matrix output truly depends on the old catalyst coordinate.
        self.assertTrue(any(row[0] for row in powered.restriction.maps[0].rows[1:]))
        step = compile_absorption(powered, naive_scheme(F2, 2, 2, 2), TensorScheme.from_tensor(certificate.S))
        self.assertEqual(step.catalyst_terms, 1)
        self.assertEqual(step.auxiliary_copies, 8*9+1)
        self.assertEqual((step.matrix_scheme.n, step.matrix_scheme.terms), (4, 73))
        self.assertTrue(step.restriction.verify())
        self.assertTrue(step.matrix_scheme.verify())

    def test_negative_map_control_and_invalid_scheme_rejection(self):
        certificate = coordinate_catalyst()
        Ds = TensorScheme.from_tensor(certificate.D)
        powered = compile_powered_catalyst(certificate, Ds)
        mapping = powered.restriction.maps[0]
        rows = list(mapping.rows)
        rows[0] = (0,)*mapping.input_size
        bad = replace(powered.restriction, maps=(LinearMap(F2, mapping.input_size, rows),)+powered.restriction.maps[1:])
        self.assertFalse(bad.verify())
        with self.assertRaisesRegex(ValueError, "coefficient equality"):
            bad.require_valid()
        with self.assertRaises(ValueError):
            compile_absorption(replace(powered, restriction=bad), naive_scheme(F2, 1, 1, 1), TensorScheme.from_tensor(certificate.S))
        zero_seed = replace(naive_scheme(F2, 1, 1, 1), c=((0,),))
        with self.assertRaises(ValueError):
            compile_absorption(powered, zero_seed, TensorScheme.from_tensor(certificate.S))
        with self.assertRaises(ValueError):
            compile_powered_catalyst(certificate, TensorScheme.from_tensor(certificate.S))
        for exponent in (0, -1, True, 1.0):
            with self.assertRaises(ValueError):
                compile_powered_catalyst(certificate, Ds, exponent)

    def test_dense_guards_fire_before_power_and_repeat_allocations(self):
        certificate = coordinate_catalyst()
        Ds = TensorScheme.from_tensor(certificate.D)
        with patch.object(FiniteTensor, "direct_sum") as allocations:
            with self.assertRaisesRegex(CompilerResourceLimit, "incomplete"):
                compile_powered_catalyst(certificate, Ds, 4)
            allocations.assert_not_called()
        with patch("reference.catalyst_compiler.naive_scheme") as seeds:
            with patch.object(FiniteTensor, "direct_sum") as allocations:
                with self.assertRaisesRegex(CompilerResourceLimit, "incomplete"):
                    compile_powered_catalyst(certificate, Ds, 10**9)
                seeds.assert_not_called()
                allocations.assert_not_called()
        powered = compile_powered_catalyst(certificate, Ds)
        with patch.object(FiniteTensor, "direct_sum") as allocations:
            with self.assertRaises(CompilerResourceLimit):
                compile_absorption(powered, naive_scheme(F2, 4, 4, 4), TensorScheme.from_tensor(certificate.S))
            allocations.assert_not_called()
        with self.assertRaises(CompilerResourceLimit):
            coordinate_catalyst(limits=CompilerLimits(max_tensor_coefficients=100))
        with self.assertRaises(CompilerResourceLimit):
            compile_absorption(powered, naive_scheme(F2, 1, 1, 1), TensorScheme.from_tensor(certificate.S),
                               limits=CompilerLimits(max_scheme_scalars=1))


class GainAbsorptionTests(unittest.TestCase):
    def test_gain_retention_recovers_classical_lengths_and_actual_arrays(self):
        certificate = coordinate_catalyst()
        D, S = TensorScheme.from_tensor(certificate.D), TensorScheme.from_tensor(certificate.S)
        seed = naive_scheme(F2, 1, 1, 1)
        for size, terms, gains in ((2, 8, 1), (4, 64, 8)):
            result = compile_gain_absorption(certificate, seed, S, D)
            self.assertEqual((result.matrix_scheme.n, result.matrix_scheme.terms), (size, terms))
            self.assertEqual(result.gain_removed, gains)
            self.assertEqual(result.raw_gain_scheme.terms, terms+gains)
            self.assertEqual(result.raw_gain_scheme.target,
                             FiniteTensor.direct_sum((result.matrix_scheme.as_tensor_scheme().target,)+(FiniteTensor.unit(F2),)*gains))
            self.assertTrue(result.restriction.verify())
            self.assertTrue(result.raw_gain_scheme.verify())
            left = [[(i+j) % 2 for j in range(size)] for i in range(size)]
            right = [[int(i == j) for j in range(size)] for i in range(size)]
            self.assertEqual(result.matrix_scheme.multiply(left, right), naive_multiply(F2, left, right))
            seed = result.matrix_scheme

    def test_cross_mixed_catalyst_cost_is_paid_once_without_promotion(self):
        certificate = cross_mixed_catalyst()
        D, S = TensorScheme.from_tensor(certificate.D), TensorScheme.from_tensor(certificate.S)
        result = compile_gain_absorption(certificate, naive_scheme(F2, 2, 2, 2), S, D)
        self.assertEqual((result.catalyst_terms, result.gain_removed), (1, 8))
        self.assertEqual((result.matrix_scheme.n, result.matrix_scheme.terms), (4, 65))
        self.assertEqual(result.raw_gain_scheme.terms, 73)
        self.assertEqual(result.restriction.source.shape, (73, 73, 73))
        self.assertTrue(result.matrix_scheme.verify())

    def test_two_gains_in_characteristic_two_are_both_eliminated(self):
        certificate = coordinate_catalyst(m=2)
        D, S = TensorScheme.from_tensor(certificate.D), TensorScheme.from_tensor(certificate.S)
        result = compile_gain_absorption(certificate, naive_scheme(F2, 1, 1, 1), S, D)
        self.assertEqual(F2.embed(result.gain_removed), 0)
        self.assertEqual(result.gain_removed, 2)
        self.assertEqual(result.raw_gain_scheme.shape, (6, 6, 6))
        self.assertEqual(result.raw_gain_scheme.terms, 10)
        self.assertEqual(result.matrix_scheme.terms, 8)
        self.assertTrue(result.matrix_scheme.verify())

    def test_nontrivial_extension_auxiliary_does_not_multiply_catalyst_cost(self):
        D = FiniteTensor(F4, (1, 2, 1), (2, 3))
        S = FiniteTensor(F4, (2, 1, 2), (2, 0, 0, 3))
        certificate = coordinate_catalyst(D=D, S=S)
        Ds, Ss = TensorScheme.from_tensor(D), TensorScheme.from_tensor(S)
        result = compile_gain_absorption(certificate, naive_scheme(F4, 1, 1, 1), Ss, Ds)
        self.assertEqual((Ds.terms, Ss.terms), (2, 2))
        self.assertEqual(result.raw_gain_scheme.terms, 2+9*2)
        self.assertEqual(result.matrix_scheme.terms, 2+(9*2-1))
        self.assertLess(result.matrix_scheme.terms, 2*(9+2))
        self.assertTrue(result.matrix_scheme.verify())
        self.assertEqual(result.descend().terms, 4*19)
        self.assertTrue(result.descend().verify())

    def test_actual_rank_seven_lp_catalyst_retains_gain(self):
        artifact = Path(__file__).with_name("research") / "unseeded-f2-rank7.json"
        data = json.loads(artifact.read_text())
        M = matrix_multiplication_tensor(F2, 2, 2, 2)
        scheme = TensorScheme(M, data["a"], data["b"], data["c"])
        system = WitnessedStateSystem(2, 8, (("unit", FiniteTensor.unit(F2)), ("M", M)),
                                      (rank_upper_constraint("unit", "M", scheme), DetectorConstraint("unit", "M")))
        dual = find_nonnegative_dual(system)
        self.assertEqual(dual.status, "found")
        certificate = assemble_catalyst(clear_rational_dual(system, dual.coefficients))
        result = compile_gain_absorption(certificate, naive_scheme(F2, 1, 1, 1),
                                         TensorScheme.from_tensor(certificate.S), TensorScheme.from_tensor(certificate.D))
        self.assertEqual(result.catalyst_terms, 7)
        self.assertEqual(result.raw_gain_scheme.terms, 15)
        self.assertEqual(result.matrix_scheme.terms, 14)
        from .extraction import simplify_scheme
        compact = simplify_scheme(result.matrix_scheme.as_tensor_scheme())
        self.assertEqual(compact.terms, 7)
        self.assertTrue(compact.verify())

    def test_bad_gain_maps_bad_decompositions_and_dense_limits_are_rejected(self):
        certificate = coordinate_catalyst()
        D, S = TensorScheme.from_tensor(certificate.D), TensorScheme.from_tensor(certificate.S)
        mapping = certificate.maps[0]
        rows = list(mapping.rows)
        rows[0] = (0,)*mapping.input_size
        bad = replace(certificate, maps=(LinearMap(F2, mapping.input_size, rows),)+certificate.maps[1:])
        with self.assertRaises(ValueError):
            compile_gain_absorption(bad, naive_scheme(F2, 1, 1, 1), S, D)
        with self.assertRaises(ValueError):
            compile_gain_absorption(certificate, naive_scheme(F2, 1, 1, 1), replace(S, c=((0,),)), D)
        with patch.object(FiniteTensor, "direct_sum") as allocations:
            with self.assertRaisesRegex(CompilerResourceLimit, "incomplete"):
                compile_gain_absorption(certificate, naive_scheme(F2, 4, 4, 4), S, D)
            allocations.assert_not_called()
        with self.assertRaises(CompilerResourceLimit):
            compile_gain_absorption(certificate, naive_scheme(F2, 1, 1, 1), S, D,
                                     limits=CompilerLimits(max_scheme_scalars=64))


class PoweredGainTests(unittest.TestCase):
    def test_positive_power_b2_retains_eleven_independent_units(self):
        certificate = coordinate_catalyst()
        D, S = TensorScheme.from_tensor(certificate.D), TensorScheme.from_tensor(certificate.S)
        for exponent in (1, 2):
            powered = compile_gain_powered_catalyst(certificate, D, exponent)
            self.assertEqual(powered.exponent, exponent)
            self.assertEqual(powered.matrix_size, 2**exponent)
            self.assertEqual(powered.copies, 9**exponent)
            self.assertEqual(powered.gain, sum(9**(exponent-1-i)*2**i for i in range(exponent)))
            self.assertTrue(powered.certificate.verify())
            self.assertTrue(powered.restriction.verify())
        self.assertEqual(powered.gain, 11)
        self.assertEqual(powered.restriction.source.shape, (81, 81, 81))
        self.assertEqual(powered.restriction.target.shape, (27, 27, 27))
        result = compile_gain_absorption(powered, naive_scheme(F2, 1, 1, 1), S, powered.catalyst_scheme)
        self.assertEqual((result.matrix_scheme.n, result.matrix_scheme.terms), (4, 70))
        self.assertEqual(result.gain_removed, 11)
        self.assertTrue(result.matrix_scheme.verify())

    def test_cross_mixed_nonempty_powered_gain_catalyst_and_absorption(self):
        certificate = cross_mixed_catalyst()
        powered = compile_gain_powered_catalyst(certificate, TensorScheme.from_tensor(certificate.D), 2)
        self.assertEqual(powered.catalyst_scheme.terms, 17)
        self.assertEqual(powered.catalyst.shape, (13, 13, 13))
        self.assertEqual(powered.gain, 11)
        self.assertTrue(powered.certificate.verify())
        result = compile_gain_absorption(powered, naive_scheme(F2, 1, 1, 1),
                                         TensorScheme.from_tensor(certificate.S), powered.catalyst_scheme)
        self.assertEqual(result.raw_gain_scheme.terms, 98)
        self.assertEqual(result.matrix_scheme.terms, 87)
        self.assertTrue(result.matrix_scheme.verify())

    def test_extension_coefficients_and_characteristic_two_gain_counts(self):
        S = FiniteTensor(F4, (1, 1, 1), (2,))
        D = FiniteTensor(F4, (1, 1, 1), (3,))
        certificate = coordinate_catalyst(D=D, S=S, m=2)
        powered = compile_gain_powered_catalyst(certificate, TensorScheme.from_tensor(D), 2)
        self.assertEqual(powered.gain, 2*(10+2))
        self.assertEqual(F4.embed(powered.gain), 0)
        self.assertEqual(powered.catalyst_scheme.terms, 8+10)
        self.assertTrue(powered.certificate.verify())
        result = compile_gain_absorption(powered, naive_scheme(F4, 1, 1, 1),
                                         TensorScheme.from_tensor(S), powered.catalyst_scheme)
        self.assertEqual(result.matrix_scheme.terms, 18+100-24)
        self.assertTrue(result.matrix_scheme.verify())
        self.assertTrue(result.descend().verify())
        self.assertEqual(result.descend().terms, 4*94)

    def test_diagonal_subrank_selector_uses_matrix_diagonal_on_all_legs(self):
        matrix = matrix_multiplication_tensor(F4, 3, 3, 3)
        diagonal = LinearMap.selection(F4, 9, (0, 4, 8))
        target = FiniteTensor.direct_sum((FiniteTensor.unit(F4),)*3)
        self.assertEqual(matrix.restrict(diagonal, diagonal, diagonal), target)
        wrong = LinearMap.selection(F4, 9, (0, 1, 2))
        self.assertNotEqual(matrix.restrict(diagonal, wrong, diagonal), target)

    def test_powered_positive_maps_preserve_unequal_catalyst_axis_widths(self):
        D = FiniteTensor(F4, (1, 2, 1), (2, 3))
        S = FiniteTensor(F4, (1, 1, 1), (2,))
        certificate = coordinate_catalyst(D=D, S=S)
        powered = compile_gain_powered_catalyst(certificate, TensorScheme.from_tensor(D), 2)
        self.assertEqual(powered.catalyst.shape, (13, 26, 13))
        self.assertEqual(powered.catalyst_scheme.terms, 34)
        self.assertEqual(powered.source.shape, (94, 107, 94))
        self.assertEqual(powered.target.shape, (40, 53, 40))
        self.assertTrue(powered.restriction.verify())
        result = compile_gain_absorption(powered, naive_scheme(F4, 1, 1, 1),
                                         TensorScheme.from_tensor(S), powered.catalyst_scheme)
        self.assertEqual(result.matrix_scheme.terms, 34+81-11)
        self.assertTrue(result.matrix_scheme.verify())

    def test_bad_maps_and_huge_exponents_rejected_before_power_allocations(self):
        certificate = coordinate_catalyst()
        D = TensorScheme.from_tensor(certificate.D)
        with patch("reference.catalyst_compiler.naive_scheme") as seeds:
            with patch.object(FiniteTensor, "direct_sum") as allocations:
                with self.assertRaisesRegex(CompilerResourceLimit, "incomplete"):
                    compile_gain_powered_catalyst(certificate, D, 10**9)
                seeds.assert_not_called()
                allocations.assert_not_called()
        powered = compile_gain_powered_catalyst(certificate, D, 2)
        mapping = powered.certificate.maps[0]
        rows = list(mapping.rows)
        rows[0] = (0,)*mapping.input_size
        bad = replace(powered.certificate, maps=(LinearMap(F2, mapping.input_size, rows),)+powered.certificate.maps[1:])
        self.assertFalse(bad.verify())
        with self.assertRaises(ValueError):
            compile_gain_absorption(replace(powered, certificate=bad), naive_scheme(F2, 1, 1, 1),
                                     TensorScheme.from_tensor(certificate.S), powered.catalyst_scheme)
        for exponent in (0, -1, True, 1.0):
            with self.assertRaises(ValueError):
                compile_gain_powered_catalyst(certificate, D, exponent)


if __name__ == "__main__":
    unittest.main()
