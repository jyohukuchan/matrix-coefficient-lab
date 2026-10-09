"""Shared-dot quotient maps and full formal border certificates over base fields."""

from dataclasses import replace
from itertools import product
import unittest
from unittest.mock import patch

from .catalyst_compiler import CompilerLimits, CompilerResourceLimit
from .fields import FiniteField
from .maps import LinearMap
from .polynomials import Polynomial
from .shared_dot_orders import (export_shared_dot_quotient_order, export_shared_dot_subrank_order,
                                recover_shared_dot_power, require_shared_dot_identity,
                                shared_dot_degeneration, shared_dot_scheme, shared_dot_tensor,
                                shared_dot_w_degeneration)
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor
from .witnessed_constraints import StateAssemblyLimit, StateAssemblyLimits


FIELDS = (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)), FiniteField(5), FiniteField(11))


class SharedDotOrderTests(unittest.TestCase):
    def test_two_term_w_family_is_formal_in_every_characteristic(self):
        for field in FIELDS:
            degeneration = shared_dot_w_degeneration(field)
            self.assertEqual((degeneration.terms, degeneration.leading, degeneration.degree_bound), (2, 1, 3))
            self.assertEqual(degeneration.target, shared_dot_tensor(field, 1))
            degeneration.require_valid()
            self.assertEqual(degeneration.tensor_polynomial(0, 0, 0), Polynomial.monomial(field, 3))
            self.assertEqual(degeneration.tensor_polynomial(1, 1, 1), Polynomial(field))
        with self.assertRaises(CompilerResourceLimit):
            shared_dot_w_degeneration(FIELDS[0], limits=CompilerLimits(max_terms=1))

    def test_actual_low_term_schemes_and_field_dependent_cancellations(self):
        for field in FIELDS:
            for q in (1, 2, 3, 4):
                with self.subTest(field=field, q=q):
                    scheme = shared_dot_scheme(field, q)
                    expected = 2*q+1 if field.p == 2 else 2*q+int(field.embed(q) != 0)
                    self.assertEqual(scheme.terms, expected)
                    self.assertEqual(scheme.target, shared_dot_tensor(field, q))
                    scheme.require_exact()
                    left = tuple((i+1) % field.order for i in range(q+1))
                    right = tuple((2*i+3) % field.order for i in range(q+1))
                    self.assertEqual(scheme.apply(left, right), scheme.target.contract(left, right))
        self.assertEqual(shared_dot_scheme(FIELDS[1], 3).terms, 6)
        self.assertEqual(shared_dot_scheme(FIELDS[3], 3).terms, 7)
        self.assertEqual(shared_dot_scheme(FIELDS[0], 4).terms, 9)
        self.assertEqual(shared_dot_scheme(FIELDS[2], 4).terms, 9)

    def test_characteristic_two_formula_cannot_be_reused_in_odd_characteristic(self):
        for q in (1, 2, 3, 4):
            binary = shared_dot_scheme(FIELDS[0], q)
            extension = TensorScheme(shared_dot_tensor(FIELDS[2], q), binary.a, binary.b, binary.c)
            extension.require_exact()
            odd = TensorScheme(shared_dot_tensor(FIELDS[1], q), binary.a, binary.b, binary.c)
            self.assertFalse(odd.verify())
            broken = TensorScheme(binary.target, binary.a[:-1], binary.b[:-1], binary.c[:-1])
            self.assertFalse(broken.verify())
        characteristic_three = shared_dot_scheme(FIELDS[1], 3)
        wrong_scalar_field = TensorScheme(shared_dot_tensor(FIELDS[3], 3),
                                          characteristic_three.a, characteristic_three.b, characteristic_three.c)
        self.assertFalse(wrong_scalar_field.verify())

    def test_actual_quotients_all_coordinates_over_all_control_fields(self):
        for field in FIELDS:
            for q in (1, 2, 3, 4):
                with self.subTest(field=field, q=q):
                    exported = export_shared_dot_quotient_order(field, q, "mixed", "shared")
                    exported.require_valid()
                    source, target = exported.registry["mixed"], exported.registry["shared"]
                    self.assertEqual(source.shape, (2*q+1,)*3)
                    self.assertEqual(target.shape, (q+1,)*3)
                    self.assertEqual(sum(bool(x) for x in target.coefficients), 3*q)
                    self.assertEqual(source.restrict(*exported.row.maps), target)
                    self.assertEqual(exported.row.positive, ("mixed",))
                    self.assertEqual(exported.row.negative, ("shared",))

    def test_quotient_missing_summand_and_wrong_coordinate_fail(self):
        field, q = FIELDS[2], 3
        exported = export_shared_dot_quotient_order(field, q, "mixed", "shared")
        rows = list(exported.row.maps[0].rows)
        rows[1] = tuple(int(i == 1) for i in range(2*q+1))
        missing = LinearMap(field, 2*q+1, rows)
        self.assertFalse(replace(exported, row=replace(exported.row,
            maps=(missing,)+exported.row.maps[1:])).verify())
        rows = list(exported.row.maps[2].rows)
        rows[0] = tuple(int(i == 2*q-1) for i in range(2*q+1))
        wrong = LinearMap(field, 2*q+1, rows)
        self.assertFalse(replace(exported, row=replace(exported.row,
            maps=exported.row.maps[:2]+(wrong,))).verify())

    def test_q_plus_two_terms_have_full_formal_degree_three_identity(self):
        for field in FIELDS:
            for q in (1, 2, 3, 4):
                degeneration = shared_dot_degeneration(field, q)
                self.assertEqual((degeneration.terms, degeneration.leading, degeneration.degree_bound), (q+2, 3, 6))
                self.assertEqual(degeneration.recovery_node_count, 4)
                self.assertEqual(degeneration.target, shared_dot_tensor(field, q))
                require_shared_dot_identity(degeneration)
                for coordinate in product(range(q+1), repeat=3):
                    polynomial = degeneration.tensor_polynomial(*coordinate)
                    self.assertEqual(tuple(polynomial.coefficient(i) for i in range(3)), (0, 0, 0))
                    self.assertEqual(polynomial.coefficient(3), degeneration.target.coefficient(*coordinate))
                    self.assertEqual(polynomial.coefficient(5), 0)
                    self.assertLessEqual(polynomial.degree, 6)
                # Noise is retained and explicitly verified, not assumed zero.
                self.assertEqual(degeneration.tensor_polynomial(0, 1, 1), Polynomial(field, (0, 0, 0, 1, field.neg(1))))
                self.assertEqual(degeneration.tensor_polynomial(1, 1, 1), Polynomial(field, (0, 0, 0, 0, 1, 0, field.neg(1))))
                if q >= 2:
                    self.assertEqual(degeneration.tensor_polynomial(0, 1, 2), Polynomial.monomial(field, 4, field.neg(1)))

    def test_characteristic_dividing_q_uses_integer_embedding(self):
        for field, q in ((FIELDS[0], 2), (FIELDS[0], 4), (FIELDS[1], 3), (FIELDS[2], 2), (FIELDS[2], 4)):
            degeneration = shared_dot_degeneration(field, q)
            self.assertEqual(field.embed(q), 0)
            self.assertEqual(degeneration.a[-1][0], Polynomial.constant(field, 1))
            degeneration.require_valid()
        extension = shared_dot_degeneration(FIELDS[2], 3)
        self.assertEqual(extension.a[-1][0].coefficient(1), 1)
        self.assertNotEqual(FIELDS[2].embed(3), 3)  # 3 encodes an extension element

    def test_corrupt_low_coefficients_and_noise_are_rejected_formally(self):
        field = FIELDS[3]
        degeneration = shared_dot_degeneration(field, 2)
        a = list(degeneration.a)
        a[-1] = (a[-1][0]+Polynomial.monomial(field, 2),)+a[-1][1:]
        bad_low = replace(degeneration, a=a)
        self.assertFalse(bad_low.verify())
        with self.assertRaises(ValueError):
            require_shared_dot_identity(bad_low)
        a[-1] = (degeneration.a[-1][0]+Polynomial.monomial(field, 7),)+degeneration.a[-1][1:]
        bad_noise = replace(degeneration, a=a)
        self.assertTrue(bad_noise.verify())  # valid generic leading coefficient, wrong advertised full identity
        with self.assertRaisesRegex(ValueError, "full formal identity"):
            require_shared_dot_identity(bad_noise)
        field = FIELDS[0]
        degeneration = shared_dot_degeneration(field, 2)
        a = list(degeneration.a)
        a[-1] = (a[-1][0]+Polynomial(field, (0, 1, 1)),)+a[-1][1:]
        forged = replace(degeneration, a=a)
        for node in range(field.order):
            self.assertEqual(forged.evaluate_at(node), degeneration.evaluate_at(node))
        with self.assertRaises(ValueError):
            require_shared_dot_identity(forged)

    def test_recovery_requires_available_nonzero_nodes_without_extension(self):
        for field in (FIELDS[0], FIELDS[1], FIELDS[2]):
            with patch("reference.shared_dot_orders.shared_dot_degeneration") as constructor:
                with self.assertRaises(ValueError):
                    recover_shared_dot_power(field, 3)
                constructor.assert_not_called()
        with self.assertRaises(ValueError):
            recover_shared_dot_power(FIELDS[4], 3, nodes=(1, 2, 3))
        for nodes in ((0, 1, 2, 3), (1, 1, 2, 3)):
            with self.assertRaises(ValueError):
                recover_shared_dot_power(FIELDS[4], 3, nodes=nodes)

    def test_four_node_exact_recovery_over_f5_and_f11(self):
        for field in (FIELDS[3], FIELDS[4]):
            for q in (1, 2, 3, 4):
                recovery = recover_shared_dot_power(field, q)
                self.assertEqual(len(recovery.nodes), 4)
                self.assertEqual(recovery.scheme.terms, 4*(q+2))
                self.assertEqual(recovery.scheme.field, field)
                self.assertEqual(recovery.scheme.target, shared_dot_tensor(field, q))
                recovery.scheme.require_exact()
                left = tuple((i+2) % field.order for i in range(q+1))
                right = tuple((2*i+1) % field.order for i in range(q+1))
                self.assertEqual(recovery.scheme.apply(left, right), recovery.scheme.target.contract(left, right))

    def test_q3_square_recovery_and_complete_powered_polynomial(self):
        field, q = FIELDS[4], 3
        base = shared_dot_degeneration(field, q)
        recovery = recover_shared_dot_power(field, q, 2)
        self.assertEqual((recovery.degeneration.terms, recovery.degeneration.leading, recovery.degeneration.degree_bound), (25, 6, 12))
        self.assertEqual(len(recovery.nodes), 7)
        self.assertEqual(recovery.scheme.terms, 175)
        self.assertEqual(recovery.scheme.target.shape, (16,)*3)
        self.assertEqual(recovery.scheme.target, shared_dot_tensor(field, q).tensor_power(2))
        recovery.scheme.require_exact()
        for coordinate in product(range(16), repeat=3):
            pairs = tuple(divmod(index, 4) for index in coordinate)
            expected = (base.tensor_polynomial(*(pair[0] for pair in pairs))
                        *base.tensor_polynomial(*(pair[1] for pair in pairs)))
            self.assertEqual(recovery.degeneration.tensor_polynomial(*coordinate), expected)

    def test_actual_three_unit_lower_rows_for_q_at_least_three(self):
        for field in FIELDS:
            for q in (3, 4):
                exported = export_shared_dot_subrank_order(field, q, "shared")
                self.assertEqual(exported.row.negative, ("unit",)*3)
                exported.require_valid()
                diagonal = exported.registry["shared"].restrict(*exported.row.maps)
                self.assertEqual(diagonal, FiniteTensor.from_function(field, (3,)*3,
                                                                      lambda i, j, k: int(i == j == k)))
        with self.assertRaises(ValueError):
            export_shared_dot_subrank_order(FIELDS[0], 2, "shared")

    def test_resource_caps_precede_polynomial_tensor_and_power_allocations(self):
        field = FIELDS[4]
        degeneration = shared_dot_degeneration(field, 2)
        with patch("reference.shared_dot_orders.TensorDegeneration.tensor_polynomial", side_effect=AssertionError("degree verification")):
            with self.assertRaises(CompilerResourceLimit):
                require_shared_dot_identity(degeneration, limits=CompilerLimits(max_scheme_scalars=1))
        with patch("reference.shared_dot_orders.FiniteTensor.from_function", side_effect=AssertionError("allocation")):
            for operation in (
                lambda: shared_dot_tensor(field, 10**9),
                lambda: export_shared_dot_quotient_order(field, 2, "s", "c", limits=StateAssemblyLimits(max_map_entries=1)),
                lambda: shared_dot_degeneration(field, 1, limits=CompilerLimits(max_scheme_scalars=181)),
                lambda: shared_dot_degeneration(field, 10**9),
                lambda: shared_dot_scheme(field, 10**9),
                lambda: shared_dot_scheme(field, 3, limits=CompilerLimits(max_terms=6)),
                lambda: recover_shared_dot_power(field, 3, 10**9),
                lambda: recover_shared_dot_power(field, 3, 2, limits=CompilerLimits(max_terms=100)),
            ):
                with self.assertRaises((CompilerResourceLimit, StateAssemblyLimit)):
                    operation()
        with patch("reference.shared_dot_orders.shared_dot_degeneration", side_effect=AssertionError("base allocation")):
            recovery = recover_shared_dot_power(field, 10**9, 0)
        self.assertEqual(recovery.scheme.target, FiniteTensor.unit(field))
        self.assertEqual(recovery.scheme.terms, 1)
        self.assertEqual(shared_dot_scheme(FIELDS[1], 3, limits=CompilerLimits(max_terms=6)).terms, 6)

    def test_invalid_inputs_and_tensor_target(self):
        for q in (0, -1, True, 1.0):
            for function in (shared_dot_tensor, shared_dot_degeneration, shared_dot_scheme, recover_shared_dot_power):
                with self.assertRaises(ValueError):
                    function(FIELDS[0], q)
        with self.assertRaises(ValueError):
            shared_dot_tensor("F2", 2)
        for exponent in (-1, True, 1.0):
            with self.assertRaises(ValueError):
                recover_shared_dot_power(FIELDS[4], 2, exponent)
        for function in (shared_dot_tensor, shared_dot_degeneration):
            with self.assertRaises(ValueError):
                function(FIELDS[4], 2, limits=None)
        with self.assertRaises(ValueError):
            export_shared_dot_quotient_order(FIELDS[0], 2, "same", "same")
        bad = replace(shared_dot_degeneration(FIELDS[0], 2), target=FiniteTensor.zero(FIELDS[0], (3,)*3))
        with self.assertRaises(ValueError):
            require_shared_dot_identity(bad)


if __name__ == "__main__":
    unittest.main()
