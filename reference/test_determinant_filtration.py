"""Source reconstruction, formal determinant filtration, and ordinary rows."""

from dataclasses import replace
from itertools import product
import unittest
from unittest.mock import patch

from .convolution import convolution_scheme
from .determinant_filtration import (DeterminantFiltration, adapted_coordinates,
                                     adapted_vectors)
from .fields import FiniteField
from .maps import LinearMap
from .polynomials import Polynomial
from .sectors import convolution_tensor
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor
from .witnessed_constraints import (StateAssemblyLimit, StateAssemblyLimits,
                                    WitnessedStateSystem, find_nonnegative_dual)


F2 = FiniteField(2)
F3 = FiniteField(3)
F4 = FiniteField(2, (1, 1, 1))
F9 = FiniteField(3, (1, 0, 1))


def source_scheme_small(field):
    """Actual six-term evaluation source C(2,2) tensor dot2, reindexed by tag."""
    convolution = convolution_scheme(field, 2, 2)
    dot = TensorScheme.from_tensor(FiniteTensor.from_function(field, (1, 2, 2), lambda i, j, k: int(j == k)))
    generic = convolution.tensor_product(dot)
    return generic.restrict(LinearMap.identity(field, 2), LinearMap.selection(field, 4, (0, 2, 1, 3)),
                            LinearMap.selection(field, 6, (0, 2, 4, 1, 3, 5)))


class DeterminantFiltrationTests(unittest.TestCase):
    def test_explicit_basis_vectors_and_inverse_coordinates(self):
        f = F3
        self.assertEqual(adapted_vectors(f, 1).rows,
                         ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 0, 1), (0, 2, 1, 0)))
        self.assertEqual(adapted_coordinates(f, 1).rows,
                         ((1, 0, 0, 0), (0, 1, 1, 0), (0, 0, 0, 1), (0, 0, 1, 0)))
        self.assertEqual(adapted_vectors(F2, 1).rows[-1], (0, 1, 1, 0))
        for field in (F2, F3, F4, F9):
            for degree in range(4):
                basis = adapted_vectors(field, degree)
                columns = LinearMap(field, basis.output_size,
                                    tuple(tuple(row[i] for row in basis.rows) for i in range(basis.input_size)))
                coordinates = adapted_coordinates(field, degree)
                identity = LinearMap.identity(field, 2*(degree+1))
                self.assertEqual(coordinates.compose(columns), identity)
                self.assertEqual(columns.compose(coordinates), identity)

    def test_grid_source_basis_change_and_formal_polynomial_every_coordinate(self):
        for field, d, e in product((F2, F3, F4, F9), range(3), range(3)):
            filtration = DeterminantFiltration(field, d, e)
            self.assertTrue(filtration.verify())
            self.assertEqual(filtration.source.shape[0], d+1)
            self.assertEqual(filtration.source.restrict(*filtration.basis_maps), filtration.adapted)
            self.assertEqual(sum(filtration.cross.coefficients), d)
            for coordinate in product(*(range(size) for size in filtration.shape)):
                self.assertEqual(filtration.polynomial_coefficient(*coordinate),
                                 Polynomial(field, (0, filtration.graded.coefficient(*coordinate),
                                                     filtration.cross.coefficient(*coordinate))))

    def test_two_graded_convolution_branches_share_exact_same_first_leg(self):
        for d, e in product(range(3), range(3)):
            filtration = DeterminantFiltration(F4, d, e)
            n = d+e
            quotient = filtration.graded.restrict(LinearMap.identity(F4, d+1),
                         LinearMap.selection(F4, filtration.shape[1], range(e+2)),
                         LinearMap.selection(F4, filtration.shape[2], range(n+2)))
            kernel = filtration.graded.restrict(LinearMap.identity(F4, d+1),
                       LinearMap.selection(F4, filtration.shape[1], range(e+2, 2*(e+1))),
                       LinearMap.selection(F4, filtration.shape[2], range(n+2, 2*(n+1))))
            self.assertEqual(quotient, convolution_tensor(F4, d+1, e+2))
            self.assertEqual(kernel, convolution_tensor(F4, d+1, e))
            self.assertEqual(quotient.shape[0], kernel.shape[0])

    def test_real_evaluation_source_scheme_generates_and_recovers_graded(self):
        for field in (F3, F4, F9):
            filtration = DeterminantFiltration(field, 1, 1)
            source = source_scheme_small(field)
            self.assertEqual(source.target, filtration.source)
            self.assertEqual(source.terms, 6)
            if field.degree > 1:
                source = source.rescale_terms(field.p, field.p+1)
            degeneration = filtration.generate_degeneration(source)
            self.assertEqual((degeneration.leading, degeneration.degree_bound, degeneration.terms), (1, 2, 6))
            self.assertTrue(degeneration.verify())
            recovered = degeneration.recover_scheme()
            self.assertEqual(recovered.terms, 12)
            self.assertEqual(recovered.target, filtration.graded)
            left = tuple((i+1) % field.order for i in range(2))
            right = tuple((j+2) % field.order for j in range(4))
            self.assertEqual(recovered.apply(left, right), filtration.graded.contract(left, right))
            for coordinate in product(range(2), range(4), range(6)):
                self.assertEqual(degeneration.tensor_polynomial(*coordinate), filtration.polynomial_coefficient(*coordinate))

    def test_interpolated_ordinary_rows_and_extra_nodes(self):
        for field, nodes in ((F3, None), (F4, None), (F9, (1, 3, 7))):
            filtration = DeterminantFiltration(field, 1, 1)
            export = filtration.export_order("source", "graded", nodes=nodes)
            self.assertTrue(export.verify())
            source = FiniteTensor.direct_sum((filtration.source,)*len(export.nodes))
            self.assertEqual(source.restrict(*export.row.maps), filtration.graded)
            self.assertEqual(export.row.positive, ("source",)*len(export.nodes))
            self.assertEqual(export.row.negative, ("graded",))
            system = WitnessedStateSystem(2, 5, (("unit", FiniteTensor.unit(field)),)+export.inventory, (export.row,))
            system.require_valid()
            self.assertEqual(find_nonnegative_dual(system).status, "finite_exhausted")

    def test_wrong_output_dual_wrong_kernel_sign_and_missing_cross_are_rejected(self):
        class WrongDual(DeterminantFiltration):
            @property
            def basis_maps(self):
                first, middle, _ = super().basis_maps
                return first, middle, adapted_vectors(self.field, self.d+self.e)
        class WrongSign(DeterminantFiltration):
            @property
            def basis_maps(self):
                first, middle, last = super().basis_maps
                rows = tuple(tuple(1 if value == self.field.neg(1) else value for value in row)
                             for row in middle.rows)
                return first, LinearMap(self.field, middle.input_size, rows), last
        class MissingCross(DeterminantFiltration):
            @property
            def adapted(self):
                return self.graded
        for cls in (WrongDual, WrongSign, MissingCross):
            wrong = cls(F3, 1, 1)
            self.assertFalse(wrong.verify())
            with self.assertRaisesRegex(ValueError, "determinant filtration"):
                wrong.require_valid()

    def test_missing_interpolation_normalization_is_detected(self):
        filtration = DeterminantFiltration(F4, 1, 1)
        export = filtration.export_order("source", "graded", nodes=(1, 2))
        old = export.row.maps[2]
        width = filtration.shape[2]
        rows = tuple(tuple(F4.mul(value, export.nodes[j//width]) for j, value in enumerate(row)) for row in old.rows)
        bad = replace(export, row=replace(export.row, maps=export.row.maps[:2]+(LinearMap(F4, old.input_size, rows),)))
        self.assertFalse(bad.verify())

    def test_bad_source_nodes_parameters_and_resource_preflights(self):
        filtration = DeterminantFiltration(F4, 1, 1)
        for nodes in ((1,), (0, 1), (1, 1), (1, 4)):
            with self.assertRaises(ValueError):
                filtration.export_order("source", "graded", nodes=nodes)
        with self.assertRaises(ValueError):
            DeterminantFiltration(F2, 1, 1).export_order("source", "graded")
        source = source_scheme_small(F4)
        bad = replace(source, c=((0,)*6,)*source.terms)
        for scheme in (bad, TensorScheme.from_tensor(filtration.graded)):
            with self.assertRaises(ValueError):
                filtration.generate_degeneration(scheme)
        for d, e in ((-1, 1), (1, True), (1.0, 1)):
            with self.assertRaises(ValueError):
                DeterminantFiltration(F4, d, e)
        with patch.object(FiniteTensor, "from_function") as allocation:
            with self.assertRaisesRegex(StateAssemblyLimit, "incomplete"):
                DeterminantFiltration(F4, 10**9, 10**9)
            allocation.assert_not_called()
        with self.assertRaises(StateAssemblyLimit):
            DeterminantFiltration(F4, 1, 1, StateAssemblyLimits(max_map_entries=1))
        with self.assertRaises(StateAssemblyLimit):
            filtration.export_order("source", "graded", limits=StateAssemblyLimits(max_blocks=1))
        with self.assertRaises(StateAssemblyLimit):
            filtration.export_order("source", "graded", limits=StateAssemblyLimits(max_tensor_entries=100))


if __name__ == "__main__":
    unittest.main()
