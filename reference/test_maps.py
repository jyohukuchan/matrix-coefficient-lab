"""Independent contraction, composition, direct sums, and shared branch tags."""

from dataclasses import replace
from itertools import permutations, product
from random import Random
import unittest

from . import (FiniteField, FiniteTensor, LinearMap, TensorScheme, TensorDegeneration,
               ThreeSectorConstruction, shared_first_tensor, tag_shared_scheme)


class MapTests(unittest.TestCase):
    def test_general_local_maps_match_full_independent_contraction(self):
        random = Random(12321)
        for field in (FiniteField(2), FiniteField(5), FiniteField(2, (1, 1, 1))):
            target = FiniteTensor(field, (2, 3, 2), tuple(random.randrange(field.order) for _ in range(12)))
            maps = tuple(LinearMap(field, n, tuple(tuple(random.randrange(field.order) for _ in range(n))
                                                   for _ in range(out))) for n, out in zip(target.shape, (3, 2, 4)))
            expected = FiniteTensor.from_function(field, (3, 2, 4), lambda x, y, z:
                       field.sum(field.mul(target.coefficient(i, j, k),
                                           field.mul(maps[0].rows[x][i], field.mul(maps[1].rows[y][j], maps[2].rows[z][k])))
                                 for i, j, k in product(range(2), range(3), range(2))))
            self.assertEqual(target.restrict(*maps), expected)
            scheme = TensorScheme.from_tensor(target).restrict(*maps)
            self.assertEqual(scheme.target, expected)
            self.assertTrue(scheme.verify())
            left, right = (1, 0, 1), (1, 1)
            self.assertEqual(scheme.apply(left, right), expected.contract(left, right))

    def test_map_composition_and_coordinate_selection(self):
        field = FiniteField(5)
        first = LinearMap(field, 3, ((1, 2, 0), (0, 1, 4)))
        second = LinearMap(field, 2, ((2, 1), (0, 3), (4, 0)))
        for vector in product(range(5), repeat=3):
            self.assertEqual(second.compose(first).apply(vector), second.apply(first.apply(vector)))
        selection = LinearMap.selection(field, 3, (2, 0, 2))
        self.assertEqual(selection.apply((1, 2, 3)), (3, 1, 3))
        self.assertEqual(first.compose(LinearMap.identity(field, 3)), first)
        self.assertEqual(LinearMap.identity(field, 2).compose(first), first)

    def test_full_direct_sum_has_independent_first_inputs_and_no_cross_terms(self):
        field = FiniteField(5)
        first = FiniteTensor(field, (1, 2, 1), (2, 3))
        empty = FiniteTensor.zero(field, (0, 1, 2))
        last = FiniteTensor(field, (2, 1, 2), (1, 4, 2, 0))
        tensor = FiniteTensor.direct_sum((first, empty, last))
        self.assertEqual(tensor.shape, (3, 4, 5))
        for i, j, k in product(range(3), range(4), range(5)):
            expected = (first.coefficient(i, j, k) if i < 1 and j < 2 and k < 1
                        else last.coefficient(i-1, j-3, k-3) if i >= 1 and j >= 3 and k >= 3 else 0)
            self.assertEqual(tensor.coefficient(i, j, k), expected)
        scheme = TensorScheme.direct_sum(TensorScheme.from_tensor(t) for t in (first, empty, last))
        self.assertEqual(scheme.target, tensor)
        self.assertEqual(scheme.terms, 4)
        self.assertTrue(scheme.verify())
        self.assertEqual(scheme.apply((1, 2, 3), (4, 1, 2, 3)), tensor.contract((1, 2, 3), (4, 1, 2, 3)))

    def test_axis_permutations_and_constant_maps_on_polynomial_families(self):
        field = FiniteField(3)
        target = FiniteTensor.from_function(field, (2, 3, 2), lambda i, j, k: (i+2*j+k) % 3)
        scheme = TensorScheme.from_tensor(target)
        for order in permutations(range(3)):
            permuted = scheme.permute_axes(iter(order))
            self.assertTrue(permuted.verify())
            inverse = tuple(order.index(i) for i in range(3))
            self.assertEqual(permuted.permute_axes(inverse).target, target)
            for coordinate in product(*(range(d) for d in permuted.shape)):
                old = tuple(coordinate[inverse[i]] for i in range(3))
                self.assertEqual(permuted.target.coefficient(*coordinate), target.coefficient(*old))
        maps = (LinearMap(field, 2, ((1, 1),)), LinearMap.selection(field, 3, (2, 0)), LinearMap.identity(field, 2))
        polynomial = TensorDegeneration.from_scheme(scheme, leading=2).restrict(*maps)
        self.assertTrue(polynomial.verify())
        self.assertEqual(polynomial.recover_scheme().target, target.restrict(*maps))

    def test_three_sector_tags_preserve_shared_input_and_can_be_forgotten(self):
        c = ThreeSectorConstruction(FiniteField(5), 2, 1)
        original = TensorScheme.from_tensor(c.retained)
        tagged = c.tag_scheme(original)
        self.assertEqual(tagged.terms, original.terms)
        self.assertEqual(tagged.shape, (2, 12, 15))
        self.assertEqual(tagged.target, shared_first_tensor(c.branches))
        self.assertNotEqual(tagged.shape[0], FiniteTensor.direct_sum(c.branches).shape[0])
        y_labels, z_labels = c.side_labels
        forgotten = tagged.restrict(LinearMap.identity(c.field, 2),
                                    LinearMap.selection(c.field, 12, (h*4+j for j, h in enumerate(y_labels))),
                                    LinearMap.selection(c.field, 15, (h*5+k for k, h in enumerate(z_labels))))
        self.assertEqual(forgotten.target, c.retained)
        self.assertTrue(forgotten.verify())
        for i, j, k in product(range(2), range(12), range(15)):
            h, y = divmod(j, 4)
            h2, z = divmod(k, 5)
            self.assertEqual(tagged.target.coefficient(i, j, k), c.branches[h].coefficient(i, y, z) if h == h2 else 0)

    def test_invalid_tags_detect_hidden_characteristic_two_overlap(self):
        field = FiniteField(2)
        branch = FiniteTensor.unit(field)
        zero = TensorScheme.from_tensor(FiniteTensor.zero(field, (1, 1, 1)))
        with self.assertRaisesRegex(ValueError, "determine"):
            tag_shared_scheme(zero, (branch, branch), (0,), (0,))
        c = ThreeSectorConstruction(field, 2, 1)
        for labels in (((0,)*4, c.side_labels[1]), (c.side_labels[0], (0,)*5), ((0, 1), c.side_labels[1])):
            with self.assertRaises(ValueError):
                tag_shared_scheme(TensorScheme.from_tensor(c.retained), c.branches, *labels)
        with self.assertRaises(ValueError):
            c.tag_scheme(TensorScheme.from_tensor(c.source))

    def test_empty_axes_and_zero_maps(self):
        field = FiniteField(2)
        tensor = FiniteTensor.zero(field, (0, 2, 3))
        maps = (LinearMap(field, 0, ((), ())), LinearMap(field, 2, ()), LinearMap.identity(field, 3))
        scheme = TensorScheme.from_tensor(tensor).restrict(*maps)
        self.assertEqual(scheme.shape, (2, 0, 3))
        self.assertTrue(scheme.verify())
        self.assertEqual(shared_first_tensor((tensor, tensor)).shape, (0, 4, 6))
        for a, h in ((0, 0), (0, 2), (1, 0), (2, 0)):
            c = ThreeSectorConstruction(field, a, h)
            self.assertTrue(c.tag_scheme().verify())

    def test_bad_maps_fields_and_certificates_are_rejected(self):
        field = FiniteField(3)
        tensor = FiniteTensor.unit(field)
        identity = LinearMap.identity(field, 1)
        for maps in ((identity, identity), (identity, identity, LinearMap.identity(field, 2)),
                     (identity, identity, LinearMap.identity(FiniteField(2), 1))):
            with self.assertRaises(ValueError):
                tensor.restrict(*maps)
        for input_size, rows in ((-1, ()), (True, ()), (2, ((1,),)), (1, ((3,),))):
            with self.assertRaises(ValueError):
                LinearMap(field, input_size, rows)
        for indices in ((True,), (-1,), (1,)):
            with self.assertRaises(ValueError):
                LinearMap.selection(field, 1, indices)
        for tensors in ((), (tensor, FiniteTensor.unit(FiniteField(2)))):
            with self.assertRaises(ValueError):
                FiniteTensor.direct_sum(tensors)
        bad = replace(TensorScheme.from_tensor(tensor), c=((2,),))
        with self.assertRaises(ValueError):
            bad.restrict(identity, identity, LinearMap(field, 1, ((0,),)))
        with self.assertRaises(ValueError):
            TensorScheme.direct_sum((bad,))
        for order in ((0, 0, 2), (0, 1, True), (0, 1)):
            with self.assertRaises(ValueError):
                tensor.permute_axes(order)


if __name__ == "__main__":
    unittest.main()
