"""Compare independent dense operations, replay, forgery rejection and caps."""

from copy import deepcopy
from itertools import product
import json
import unittest
from unittest.mock import patch

from .fields import FiniteField
from .maps import LinearMap
from .schemes import naive_scheme, strassen_scheme
from .structured_certificates import CertificateGraph, GraphLimits, GraphResourceLimit
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor


class StructuredCertificateTests(unittest.TestCase):
    def test_asymmetric_tensor_operations_match_dense_every_coordinate(self):
        for f in (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1))):
            g = CertificateGraph(f)
            a = FiniteTensor.from_function(f, (2, 1, 3), lambda i, j, k: (i+k) % f.order)
            b = FiniteTensor.from_function(f, (1, 2, 1), lambda i, j, k: j)
            x, y = g.tensor(a), g.tensor(b)
            cases = ((g.tensor_product(x, y), a.tensor_product(b)),
                     (g.direct_sum((x, y, x)), FiniteTensor.direct_sum((a, b, a))),
                     (g.repeat(x, 3), FiniteTensor.direct_sum((a, a, a))),
                     (g.permute(x, (2, 0, 1)), a.permute_axes((2, 0, 1))),
                     (g.power(x, 0), FiniteTensor.unit(f)),
                     (g.power(x, 2), a.tensor_power(2)))
            for node, expected in cases:
                self.assertEqual(g.materialize_tensor(node), expected)

    def test_scheme_compositions_match_dense_families(self):
        f = FiniteField(3)
        g = CertificateGraph(f)
        a = TensorScheme.from_tensor(FiniteTensor(f, (2, 1, 2), (1, 0, 0, 2)))
        b = TensorScheme.from_tensor(FiniteTensor(f, (1, 2, 1), (1, 2)))
        x, y = g.scheme(a), g.scheme(b)
        cases = ((g.scheme_product(x, y), a.tensor_product(b)),
                 (g.direct_sum((x, y, x)), TensorScheme.direct_sum((a, b, a))),
                 (g.repeat(x, 2), TensorScheme.direct_sum((a, a))),
                 (g.permute(x, (2, 0, 1)), a.permute_axes((2, 0, 1))),
                 (g.power(x, 2), a.tensor_power(2)))
        for node, expected in cases:
            self.assertEqual(g.materialize_scheme(node), expected)

    def test_map_composition_product_and_local_substitution(self):
        f = FiniteField(3)
        g = CertificateGraph(f)
        a = LinearMap(f, 2, ((1, 2), (0, 1)))
        b = LinearMap(f, 3, ((1, 0, 1), (2, 1, 0)))
        x, y = g.mapping(a), g.mapping(b)
        self.assertEqual(g.materialize_map(g.compose(x, y)), a.compose(b))
        tensor_map = g.map_product(x, g.identity(2))
        expected = LinearMap(f, 4, tuple(tuple(a.rows[o//2][i//2] if o % 2 == i % 2 else 0 for i in range(4)) for o in range(4)))
        self.assertEqual(g.materialize_map(tensor_map), expected)
        supplied = strassen_scheme(f).as_tensor_scheme()
        node = g.scheme(supplied)
        maps = (g.selection(4, (0, 3)), g.identity(4), g.selection(4, (3, 1, 0)))
        dense_maps = tuple(g.materialize_map(i) for i in maps)
        self.assertEqual(g.materialize_scheme(g.restrict(node, maps)), supplied.restrict(*dense_maps))

    def test_extension_powers_descend_once(self):
        f = FiniteField(2, (1, 1, 1))
        supplied = naive_scheme(f, 2, 2, 2).as_tensor_scheme().rescale_terms(2, 3)
        g = CertificateGraph(f)
        node = g.descend(g.power(g.scheme(supplied), 2))
        expected = supplied.tensor_power(2).descend()
        self.assertEqual(g.node(node).terms, 4*supplied.terms**2)
        self.assertEqual(g.materialize_scheme(node), expected)
        # A non-base-valued target cannot be relabelled as unchanged descent.
        bad = g.scheme(TensorScheme.from_tensor(FiniteTensor(f, (1, 1, 1), (2,))))
        with self.assertRaisesRegex(ValueError, "base-valued"):
            g.descend(bad)

    def test_checked_witness_and_invalid_identity(self):
        f = FiniteField(3)
        g = CertificateGraph(f)
        source = g.tensor(FiniteTensor(f, (2, 2, 2), (1, 0, 0, 0, 0, 0, 0, 1)))
        target = g.tensor(FiniteTensor.unit(f))
        maps = (g.selection(2, (0,)),)*3
        witness = g.witness(source, target, maps, "coordinate-selection")
        self.assertEqual(g.materialize_tensor(witness), FiniteTensor.unit(f))
        wrong = g.tensor(FiniteTensor(f, (1, 1, 1), (2,)))
        with self.assertRaisesRegex(ValueError, "exact coefficient"):
            g.witness(source, wrong, maps, "false-label")

    def test_large_power_queries_and_replay_without_large_allocation(self):
        f = FiniteField(2)
        g = CertificateGraph(f)
        node = g.power(g.scheme(strassen_scheme(f).as_tensor_scheme()), 25)
        self.assertEqual((g.node(node).terms, g.node(node).shape), (7**25, (4**25,)*3))
        self.assertEqual(g.scheme_coefficient(node, 0, 0, 0), 1)
        self.assertEqual(g.tensor_coefficient(g.node(node).target, 0, 0, 0), 1)
        payload = json.loads(json.dumps(g.to_dict(node)))
        original = FiniteTensor.from_function
        def small_only(field, shape, callback):
            self.assertLessEqual(shape[0]*shape[1]*shape[2], 64)
            return original(field, shape, callback)
        with patch.object(FiniteTensor, "from_function", side_effect=small_only):
            h, root = CertificateGraph.from_dict(payload)
            self.assertEqual(h.scheme_coefficient(root, 0, 0, 0), 1)
        self.assertLess(len(h.nodes), 40)
        with self.assertRaises(GraphResourceLimit):
            h.materialize_scheme(root)

    def test_analytic_matrix_node_and_corrupt_serialized_claims(self):
        g = CertificateGraph(FiniteField(2))
        root = g.matrix(2**25)
        self.assertEqual(g.tensor_coefficient(root, 0, 0, 0), 1)
        self.assertEqual(g.tensor_coefficient(root, 0, 1, 0), 0)
        payload = json.loads(json.dumps(g.to_dict(root)))
        h, r = CertificateGraph.from_dict(payload)
        self.assertEqual(h.node(r).shape, (2**50,)*3)
        for mutate in (lambda p: p['nodes'][0].update(op='rank_oracle'),
                       lambda p: p['nodes'][0].update(args=[0]),
                       lambda p: p['nodes'][0].update(terms=1),
                       lambda p: p.update(root=True)):
            bad = deepcopy(payload)
            mutate(bad)
            with self.assertRaises(ValueError):
                CertificateGraph.from_dict(bad)

    def test_corrupt_leaf_cannot_survive_replay(self):
        g = CertificateGraph(FiniteField(2))
        root = g.scheme(strassen_scheme(g.field).as_tensor_scheme())
        payload = json.loads(json.dumps(g.to_dict(root)))
        payload['nodes'][root]['data'][2][0] = [0]*4
        with self.assertRaisesRegex(ValueError, "invalid tensor coefficient"):
            CertificateGraph.from_dict(payload)

    def test_zero_axes_zero_repetition_and_invalid_coordinates(self):
        g = CertificateGraph(FiniteField(2))
        empty = FiniteTensor.zero(g.field, (0, 2, 3))
        node = g.scheme(TensorScheme.from_tensor(empty))
        self.assertEqual(g.materialize_scheme(g.repeat(node, 2)).shape, (0, 4, 6))
        self.assertEqual(g.materialize_scheme(g.repeat(node, 0)).shape, (0, 0, 0))
        self.assertEqual(g.materialize_tensor(g.power(g.tensor(empty), 2)), empty.tensor_power(2))
        with self.assertRaises(ValueError):
            g.scheme_coefficient(node, 0, 0, 0)
        with self.assertRaises(ValueError):
            g.tensor_coefficient(g.matrix(2), True, 0, 0)
        with self.assertRaises(ValueError):
            g.power(node, True)

    def test_query_and_allocation_caps_are_incomplete(self):
        g = CertificateGraph(FiniteField(2), GraphLimits(max_query_work=1))
        node = g.tensor_product(g.matrix(2), g.matrix(2))
        with self.assertRaisesRegex(GraphResourceLimit, "query incomplete"):
            g.tensor_coefficient(node, 0, 0, 0)
        g = CertificateGraph(FiniteField(2), GraphLimits(max_tensor_entries=1))
        large = g.matrix(2)
        with patch.object(FiniteTensor, "from_function") as allocate:
            with self.assertRaises(GraphResourceLimit):
                g.materialize_tensor(large)
            allocate.assert_not_called()
        with self.assertRaises(GraphResourceLimit):
            CertificateGraph.from_dict(g.to_dict(large), GraphLimits(max_nodes=0))

    def test_support_counts_match_expanded_coefficients_and_scale_symbolically(self):
        g = CertificateGraph(FiniteField(2))
        small = g.scheme(strassen_scheme(g.field).as_tensor_scheme())
        root = g.direct_sum((g.power(small, 2), g.repeat(small, 3)))
        expanded = g.materialize_scheme(root)
        counts = g.scheme_support_summary(root)
        expected = tuple(sum(bool(v) for row in family for v in row)
                         for family in (expanded.a, expanded.b, expanded.c))
        self.assertEqual(counts['coefficient_nonzeros'], expected)
        large = g.power(small, 25)
        base = g.scheme_support_summary(small)
        self.assertEqual(g.scheme_support_summary(large)['coefficient_nonzeros'],
                         tuple(x**25 for x in base['coefficient_nonzeros']))
        mapped = g.restrict(small, (g.identity(4),)*3)
        self.assertEqual(g.scheme_support_summary(mapped)['status'], 'unknown_after_linear_substitution')


if __name__ == '__main__':
    unittest.main()
