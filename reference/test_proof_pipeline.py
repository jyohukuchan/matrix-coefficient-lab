"""Finite integration must report its control dependence and failed searches."""

from dataclasses import replace
from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest

from .proof_pipeline import (ConstraintInventory, PipelineConfig, generate_proof_inventory,
                             run_pipeline)
from .fields import FiniteField
from .structured_certificates import CertificateGraph, GraphLimits
from .tensors import FiniteTensor
from .witnessed_constraints import StateAssemblyLimits


class ProofPipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = PipelineConfig()
        cls.report = run_pipeline(cls.config)

    def test_end_to_end_known_controls_and_gain_compilation(self):
        report = self.report
        self.assertEqual(report['status'], 'found')
        self.assertEqual([(x['size'], x['terms']) for x in report['known_reconstruction']], [(2, 7), (4, 49)])
        self.assertEqual(report['gain_compilation']['status'], 'verified')
        self.assertFalse(report['new_useful_certificate'])
        graph, root = CertificateGraph.from_dict(report['compiled_scheme_dag'])
        self.assertTrue(graph.materialize_scheme(root).verify())

    def test_attribution_does_not_credit_unused_proof_rows(self):
        used = self.report['search']['used_rows']
        self.assertEqual({x['family'] for x in used}, {'known_control', 'detector'})
        analyses = self.report['family_analysis']
        self.assertEqual(analyses['without_known_control']['status'], 'finite_exhausted')
        self.assertEqual(analyses['controls_only']['status'], 'found')
        for family in ('sector', 'exact_type', 'fourier', 'determinant'):
            self.assertEqual(analyses['without_'+family]['status'], 'found')

    def test_tensor_interning_connects_type_and_fourier(self):
        book = generate_proof_inventory(self.config)
        type_row = book.rows[book.families.index('exact_type')]
        fourier_row = book.rows[book.families.index('fourier')]
        self.assertEqual(set(type_row.negative), set(fourier_row.positive))
        self.assertEqual(len(fourier_row.positive), 20)
        book.system(2, 8).require_valid()
        self.assertIsNone(book.witness_nodes[-1])

    def test_k5_and_proof_only_are_finite_failures_not_global_obstructions(self):
        for config in (replace(self.config, k=5), replace(self.config, include_known_control=False)):
            report = run_pipeline(replace(config, analyze_families=False))
            self.assertEqual(report['status'], 'finite_exhausted')
            self.assertNotIn('gain_compilation', report)
            self.assertFalse(report['new_useful_certificate'])

    def test_caps_and_unreached_target_are_reported_honestly(self):
        capped = run_pipeline(replace(self.config, max_subsets=0, analyze_families=False))
        self.assertEqual(capped['status'], 'resource_cap')
        tiny = run_pipeline(self.config, graph_limits=GraphLimits(max_nodes=1))
        self.assertEqual(tiny['status'], 'resource_cap')
        tiny = run_pipeline(self.config, limits=StateAssemblyLimits(max_tensor_entries=1))
        self.assertEqual(tiny['status'], 'resource_cap')
        target = run_pipeline(replace(self.config, tau=Fraction(12, 5), analyze_families=False))
        self.assertEqual(target['status'], 'found')
        self.assertEqual(target['gain_compilation']['status'], 'resource_cap')
        self.assertNotIn('compiled_scheme_dag', target)

    def test_saved_dag_replays_actual_49_term_control(self):
        with tempfile.TemporaryDirectory() as directory:
            report = run_pipeline(replace(self.config, analyze_families=False), output=directory)
            saved = json.loads((Path(directory)/'report.json').read_text())
            self.assertEqual(saved['status'], report['status'])
            graph, root = CertificateGraph.from_dict(json.loads((Path(directory)/'construction-dag.json').read_text()))
            scheme = graph.materialize_scheme(root)
            self.assertEqual((scheme.shape, scheme.terms), ((16, 16, 16), 49))
            self.assertTrue(scheme.verify())

    def test_inventory_names_and_configuration_are_checked(self):
        book = ConstraintInventory(FiniteField(11))
        self.assertEqual(book.register('scalar-alias', FiniteTensor.unit(book.field)), 'unit')
        with self.assertRaises(ValueError):
            book.register('unit', FiniteTensor(book.field, (1, 1, 1), (3,)))
        for kwargs in ({'d': 3}, {'tau': 2.4}, {'max_subsets': True}, {'type_counts': (1, True)}):
            with self.assertRaises(ValueError):
                PipelineConfig(**kwargs)


if __name__ == '__main__':
    unittest.main()
