"""Bounded proof-row generation, exact dual search, attribution and replay.

python3 -m reference.proof_pipeline --output /path/to/artifacts
The known-rank control is supplied independently of the proof rows. Its
success never implies that proof-derived rows discovered a new algorithm.
"""

import argparse
from dataclasses import asdict, dataclass, replace
from fractions import Fraction
import json
from pathlib import Path

from .catalyst_compiler import compile_absorption, compile_powered_catalyst
from .constraint_generators import WitnessedOrderExport, export_sector_order
from .determinant_filtration import DeterminantFiltration
from .extraction import simplify_scheme
from .fields import FiniteField
from .finite_type_constraints import export_exact_type_order, export_finite_separation_order
from .gain_planning import compile_gain_target
from .resource_planning import forecast_verified_inputs
from .schemes import naive_multiply, naive_scheme, strassen_scheme
from .sectors import ThreeSectorConstruction
from .separation import FiniteSeparation
from .structured_certificates import CertificateGraph, GraphLimits, GraphResourceLimit
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, matrix_multiplication_tensor
from .witnessed_constraints import (DEFAULT_LIMITS, DetectorConstraint, StateAssemblyLimit,
                                    WitnessedStateSystem, assemble_catalyst,
                                    clear_rational_dual, find_nonnegative_dual, rank_upper_constraint)


@dataclass(frozen=True)
class PipelineConfig:
    """Small reproducible controls; not a useful d=2,k=5 certificate."""

    d: int = 2
    k: int = 8
    sector_a: int = 2
    sector_h: int = 1
    sector_exponent: int = 1
    type_counts: tuple[int, int] = (1, 1)
    determinant_d: int = 1
    determinant_e: int = 1
    max_subsets: int = 10_000
    max_b: int = 4
    max_iterations: int = 4
    tau: Fraction = Fraction(4)
    include_known_control: bool = True
    analyze_families: bool = True

    def __post_init__(self):
        if self.d != 2 or type(self.d) is not int:
            raise ValueError("this finite control pipeline currently supports d=2")
        if type(self.k) is not int or self.k < 1:
            raise ValueError("k must be a positive integer")
        if type(self.include_known_control) is not bool or type(self.analyze_families) is not bool:
            raise ValueError("control and analysis flags must be booleans")
        if not isinstance(self.tau, Fraction) or self.tau <= 0:
            raise ValueError("tau must be a positive Fraction")
        for name in ("sector_a", "sector_h", "sector_exponent", "determinant_d",
                     "determinant_e", "max_subsets", "max_b", "max_iterations"):
            value = getattr(self, name)
            if type(value) is not int or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
        if self.sector_a == 0 or self.sector_h == 0:
            raise ValueError("sector export requires positive a,h")
        counts = tuple(self.type_counts)
        if len(counts) != 2 or any(type(x) is not int or x < 0 for x in counts):
            raise ValueError("provide two nonnegative exact-type counts")
        object.__setattr__(self, "type_counts", counts)


class ConstraintInventory:
    """Intern equal tensors across families and retain every ordinary witness."""

    def __init__(self, field, limits=DEFAULT_LIMITS, graph_limits=GraphLimits()):
        self.field, self.limits = field, limits
        self.graph = CertificateGraph(field, graph_limits)
        self.inventory = {"unit": FiniteTensor.unit(field)}
        self._canonical = {self.inventory["unit"]: "unit"}
        self.tensor_nodes = {"unit": self.graph.tensor(self.inventory["unit"])}
        self.rows, self.families, self.witness_nodes = [], [], []

    def register(self, name, tensor):
        if not isinstance(name, str) or not name or not isinstance(tensor, FiniteTensor) or tensor.field != self.field:
            raise ValueError("register a named tensor over the inventory field")
        if name in self.inventory and self.inventory[name] != tensor:
            raise ValueError("tensor name already denotes a different tensor")
        if tensor in self._canonical:
            return self._canonical[tensor]
        self.limits.tensor(tensor.shape)
        self.inventory[name] = tensor
        self._canonical[tensor] = name
        self.tensor_nodes[name] = self.graph.tensor(tensor)
        return name

    def add_export(self, export, family):
        if not isinstance(export, WitnessedOrderExport):
            raise ValueError("supply an exact ordinary export")
        export.require_valid(self.limits)
        renaming = {key: self.register(key, tensor) for key, tensor in export.inventory}
        row = replace(export.row, positive=tuple(renaming[key] for key in export.row.positive),
                      negative=tuple(renaming[key] for key in export.row.negative))
        def graph_sum(keys):
            return (self.graph.direct_sum(self.tensor_nodes[key] for key in keys) if keys
                    else self.graph.tensor(FiniteTensor.zero(self.field, (0, 0, 0))))
        source, target = graph_sum(row.positive), graph_sum(row.negative)
        maps = tuple(self.graph.mapping(mapping) for mapping in row.maps)
        node = self.graph.witness(source, target, maps, row.label)
        self.rows.append(row)
        self.families.append(family)
        self.witness_nodes.append(node)
        return row

    def add_detector(self, test_key, d):
        test = self.inventory[test_key]
        shape = tuple(d*d*x for x in test.shape)
        self.limits.tensor(shape)
        product_key = self.register(f"M{d}*{test_key}", matrix_multiplication_tensor(self.field, d, d, d).tensor_product(test))
        self.rows.append(DetectorConstraint(test_key, product_key))
        self.families.append("detector")
        self.witness_nodes.append(None)  # A hypothetical state row is not a restriction.

    def system(self, d, k, indices=None):
        rows = tuple(self.rows if indices is None else (self.rows[i] for i in indices))
        return WitnessedStateSystem(d, k, tuple(self.inventory.items()), rows)


def generate_proof_inventory(config=PipelineConfig(), *, limits=DEFAULT_LIMITS,
                             graph_limits=GraphLimits()):
    # F_11 contains a primitive tenth root for the two-word Fourier control.
    field = FiniteField(11)
    book = ConstraintInventory(field, limits, graph_limits)
    sector = ThreeSectorConstruction(field, config.sector_a, config.sector_h)
    book.add_export(export_sector_order(sector, "sector-source", "sector-retained",
                                       exponent=config.sector_exponent, limits=limits), "sector")
    # Supplied scalar branches exercise exact word ordering and shared inputs.
    # They are finite test inputs, not the missing spectral existence argument.
    branches = (FiniteTensor.unit(field), FiniteTensor(field, (1, 1, 1), (3,)))
    typed = export_exact_type_order(branches, config.type_counts, "type-power", "type-words", limits=limits)
    book.add_export(typed.restriction, "exact_type")
    separation = FiniteSeparation(typed.branches)
    book.add_export(export_finite_separation_order(separation, "Fourier-shared", "Fourier-direct", limits=limits), "fourier")
    # The Fourier source is canonicalized to the exact-type target, connecting
    # the two actual rows even though their exported IDs differ.
    filtration = DeterminantFiltration(field, config.determinant_d, config.determinant_e, limits)
    book.add_export(filtration.export_order("determinant-source", "determinant-graded", limits=limits), "determinant")
    if config.include_known_control:
        supplied = strassen_scheme(field).as_tensor_scheme()
        key = book.register("known-M2", supplied.target)
        upper = rank_upper_constraint("unit", key, supplied)
        book.add_export(WitnessedOrderExport(upper, (("unit", book.inventory["unit"]), (key, supplied.target))), "known_control")
    book.add_detector("unit", config.d)
    return book


def _search_summary(book, config, indices):
    system = book.system(config.d, config.k, indices)
    result = find_nonnegative_dual(system, max_subsets=config.max_subsets, limits=book.limits)
    summary = {"status": result.status, "checked_subsets": result.checked_subsets, "reason": result.reason}
    if result.status == "found":
        integer = clear_rational_dual(system, result.coefficients)
        summary["gain"] = integer.m
        summary["used_rows"] = [{"index": original, "family": book.families[original], "weight": weight}
                                for original, weight in zip(indices, integer.weights) if weight]
    return result, system, summary


def _known_reconstruction(certificate, S, D, graph):
    powered = compile_powered_catalyst(certificate, D)
    seed, outputs = naive_scheme(certificate.field, 1, 1, 1), []
    for _ in range(2):
        step = compile_absorption(powered, seed, S)
        raw = step.matrix_scheme
        seed = simplify_scheme(raw.as_tensor_scheme()).to_matrix_scheme(raw.n, raw.n, raw.n)
        seed.require_exact()
        left = [[(i+j) % certificate.field.p for j in range(seed.n)] for i in range(seed.n)]
        right = [[int(i == j) for j in range(seed.n)] for i in range(seed.n)]
        if seed.multiply(left, right) != naive_multiply(certificate.field, left, right):
            raise ValueError("reconstructed scheme fails independent matrix evaluation")
        root = graph.scheme(seed.as_tensor_scheme())
        outputs.append({"size": seed.n, "raw_terms": raw.terms, "terms": seed.terms,
                        "exact_coefficients": True, "graph_node": root,
                        "support": graph.scheme_support_summary(root)})
    return outputs


def run_pipeline(config=PipelineConfig(), *, limits=DEFAULT_LIMITS,
                 graph_limits=GraphLimits(), output=None):
    report = {"status": "incomplete", "config": {**asdict(config), "tau": str(config.tau)},
              "new_useful_certificate": False}
    graph, root = None, None
    try:
        book = generate_proof_inventory(config, limits=limits, graph_limits=graph_limits)
        graph = book.graph
        all_indices = tuple(range(len(book.rows)))
        result, system, summary = _search_summary(book, config, all_indices)
        report.update(status=result.status, search=summary, inventory_size=len(book.inventory),
                      rows=[{"index": i, "family": family,
                             "vector": dict(vector), "witness_node": witness}
                            for i, (family, vector, witness) in enumerate(zip(book.families, system.sparse_vectors, book.witness_nodes))])
        if config.analyze_families:
            analyses = {}
            for family in dict.fromkeys(book.families):
                if family != "detector":
                    indices = tuple(i for i in all_indices if book.families[i] != family)
                    analyses["without_"+family] = _search_summary(book, config, indices)[2]
            control_indices = tuple(i for i in all_indices if book.families[i] in ("known_control", "detector"))
            analyses["controls_only"] = _search_summary(book, config, control_indices)[2]
            report["family_analysis"] = analyses
        if result.status == "found":
            integer = clear_rational_dual(system, result.coefficients)
            certificate = assemble_catalyst(integer, limits=limits)
            S, D = TensorScheme.from_tensor(certificate.S), TensorScheme.from_tensor(certificate.D)
            forecast = forecast_verified_inputs(certificate, S, D, config.tau,
                                                max_b=config.max_b, max_iterations=config.max_iterations)
            report["cost_forecast"] = forecast.to_dict()
            compiled = compile_gain_target(certificate, S, D, config.tau,
                                           max_b=config.max_b, max_iterations=config.max_iterations)
            report["gain_compilation"] = {"status": compiled.status, "reason": compiled.reason,
                                           "size": compiled.plan.predicted_size,
                                           "terms": compiled.plan.predicted_terms}
            if compiled.prime_scheme is not None:
                report["gain_compilation"]["exact_coefficients"] = True
                # Separate field graph: extension descent is performed once.
                final_graph = CertificateGraph(compiled.prime_scheme.field, graph_limits)
                final_root = final_graph.scheme(compiled.prime_scheme.as_tensor_scheme())
                report["compiled_scheme_dag"] = final_graph.to_dict(final_root)
            if config.k == 8 and config.include_known_control:
                report["known_reconstruction"] = _known_reconstruction(certificate, S, D, graph)
                root = report["known_reconstruction"][-1]["graph_node"]
        if root is None:
            root = book.witness_nodes[0]
        report["graph_nodes"] = len(graph.nodes)
        report["interpretation"] = ("Known input control reconstructed; proof-row contribution is reported separately. "
                                    "Finite exhaustion and resource caps do not exclude other row families or catalysts.")
    except (StateAssemblyLimit, GraphResourceLimit) as error:
        report.update(status="resource_cap", reason=str(error))
    if output is not None:
        destination = Path(output)
        destination.mkdir(parents=True, exist_ok=True)
        (destination / "report.json").write_text(json.dumps(report, indent=2)+"\n")
        if graph is not None and root is not None:
            (destination / "construction-dag.json").write_text(json.dumps(graph.to_dict(root), indent=2)+"\n")
        else:
            # An earlier successful run must not masquerade as this capped run.
            (destination / "construction-dag.json").unlink(missing_ok=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--k", type=int, default=8)
    parser.add_argument("--tau", type=Fraction, default=Fraction(4))
    parser.add_argument("--max-subsets", type=int, default=10_000)
    parser.add_argument("--max-b", type=int, default=4)
    parser.add_argument("--max-iterations", type=int, default=4)
    parser.add_argument("--proof-only", action="store_true")
    parser.add_argument("--skip-family-analysis", action="store_true")
    args = parser.parse_args()
    config = PipelineConfig(k=args.k, tau=args.tau, max_subsets=args.max_subsets,
                            max_b=args.max_b, max_iterations=args.max_iterations,
                            include_known_control=not args.proof_only,
                            analyze_families=not args.skip_family_analysis)
    print(json.dumps(run_pipeline(config, output=args.output), indent=2))


if __name__ == "__main__":
    main()
