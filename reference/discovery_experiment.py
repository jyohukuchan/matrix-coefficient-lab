"""Run and replay capped finite-row and direct-decomposition experiments.

python3 -m reference.discovery_experiment --field 11 --output /tmp/discovery
Optional SciPy proposes LP solutions; optional Z3 proposes F2/F4 coefficients.
No numeric success or timeout is accepted as a mathematical conclusion.
"""

import argparse
from dataclasses import asdict
from fractions import Fraction
import json
from pathlib import Path
from time import monotonic
from math import isfinite

from .catalyst_search_experiment import bounded_restriction_sat, diagonal_tensor
from .discovery_inventory import DiscoveryBudget, generate_discovery_inventory, star_tensor
from .fields import FiniteField
from .finite_state_solver import RationalStateCertificate, solve_finite_states
from .flip_search import FlipBudget, search_flips
from .resource_planning import forecast_verified_inputs
from .schemes import naive_scheme
from .structured_certificates import CertificateGraph
from .tensor_schemes import TensorScheme
from .tensors import matrix_multiplication_tensor
from .witnessed_constraints import (DetectorConstraint, IntegerStateCertificate, StateAssemblyLimit,
                                    WitnessedOrder, WitnessedStateSystem, assemble_catalyst)


def serialize_finite(book, system, result):
    """Store the actual representatives/maps, not only row vectors or flags."""
    if not isinstance(system, WitnessedStateSystem) or system.registry != book.inventory:
        raise ValueError('serialize the actual inventory of the supplied constraint book')
    system.require_valid()
    for certificate in (result.state, result.dual):
        if certificate is not None:
            if certificate.system != system:
                raise ValueError('the result certificate belongs to a different finite system')
            certificate.require_valid()
    rows = []
    for row in system.constraints:
        matches = [witness for original, witness in zip(book.rows, book.witness_nodes) if original == row]
        if not matches:
            raise ValueError('finite row has no registered witness in this constraint book')
        witness = matches[0]
        if isinstance(row, WitnessedOrder):
            rows.append({'positive': row.positive, 'negative': row.negative,
                         'maps': book.graph.node(witness).args[2:], 'label': row.label})
        else:
            rows.append({'test': row.test, 'product': row.product})
    data = {'graph': book.graph.to_dict(book.tensor_nodes['unit']),
            'inventory': book.tensor_nodes, 'd': system.d, 'k': system.k, 'rows': rows}
    if result.state is not None:
        data['state'] = [str(v) for v in result.state.values]
    if result.dual is not None:
        data['dual'] = {'weights': result.dual.weights, 'gain': result.dual.m}
    return data


def replay_finite(data):
    graph, _ = CertificateGraph.from_dict(data['graph'])
    inventory = tuple((key, graph.materialize_tensor(node)) for key, node in data['inventory'].items())
    rows = tuple(WitnessedOrder(tuple(row['positive']), tuple(row['negative']),
                                tuple(graph.materialize_map(node) for node in row['maps']), row['label'])
                 if 'maps' in row else DetectorConstraint(row['test'], row['product']) for row in data['rows'])
    system = WitnessedStateSystem(data['d'], data['k'], inventory, rows)
    system.require_valid()
    checked = []
    if 'state' in data:
        state = RationalStateCertificate(system, tuple(Fraction(v) for v in data['state']))
        state.require_valid()
        checked.append('finite_feasible')
    if 'dual' in data:
        dual = IntegerStateCertificate(system, tuple(data['dual']['weights']), data['dual']['gain'])
        dual.require_valid()
        checked.append('negative_unit_dual')
    if len(checked) > 1:
        raise ValueError('a normalized feasible state and a negative-unit dual cannot coexist')
    return system, tuple(checked)


def _flip_case(scheme, goal, seconds, steps, seed):
    result = search_flips(scheme, goal_terms=goal, seed=seed,
                          budget=FlipBudget(timeout_seconds=seconds, max_steps=steps, probes_per_step=8))
    summary = {key: value for key, value in asdict(result).items() if key != 'scheme'}
    dag = None
    if result.scheme is not None:
        graph = CertificateGraph(scheme.field)
        root = graph.scheme(result.scheme)
        dag = graph.to_dict(root)
        summary['exact_coefficients'] = True
    summary['supplied_seed'] = 'coordinate-pair decomposition; no seven-term matrix scheme supplied'
    return summary, dag


def run_experiment(field=FiniteField(11), *, k=5, include_known_control=False,
                   contexts=('M2', 'C22', 'star'), budget=DiscoveryBudget(), lp_seconds=10.0,
                   flip_seconds=5.0, core_seconds=30.0, flip_steps=10_000,
                   seeds=(0,), sat_seconds=0, solver_path=None, output=None):
    if not isinstance(field, FiniteField) or not isinstance(budget, DiscoveryBudget):
        raise ValueError('supply an explicit finite field and discovery budget')
    contexts = tuple(contexts)
    if type(k) is not int or k < 1 or type(flip_steps) is not int or flip_steps < 0:
        raise ValueError('k positive and flip steps nonnegative are required')
    seeds = tuple(seeds)
    if not seeds or any(type(seed) is not int for seed in seeds) or len(set(seeds)) != len(seeds):
        raise ValueError('supply distinct integer search seeds')
    if any(not isinstance(v, (int, float)) or isinstance(v, bool) or not isfinite(v) or v <= 0
           for v in (lp_seconds, flip_seconds, core_seconds)):
        raise ValueError('proposal and flip time budgets must be positive finite numbers')
    if type(sat_seconds) is not int or sat_seconds < 0:
        raise ValueError('SAT seconds must be a nonnegative integer')
    start = monotonic()
    report = {'field': {'p': field.p, 'modulus': field.modulus}, 'd': 2, 'k': k,
              'include_known_control': include_known_control, 'contexts': list(contexts),
              'budget': asdict(budget), 'new_useful_certificate': False,
              'timing_scope': 'wall timings include construction and exact verification; native solver limits are advisory'}
    artifacts = {}
    try:
        book, skipped = generate_discovery_inventory(field, include_known_control=include_known_control,
                                                    contexts=contexts, budget=budget)
        built = monotonic()
        system = book.system(2, k)
        result = solve_finite_states(system, timeout_seconds=lp_seconds)
        finite = {'status': result.status, 'reason': result.reason, 'atoms': len(book.inventory),
                  'rows': len(book.rows), 'construction_seconds': built-start,
                  'proposal_and_verification_seconds': monotonic()-built,
                  'numeric_statuses': result.numeric_statuses, 'skipped': skipped,
                  'family_counts': {family: book.families.count(family) for family in dict.fromkeys(book.families)}}
        if result.state is not None:
            finite['exact_state'] = {key: str(value) for key, value in result.state.assignment.items()}
            finite['interpretation'] = 'No negative-unit dual for these selected rows; no global catalyst exclusion.'
        if result.dual is not None:
            finite['gain'] = result.dual.m
            finite['used_rows'] = [{'index': i, 'family': book.families[i], 'weight': weight}
                                   for i, weight in enumerate(result.dual.weights) if weight]
            try:
                certificate = assemble_catalyst(result.dual)
                report['new_useful_certificate'] = k == 5
                finite['actual_catalyst_verified'] = True
                forecast = forecast_verified_inputs(certificate, TensorScheme.from_tensor(certificate.S),
                    TensorScheme.from_tensor(certificate.D), Fraction(12, 5), max_b=8, max_iterations=8)
                finite['cost_forecast'] = forecast.to_dict()
            except StateAssemblyLimit as error:
                finite['assembly_status'] = 'resource_cap'
                finite['assembly_reason'] = str(error)
        report['finite'] = finite
        artifacts['finite-replay.json'] = serialize_finite(book, system, result)
    except StateAssemblyLimit as error:
        report['finite'] = {'status': 'resource_cap', 'reason': str(error)}

    direct = {}
    if field.order <= 16:
        matrix = naive_scheme(field, 2, 2, 2).as_tensor_scheme()
        core = matrix.tensor_product(TensorScheme.from_tensor(star_tensor(field)))
        for seed in seeds:
            for name, supplied, goal, seconds in (('matrix-control', matrix, 7, flip_seconds),
                                                  ('star-product', core, 19, core_seconds)):
                label = f'{name}-seed{seed}'
                direct[label], dag = _flip_case(supplied, goal, seconds, flip_steps, seed)
                if dag is not None:
                    artifacts[label+'.json'] = dag
    else:
        direct['status'] = 'unavailable'
        direct['reason'] = 'the local flip engine supports field orders at most 16'
    report['direct_seeded'] = direct
    if sat_seconds:
        if field not in (FiniteField(2), FiniteField(2, (1, 1, 1))):
            report['direct_unseeded'] = {'status': 'unavailable', 'reason': 'SAT encoding supports F2 and F4'}
        else:
            target = matrix_multiplication_tensor(field, 2, 2, 2)
            sat_start = monotonic()
            result = bounded_restriction_sat(diagonal_tensor(field, 7), target,
                                             identical_blocks=7, timeout_ms=1000*sat_seconds,
                                             solver_path=solver_path)
            report['direct_unseeded'] = {'status': result.status, 'reason': result.reason,
                                          'target': 'M2', 'goal_terms': 7,
                                          'elapsed_seconds': monotonic()-sat_start}
            if result.maps is not None:
                scheme = TensorScheme(target, *(tuple(tuple(mapping.rows[out][term] for out in range(mapping.output_size))
                                                       for term in range(7)) for mapping in result.maps))
                scheme.require_exact()
                graph = CertificateGraph(field)
                artifacts['unseeded-matrix-control.json'] = graph.to_dict(graph.scheme(scheme))
                report['direct_unseeded']['exact_coefficients'] = True
    report['elapsed_seconds'] = monotonic()-start
    report['interpretation'] = ('These searches solve different finite problems. They are baselines, not a demonstrated '
                                'performance advantage or a complete spectral construction. Rank 19 for M2 tensor star '
                                'would not by itself supply a restriction from five stars.')
    if output is not None:
        destination = Path(output)
        destination.mkdir(parents=True, exist_ok=True)
        artifacts['report.json'] = report
        expected = set(artifacts)
        # Clear only files managed by this runner from an earlier capped run.
        for old in destination.glob('*.json'):
            if (old.name in ('finite-replay.json', 'unseeded-matrix-control.json')
                    or old.name.startswith(('matrix-control-seed', 'star-product-seed'))) and old.name not in expected:
                old.unlink()
        for name, data in artifacts.items():
            (destination/name).write_text(json.dumps(data, indent=2)+'\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--field', type=int, choices=(2, 3, 4, 11), default=11)
    parser.add_argument('--k', type=int, default=5)
    parser.add_argument('--known-control', action='store_true')
    parser.add_argument('--no-contexts', action='store_true')
    parser.add_argument('--lp-seconds', type=float, default=10)
    parser.add_argument('--flip-seconds', type=float, default=5)
    parser.add_argument('--core-seconds', type=float, default=30)
    parser.add_argument('--flip-steps', type=int, default=10_000)
    parser.add_argument('--seeds', type=int, nargs='+', default=[0])
    parser.add_argument('--sat-seconds', type=int, default=0)
    parser.add_argument('--solver-path')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--replay', type=Path)
    args = parser.parse_args()
    if args.replay is not None:
        system, checked = replay_finite(json.loads(args.replay.read_text()))
        print(json.dumps({'checked': checked, 'rows': len(system.constraints), 'atoms': len(system.inventory)}))
        return
    field = FiniteField(2, (1, 1, 1)) if args.field == 4 else FiniteField(args.field)
    report = run_experiment(field, k=args.k, include_known_control=args.known_control,
                            contexts=() if args.no_contexts else ('M2', 'C22', 'star'),
                            lp_seconds=args.lp_seconds, flip_seconds=args.flip_seconds,
                            core_seconds=args.core_seconds, flip_steps=args.flip_steps, seeds=args.seeds,
                            sat_seconds=args.sat_seconds, solver_path=args.solver_path, output=args.output)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
