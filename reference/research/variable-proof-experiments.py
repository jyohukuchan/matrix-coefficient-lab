"""Bounded non-scalar proof-row experiments, with literal coefficient replay.

Run from the repository root. No supplied matrix rank-seven row is included.
Every accepted ordinary row is checked on its full coefficient array.
The field is never extended. A feasible result means an exact rational state
with unit value 1 and M2 value 5; a dual means an exact negative-unit balance.
Coordinate search exhaustion concerns only its selected coordinate family.
Construction windows are checked between operations. Full coefficient
validation and saved-certificate replay are timed separately and can extend
beyond a window; native linear-program proposals have a ten-second budget.
"""

import argparse
from collections import Counter
from dataclasses import asdict
from fractions import Fraction
import json
from math import prod
from pathlib import Path
import sys
import tempfile
from time import monotonic

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from reference.constraint_generators import (WitnessedOrderExport,
    export_sector_order, lift_witnessed_order)
from reference.constructions import primitive_root
from reference.coordinate_subrank import CoordinateSearchBudget, find_coordinate_subrank
from reference.determinant_filtration import DeterminantFiltration
from reference.discovery_inventory import commute_product_export, scalar_lower_export
from reference.fields import FiniteField
from reference.finite_state_solver import RationalStateCertificate, solve_finite_states
from reference.finite_type_constraints import (export_exact_type_order,
    export_finite_separation_order)
from reference.maps import LinearMap
from reference.proof_pipeline import ConstraintInventory
from reference.sectors import ThreeSectorConstruction, convolution_tensor
from reference.separation import FiniteSeparation
from reference.sector_matrix import (export_matrix_column_sum_order,
    export_retained_matrix_order, export_retained_two_copy_matrix_order,
    export_star_matrix_order)
from reference.structured_certificates import GraphLimits, GraphResourceLimit
from reference.subrank_constraints import (export_convolution_subrank_order,
    export_matrix_subrank_order, export_projective_convolution_subrank_order,
    export_selected_subrank_order, export_subrank_product_order)
from reference.tensor_schemes import TensorScheme
from reference.tensors import FiniteTensor, shared_first_tensor, matrix_multiplication_tensor
from reference.witnessed_constraints import (DetectorConstraint,
    IntegerStateCertificate, StateAssemblyLimit, StateAssemblyLimits,
    WitnessedOrder, WitnessedStateSystem, rank_upper_constraint)


LIMITS = StateAssemblyLimits(max_tensor_entries=2_000_000,
    max_map_entries=1_000_000, max_constraint_entries=1_000_000)
FIELDS = {4: FiniteField(2, (1, 1, 1)), 11: FiniteField(11),
    16: FiniteField(2, (1, 1, 0, 0, 1)), 31: FiniteField(31)}


def literal_export(book, result, limits=LIMITS):
    system = book.system(2, 5)
    system.require_valid(limits)
    for certificate in (result.state, result.dual):
        if certificate is not None:
            if certificate.system != system:
                raise ValueError('certificate belongs to a different exact finite system')
            # The immutable system has just been fully checked. Check the
            # rational inequalities/integer balance without replaying all maps.
            if isinstance(certificate, RationalStateCertificate):
                if not certificate.verify_inequalities():
                    raise ValueError('state fails exact normalization or inequality')
            elif not certificate.verify_balance():
                raise ValueError('dual fails exact integer negative-unit balance')
    rows = []
    for row, family in zip(book.rows, book.families):
        if isinstance(row, WitnessedOrder):
            rows.append(dict(positive=row.positive, negative=row.negative,
                maps=[dict(input_size=m.input_size, rows=m.rows) for m in row.maps],
                label=row.label, family=family))
        else:
            rows.append(dict(test=row.test, product=row.product, family=family))
    data = dict(format='literal-witnessed-system-v1', field=dict(p=book.field.p,
        modulus=book.field.modulus), limits=asdict(limits), d=2, k=5,
        inventory={key: dict(shape=t.shape, coefficients=t.coefficients)
                   for key, t in book.inventory.items()}, rows=rows)
    if result.state is not None:
        data['state'] = {key: str(value) for key, value in result.state.assignment.items()}
    if result.dual is not None:
        data['dual'] = dict(weights=result.dual.weights, gain=result.dual.m)
    return data


def replay(data):
    field = FiniteField(data['field']['p'], tuple(data['field']['modulus']))
    limits = StateAssemblyLimits(**data['limits'])
    registry = tuple((key, FiniteTensor(field, tuple(value['shape']),
        tuple(value['coefficients']))) for key, value in data['inventory'].items())
    rows = tuple(WitnessedOrder(tuple(r['positive']), tuple(r['negative']),
        tuple(LinearMap(field, m['input_size'], tuple(tuple(x) for x in m['rows']))
              for m in r['maps']), r['label']) if 'maps' in r else
        DetectorConstraint(r['test'], r['product']) for r in data['rows'])
    system = WitnessedStateSystem(data['d'], data['k'], registry, rows)
    system.require_valid(limits)
    checked = []
    if 'state' in data:
        state = RationalStateCertificate(system, tuple(Fraction(data['state'][key])
            for key in system.keys))
        if not state.verify_inequalities():
            raise ValueError('replayed state fails exact inequalities/normalization')
        checked.append('finite_feasible')
    if 'dual' in data:
        dual = IntegerStateCertificate(system, tuple(data['dual']['weights']), data['dual']['gain'])
        if not dual.verify_balance():
            raise ValueError('replayed dual fails exact integer balance')
        checked.append('negative_unit_dual')
    if len(checked) > 1:
        raise ValueError('inconsistent simultaneous state and dual')
    return dict(checked=checked, atoms=len(registry), rows=len(rows))


def diagnostics(data):
    """Exact finite-state slacks; no floating-point interpretation is used."""
    if 'state' not in data:
        return dict(status='no exact feasible state in artifact')
    values = {key: Fraction(value) for key, value in data['state'].items()}
    detector = {row['test']: row for row in data['rows'] if 'test' in row}
    active, sectors, fourier, lowers = Counter(), [], [], []
    for index, row in enumerate(data['rows']):
        if 'test' in row:
            slack = values[row['product']]-data['k']*values[row['test']]
        else:
            slack = sum((values[key] for key in row['positive']), Fraction())-sum(
                (values[key] for key in row['negative']), Fraction())
        if slack == 0:
            active[row['family']] += 1
        if row['family'] == 'sector':
            source, retained = row['positive'][0], row['negative'][0]
            info = dict(row=index, source=source, retained=retained,
                source_value=str(values[source]), retained_value=str(values[retained]),
                copies=len(row['positive']), interpolation_slack=str(slack))
            test = detector.get(retained)
            if test is not None:
                product_value = values[test['product']]
                info.update(product=test['product'], product_value=str(product_value),
                    detector_slack=str(product_value-data['k']*values[retained]),
                    detector_tight=product_value == data['k']*values[retained])
            else:
                info['detector_status'] = 'skipped by operational or allocation cap'
            info['branch_lowers'] = [dict(target=r['negative'][0],
                target_value=str(values[r['negative'][0]]),
                slack=str(values[retained]-values[r['negative'][0]]))
                for r in data['rows'] if r['family'] == 'branch_lower'
                and r['positive'] == [retained] and len(r['negative']) == 1]
            sectors.append(info)
        if row['family'] == 'fourier':
            fourier.append(dict(row=index, source=row['positive'][0],
                target=row['negative'][0], copies=len(row['positive']),
                source_value=str(values[row['positive'][0]]),
                target_value=str(values[row['negative'][0]]), slack=str(slack)))
        if row['family'] in ('subrank_lower', 'subrank_product'):
            lowers.append(dict(row=index, source=row['positive'][0],
                value=str(values[row['positive'][0]]),
                independent_units=len(row['negative']), slack=str(slack)))
    return dict(normalization=str(values['unit']),
        nonzero_state_coordinates=sum(value != 0 for value in values.values()),
        active_family_counts=dict(active), sectors=sectors, fourier=fourier,
        actual_subrank_lowers=lowers,
        interpretation='Slacks diagnose only these selected exact rows. Missing context/product or asymptotic relations cannot be inferred from numeric solver statuses.')


def refine(source_path, destination):
    """Add the explicitly diagnosed branch-mixing and powered lower maps."""
    start = monotonic()
    data = json.loads(source_path.read_text())
    original_check = replay(data)
    field = FiniteField(data['field']['p'], tuple(data['field']['modulus']))
    experiment = Experiment(field, start+240, strong=True)
    book = experiment.book
    for key, entry in data['inventory'].items():
        book.register(key, FiniteTensor(field, tuple(entry['shape']), tuple(entry['coefficients'])))
    for row in data['rows']:
        if 'maps' in row:
            actual = WitnessedOrder(tuple(row['positive']), tuple(row['negative']),
                tuple(LinearMap(field, m['input_size'], tuple(tuple(x) for x in m['rows']))
                      for m in row['maps']), row['label'])
        else:
            actual = DetectorConstraint(row['test'], row['product'])
        book.rows.append(actual)
        book.families.append(row['family'])
        book.witness_nodes.append(None)  # Literal replay already checked existing maps.
    matrix = matrix_multiplication_tensor(field, 2, 2, 2)
    new_lowers = []
    for a, h in ((2, 1), (3, 1), (2, 2)):
        construction = ThreeSectorConstruction(field, a, h)
        key = book._canonical.get(construction.retained)
        if key is None:
            continue
        if a >= 3:
            indices = ((0, 1, a-1), (0, construction.middle_z_start-1, construction.right_start),
                (0, construction.middle_z_start, construction.right_start+a-1))
        else:
            indices = ((0, 1), (0, construction.right_start), (0, construction.right_start+1))
        exported = export_selected_subrank_order(construction.retained, indices, key, limits=LIMITS)
        row = experiment.attempt('retained-mixed-subrank-'+key, lambda exported=exported: exported,
            'retained_subrank', dict(a=a, h=h, independent_units=len(indices[0]), selections=indices))
        if row is not None:
            new_lowers.append(exported)
            names = {k: 'refined-context-'+key+'-'+k
                for k in dict.fromkeys(row.positive+row.negative)}
            experiment.attempt('retained-subrank-context-'+key,
                lambda row=row, names=names: lift_witnessed_order(row, book.inventory,
                    matrix, names, limits=LIMITS), 'context:M2:retained_subrank')
    # The e=2 power atoms already exist in the selected family. Connect them
    # with actual product lower maps instead of imposing multiplicative states.
    source = ThreeSectorConstruction(field, 2, 1).source
    source_key = book._canonical.get(source)
    if source_key is not None:
        new_lowers.append(export_convolution_subrank_order(field, 2, 4, 2, source_key, limits=LIMITS))
    for exported in new_lowers:
        tensor = exported.inventory[0][1]
        squared = tensor.tensor_product(tensor)
        key = book._canonical.get(squared)
        if key is not None:
            experiment.attempt('powered-subrank-'+key,
                lambda exported=exported, key=key: export_subrank_product_order(
                    exported, exported, key, limits=LIMITS), 'powered_subrank')
    result = solve_finite_states(book.system(2, 5), timeout_seconds=10, limits=LIMITS)
    refined = literal_export(book, result)
    checked = replay(refined)
    destination.mkdir(parents=True, exist_ok=True)
    (destination/'finite-replay.json').write_text(json.dumps(refined)+'\n')
    (destination/'diagnostics.json').write_text(json.dumps(diagnostics(
        json.loads(json.dumps(refined))), indent=2)+'\n')
    report = dict(status=result.status, reason=result.reason, original=original_check,
        atoms=len(book.inventory), rows=len(book.rows), accepted=experiment.accepted,
        skipped=experiment.skipped, replay=checked, total_seconds=monotonic()-start)
    if result.state is not None:
        report['exact_state'] = {key: str(value) for key, value in result.state.assignment.items()}
    if result.dual is not None:
        report['gain'] = result.dual.m
        report['used_rows'] = [dict(index=i, family=book.families[i], weight=w)
            for i, w in enumerate(result.dual.weights) if w]
    (destination/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: report[key] for key in ('status', 'atoms', 'rows', 'total_seconds')}), flush=True)


def sector_pair_matrix_export(field, source_key, target_key):
    """Compose the padded LM seed square -> star square -> M2 maps."""
    construction = ThreeSectorConstruction(field, 2, 1)
    seed = shared_first_tensor(construction.branches[:2])
    selectors = ((0, 1), (0, 6, 5), (7, 0, 1))
    squared = tuple(tuple(i*size+j for i in selection for j in selection)
        for size, selection in zip(seed.shape, selectors))
    source = seed.tensor_power(2)
    selection_maps = tuple(LinearMap.selection(field, size, indices)
        for size, indices in zip(source.shape, squared))
    star = export_star_matrix_order(field, 'star-square-composition', target_key, limits=LIMITS)
    if source.restrict(*selection_maps) != dict(star.inventory)['star-square-composition']:
        raise ValueError('padded LM seed square is not the supplied actual star square')
    maps = tuple(after.compose(before) for after, before in zip(star.row.maps, selection_maps))
    row = WitnessedOrder((source_key,), (target_key,), maps, 'sector-pair-square-matrix-trace-correction')
    exported = WitnessedOrderExport(row, ((source_key, source),
        (target_key, dict(star.inventory)[target_key])), exponent=2)
    exported.require_valid(LIMITS)
    return exported


def sector_word_columns_export(field, source_key, target_key, word_index=None):
    """Compose a checked padded word rectangle with actual column-sum maps."""
    left, middle = ThreeSectorConstruction(field, 2, 1).branches[:2]
    words = (left.tensor_product(middle), middle.tensor_product(left))
    source = shared_first_tensor(words) if word_index is None else words[word_index]
    LIMITS.tensor(tuple(2*size for size in source.shape))
    selectors = (((0, 2, 1, 3), (8, 4), (10, 11)) if word_index == 1 else
        ((0, 1, 2, 3), (2, 1), (2, 7)))
    extraction = tuple(LinearMap.selection(field, size, indices)
        for size, indices in zip(source.shape, selectors))
    if source.restrict(*extraction) != matrix_multiplication_tensor(field, 2, 2, 1):
        raise ValueError('selected LM word is not the actual matrix-column tensor')
    doubled_extraction = tuple(LinearMap.selection(field, 2*size,
        tuple(copy*size+index for copy in range(2) for index in indices))
        for size, indices in zip(source.shape, selectors))
    columns = export_matrix_column_sum_order(field, 2, 2, 2, 'column-source',
        target_key, limits=LIMITS)
    maps = tuple(after.compose(before) for after, before in zip(columns.row.maps, doubled_extraction))
    target = dict(columns.inventory)[target_key]
    row = WitnessedOrder((source_key, source_key), (target_key,), maps,
        'sector-two-word-matrix-columns' if word_index is None else
        f'sector-word-{word_index}-matrix-columns')
    exported = WitnessedOrderExport(row, ((source_key, source), (target_key, target)))
    exported.require_valid(LIMITS)
    return exported


def retained_two_copy_matrix_export(construction, source_key, target_key):
    """Select two extreme first coordinates and stack two retained row products."""
    return export_retained_two_copy_matrix_order(construction, source_key,
        target_key, limits=LIMITS)


def existing_retained_square_context_export(book, matrix_key):
    """Lift the two-copy R21 row only through already registered atoms."""
    construction = ThreeSectorConstruction(book.field, 2, 1)
    retained = construction.retained
    retained_key = book._canonical.get(retained)
    square_key = book._canonical.get(retained.tensor_power(2))
    product_key = book._canonical.get(book.inventory[matrix_key].tensor_product(retained))
    if any(key is None for key in (retained_key, square_key, product_key)):
        return None
    base = retained_two_copy_matrix_export(construction, retained_key, matrix_key)
    return lift_witnessed_order(base.row, dict(base.inventory), retained,
        {retained_key: square_key, matrix_key: product_key}, limits=LIMITS)


def cross_level(source_path, destination):
    """Add actual sector/word matrix restrictions and the M2 detector."""
    start = monotonic()
    data = json.loads(source_path.read_text())
    original_check = replay(data)
    field = FiniteField(data['field']['p'], tuple(data['field']['modulus']))
    experiment = Experiment(field, float('inf'), max_rows=1000, max_atoms=400, strong=True)
    book = experiment.book
    for key, entry in data['inventory'].items():
        book.register(key, FiniteTensor(field, tuple(entry['shape']), tuple(entry['coefficients'])))
    for r in data['rows']:
        row = WitnessedOrder(tuple(r['positive']), tuple(r['negative']),
            tuple(LinearMap(field, m['input_size'], tuple(tuple(x) for x in m['rows']))
                  for m in r['maps']), r['label']) if 'maps' in r else DetectorConstraint(r['test'], r['product'])
        book.rows.append(row)
        book.families.append(r['family'])
        book.witness_nodes.append(None)
    construction = ThreeSectorConstruction(field, 2, 1)
    source_key = book._canonical.get(construction.retained.tensor_power(2))
    if source_key is None:
        raise ValueError('saved inventory has no actual R(2,1) square')
    matrix_key = book._canonical.get(matrix_multiplication_tensor(field, 2, 2, 2))
    if matrix_key is None:
        raise ValueError('saved inventory has no actual M2')
    exported = export_retained_matrix_order(construction, source_key, matrix_key, limits=LIMITS)
    if exported.row not in book.rows:
        row = experiment.attempt('retained-square-matrix', lambda: exported, 'sector_matrix',
            dict(a=2, h=1, source=source_key, target=matrix_key,
                branch_words=(('left', 'middle'), ('right', 'middle'))))
        if row is None:
            raise ValueError('matrix restriction was not accepted')
    pair_source = shared_first_tensor(construction.branches[:2]).tensor_power(2)
    pair_key = book._canonical.get(pair_source)
    if pair_key is None:
        raise ValueError('saved inventory has no actual shared LM seed square')
    pair_export = sector_pair_matrix_export(field, pair_key, matrix_key)
    if pair_export.row not in book.rows:
        row = experiment.attempt('sector-pair-square-matrix', lambda: pair_export,
            'sector_matrix', dict(source=pair_key, target=matrix_key,
                star_seed_selectors=((0, 1), (0, 6, 5), (7, 0, 1)), actual_composition=True))
        if row is None:
            raise ValueError('composed star matrix restriction was not accepted')
    word_source = dict(sector_word_columns_export(field, 'words', matrix_key).inventory)['words']
    word_key = book._canonical.get(word_source)
    if word_key is None:
        raise ValueError('saved inventory has no actual LM/ML shared-word tensor')
    column_export = sector_word_columns_export(field, word_key, matrix_key)
    if column_export.row not in book.rows:
        row = experiment.attempt('sector-word-matrix-columns', lambda: column_export,
            'matrix_columns', dict(source=word_key, copies=2, target=matrix_key,
                rectangle_selectors=((0, 1, 2, 3), (2, 1), (2, 7))))
        if row is None:
            raise ValueError('word-column matrix restriction was not accepted')
    for index in (0, 1):
        branch_tensor = dict(sector_word_columns_export(field, 'branch', matrix_key,
            index).inventory)['branch']
        branch_key = book._canonical.get(branch_tensor)
        if branch_key is None:
            raise ValueError(f'saved inventory has no actual word-{index} rectangle')
        branch_export = sector_word_columns_export(field, branch_key, matrix_key, index)
        if branch_export.row not in book.rows:
            row = experiment.attempt(f'sector-word-{index}-matrix-columns',
                lambda branch_export=branch_export: branch_export,
                'matrix_columns', dict(source=branch_key, word_index=index, copies=2,
                    actual_rectangle_isomorphism_checked=True))
            if row is None:
                raise ValueError('individual word-column matrix restriction was not accepted')
    for a, h in ((2, 1), (3, 1), (2, 2)):
        retained = ThreeSectorConstruction(field, a, h)
        retained_key = book._canonical.get(retained.retained)
        if retained_key is None:
            continue
        retained_export = retained_two_copy_matrix_export(retained, retained_key, matrix_key)
        if retained_export.row not in book.rows:
            row = experiment.attempt(f'retained-two-copy-matrix-{a}-{h}',
                lambda retained_export=retained_export: retained_export,
                'retained_matrix_rows', dict(a=a, h=h, source=retained_key, copies=2,
                    merges_second_leg_across_copies=True))
            if row is None:
                raise ValueError('retained two-copy matrix restriction was not accepted')
    square_context = existing_retained_square_context_export(book, matrix_key)
    if square_context is None:
        raise ValueError('saved inventory lacks an atom for the fixed R21 square context')
    if square_context.row not in book.rows:
        row = experiment.attempt('retained-square-existing-context', lambda: square_context,
            'context:retained:retained_matrix_rows', dict(a=2, h=1,
                existing_atoms_only=True, context='R21'))
        if row is None:
            raise ValueError('fixed retained-square context was not accepted')
    if not any(isinstance(row, DetectorConstraint) and row.test == matrix_key for row in book.rows):
        detector_start = monotonic()
        book.add_detector(matrix_key, 2)
        experiment.accepted.append(dict(operation='fill-M2-detector', family='detector',
            parameters=dict(test=matrix_key, product=book.rows[-1].product),
            row=len(book.rows)-1, seconds=monotonic()-detector_start))
    matrix_product_key = next(row.product for row in book.rows
        if isinstance(row, DetectorConstraint) and row.test == matrix_key)
    solve_start = monotonic()
    result = solve_finite_states(book.system(2, 5), timeout_seconds=10, limits=LIMITS)
    solve_seconds = monotonic()-solve_start
    final = literal_export(book, result)
    replay_start = monotonic()
    checked = replay(final)
    destination.mkdir(parents=True, exist_ok=True)
    (destination/'finite-replay.json').write_text(json.dumps(final)+'\n')
    (destination/'diagnostics.json').write_text(json.dumps(diagnostics(
        json.loads(json.dumps(final))), indent=2)+'\n')
    report = dict(status=result.status, reason=result.reason, original=original_check,
        atoms=len(book.inventory), rows=len(book.rows), accepted=experiment.accepted,
        skipped=experiment.skipped, replay=checked, replay_seconds=monotonic()-replay_start,
        proposal_and_exact_verification_seconds=solve_seconds, native_LP_seconds=10,
        row_budget=experiment.max_rows, atom_budget=experiment.max_atoms,
        family_counts=dict(Counter(book.families)),
        limits=asdict(LIMITS), total_seconds=monotonic()-start,
        interpretation='These cross-level restrictions refute the proposed scaled flattening assignment. Actual two-copy column rows connect both individual rectangle words and their shared tensor to M2. The LM seed square is related by actual composed maps to star square. LM/ML words remain distinct from LM/RM.')
    if result.state is not None:
        report['exact_state'] = {key: str(value) for key, value in result.state.assignment.items()}
        report['cross_level_values'] = dict(source=source_key, target=matrix_key,
            source_value=str(result.state.assignment[source_key]),
            target_value=str(result.state.assignment[matrix_key]), pair_source=pair_key,
            pair_source_value=str(result.state.assignment[pair_key]), word_source=word_key,
            word_source_value=str(result.state.assignment[word_key]),
            matrix_product=matrix_product_key,
            matrix_product_value=str(result.state.assignment[matrix_product_key]))
    if result.dual is not None:
        report['gain'] = result.dual.m
        report['used_rows'] = [dict(index=i, family=book.families[i], weight=w)
            for i, w in enumerate(result.dual.weights) if w]
    (destination/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: report[key] for key in ('status', 'atoms', 'rows', 'total_seconds')}), flush=True)


def adaptive_cut(source_path, destination, max_rounds=6, wall_seconds=600,
                 fill_existing_detectors=False):
    """Cut exact states with checked coordinate restrictions, then seek a dual."""
    start = monotonic()
    data = json.loads(source_path.read_text())
    original_check = replay(data)
    if 'state' not in data:
        raise ValueError('adaptive cuts require an exact feasible-state artifact')
    field = FiniteField(data['field']['p'], tuple(data['field']['modulus']))
    experiment = Experiment(field, start+wall_seconds, max_rows=1000, max_atoms=400, strong=True)
    book = experiment.book
    for key, entry in data['inventory'].items():
        book.register(key, FiniteTensor(field, tuple(entry['shape']), tuple(entry['coefficients'])))
    for r in data['rows']:
        row = WitnessedOrder(tuple(r['positive']), tuple(r['negative']),
            tuple(LinearMap(field, m['input_size'], tuple(tuple(x) for x in m['rows']))
                  for m in r['maps']), r['label']) if 'maps' in r else DetectorConstraint(r['test'], r['product'])
        book.rows.append(row)
        book.families.append(r['family'])
        book.witness_nodes.append(None)
    values = {key: Fraction(value) for key, value in data['state'].items()}
    # Projective variants repeat the affine vector on these selected q>=4
    # atoms, while checking distinct maps that genuinely use infinity.
    for key, tensor in tuple(book.inventory.items()):
        a, b, z = tensor.shape
        n = min(a, b, field.order+1)
        if n < 2 or z != a+b-1 or tensor != convolution_tensor(field, a, b):
            continue
        experiment.attempt('projective-subrank-'+key,
            lambda key=key, a=a, b=b, n=n: export_projective_convolution_subrank_order(
                field, a, b, n, key, limits=LIMITS), 'projective_subrank',
            dict(a=a, b=b, n=n, finite_nodes=n-1, infinity=True))
    setup = None
    if fill_existing_detectors:
        existing = {row.test for row in book.rows if isinstance(row, DetectorConstraint)}
        for key in tuple(book.inventory):
            if key in existing or 'M2@'+key not in book.inventory:
                continue
            if monotonic() >= experiment.deadline or len(book.rows) >= experiment.max_rows:
                break
            detector_start = monotonic()
            try:
                LIMITS.tensor(tuple(4*x for x in book.inventory[key].shape))
                book.add_detector(key, 2)
            except (StateAssemblyLimit, GraphResourceLimit) as error:
                experiment.skipped.append(dict(operation='fill-existing-detector-'+key,
                    reason=str(error)))
                continue
            experiment.accepted.append(dict(operation='fill-existing-detector-'+key,
                family='detector', parameters=dict(test=key, product=book.rows[-1].product,
                    actual_factor_swap_already_registered=True), row=len(book.rows)-1,
                seconds=monotonic()-detector_start))
            existing.add(key)
        setup_start = monotonic()
        setup_result = solve_finite_states(book.system(2, 5), timeout_seconds=10, limits=LIMITS)
        setup = dict(status=setup_result.status, reason=setup_result.reason,
            proposal_and_exact_verification_seconds=monotonic()-setup_start)
        if setup_result.state is not None:
            values = setup_result.state.assignment
    budget = CoordinateSearchBudget(max_nodes=100_000, max_support=2048, max_diagonal=64, seconds=1)
    rounds, final_result = [], setup_result if fill_existing_detectors else None
    for round_index in range(max_rounds):
        round_start, proposals, cuts = monotonic(), [], 0
        if fill_existing_detectors and setup_result.state is None:
            break
        if round_start >= experiment.deadline:
            break
        for key, tensor in tuple(book.inventory.items()):
            if monotonic() >= experiment.deadline:
                break
            value = values.get(key)
            if value is None or value < 0:
                continue
            goal = value.numerator//value.denominator+1
            if goal > min(tensor.shape):
                continue
            search_start = monotonic()
            candidate = find_coordinate_subrank(tensor, key, goal, budget=budget, limits=LIMITS)
            proposed = dict(key=key, state_value=str(value), goal=goal, status=candidate.status,
                reason=candidate.reason, checked_nodes=candidate.checked_nodes,
                seconds=monotonic()-search_start)
            proposals.append(proposed)
            if candidate.export is None:
                continue
            row = experiment.attempt(f'coordinate-cut-{round_index}-{key}',
                lambda candidate=candidate: candidate.export, 'coordinate_cut',
                dict(round=round_index, state_value=str(value), independent_units=goal,
                    selections=candidate.selections))
            if row is not None:
                cuts += 1
                proposed['accepted_row'] = len(book.rows)-1
        solve_start = monotonic()
        final_result = solve_finite_states(book.system(2, 5), timeout_seconds=10, limits=LIMITS)
        rounds.append(dict(index=round_index, cuts=cuts, proposals=proposals,
            search_seconds=solve_start-round_start, status=final_result.status,
            proposal_and_exact_verification_seconds=monotonic()-solve_start))
        print(json.dumps(dict(round=round_index, cuts=cuts, status=final_result.status,
            rows=len(book.rows))), flush=True)
        if final_result.state is None or cuts == 0:
            break
        values = final_result.state.assignment
    if final_result is None:
        final_result = solve_finite_states(book.system(2, 5), timeout_seconds=10, limits=LIMITS)
    final = literal_export(book, final_result)
    replay_start = monotonic()
    checked = replay(final)
    replay_seconds = monotonic()-replay_start
    destination.mkdir(parents=True, exist_ok=True)
    (destination/'finite-replay.json').write_text(json.dumps(final)+'\n')
    (destination/'diagnostics.json').write_text(json.dumps(diagnostics(
        json.loads(json.dumps(final))), indent=2)+'\n')
    report = dict(status=final_result.status, reason=final_result.reason,
        original=original_check, atoms=len(book.inventory), rows=len(book.rows),
        rounds=rounds, accepted=experiment.accepted, skipped=experiment.skipped,
        coordinate_budget=asdict(budget), max_rounds=max_rounds,
        row_budget=experiment.max_rows, atom_budget=experiment.max_atoms,
        fill_existing_detectors=fill_existing_detectors, setup=setup,
        search_wall_window_seconds=wall_seconds, replay=checked, replay_seconds=replay_seconds,
        total_seconds=monotonic()-start,
        interpretation='Exhaustion concerns only coordinate induced matchings; no general subrank or catalyst upper bound is claimed.')
    if final_result.state is not None:
        report['exact_state'] = {key: str(value) for key, value in final_result.state.assignment.items()}
    if final_result.dual is not None:
        report['gain'] = final_result.dual.m
        report['used_rows'] = [dict(index=i, family=book.families[i], weight=w)
            for i, w in enumerate(final_result.dual.weights) if w]
    (destination/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: report[key] for key in ('status', 'atoms', 'rows', 'total_seconds')}), flush=True)


class Experiment:
    def __init__(self, field, deadline, max_rows=700, max_atoms=300, strong=False):
        self.field, self.deadline = field, deadline
        self.book = ConstraintInventory(field, LIMITS,
            GraphLimits(max_nodes=20_000, max_map_entries=1_000_000))
        self.skipped, self.accepted = [], []
        self.max_rows, self.max_atoms = max_rows, max_atoms
        self.strong = strong
        self.products = {}
        # Separate operational cap for repeated graph materialization in contexts.
        # The construction/export hard cap remains two million coefficients.
        self.max_context_entries = 100_000 if strong else LIMITS.max_tensor_entries
        self.max_upper_terms = 32 if strong else 100

    def attempt(self, operation, build, family, parameters=None):
        start = monotonic()
        if start >= self.deadline:
            self.skipped.append(dict(operation=operation, reason='total wall cap reached'))
            return None
        if len(self.book.rows) >= self.max_rows or len(self.book.inventory) >= self.max_atoms:
            self.skipped.append(dict(operation=operation, reason='finite inventory cap reached'))
            return None
        try:
            exported = build()
            row = self.book.add_export(exported, family)
        except (ValueError, GraphResourceLimit) as error:
            self.skipped.append(dict(operation=operation, parameters=parameters,
                reason=str(error), seconds=monotonic()-start))
            return None
        self.accepted.append(dict(operation=operation, family=family,
            parameters=parameters, row=len(self.book.rows)-1, seconds=monotonic()-start))
        return row

    def selected(self, operation, source_key, indices, target_key, family='branch_lower'):
        tensor = self.book.inventory[source_key]
        maps = tuple(LinearMap.selection(self.field, size, selection)
            for size, selection in zip(tensor.shape, indices))
        target = tensor.restrict(*maps)
        return self.attempt(operation, lambda: WitnessedOrderExport(
            WitnessedOrder((source_key,), (target_key,), maps, operation),
            ((source_key, tensor), (target_key, target))), family,
            dict(source=source_key, indices=indices))

    def fourier(self, branches, label):
        # Enumerate only divisors of the ACTUAL field's multiplicative order.
        count = len(branches)
        periods = [n for n in range(5*count, self.field.order)
            if (self.field.order-1) % n == 0 and self.field.embed(n)]
        if not periods:
            self.skipped.append(dict(operation=label, reason='no available primitive-root period',
                count=count, actual_field_order=self.field.order))
            return None
        period, nodes = periods[0], (count-1)**2+1
        copies = period*nodes
        x, y, z = branches[0].shape
        source_shape = (copies*x, copies*count*y, copies*count*z)
        target_shape = (count*x, count*count*y, count*count*z)
        params = dict(count=count, branch_shape=branches[0].shape, period=period,
            nodes=nodes, copies=copies, source_shape=source_shape,
            source_entries=prod(source_shape), target_shape=target_shape)
        try:
            LIMITS.tensor(source_shape)
            LIMITS.tensor(target_shape)
            LIMITS.maps(source_shape, target_shape)
        except StateAssemblyLimit as error:
            self.skipped.append(dict(operation=label, reason=str(error), parameters=params))
            return None
        row = self.attempt(label, lambda: export_finite_separation_order(
            FiniteSeparation(branches, period=period,
                root=primitive_root(self.field, period)), label+'-shared',
            label+'-direct', limits=LIMITS), 'fourier', params)
        if row is not None and self.strong:
            for index, branch in enumerate(branches):
                self.selected(label+'-seed-branch-'+str(index), row.positive[0],
                    (tuple(range(x)), tuple(range(index*y, (index+1)*y)),
                     tuple(range(index*z, (index+1)*z))), label+'-seed-branch-'+str(index))
            dot = FiniteTensor.from_function(self.field, (1, count, count),
                lambda i, j, k: int(j == k))
            components = tuple((label+'-dot-branch-'+str(i), b.tensor_product(dot))
                for i, b in enumerate(branches))
            whole = self.book.inventory[row.negative[0]]
            component_keys = tuple(key for key, _ in components)
            maps = tuple(LinearMap.identity(self.field, size) for size in whole.shape)
            items = ((row.negative[0], whole),)+components
            for reverse in (False, True):
                self.attempt(label+'-split-'+str(reverse), lambda reverse=reverse:
                    WitnessedOrderExport(WitnessedOrder(
                        component_keys if reverse else row.negative,
                        row.negative if reverse else component_keys, maps,
                        'Fourier-direct-sum-components'), items), 'fourier_additivity')
        return row

    def seeds(self):
        for a, h, exponent in ((2, 1, 1), (2, 1, 2), (3, 1, 1), (2, 2, 1)):
            construction = ThreeSectorConstruction(self.field, a, h)
            name = f'sector-{a}-{h}-e{exponent}'
            row = self.attempt(name, lambda c=construction, n=name, e=exponent:
                export_sector_order(c, n+'-source', n+'-retained',
                    exponent=e, limits=LIMITS), 'sector', dict(a=a, h=h, exponent=exponent))
            if row is not None and exponent == 1:
                for branch in ('left', 'middle', 'right'):
                    if branch == 'left':
                        indices = (tuple(range(a)), tuple(range(h)), tuple(range(a+h-1)))
                    elif branch == 'middle':
                        indices = (tuple(range(a-1, -1, -1)),
                            tuple(range(h, h+a+h-1)),
                            tuple(range(h+a-1, h+a-1+h)))
                    else:
                        indices = (tuple(range(a)),
                            tuple(range(2*h+a-1, 3*h+a-1)),
                            tuple(range(2*h+a-1, 2*h+a-1+a+h-1)))
                    self.selected(name+'-'+branch, row.negative[0], indices,
                        name+'-'+branch+'-compact')
        for d, e in ((1, 1), (2, 1), (1, 2), (2, 2)):
            name = f'determinant-{d}-{e}'
            row = self.attempt(name, lambda d=d, e=e, n=name:
                DeterminantFiltration(self.field, d, e, LIMITS).export_order(
                    n+'-source', n+'-graded', limits=LIMITS), 'determinant', dict(d=d, e=e))
            if row is not None:
                self.selected(name+'-quotient', row.negative[0],
                    (tuple(range(d+1)), tuple(range(e+2)), tuple(range(d+e+2))),
                    name+'-quotient')
                self.selected(name+'-kernel', row.negative[0],
                    (tuple(range(d+1)), tuple(range(e+2, 2*(e+1))),
                     tuple(range(d+e+2, 2*(d+e+1)))), name+'-kernel')
        sector = ThreeSectorConstruction(self.field, 2, 1)
        actual_branches = sector.branches[:2]
        compact = convolution_tensor(self.field, 2, 1)
        # The second compact branch is an explicitly coefficient-distinct,
        # nonscalar branch, not the old scalar (1,3) testing family.
        twisted = FiniteTensor.from_function(self.field, (2, 1, 2),
            lambda i, j, k: int(i+k == 1))
        for label, branches in (('sector-pair', actual_branches),
                                ('compact-pair', (compact, twisted))):
            typed_holder = []
            def build(branches=branches, label=label):
                typed = export_exact_type_order(branches, (1, 1),
                    label+'-power', label+'-words', limits=LIMITS)
                typed_holder.append(typed)
                return typed.restriction
            row = self.attempt(label+'-type', build, 'exact_type',
                dict(counts=(1, 1), branch_shape=branches[0].shape,
                    branch_coefficients=[b.coefficients for b in branches]))
            self.fourier(branches, label+'-Fourier-seed')
            if row is not None:
                self.fourier(typed_holder[0].branches, label+'-Fourier-words')
                for index, branch in enumerate(typed_holder[0].branches):
                    x, y, z = branch.shape
                    self.selected(label+'-word-'+str(index), row.negative[0],
                        (tuple(range(x)), tuple(range(index*y, (index+1)*y)),
                         tuple(range(index*z, (index+1)*z))),
                        label+'-word-'+str(index))

    def connect(self):
        book, field = self.book, self.field
        matrix = matrix_multiplication_tensor(field, 2, 2, 2)
        book.register('M2', matrix)
        construction = ThreeSectorConstruction(field, 2, 1)
        squared_key = book._canonical.get(construction.retained.tensor_power(2))
        if squared_key is not None:
            self.attempt('retained-square-matrix', lambda: export_retained_matrix_order(
                construction, squared_key, 'M2', limits=LIMITS), 'sector_matrix',
                dict(a=2, h=1, branch_words=(('left', 'middle'), ('right', 'middle'))))
        pair_square = shared_first_tensor(construction.branches[:2]).tensor_power(2)
        pair_key = book._canonical.get(pair_square)
        if pair_key is not None:
            self.attempt('sector-pair-square-matrix', lambda: sector_pair_matrix_export(
                field, pair_key, 'M2'), 'sector_matrix', dict(actual_composition=True))
        word_tensor = dict(sector_word_columns_export(field, 'words', 'M2').inventory)['words']
        word_key = book._canonical.get(word_tensor)
        if word_key is not None:
            self.attempt('sector-word-matrix-columns', lambda: sector_word_columns_export(
                field, word_key, 'M2'), 'matrix_columns', dict(copies=2))
        for index in (0, 1):
            branch_tensor = dict(sector_word_columns_export(field, 'branch', 'M2',
                index).inventory)['branch']
            branch_key = book._canonical.get(branch_tensor)
            if branch_key is not None:
                self.attempt(f'sector-word-{index}-matrix-columns',
                    lambda branch_key=branch_key, index=index: sector_word_columns_export(
                        field, branch_key, 'M2', index), 'matrix_columns',
                    dict(word_index=index, copies=2))
        for a, h in ((2, 1), (3, 1), (2, 2)):
            retained = ThreeSectorConstruction(field, a, h)
            retained_key = book._canonical.get(retained.retained)
            if retained_key is not None:
                self.attempt(f'retained-two-copy-matrix-{a}-{h}',
                    lambda retained=retained, retained_key=retained_key:
                        retained_two_copy_matrix_export(retained, retained_key, 'M2'),
                    'retained_matrix_rows', dict(a=a, h=h, copies=2))
        detector_start = monotonic()
        book.add_detector('M2', 2)
        self.accepted.append(dict(operation='detector-M2', family='detector',
            parameters=dict(test='M2', product=book.rows[-1].product),
            row=len(book.rows)-1, seconds=monotonic()-detector_start))
        subranks = {}
        if self.strong:
            for a, b in ((2, 2), (2, 3), (3, 3)):
                book.register(f'C{a}{b}', convolution_tensor(field, a, b))
            matrix_lower = export_matrix_subrank_order(field, 2, 'M2', limits=LIMITS)
            self.attempt('subrank-M2', lambda: matrix_lower, 'subrank_lower', dict(d=2, n=2))
            subranks['M2'] = matrix_lower
            for key, tensor in tuple(book.inventory.items()):
                a, b, z = tensor.shape
                if z == a+b-1 and tensor == convolution_tensor(field, a, b):
                    n = min(a, b, field.order)
                    exported = export_convolution_subrank_order(field, a, b, n, key, limits=LIMITS)
                    self.attempt('subrank-'+key, lambda exported=exported: exported,
                        'subrank_lower', dict(a=a, b=b, n=n))
                    subranks[key] = exported
        initial_rows = tuple(book.rows)
        initial_atoms = tuple(book.inventory)
        # All source and target atoms get M2 detectors when their exact product fits.
        for key in initial_atoms:
            if key == 'M2':
                continue
            try:
                product_shape = tuple(4*x for x in book.inventory[key].shape)
                LIMITS.tensor(product_shape)
                if prod(product_shape) > self.max_context_entries:
                    raise StateAssemblyLimit('context operational coefficient cap reached')
                book.add_detector(key, 2)
                self.products[(matrix, book.inventory[key])] = book.inventory[book.rows[-1].product]
            except StateAssemblyLimit as error:
                self.skipped.append(dict(operation='detector-'+key, reason=str(error)))
        for index, row in enumerate(initial_rows):
            if not isinstance(row, WitnessedOrder):
                continue
            keys = tuple(dict.fromkeys(row.positive+row.negative))
            names = {key: f'M2-context-{index}-{key}' for key in keys}
            source_shape = tuple(4*sum(book.inventory[key].shape[axis]
                for key in row.positive) for axis in range(3))
            target_shape = tuple(4*sum(book.inventory[key].shape[axis]
                for key in row.negative) for axis in range(3))
            def build_lift(row=row, names=names, source_shape=source_shape, target_shape=target_shape):
                if max(prod(source_shape), prod(target_shape)) > self.max_context_entries:
                    raise StateAssemblyLimit('context operational coefficient cap reached')
                return lift_witnessed_order(row, book.inventory, matrix, names, limits=LIMITS)
            lifted = self.attempt('M2-context-'+str(index), build_lift,
                'context:M2', dict(seed_row=index))
            if lifted is not None:
                for old_key, new_key in zip(row.positive+row.negative,
                                           lifted.positive+lifted.negative):
                    self.products[(book.inventory[old_key], matrix)] = book.inventory[new_key]
            for side, side_keys in (('positive', row.positive), ('negative', row.negative)):
                if len(side_keys) < 2:
                    continue
                shape = tuple(sum(book.inventory[key].shape[a] for key in side_keys)
                    for a in range(3))
                name = f'sum-{index}-{side}'
                def build_sum(side_keys=side_keys, shape=shape, name=name, reverse=False):
                    LIMITS.tensor(shape)
                    LIMITS.maps(shape, shape)
                    if prod(shape) > self.max_context_entries:
                        raise StateAssemblyLimit('context operational coefficient cap reached')
                    tensor = FiniteTensor.direct_sum(book.inventory[k] for k in side_keys)
                    maps = tuple(LinearMap.identity(field, size) for size in shape)
                    items = tuple((k, book.inventory[k]) for k in dict.fromkeys(side_keys))+((name, tensor),)
                    return WitnessedOrderExport(WitnessedOrder(side_keys if reverse else (name,),
                        (name,) if reverse else side_keys, maps, 'direct-sum-identity'), items)
                self.attempt(name+'-forward', build_sum, 'additivity')
                self.attempt(name+'-reverse', lambda build_sum=build_sum:
                    build_sum(reverse=True), 'additivity')
        for key in initial_atoms:
            if key == 'M2':
                continue
            for reverse in (False, True):
                def swap(key=key, reverse=reverse):
                    shape = tuple(4*x for x in book.inventory[key].shape)
                    LIMITS.tensor(shape)
                    LIMITS.maps(shape, shape)
                    first = book.inventory[key]
                    for pair in ((first, matrix), (matrix, first)):
                        if pair not in self.products:
                            self.products[pair] = pair[0].tensor_product(pair[1])
                    source, target = self.products[(first, matrix)], self.products[(matrix, first)]
                    names = (key+'@M2', 'M2@'+key)
                    maps = tuple(LinearMap.selection(field, a*b,
                        tuple(j*a+i for i in range(a) for j in range(b))
                        if reverse else tuple(i*b+j for j in range(b) for i in range(a)))
                        for a, b in zip(first.shape, matrix.shape))
                    if reverse:
                        source, target, names = target, source, names[::-1]
                    return WitnessedOrderExport(WitnessedOrder((names[0],), (names[1],),
                        maps, 'tensor-factor-swap'), tuple(zip(names, (source, target))))
                self.attempt(f'swap-{key}-{reverse}', swap, 'factor_swap')
        for key, tensor in tuple(book.inventory.items()):
            if key == 'unit':
                continue
            # Known convolutions use their full n-unit lower row in this pass.
            if any(tensor.coefficients) and key not in subranks:
                self.attempt('lower-'+key, lambda key=key, tensor=tensor:
                    scalar_lower_export(key, tensor), 'scalar_lower')
            # Coordinate-pair schemes supply safe rank upper rows. No Strassen.
            # Preflight term counts BEFORE allocating the cubic unit direct sum.
            z = tensor.shape[2]
            terms = sum(any(tensor.coefficients[i:i+z])
                for i in range(0, len(tensor.coefficients), z)) if z else 0
            if terms > self.max_upper_terms or terms**3 > LIMITS.max_tensor_entries:
                self.skipped.append(dict(operation='upper-'+key,
                    reason='rank-upper dense/term cap', terms=terms))
                continue
            def upper(key=key, tensor=tensor):
                scheme = TensorScheme.from_tensor(tensor)
                return WitnessedOrderExport(rank_upper_constraint('unit', key, scheme),
                    (('unit', book.inventory['unit']), (key, tensor)))
            self.attempt('upper-'+key, upper, 'rank_upper')
        if self.strong:
            for key, exported in subranks.items():
                if key == 'M2':
                    continue
                for reverse in (False, True):
                    first, second = (exported, matrix_lower) if reverse else (matrix_lower, exported)
                    tensor = first.inventory[0][1].tensor_product(second.inventory[0][1])
                    if tensor not in book._canonical:
                        continue
                    product_key = book._canonical[tensor]
                    self.attempt('product-subrank-'+product_key,
                        lambda first=first, second=second, product_key=product_key:
                            export_subrank_product_order(first, second, product_key, limits=LIMITS),
                        'subrank_product', dict(first=first.row.positive[0],
                            second=second.row.positive[0], n=len(first.row.negative)*len(second.row.negative)))
        square_context = existing_retained_square_context_export(book, 'M2')
        if square_context is not None:
            self.attempt('retained-square-existing-context', lambda: square_context,
                'context:retained:retained_matrix_rows', dict(a=2, h=1,
                    existing_atoms_only=True, context='R21'))


def run(fields, destination, wall_seconds, strong=False):
    start, reports = monotonic(), []
    deadline = start+wall_seconds
    destination.mkdir(parents=True, exist_ok=True)
    for order in fields:
        field_start = monotonic()
        experiment = Experiment(FIELDS[order], deadline, strong=strong)
        experiment.seeds()
        experiment.connect()
        book = experiment.book
        built = monotonic()
        result = solve_finite_states(book.system(2, 5), timeout_seconds=10, limits=LIMITS)
        report = dict(field_order=order, field=dict(p=book.field.p, modulus=book.field.modulus),
            atoms=len(book.inventory), rows=len(book.rows), status=result.status,
            reason=result.reason, native_LP_seconds=10, limits=asdict(LIMITS),
            family_counts=dict(Counter(book.families)), accepted=experiment.accepted,
            skipped=experiment.skipped, construction_seconds=built-field_start,
            proposal_and_verification_seconds=monotonic()-built,
            supplied_rank_seven_control=False)
        report['strong_lower_and_Fourier_splitting'] = strong
        report['context_operational_entries_cap'] = experiment.max_context_entries
        report['operational_cap_scope'] = 'detector construction, context lifts and direct-sum identities; factor swaps use the hard tensor/map caps'
        report['upper_terms_cap'] = experiment.max_upper_terms
        report['auxiliary_scope'] = ('Complete convolution C(a,b) alone is excluded as a '
            'fixed k=5 auxiliary by the separately audited output-substitution/commutator '
            'inequality; its rows here connect finite states. Retained sector, determinant '
            'graded, exact-type and Fourier direct tensors remain distinct varied tests.')
        if result.state is not None:
            report['exact_state'] = {key: str(value) for key, value in result.state.assignment.items()}
        if result.dual is not None:
            report['gain'] = result.dual.m
            report['used_rows'] = [{ 'index': i, 'family': book.families[i], 'weight': w}
                for i, w in enumerate(result.dual.weights) if w]
        data = literal_export(book, result)
        replay_start = monotonic()
        report['replay'] = replay(data)
        report['replay_seconds'] = monotonic()-replay_start
        folder = destination/f'f{order}'
        folder.mkdir(exist_ok=True)
        (folder/'finite-replay.json').write_text(json.dumps(data)+'\n')
        report['total_seconds'] = monotonic()-field_start
        (folder/'report.json').write_text(json.dumps(report, indent=2)+'\n')
        reports.append(report)
        print(json.dumps({k: report[k] for k in ('field_order', 'atoms', 'rows', 'status',
            'family_counts', 'total_seconds')}), flush=True)
        if monotonic() >= deadline:
            break
    (destination/'summary.json').write_text(json.dumps(dict(reports=reports,
        total_seconds=monotonic()-start,
        interpretation='Exact feasible states exclude negative-unit duals only for these finite selected rows; no arbitrary catalyst exclusion.'), indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fields', type=int, nargs='+', choices=tuple(FIELDS), default=[11, 4, 16, 31])
    parser.add_argument('--output', type=Path, default=Path(tempfile.gettempdir())/'variable-proof-experiments')
    parser.add_argument('--wall-seconds', type=float, default=1000,
        help='Construction/search window checked between operations; replay can extend beyond it')
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--replay', type=Path)
    parser.add_argument('--strong', action='store_true',
        help='Add actual convolution/matrix subrank and Fourier-component identity rows')
    modes.add_argument('--diagnose', type=Path,
        help='Read a literal replay artifact and print exact state slacks')
    modes.add_argument('--refine', type=Path,
        help='Add diagnosed retained and powered subrank rows to an existing replay artifact')
    modes.add_argument('--cut', type=Path,
        help='Adaptively add actual coordinate subrank rows that violate the saved exact state')
    modes.add_argument('--cross-level', type=Path,
        help='Add actual sector-square/rectangle-column matrix rows and the M2 detector to a saved literal system')
    parser.add_argument('--cut-rounds', type=int, default=6)
    parser.add_argument('--fill-existing-detectors', action='store_true',
        help='Before cutting, add missing M2 detectors whose exact reverse product is already registered')
    args = parser.parse_args()
    if args.wall_seconds <= 0 or args.cut_rounds < 0:
        parser.error('--wall-seconds must be positive and --cut-rounds nonnegative')
    if args.cross_level:
        cross_level(args.cross_level, args.output)
    elif args.cut:
        adaptive_cut(args.cut, args.output, args.cut_rounds, args.wall_seconds,
            args.fill_existing_detectors)
    elif args.refine:
        refine(args.refine, args.output)
    elif args.diagnose:
        print(json.dumps(diagnostics(json.loads(args.diagnose.read_text())), indent=2))
    elif args.replay:
        print(json.dumps(replay(json.loads(args.replay.read_text()))))
    else:
        run(args.fields, args.output, args.wall_seconds, strong=args.strong)
