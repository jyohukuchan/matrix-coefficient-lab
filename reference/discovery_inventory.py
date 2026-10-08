"""Connected, capped finite experiments; every ordinary row has actual maps.

These are finite row families, not the proof's missing spectral witness.
An exact feasible state rules out a dual only for the selected rows.
"""

from dataclasses import dataclass
from math import prod

from .constraint_generators import (WitnessedOrderExport, export_sector_order,
                                    lift_witnessed_order)
from .determinant_filtration import DeterminantFiltration
from .fields import FiniteField
from .maps import LinearMap
from .projective_convolution import projective_convolution_scheme
from .proof_pipeline import ConstraintInventory, PipelineConfig, generate_proof_inventory
from .schemes import naive_scheme, strassen_scheme
from .sectors import ThreeSectorConstruction, convolution_tensor
from .structured_certificates import GraphResourceLimit
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, matrix_multiplication_tensor
from .witnessed_constraints import (DEFAULT_LIMITS, StateAssemblyLimit,
                                    WitnessedOrder, rank_upper_constraint)


@dataclass(frozen=True)
class DiscoveryBudget:
    max_rows: int = 256
    max_atoms: int = 128
    max_upper_terms: int = 64
    max_lift_entries: int = 100_000

    def __post_init__(self):
        if any(type(v) is not int or v < 1 for v in
               (self.max_rows, self.max_atoms, self.max_upper_terms, self.max_lift_entries)):
            raise ValueError('discovery caps must be positive integers')


def star_tensor(field):
    """The concise symmetric singular (2,3,3) pencil."""
    support = {(0, 0, 1), (0, 1, 0), (1, 0, 2), (1, 2, 0)}
    return FiniteTensor.from_function(field, (2, 3, 3), lambda i, j, k: int((i, j, k) in support))


def scalar_lower_export(key, tensor, unit_key='unit'):
    """Extract one unit from an actual nonzero coefficient, with its inverse."""
    index = next((i for i, value in enumerate(tensor.coefficients) if value), None)
    if index is None:
        raise ValueError('a zero tensor has no scalar lower witness')
    z, y = tensor.shape[2], tensor.shape[1]
    i, remainder = divmod(index, y*z)
    j, k = divmod(remainder, z)
    scalar = tensor.field.inv(tensor.coefficients[index])
    maps = tuple(LinearMap(tensor.field, size,
                           (tuple((scalar if axis == 2 else 1) if x == coordinate else 0
                                  for x in range(size)),))
                 for axis, (size, coordinate) in enumerate(zip(tensor.shape, (i, j, k))))
    items = ((key, tensor),) if key == unit_key else ((key, tensor), (unit_key, FiniteTensor.unit(tensor.field)))
    if key == unit_key and tensor != FiniteTensor.unit(tensor.field):
        raise ValueError('the unit key must name the actual unit')
    return WitnessedOrderExport(WitnessedOrder((key,), (unit_key,), maps, 'nonzero-coefficient-unit'), items)


def commute_product_export(first_key, first, second_key, second, *, reverse=False):
    """An actual coordinate permutation, not equality of product orderings."""
    source, target = first.tensor_product(second), second.tensor_product(first)
    names = f'{first_key}@{second_key}', f'{second_key}@{first_key}'
    maps = tuple(LinearMap.selection(first.field, a*b, tuple(i*b+j for j in range(b) for i in range(a)))
                 for a, b in zip(first.shape, second.shape))
    if reverse:
        maps = tuple(LinearMap.selection(first.field, a*b, tuple(j*a+i for i in range(a) for j in range(b)))
                     for a, b in zip(first.shape, second.shape))
        source, target = target, source
        names = names[::-1]
    return WitnessedOrderExport(WitnessedOrder((names[0],), (names[1],), maps, 'tensor-factor-swap'),
                                tuple(zip(names, (source, target))))


def generate_discovery_inventory(field=FiniteField(11), *, include_known_control=False,
                                 contexts=('M2', 'C22', 'star'), budget=DiscoveryBudget(),
                                 limits=DEFAULT_LIMITS):
    if not isinstance(field, FiniteField) or type(include_known_control) is not bool:
        raise ValueError('supply an explicit finite field and a boolean control flag')
    contexts = tuple(contexts)
    if len(set(contexts)) != len(contexts) or any(x not in ('M2', 'C22', 'star') for x in contexts):
        raise ValueError('choose distinct supported tensor contexts')
    skipped = []
    if field == FiniteField(11):
        book = generate_proof_inventory(PipelineConfig(include_known_control=False), limits=limits)
    else:
        book = ConstraintInventory(field, limits)
        try:
            book.add_export(export_sector_order(ThreeSectorConstruction(field, 2, 1),
                                                 'sector-source', 'sector-retained', limits=limits), 'sector')
        except (ValueError, StateAssemblyLimit) as error:
            skipped.append({'operation': 'sector', 'reason': str(error)})
        # The available exact-type/Fourier control needs a primitive tenth root.
        skipped.append({'operation': 'exact_type_fourier',
                        'reason': 'this control is generated over F11; no field extension is implicit'})
        try:
            book.add_export(DeterminantFiltration(field, 1, 1, limits).export_order(
                'determinant-source', 'determinant-graded', limits=limits), 'determinant')
        except (ValueError, StateAssemblyLimit) as error:
            skipped.append({'operation': 'determinant', 'reason': str(error)})
        book.add_detector('unit', 2)

    if len(book.rows) > budget.max_rows or len(book.inventory) > budget.max_atoms:
        raise StateAssemblyLimit('seed inventory exceeds discovery budget')

    def add(export, family):
        if len(book.rows) >= budget.max_rows:
            skipped.append({'operation': family, 'reason': 'row cap reached'})
            return
        new_atoms = len({tensor for _, tensor in export.inventory if tensor not in book._canonical})
        if len(book.inventory)+new_atoms > budget.max_atoms:
            skipped.append({'operation': family, 'reason': 'atom cap reached'})
            return
        try:
            book.add_export(export, family)
        except (StateAssemblyLimit, GraphResourceLimit) as error:
            skipped.append({'operation': family, 'reason': str(error)})

    def register(name, tensor):
        if tensor not in book._canonical and len(book.inventory) >= budget.max_atoms:
            raise StateAssemblyLimit('mandatory discovery seed exceeds atom cap')
        return book.register(name, tensor)

    schemes = {}
    matrix = strassen_scheme(field) if include_known_control else naive_scheme(field, 2, 2, 2)
    matrix = matrix.as_tensor_scheme()
    matrix_key = register('M2', matrix.target)
    schemes[matrix_key] = matrix
    if include_known_control:
        add(WitnessedOrderExport(rank_upper_constraint('unit', matrix_key, matrix),
                                 (('unit', book.inventory['unit']), (matrix_key, matrix.target))), 'known_control')
    for a, b in ((2, 2), (2, 3), (3, 3)):
        target = convolution_tensor(field, a, b)
        try:
            scheme = projective_convolution_scheme(field, a, b)
        except ValueError as error:
            skipped.append({'operation': f'projective-C{a}{b}', 'reason': str(error),
                            'fallback': 'checked coordinate-pair scheme'})
            scheme = TensorScheme.from_tensor(target)
        key = register(f'C{a}{b}', target)
        schemes[key] = scheme
    star = star_tensor(field)
    star_key = register('star', star)
    schemes[star_key] = TensorScheme.from_tensor(star)
    seed_rows = tuple(zip(book.rows, book.families))
    context_tensors = {'M2': matrix.target, 'C22': convolution_tensor(field, 2, 2), 'star': star}
    for index, (row, family) in enumerate(seed_rows):
        if not isinstance(row, WitnessedOrder):
            continue
        # Register whole direct sums as atoms and enforce both identity directions.
        for side, keys in (('positive', row.positive), ('negative', row.negative)):
            if len(keys) < 2:
                continue
            shape = tuple(sum(book.inventory[key].shape[axis] for key in keys) for axis in range(3))
            try:
                limits.tensor(shape)
                direct = FiniteTensor.direct_sum(book.inventory[key] for key in keys)
                name = f'sum-{index}-{side}'
                maps = tuple(LinearMap.identity(field, size) for size in shape)
                items = tuple((key, book.inventory[key]) for key in dict.fromkeys(keys))+((name, direct),)
                add(WitnessedOrderExport(WitnessedOrder((name,), keys, maps, 'direct-sum-identity'), items), 'additivity')
                add(WitnessedOrderExport(WitnessedOrder(keys, (name,), maps, 'direct-sum-identity'), items), 'additivity')
            except StateAssemblyLimit as error:
                skipped.append({'operation': f'sum-{index}-{side}', 'reason': str(error)})
        for context in contexts:
            context_tensor = context_tensors[context]
            shape = tuple(sum(book.inventory[key].shape[axis] for key in row.positive)*context_tensor.shape[axis]
                          for axis in range(3))
            if prod(shape) > budget.max_lift_entries:
                skipped.append({'operation': f'lift-{index}-{context}', 'reason': 'lift construction cap reached'})
                continue
            names = {key: f'lift-{index}-{context}-{key}' for key in dict.fromkeys(row.positive+row.negative)}
            try:
                add(lift_witnessed_order(row, book.inventory, context_tensor, names, limits=limits),
                    f'context:{context}:{family}')
            except StateAssemblyLimit as error:
                skipped.append({'operation': f'lift-{index}-{context}', 'reason': str(error)})

    # Context lifts use T tensor M2, while detectors use M2 tensor T. These
    # coefficient arrays are usually different; explicit swaps connect them.
    if 'M2' in contexts:
        keys = tuple(dict.fromkeys(key for row, _ in seed_rows if isinstance(row, WitnessedOrder)
                                  for key in row.positive+row.negative))
        for key in keys:
            for reverse in (False, True):
                try:
                    limits.tensor(tuple(a*b for a, b in zip(book.inventory[key].shape, matrix.target.shape)))
                    add(commute_product_export(key, book.inventory[key], 'M2-context', matrix.target,
                                               reverse=reverse), 'factor_swap')
                except StateAssemblyLimit as error:
                    skipped.append({'operation': f'swap-{key}', 'reason': str(error)})

    # Pin detectors on all initial proof atoms and selected non-matrix tensors.
    tests = tuple(dict.fromkeys(tuple(key for key in book.inventory if not key.startswith(('lift-', 'sum-')) and '@' not in key)
                               +(star_key,)))
    existing = {getattr(row, 'test', None) for row in book.rows}
    for key in tests:
        if key in existing or key == matrix_key:
            continue
        if len(book.rows) >= budget.max_rows or len(book.inventory) >= budget.max_atoms:
            skipped.append({'operation': f'detector-{key}', 'reason': 'inventory cap reached'})
            continue
        try:
            book.add_detector(key, 2)
        except (StateAssemblyLimit, GraphResourceLimit) as error:
            skipped.append({'operation': f'detector-{key}', 'reason': str(error)})
            continue
        product_key = book._canonical[matrix.target.tensor_product(book.inventory[key])]
        if key in schemes:
            schemes[product_key] = matrix.tensor_product(schemes[key])

    for key, tensor in tuple(book.inventory.items()):
        if key == 'unit':
            continue
        if any(tensor.coefficients):
            add(scalar_lower_export(key, tensor), 'scalar_lower')
        scheme = schemes.get(key) or TensorScheme.from_tensor(tensor)
        if scheme.terms > budget.max_upper_terms:
            skipped.append({'operation': f'upper-{key}', 'reason': 'supplied decomposition term cap reached'})
            continue
        add(WitnessedOrderExport(rank_upper_constraint('unit', key, scheme),
                                 (('unit', book.inventory['unit']), (key, tensor))), 'rank_upper')
    return book, skipped
