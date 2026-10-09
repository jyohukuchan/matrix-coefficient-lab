"""Sparse integral Koszul certificate for Unit direct_sum (M2 tensor Star).

Every replay operation is invertible over the integers and hence every field.
The resulting rank1054 / rank-one-factor70 proves border rank at least16.
The arbitrary-catalyst transfer is separate mathematical research, not a
consequence of finite-field sampling or an ordinary normalized-state row.
"""
from dataclasses import dataclass
from itertools import combinations
import json
from pathlib import Path


class StarKoszulLimit(ValueError):
    """A declared resource cap prevents completing the integral replay."""


@dataclass(frozen=True)
class StarKoszulLimits:
    max_nonzero_entries: int = 10000
    max_coordinates: int = 4096
    max_component_entries: int = 144
    max_components: int = 2000
    max_operations: int = 10000
    max_integer_bits: int = 4096
    max_certificate_bytes: int = 100_000

    def __post_init__(self):
        if any(type(v) is not int or v < 0 for v in self.__dict__.values()):
            raise ValueError('Koszul replay caps must be nonnegative integers')


DEFAULT_LIMITS = StarKoszulLimits()


def _preflight(limits):
    if not isinstance(limits, StarKoszulLimits):
        raise ValueError('supply explicit star Koszul replay limits')
    if (limits.max_nonzero_entries < 2310 or limits.max_coordinates < 3276
            or limits.max_integer_bits < 1):
        raise StarKoszulLimit('sparse matrix, coordinate or integer cap reached')


def star_integer_target_support():
    """Reconstruct the literal (9,13,13) target from support triples."""
    support = {(0, 0, 0)}
    for i in range(2):
        for j in range(2):
            for h in range(2):
                for a in range(2):
                    for b, c in ((0, a+1), (a+1, 0)):
                        support.add((1+(2*i+j)*2+a, 1+(2*j+h)*3+b, 1+(2*i+h)*3+c))
    if len(support) != 33:
        raise ValueError('the independently reconstructed target has incorrect support')
    return tuple(sorted(support))


def star_integer_koszul(*, limits=DEFAULT_LIMITS):
    """Sparse signed matrix; row=wedge5*13+z, column=wedge4*13+y."""
    _preflight(limits)
    domain = tuple(combinations(range(9), 4))
    codomain = {subset: i for i, subset in enumerate(combinations(range(9), 5))}
    matrix = {}
    for x, y, z in star_integer_target_support():
        for column, subset in enumerate(domain):
            if x in subset:
                continue
            row = codomain[tuple(sorted((*subset, x)))]
            coordinate = (row*13+z, column*13+y)
            if coordinate in matrix:
                raise ValueError('unexpected overlapping integer Koszul entries')
            matrix[coordinate] = -1 if sum(v < x for v in subset) % 2 else 1
    if len(matrix) != 2310:
        raise ValueError('the signed sparse matrix has incorrect support size')
    return matrix


def _certificate(certificate, limits):
    if certificate is None:
        path = Path(__file__).with_name('star_koszul.json')
        if path.stat().st_size > limits.max_certificate_bytes:
            raise StarKoszulLimit('certificate byte cap reached')
        certificate = json.loads(path.read_text())
    if (not isinstance(certificate, dict)
            or certificate.get('format') != 'star-full-koszul-unimodular-components-v1'
            or certificate.get('koszul_degree') != 4
            or certificate.get('rank_one_factor') != 70
            or certificate.get('unit_rank') != 1054
            or certificate.get('matrix_shape') != [1638, 1638]):
        raise ValueError('unsupported star certificate format or asserted parameters')
    components = certificate.get('components')
    if not isinstance(components, (list, tuple)):
        raise ValueError('supply support-component certificates')
    if len(components) > limits.max_components:
        raise StarKoszulLimit('component count cap reached')
    operation_count = coordinate_count = 0
    for component in components:
        if not isinstance(component, dict):
            raise ValueError('invalid support-component certificate')
        rows, columns, operations = (component.get(k) for k in ('rows', 'columns', 'operations'))
        if any(not isinstance(v, (list, tuple)) for v in (rows, columns, operations)):
            raise ValueError('component indices and operations must be sequences')
        if len(rows)*len(columns) > limits.max_component_entries:
            raise StarKoszulLimit('dense component cap reached')
        coordinate_count += len(rows)+len(columns)
        operation_count += len(operations)
    if coordinate_count > limits.max_coordinates or operation_count > limits.max_operations:
        raise StarKoszulLimit('certificate coordinate or operation cap reached')
    return components


def _replay(block, operations, limits):
    nr, nc = len(block), len(block[0]) if block else 0
    largest_bits = 1
    for operation in operations:
        if not isinstance(operation, (list, tuple)) or len(operation) not in (3, 4):
            raise ValueError('invalid elementary operation')
        kind, i, j = operation[:3]
        dimension = nr if kind in ('row_swap', 'row_add') else nc
        if any(type(v) is not int or not 0 <= v < dimension for v in (i, j)) or i == j:
            raise ValueError('elementary indices must be distinct component coordinates')
        if kind in ('row_swap', 'column_swap') and len(operation) == 3:
            if kind == 'row_swap':
                block[i], block[j] = block[j], block[i]
            else:
                for row in block:
                    row[i], row[j] = row[j], row[i]
        elif kind in ('row_add', 'column_add') and len(operation) == 4:
            factor = operation[3]
            if type(factor) is not int:
                raise ValueError('elementary addition needs an integer factor')
            if abs(factor).bit_length() > limits.max_integer_bits:
                raise StarKoszulLimit('integer factor bit cap reached')
            if kind == 'row_add':
                values = [a+factor*b for a, b in zip(block[i], block[j])]
                locations = ((i, column) for column in range(nc))
            else:
                values = [row[i]+factor*row[j] for row in block]
                locations = ((row, i) for row in range(nr))
            for (r, c), value in zip(locations, values):
                bits = abs(value).bit_length()
                if bits > limits.max_integer_bits:
                    raise StarKoszulLimit('integer replay bit cap reached')
                largest_bits = max(largest_bits, bits)
                block[r][c] = value
        else:
            raise ValueError('operation is not an elementary unimodular operation')
    return largest_bits


@dataclass(frozen=True)
class StarKoszulControl:
    unit_rank: int
    rank_one_factor: int
    component_count: int
    operation_count: int
    largest_integer_bits: int

    @property
    def rank_lower(self):
        return self.unit_rank

    @property
    def border_rank_lower(self):
        return (self.unit_rank+self.rank_one_factor-1)//self.rank_one_factor


def verify_star_integer_components(*, certificate=None, limits=DEFAULT_LIMITS):
    """Verify disjoint matrix blocks and every integer unit-minor operation."""
    _preflight(limits)
    components = _certificate(certificate, limits)
    matrix = star_integer_koszul(limits=limits)
    row_owner, column_owner = {}, {}
    for index, component in enumerate(components):
        for coordinates, owner in ((component['rows'], row_owner), (component['columns'], column_owner)):
            for coordinate in coordinates:
                if type(coordinate) is not int or not 0 <= coordinate < 1638 or coordinate in owner:
                    raise ValueError('component indices must be distinct disjoint matrix coordinates')
                owner[coordinate] = index
    for row, column in matrix:
        if row not in row_owner or column not in column_owner or row_owner[row] != column_owner[column]:
            raise ValueError('component partition does not cover the full signed matrix support')
    rank = operations = 0
    largest_bits = 1
    for component in components:
        rows, columns = component['rows'], component['columns']
        r = component.get('unit_rank')
        if type(r) is not int or not 0 <= r <= min(len(rows), len(columns)):
            raise ValueError('invalid component unit rank')
        block = [[matrix.get((i, j), 0) for j in columns] for i in rows]
        largest_bits = max(largest_bits, _replay(block, component['operations'], limits))
        determinant = 1
        for i in range(r):
            if block[i][i] not in (-1, 1) or any(block[i][j] for j in range(i)):
                raise ValueError('component minor is not triangular with unit diagonal')
            determinant *= block[i][i]
        if component.get('determinant') not in (-1, 1) or determinant != component['determinant']:
            raise ValueError('component minor has the wrong unit determinant')
        rank += r
        operations += len(component['operations'])
    if rank != 1054:
        raise ValueError('component unit ranks do not total the asserted universal lower bound')
    return StarKoszulControl(rank, 70, len(components), operations, largest_bits)


verify_star_integer_koszul = verify_star_integer_components
