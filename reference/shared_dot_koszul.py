"""A small integral minor certifying an all-characteristic C3 obstruction.

Only the binary first map and original minor indices are stored. Replay
rebuilds the target and signed exterior matrix over the integers, then proves
the selected minor unimodular without modular sampling or external packages.
The arbitrary-catalyst transfer is the argument in research/shared-dot-catalysts.md.
"""

from dataclasses import dataclass
from itertools import combinations
import json
from pathlib import Path

from .fields import FiniteField


class SharedDotKoszulLimit(ValueError):
    """A declared cap prevented completion of the certificate replay."""


@dataclass(frozen=True)
class SharedDotKoszulLimits:
    max_matrix_entries: int = 100_000
    max_integer_bits: int = 4096
    max_operations: int = 4096

    def __post_init__(self):
        if any(type(v) is not int or v < 0 for v in (self.max_matrix_entries, self.max_integer_bits, self.max_operations)):
            raise ValueError("Koszul replay caps must be nonnegative integers")


DEFAULT_LIMITS = SharedDotKoszulLimits()


def _preflight(limits, matrix_size=170):
    if not isinstance(limits, SharedDotKoszulLimits):
        raise ValueError("supply explicit integral Koszul replay limits")
    if matrix_size**2 > limits.max_matrix_entries or limits.max_integer_bits < 1:
        raise SharedDotKoszulLimit("integral exterior matrix or integer cap reached")


def _certificate(data=None):
    if data is None:
        data = json.loads(Path(__file__).with_name('shared_dot_three_koszul.json').read_text())
    if not isinstance(data, dict) or data.get('format') != 'shared-dot-three-unimodular-minor-v1':
        raise ValueError("unsupported integral shared-dot certificate format")
    first = tuple(tuple(row) for row in data.get('first_map', ()))
    if len(first) != 5 or any(len(row) != 17 or any(type(v) is not int or v not in (0, 1) for v in row) for row in first):
        raise ValueError("the stored first map must be a binary 5 by 17 matrix")
    rows, columns = tuple(data.get('minor_rows', ())), tuple(data.get('minor_columns', ()))
    if (len(rows) != 153 or len(columns) != 153 or len(set(rows)) != 153 or len(set(columns)) != 153
            or any(type(i) is not int or not 0 <= i < 170 for i in rows+columns)):
        raise ValueError("the stored minor must select 153 distinct rows and columns")
    if (data.get('koszul_degree') != 2 or data.get('rank_one_koszul_factor') != 6
            or data.get('all_characteristic_rank_lower') != 153 or data.get('minor_determinant') != 1):
        raise ValueError("unsupported asserted degree, rank factor, or determinant")
    return first, rows, columns


def shared_dot_three_integer_koszul(*, certificate=None, limits=DEFAULT_LIMITS):
    """Reconstruct Unit + M2*C3 from support triples and its signed p2 matrix."""
    _preflight(limits)
    first, _, _ = _certificate(certificate)
    return _integer_koszul(first, 3)


def _integer_koszul(first, q):
    size, n = q+1, 1+4*(q+1)
    support = [(0, 0, 0)]
    shared = tuple(triple for i in range(1, q+1) for triple in ((0, i, i), (i, 0, i), (i, i, 0)))
    for i in range(2):
        for k in range(2):
            for j in range(2):
                support.extend((1+(2*i+k)*size+x, 1+(2*k+j)*size+y, 1+(2*i+j)*size+z) for x, y, z in shared)
    if len(set(support)) != 1+24*q:
        raise ValueError("integer target reconstruction has the wrong support")
    slices = [[[0]*n for _ in range(n)] for _ in range(5)]
    for x, y, z in support:
        for axis in range(5):
            slices[axis][y][z] += first[axis][x]
    domain = tuple(combinations(range(5), 2))
    codomain = tuple(combinations(range(5), 3))
    positions = {subset: i for i, subset in enumerate(codomain)}
    matrix = [[0]*(10*n) for _ in range(10*n)]
    for column, subset in enumerate(domain):
        for axis in range(5):
            if axis in subset:
                continue
            row = positions[tuple(sorted((*subset, axis)))]
            sign = -1 if sum(i < axis for i in subset) % 2 else 1
            for y in range(n):
                for z in range(n):
                    matrix[row*n+z][column*n+y] = sign*slices[axis][y][z]
    return tuple(tuple(row) for row in matrix)


@dataclass(frozen=True)
class SharedDotKoszulControl:
    minor_size: int
    minor_determinant: int
    rank_one_factor: int
    largest_integer_bits: int

    @property
    def border_rank_lower(self):
        return (self.minor_size+self.rank_one_factor-1)//self.rank_one_factor


def verify_shared_dot_three_integer_minor(*, certificate=None, limits=DEFAULT_LIMITS):
    """Ordered integral row elimination checks every pivot is +/-1 and det=1."""
    _preflight(limits)
    _, rows, columns = _certificate(certificate)
    matrix = shared_dot_three_integer_koszul(certificate=certificate, limits=limits)
    minor = [[matrix[i][j] for j in columns] for i in rows]
    determinant, largest_bits = 1, 1
    for stage in range(153):
        pivot = minor[stage][stage]
        if pivot not in (-1, 1):
            raise ValueError(f"integral minor lost its unit pivot at {stage}")
        determinant *= pivot
        for i in range(stage+1, 153):
            factor = minor[i][stage]*pivot
            if factor:
                for j in range(stage, 153):
                    value = minor[i][j]-factor*minor[stage][j]
                    bits = abs(value).bit_length()
                    if bits > limits.max_integer_bits:
                        raise SharedDotKoszulLimit("integral elimination bit cap reached")
                    largest_bits = max(largest_bits, bits)
                    minor[i][j] = value
    if determinant != 1:
        raise ValueError("the selected integral minor does not have determinant one")
    return SharedDotKoszulControl(153, determinant, 6, largest_bits)


def shared_dot_three_field_koszul(field, *, limits=DEFAULT_LIMITS):
    """Embed integer coefficients, including signs, into the represented field."""
    if not isinstance(field, FiniteField):
        raise ValueError("supply an explicitly represented finite field")
    return tuple(tuple(field.embed(value) for value in row) for row in
                 shared_dot_three_integer_koszul(limits=limits))


def _two_certificate(certificate, limits):
    if certificate is None:
        certificate = json.loads(Path(__file__).with_name('shared_dot_two_koszul.json').read_text())
    if not isinstance(certificate, dict) or certificate.get('format') != 'shared-dot-two-unimodular-transform-v1':
        raise ValueError("unsupported integral two-dot certificate format")
    first = tuple(tuple(row) for row in certificate.get('first_map', ()))
    if len(first) != 5 or any(len(row) != 13 or any(type(v) is not int or v not in (0, 1) for v in row) for row in first):
        raise ValueError("the stored C2 first map must be binary 5 by 13")
    operations = certificate.get('operations', ())
    if len(operations) > limits.max_operations:
        raise SharedDotKoszulLimit("unimodular operation cap reached")
    if (certificate.get('unit_rank') != 124 or certificate.get('rank_one_factor') != 6
            or certificate.get('determinant') != -1):
        raise ValueError("unsupported two-dot rank or determinant assertion")
    return first, operations


def shared_dot_two_integer_koszul(*, certificate=None, limits=DEFAULT_LIMITS):
    _preflight(limits, 130)
    first, _ = _two_certificate(certificate, limits)
    return _integer_koszul(first, 2)


def verify_shared_dot_two_integer_transform(*, certificate=None, limits=DEFAULT_LIMITS):
    """Replay only unimodular integer operations, then inspect a 124-unit block."""
    _preflight(limits, 130)
    first, operations = _two_certificate(certificate, limits)
    matrix = [list(row) for row in _integer_koszul(first, 2)]
    largest_bits = 1
    for operation in operations:
        if not isinstance(operation, (list, tuple)) or len(operation) not in (3, 4):
            raise ValueError("invalid integral elementary operation")
        kind, i, j = operation[:3]
        if any(type(v) is not int or not 0 <= v < 130 for v in (i, j)) or i == j:
            raise ValueError("elementary operation indices must be distinct matrix coordinates")
        if kind in ('row_swap', 'column_swap') and len(operation) == 3:
            if kind == 'row_swap':
                matrix[i], matrix[j] = matrix[j], matrix[i]
            else:
                for row in matrix:
                    row[i], row[j] = row[j], row[i]
        elif kind in ('row_add', 'column_add') and len(operation) == 4:
            factor = operation[3]
            if type(factor) is not int:
                raise ValueError("elementary addition needs an integer factor")
            if abs(factor).bit_length() > limits.max_integer_bits:
                raise SharedDotKoszulLimit("unimodular factor bit cap reached")
            if kind == 'row_add':
                values = [a+factor*b for a, b in zip(matrix[i], matrix[j])]
                locations = ((i, column) for column in range(130))
            else:
                values = [row[i]+factor*row[j] for row in matrix]
                locations = ((row, i) for row in range(130))
            for (r, c), value in zip(locations, values):
                bits = abs(value).bit_length()
                if bits > limits.max_integer_bits:
                    raise SharedDotKoszulLimit("unimodular replay integer bit cap reached")
                largest_bits = max(largest_bits, bits)
                matrix[r][c] = value
        else:
            raise ValueError("operation is not an elementary unimodular row or column operation")
    determinant = 1
    for i in range(124):
        if matrix[i][i] not in (-1, 1) or any(matrix[i][j] for j in range(i)):
            raise ValueError("transformed minor is not triangular with unit diagonal")
        determinant *= matrix[i][i]
    if determinant != -1:
        raise ValueError("transformed minor has the wrong unit determinant")
    return SharedDotKoszulControl(124, determinant, 6, largest_bits)
