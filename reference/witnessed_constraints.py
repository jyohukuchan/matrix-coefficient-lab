"""Exact finite state certificates assembled into actual catalytic maps.

Ordinary constraints retain source/target tensor blocks AND local maps.
Detector constraints are hypothetical state inequalities, not restriction
witnesses. An integer negative-unit balance joins them without cancelling
any tensor summand. This implements the finite construction described in
research/witnessed-lp.md; it does not find a useful 9/4 catalyst by itself.
"""

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from functools import cached_property
from itertools import combinations
from math import gcd, lcm, prod

from .catalysts import CatalyticCertificate
from .maps import LinearMap, validate_local_maps
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, matrix_multiplication_tensor


class StateAssemblyLimit(ValueError):
    """Incomplete construction because a declared dense-resource cap was hit."""


@dataclass(frozen=True)
class StateAssemblyLimits:
    max_blocks: int = 4096
    max_tensor_entries: int = 1_000_000
    max_map_entries: int = 1_000_000
    max_constraint_entries: int = 1_000_000

    def __post_init__(self):
        for value in (self.max_blocks, self.max_tensor_entries, self.max_map_entries,
                      self.max_constraint_entries):
            if type(value) is not int or value < 0:
                raise ValueError("assembly caps must be nonnegative integers")

    def tensor(self, shape):
        if prod(shape) > self.max_tensor_entries:
            raise StateAssemblyLimit("dense tensor cap reached; construction incomplete")

    def maps(self, source_shape, target_shape):
        if sum(a*b for a, b in zip(source_shape, target_shape)) > self.max_map_entries:
            raise StateAssemblyLimit("dense map cap reached; construction incomplete")

    def constraints(self, rows, coordinates):
        if rows*coordinates > self.max_constraint_entries:
            raise StateAssemblyLimit("dense constraint cap reached; search incomplete")


DEFAULT_LIMITS = StateAssemblyLimits()


@dataclass(frozen=True)
class WitnessedOrder:
    """vector=count(positive)-count(negative), maps DS(positive)->DS(negative).

    This covers lower, certified rank upper bound, order, and both additive
    relation directions. The label is descriptive; it carries no authority.
    """

    positive: tuple[str, ...]
    negative: tuple[str, ...]
    maps: tuple[LinearMap, LinearMap, LinearMap]
    label: str = "order"

    def __post_init__(self):
        object.__setattr__(self, "positive", tuple(self.positive))
        object.__setattr__(self, "negative", tuple(self.negative))
        object.__setattr__(self, "maps", tuple(self.maps))


@dataclass(frozen=True)
class DetectorConstraint:
    """vector=e_product-k*e_test, requiring product == M_d tensor test."""

    test: str
    product: str


def _sum_shape(tensors):
    return tuple(sum(t.shape[axis] for t in tensors) for axis in range(3))


def _direct_sum(field, tensors, limits):
    tensors = tuple(tensors)
    if len(tensors) > limits.max_blocks:
        raise StateAssemblyLimit("block cap reached; construction incomplete")
    limits.tensor(_sum_shape(tensors))
    return FiniteTensor.direct_sum(tensors) if tensors else FiniteTensor.zero(field, (0, 0, 0))


@dataclass(frozen=True)
class WitnessedStateSystem:
    """Fixed exact tensor representatives and a finite witnessed row family."""

    d: int
    k: int
    inventory: tuple[tuple[str, FiniteTensor], ...]
    constraints: tuple[WitnessedOrder | DetectorConstraint, ...]
    unit_key: str = "unit"

    def __post_init__(self):
        if type(self.d) is not int or self.d < 2 or type(self.k) is not int or self.k < 1:
            raise ValueError("d>=2 and k>=1 must be integers")
        inventory = tuple(tuple(item) for item in self.inventory)
        if not inventory or any(len(item) != 2 or not isinstance(item[0], str) or
                                not isinstance(item[1], FiniteTensor) for item in inventory):
            raise ValueError("inventory must contain named exact tensors")
        if len({item[0] for item in inventory}) != len(inventory):
            raise ValueError("tensor IDs must be unique")
        field = inventory[0][1].field
        if any(t.field != field for _, t in inventory):
            raise ValueError("registry tensors must use one fixed field representation")
        object.__setattr__(self, "inventory", inventory)
        object.__setattr__(self, "constraints", tuple(self.constraints))
        if self.registry.get(self.unit_key) != FiniteTensor.unit(field):
            raise ValueError("unit ID must denote the exact scalar unit")
        for row in self.constraints:
            if isinstance(row, WitnessedOrder):
                keys = row.positive + row.negative
            elif isinstance(row, DetectorConstraint):
                keys = (row.test, row.product)
            else:
                raise ValueError("unknown state constraint record")
            if any(key not in self.registry for key in keys):
                raise ValueError("constraint references an unregistered tensor ID")

    @property
    def field(self):
        return self.inventory[0][1].field

    @property
    def registry(self):
        # Return a fresh dict so callers cannot mutate cached representatives.
        return dict(self.inventory)

    @property
    def keys(self):
        return tuple(key for key, _ in self.inventory)

    @cached_property
    def sparse_vectors(self):
        vectors = []
        for row in self.constraints:
            if isinstance(row, WitnessedOrder):
                counts = Counter(row.positive)
                counts.subtract(row.negative)
            else:
                counts = Counter({row.product: 1})
                counts[row.test] -= self.k
            vectors.append(tuple((key, count) for key, count in counts.items() if count))
        return tuple(vectors)

    def dense_vectors(self, limits=DEFAULT_LIMITS, keys=None):
        keys = self.keys if keys is None else tuple(keys)
        limits.constraints(len(self.constraints), len(keys))
        result = []
        for row in self.sparse_vectors:
            counts = dict(row)
            result.append(tuple(counts.get(key, 0) for key in keys))
        return tuple(result)

    @property
    def vectors(self):
        """Dense inspection under default caps; certificate balance stays sparse."""
        return self.dense_vectors()

    def require_valid(self, limits=DEFAULT_LIMITS):
        """Validate EVERY ordinary coefficient map and detector identity."""
        registry = self.registry
        limits.tensor((self.d**2,)*3)
        matrix = matrix_multiplication_tensor(self.field, self.d, self.d, self.d)
        for row in self.constraints:
            if isinstance(row, WitnessedOrder):
                source = _direct_sum(self.field, (registry[key] for key in row.positive), limits)
                target = _direct_sum(self.field, (registry[key] for key in row.negative), limits)
                maps = validate_local_maps(self.field, source.shape, row.maps)
                if tuple(m.output_size for m in maps) != target.shape:
                    raise ValueError("ordinary constraint map target dimensions are wrong")
                limits.maps(source.shape, target.shape)
                if source.restrict(*maps) != target:
                    raise ValueError("ordinary constraint lacks an exact coefficient witness")
            else:
                test = registry[row.test]
                limits.tensor(tuple(a*b for a, b in zip(matrix.shape, test.shape)))
                if registry[row.product] != matrix.tensor_product(test):
                    raise ValueError("detector product must be the exact registered M_d tensor test")


@dataclass(frozen=True)
class IntegerStateCertificate:
    """sum weights[i]*vectors[i] == -m*e_unit over the INTEGERS."""

    system: WitnessedStateSystem
    weights: tuple[int, ...]
    m: int

    def __post_init__(self):
        if not isinstance(self.system, WitnessedStateSystem):
            raise ValueError("certificate requires a witnessed state system")
        weights = tuple(self.weights)
        if len(weights) != len(self.system.constraints) or any(type(n) is not int or n < 0 for n in weights):
            raise ValueError("certificate weights must be nonnegative integers, one per constraint")
        if type(self.m) is not int or self.m < 1:
            raise ValueError("scalar gain must be a positive integer")
        object.__setattr__(self, "weights", weights)

    def verify_balance(self):
        counts = Counter()
        for n, vector in zip(self.weights, self.system.sparse_vectors):
            if n:
                for key, value in vector:
                    counts[key] += n*value
        counts[self.system.unit_key] += self.m
        return not any(counts.values())

    def require_valid(self, limits=DEFAULT_LIMITS):
        if not self.verify_balance():
            raise ValueError("negative-unit balance fails over the integers")
        self.system.require_valid(limits)


def clear_rational_dual(system, coefficients):
    """Check an exact rational dual, clear denominators, and reduce its gcd.

    A solver's floating-point status is never accepted. Ordinary map witnesses
    are checked by certificate.require_valid() or assemble_catalyst().
    """
    coefficients = tuple(coefficients)
    if len(coefficients) != len(system.constraints) or any(
            type(q) not in (int, Fraction) or q < 0 for q in coefficients):
        raise ValueError("dual coefficients must be nonnegative exact rationals")
    coefficients = tuple(Fraction(q) for q in coefficients)
    m = lcm(*(q.denominator for q in coefficients)) if coefficients else 1
    weights = tuple(int(m*q) for q in coefficients)
    divisor = gcd(m, *weights)
    certificate = IntegerStateCertificate(system, tuple(n//divisor for n in weights), m//divisor)
    if not certificate.verify_balance():
        raise ValueError("rational dual does not equal negative unit exactly")
    return certificate


@dataclass(frozen=True)
class DualSearchResult:
    status: str
    checked_subsets: int
    coefficients: tuple[Fraction, ...] | None
    reason: str


def _independent_solution(columns, target):
    """Unique solution to an independent-column exact rational system, or None."""
    width = len(columns)
    rows = [[Fraction(column[j]) for column in columns]+[Fraction(value)]
            for j, value in enumerate(target)]
    pivot_row = 0
    for column in range(width):
        pivot = next((i for i in range(pivot_row, len(rows)) if rows[i][column]), None)
        if pivot is None:
            return None
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [value/scale for value in rows[pivot_row]]
        for i in range(len(rows)):
            if i != pivot_row and rows[i][column]:
                scale = rows[i][column]
                rows[i] = [a-scale*b for a, b in zip(rows[i], rows[pivot_row])]
        pivot_row += 1
    if any(not any(row[:width]) and row[-1] for row in rows):
        return None
    return tuple(rows[i][-1] for i in range(width))


def find_nonnegative_dual(system, *, max_subsets=10_000, max_support=None,
                         limits=DEFAULT_LIMITS):
    """Small exact finite-cone search by independent supports.

    This is complete for the supplied FINITE row family when all supports up
    to min(number of rows, number of coordinates) are examined. A minimal
    nonnegative representation has independent support. It is not a search
    over all tensors, extensions, or restriction witnesses.
    """
    if type(max_subsets) is not int or max_subsets < 0:
        raise ValueError("subset cap must be a nonnegative integer")
    if not isinstance(system, WitnessedStateSystem):
        raise ValueError("dual search requires a witnessed state system")
    if max_support is not None and (type(max_support) is not int or max_support < 0):
        raise ValueError("support cap must be a nonnegative integer")
    if system.constraints and (max_subsets == 0 or max_support == 0):
        return DualSearchResult("resource_cap", 0, None, "search cap reached; finite-family search incomplete")
    limits.constraints(len(system.constraints), 1)
    used = {system.unit_key}
    used.update(key for row in system.sparse_vectors for key, _ in row)
    keys = tuple(key for key in system.keys if key in used)
    full_support = min(len(system.constraints), len(keys))
    if max_support is None:
        max_support = full_support
    support = min(max_support, full_support)
    if full_support and (max_subsets == 0 or support == 0):
        return DualSearchResult("resource_cap", 0, None, "search cap reached; finite-family search incomplete")
    vectors = system.dense_vectors(limits, keys)
    system.require_valid(limits)
    target = tuple(-int(key == system.unit_key) for key in keys)
    checked = 0
    for size in range(1, support+1):
        for selected in combinations(range(len(vectors)), size):
            if checked >= max_subsets:
                return DualSearchResult("resource_cap", checked, None, "subset cap reached; finite-family search incomplete")
            checked += 1
            solution = _independent_solution(tuple(vectors[i] for i in selected), target)
            if solution is not None and all(q >= 0 for q in solution):
                coefficients = [Fraction(0)]*len(vectors)
                for i, q in zip(selected, solution):
                    coefficients[i] = q
                clear_rational_dual(system, coefficients)
                return DualSearchResult("found", checked, tuple(coefficients), "exact rational negative-unit certificate")
    if support < full_support:
        return DualSearchResult("resource_cap", checked, None, "support cap reached; finite-family search incomplete")
    return DualSearchResult("finite_exhausted", checked, None, "no dual in the supplied finite row family; no global impossibility claim")


def _block_diagonal(field, groups, limits):
    """Each group is a checked tuple of three local maps."""
    shapes = tuple(tuple(m.input_size for m in group) for group in groups)
    outputs = tuple(tuple(m.output_size for m in group) for group in groups)
    source = tuple(sum(shape[axis] for shape in shapes) for axis in range(3))
    target = tuple(sum(shape[axis] for shape in outputs) for axis in range(3))
    limits.maps(source, target)
    result = []
    for axis in range(3):
        rows, offset = [], 0
        for group in groups:
            mapping = group[axis]
            rows.extend((0,)*offset+row+(0,)*(source[axis]-offset-mapping.input_size) for row in mapping.rows)
            offset += mapping.input_size
        result.append(LinearMap(field, source[axis], tuple(rows)))
    return tuple(result)


def _occurrence_permutation(field, registry, source_keys, target_keys, limits):
    if Counter(source_keys) != Counter(target_keys):
        raise ValueError("tensor-ID occurrence multisets do not balance")
    available = defaultdict(deque)
    offsets = [0, 0, 0]
    for key in source_keys:
        available[key].append(tuple(offsets))
        offsets = [a+b for a, b in zip(offsets, registry[key].shape)]
    indices = [[], [], []]
    for key in target_keys:
        starts = available[key].popleft()
        for axis in range(3):
            indices[axis].extend(range(starts[axis], starts[axis]+registry[key].shape[axis]))
    limits.maps(tuple(offsets), tuple(len(selected) for selected in indices))
    return tuple(LinearMap.selection(field, size, selected) for size, selected in zip(offsets, indices))


def assemble_catalyst(certificate, *, limits=DEFAULT_LIMITS):
    """Turn a fully witnessed integer dual into three ACTUAL catalytic maps."""
    certificate.require_valid(limits)
    system, registry = certificate.system, certificate.system.registry
    ordinary = tuple((row, n) for row, n in zip(system.constraints, certificate.weights)
                     if isinstance(row, WitnessedOrder) and n)
    detectors = tuple((row, n) for row, n in zip(system.constraints, certificate.weights)
                     if isinstance(row, DetectorConstraint) and n)
    d_count = sum(n*len(row.positive) for row, n in ordinary)
    n_count = sum(n*len(row.negative) for row, n in ordinary)
    s_count = sum(n for _, n in detectors)
    if max(d_count+system.k*s_count, n_count+system.k*s_count,
           d_count+certificate.m+s_count, sum(n for _, n in ordinary)) > limits.max_blocks:
        raise StateAssemblyLimit("expanded block cap reached; construction incomplete")
    D_keys = tuple(key for row, n in ordinary for _ in range(n) for key in row.positive)
    N_keys = tuple(key for row, n in ordinary for _ in range(n) for key in row.negative)
    S_keys = tuple(row.test for row, n in detectors for _ in range(n))
    product_keys = tuple(row.product for row, n in detectors for _ in range(n))
    D = _direct_sum(system.field, (registry[key] for key in D_keys), limits)
    S = _direct_sum(system.field, (registry[key] for key in S_keys), limits)
    if not any(S.coefficients):
        raise ValueError("balanced certificate must produce a nonzero auxiliary tensor")
    source_shape = tuple(a+system.k*b for a, b in zip(D.shape, S.shape))
    middle_shape = _sum_shape(tuple(registry[key] for key in N_keys+S_keys*system.k))
    target_shape = tuple(a+certificate.m+system.d**2*b for a, b in zip(D.shape, S.shape))
    for shape in (source_shape, middle_shape, target_shape):
        limits.tensor(shape)
    limits.maps(source_shape, middle_shape)
    limits.maps(middle_shape, target_shape)
    limits.maps(source_shape, target_shape)
    groups = tuple(row.maps for row, n in ordinary for _ in range(n))
    identity = tuple(LinearMap.identity(system.field, system.k*size) for size in S.shape)
    maps = _block_diagonal(system.field, groups+(identity,), limits)
    target_keys = D_keys+(system.unit_key,)*certificate.m+product_keys
    reordered = _occurrence_permutation(system.field, registry, N_keys+S_keys*system.k, target_keys, limits)
    maps = tuple(p.compose(previous) for p, previous in zip(reordered, maps))
    # Pack the product direct-sum blocks into M_d tensor S. Axis offsets
    # are independent; unequal block dimensions are deliberately supported.
    packed = []
    for axis, size in enumerate(target_shape):
        prefix = D.shape[axis]+certificate.m
        selected = list(range(prefix))
        blocks = []
        product_offset = prefix
        for key in S_keys:
            width = registry[key].shape[axis]
            blocks.append((product_offset, width))
            product_offset += system.d**2*width
        for i in range(system.d**2):
            for offset, width in blocks:
                selected.extend(range(offset+i*width, offset+(i+1)*width))
        if len(selected) != size:
            raise ValueError("product distribution map has the wrong axis size")
        packed.append(LinearMap.selection(system.field, size, selected))
    maps = tuple(p.compose(previous) for p, previous in zip(packed, maps))
    result = CatalyticCertificate(system.d, system.k, certificate.m, D, S, maps)
    result.require_valid()
    return result


def rank_upper_constraint(unit_key, target_key, scheme):
    """Witness an upper row using supplied exact coefficient families."""
    if not isinstance(scheme, TensorScheme):
        raise ValueError("upper relation requires a tensor scheme")
    scheme.require_exact()
    maps = tuple(LinearMap(scheme.field, scheme.terms,
                           tuple(tuple(vector[i] for vector in family) for i in range(size)))
                 for family, size in zip((scheme.a, scheme.b, scheme.c), scheme.shape))
    return WitnessedOrder((unit_key,)*scheme.terms, (target_key,), maps, "upper")
