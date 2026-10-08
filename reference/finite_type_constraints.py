"""Actual exact-type and finite Fourier restrictions for witnessed LP rows.

These retain finite interpolation/Fourier overheads. They do not assert an
entropy limit or multiply a merely additive state like a character.
"""

from dataclasses import dataclass
from itertools import product
from math import comb, prod

from .constructions import constant_weights, nonzero_nodes
from .constraint_generators import WitnessedOrderExport
from .maps import LinearMap
from .separation import FiniteSeparation
from .tensors import FiniteTensor, _tensor_family, shared_first_tensor
from .witnessed_constraints import DEFAULT_LIMITS, StateAssemblyLimit, WitnessedOrder


def _names(source_key, target_key):
    if any(not isinstance(key, str) or not key for key in (source_key, target_key)) or source_key == target_key:
        raise ValueError("source and target IDs must be distinct nonempty strings")


def _power_shape(shape, exponent, limits):
    # Avoid constructing huge integers, let alone dense powers.
    lower_bits = prod(shape).bit_length()-1
    if exponent*lower_bits >= limits.max_tensor_entries.bit_length():
        raise StateAssemblyLimit("powered tensor cap reached; construction incomplete")
    axis_bits = max(shape).bit_length()-1
    if exponent*axis_bits >= max(limits.max_map_entries, limits.max_tensor_entries).bit_length():
        raise StateAssemblyLimit("powered axis cap reached; construction incomplete")
    result = tuple(size**exponent for size in shape)
    limits.tensor(result)
    return result


@dataclass(frozen=True)
class ExactTypeExport:
    restriction: WitnessedOrderExport
    counts: tuple[int, ...]
    words: tuple[tuple[int, ...], ...]
    branches: tuple[FiniteTensor, ...]

    def verify(self, limits=DEFAULT_LIMITS):
        return self.restriction.verify(limits)


def export_exact_type_order(branches, counts, source_key, target_key, *, limits=DEFAULT_LIMITS):
    """Restrict shared(B)^n to the shared family of all words with these counts.

    All three branch axes are uniform. Word positions retain their tensor
    coordinate order; no unverified identification of permutations occurs.
    """
    branches = _tensor_family(branches, uniform=True)
    counts = tuple(counts)
    _names(source_key, target_key)
    if len(counts) != len(branches) or any(type(c) is not int or c < 0 for c in counts):
        raise ValueError("one nonnegative integer count is required per branch")
    n, M, placed = sum(counts), 1, 0
    if n > limits.max_blocks:
        raise StateAssemblyLimit("word length cap reached; construction incomplete")
    for count in counts:
        M *= comb(placed+count, count)
        placed += count
        if M > limits.max_blocks:
            raise StateAssemblyLimit("exact-word cap reached; construction incomplete")
    # Empty tensor axes make dense coefficient/map caps vacuous. Word labels
    # still occupy n*M positions and require that many product steps.
    if n*M > limits.max_constraint_entries:
        raise StateAssemblyLimit("exact-word position cap reached; construction incomplete")
    f, s = branches[0].field, len(branches)
    x, y, z = branches[0].shape
    source_seed_shape = (x, s*y, s*z)
    source_shape = _power_shape(source_seed_shape, n, limits)
    branch_shape = _power_shape((x, y, z), n, limits)
    target_shape = (branch_shape[0], M*branch_shape[1], M*branch_shape[2])
    limits.tensor(source_seed_shape)
    limits.tensor(target_shape)
    limits.maps(source_shape, target_shape)

    words, stack = [], [((), counts)]
    active = tuple(label for label, count in enumerate(counts) if count)
    while stack:
        prefix, left = stack.pop()
        if len(prefix) == n:
            words.append(prefix)
            continue
        for label in reversed(active):
            count = left[label]
            if count:
                next_left = list(left)
                next_left[label] -= 1
                stack.append((prefix+(label,), tuple(next_left)))
    if len(words) != M:
        raise ValueError("exact-word enumeration disagrees with the integer count")
    source = shared_first_tensor(branches).tensor_power(n)
    products = []
    for word in words:
        tensor = FiniteTensor.unit(f)
        for label in word:
            tensor = tensor.tensor_product(branches[label])
        products.append(tensor)
    target = shared_first_tensor(products)
    selections = []
    for axis, size in enumerate((y, z), start=1):
        indices = []
        for word in words:
            for digits in product(range(size), repeat=n):
                old = 0
                for label, digit in zip(word, digits):
                    old = old*(s*size)+label*size+digit
                indices.append(old)
        selections.append(LinearMap.selection(f, source_shape[axis], indices))
    maps = (LinearMap.identity(f, source_shape[0]), *selections)
    row = WitnessedOrder((source_key,), (target_key,), maps, "exact-type-power")
    exported = WitnessedOrderExport(row, ((source_key, source), (target_key, target)), exponent=n)
    exported.require_valid(limits)
    return ExactTypeExport(exported, counts, tuple(words), tuple(products))


def export_finite_separation_order(separation, source_key, target_key, *, nodes=None,
                                   limits=DEFAULT_LIMITS):
    """Export L*J shared source copies -> sum_h B_h tensor dot_M.

    L is the supplied Fourier period; J interpolation nodes are retained.
    Actual unshifted square maps may use negative powers at nonzero nodes.
    """
    if not isinstance(separation, FiniteSeparation):
        raise ValueError("supply an actual finite Fourier separation")
    _names(source_key, target_key)
    count = separation.normalized_degree_bound+1
    if count > limits.max_blocks or separation.period > limits.max_blocks:
        raise StateAssemblyLimit("finite separation block cap reached; construction incomplete")
    f = separation.field
    nodes = nonzero_nodes(f, count) if nodes is None else tuple(nodes)
    if len(nodes) < count:
        raise ValueError("too few distinct nonzero square-filter interpolation nodes")
    copies = separation.period*len(nodes)
    if copies > limits.max_blocks:
        raise StateAssemblyLimit("finite separation block cap reached; construction incomplete")
    source_shape = tuple(copies*size for size in (separation.branch_shape[0],
                         separation.count*separation.branch_shape[1],
                         separation.count*separation.branch_shape[2]))
    limits.tensor(source_shape)
    limits.tensor(separation.shape)
    limits.maps(source_shape, separation.shape)
    weights = constant_weights(f, nodes)
    separation.require_valid()
    weighted = tuple(separation.weighted_maps(node) for node in nodes)
    maps = []
    for axis in range(3):
        rows = tuple(tuple(f.mul(weights[h] if axis == 2 else 1, value)
                           for h, family in enumerate(weighted) for value in family[axis].rows[i])
                     for i in range(separation.shape[axis]))
        mapping = LinearMap(f, source_shape[axis], rows)
        if axis == 0:
            x, M = separation.branch_shape[0], separation.count
            mapping = LinearMap.selection(f, separation.shape[0],
                       (i*M+h for h in range(M) for i in range(x))).compose(mapping)
        maps.append(mapping)
    source, target = separation.source, separation.direct_sum_target
    row = WitnessedOrder((source_key,)*copies, (target_key,), tuple(maps), "finite-Fourier-square-separation")
    exported = WitnessedOrderExport(row, ((source_key, source), (target_key, target)), nodes, 1)
    exported.require_valid(limits)
    return exported
