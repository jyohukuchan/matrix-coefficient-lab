"""Finite ordinary restriction rows from proof sectors and tensor contexts.

Sector export retains the interpolation-copy overhead; it does not substitute
the asymptotic character inequality for a finite restriction. No rank
decomposition, detector assertion, or near-9/4 catalyst is inferred here.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from math import prod

from .constructions import constant_weights, nonzero_nodes
from .maps import LinearMap, validate_local_maps
from .sectors import ThreeSectorConstruction
from .tensors import FiniteTensor
from .witnessed_constraints import (DEFAULT_LIMITS, StateAssemblyLimit,
                                    WitnessedOrder)


def _registry(inventory):
    items = tuple(inventory.items()) if isinstance(inventory, Mapping) else tuple(inventory)
    if not items or any(len(item) != 2 or not isinstance(item[0], str) or not item[0]
                        or not isinstance(item[1], FiniteTensor) for item in items):
        raise ValueError("inventory must name exact finite tensors")
    if len({key for key, _ in items}) != len(items):
        raise ValueError("inventory tensor IDs must be unique")
    field = items[0][1].field
    if any(tensor.field != field for _, tensor in items):
        raise ValueError("inventory must use one fixed field")
    return dict(items)


def _shape(registry, keys):
    return tuple(sum(registry[key].shape[axis] for key in keys) for axis in range(3))


def _direct(registry, keys, field):
    return FiniteTensor.direct_sum(registry[key] for key in keys) if keys else FiniteTensor.zero(field, (0, 0, 0))


def _check_row(row, registry, limits):
    if not isinstance(row, WitnessedOrder):
        raise ValueError("only ordinary witnessed restrictions may be lifted")
    if any(key not in registry for key in row.positive+row.negative):
        raise ValueError("ordinary row references an unregistered tensor")
    if max(len(row.positive), len(row.negative)) > limits.max_blocks:
        raise StateAssemblyLimit("block cap reached; construction incomplete")
    field = next(iter(registry.values())).field
    source_shape, target_shape = _shape(registry, row.positive), _shape(registry, row.negative)
    limits.tensor(source_shape)
    limits.tensor(target_shape)
    limits.maps(source_shape, target_shape)
    maps = validate_local_maps(field, source_shape, row.maps)
    if tuple(mapping.output_size for mapping in maps) != target_shape:
        raise ValueError("ordinary row output axes have the wrong sizes")
    source, target = _direct(registry, row.positive, field), _direct(registry, row.negative, field)
    if source.restrict(*maps) != target:
        raise ValueError("generated ordinary row fails full coefficient equality")


@dataclass(frozen=True)
class WitnessedOrderExport:
    """One actual restriction row and its exact named tensor representatives."""

    row: WitnessedOrder
    inventory: tuple[tuple[str, FiniteTensor], ...]
    nodes: tuple[int, ...] = ()
    exponent: int | None = None
    key_pairs: tuple[tuple[str, str], ...] = ()

    def __post_init__(self):
        if not isinstance(self.row, WitnessedOrder):
            raise ValueError("export must contain an ordinary witnessed row")
        if self.exponent is not None and (type(self.exponent) is not int or self.exponent < 0):
            raise ValueError("export exponent must be nonnegative or omitted")
        object.__setattr__(self, "inventory", tuple(tuple(item) for item in self.inventory))
        object.__setattr__(self, "nodes", tuple(self.nodes))
        object.__setattr__(self, "key_pairs", tuple(tuple(pair) for pair in self.key_pairs))

    @property
    def registry(self):
        return dict(self.inventory)

    def require_valid(self, limits=DEFAULT_LIMITS):
        _check_row(self.row, _registry(self.inventory), limits)

    def verify(self, limits=DEFAULT_LIMITS):
        try:
            self.require_valid(limits)
            return True
        except StateAssemblyLimit:
            raise
        except ValueError:
            return False


def export_sector_order(construction, source_key, retained_key, *, exponent=1,
                        nodes=None, limits=DEFAULT_LIMITS):
    """Export (e+1)*source^e -> retained^e using the actual proof diagonals.

    The finite implementation supports positive a,h, including e=0 (unit).
    Even if the actual erased tensor vanishes, e+1 distinct nonzero nodes
    are required. All powers share the same parameter. Empty-axis sector
    boundaries are deliberately excluded rather than silently changing this
    degree contract.
    """
    if not isinstance(construction, ThreeSectorConstruction) or construction.a < 1 or construction.h < 1:
        raise ValueError("sector row export requires positive a,h")
    if type(exponent) is not int or exponent < 0:
        raise ValueError("sector exponent must be a nonnegative integer")
    if any(not isinstance(key, str) or not key for key in (source_key, retained_key)) or source_key == retained_key:
        raise ValueError("source and retained IDs must be distinct nonempty strings")
    # Formal validation materializes the seed's local polynomial matrices.
    limits.tensor(construction.shape)
    limits.maps(construction.shape, construction.shape)
    # Cheap lower bound before powers or exponent-length digit loops.
    seed_bits = prod(construction.shape).bit_length()-1
    if exponent*seed_bits >= limits.max_tensor_entries.bit_length():
        raise StateAssemblyLimit("powered tensor cap reached; construction incomplete")
    count = exponent+1
    if count > limits.max_blocks:
        raise StateAssemblyLimit("interpolation block cap reached; construction incomplete")
    f = construction.field
    nodes = nonzero_nodes(f, count) if nodes is None else tuple(nodes)
    if len(nodes) > limits.max_blocks:
        raise StateAssemblyLimit("interpolation block cap reached; construction incomplete")
    if len(nodes) < count:
        raise ValueError("sector power requires at least exponent+1 nonzero nodes")
    powered_shape = tuple(size**exponent for size in construction.shape)
    source_shape = tuple(len(nodes)*size for size in powered_shape)
    limits.tensor(source_shape)
    limits.tensor(powered_shape)
    limits.maps(source_shape, powered_shape)
    weights = constant_weights(f, nodes)  # only after copy/map allocation preflight
    construction.require_valid()  # formal P=tR+t^2E, not sampled evaluations
    source = construction.source.tensor_power(exponent)
    retained = construction.retained.tensor_power(exponent)
    maps = []
    for axis, (dimension, seed_size, diagonal) in enumerate(zip(powered_shape, construction.shape,
                                                             construction.diagonal_maps)):
        rows = []
        for coordinate in range(dimension):
            values = []
            for node, weight in zip(nodes, weights):
                digits, value = coordinate, 1
                for _ in range(exponent):
                    digits, digit = divmod(digits, seed_size)
                    value = f.mul(value, diagonal[digit].evaluate(node))
                if axis == 2:
                    value = f.mul(value, f.div(weight, f.pow(node, exponent)))
                values.append(value)
            rows.append(tuple(value if j == coordinate else 0
                              for value in values for j in range(dimension)))
        maps.append(LinearMap(f, len(nodes)*dimension, rows))
    row = WitnessedOrder((source_key,)*len(nodes), (retained_key,), tuple(maps), "proof-sector-interpolation")
    result = WitnessedOrderExport(row, ((source_key, source), (retained_key, retained)), nodes, exponent)
    result.require_valid(limits)
    return result


def lift_witnessed_order(row, inventory, context_tensor, key_map, *, limits=DEFAULT_LIMITS):
    """Tensor an actual ordinary row with context and distribute its blocks.

    key_map assigns a fresh unique ID to every used old ID. The returned
    inventory contains these exact tensor products, preserving repeated
    occurrences in positive and negative. Kronecker identity maps use
    old-coordinate first, context-coordinate second on each leg.
    """
    registry = _registry(inventory)
    if not isinstance(context_tensor, FiniteTensor) or context_tensor.field != next(iter(registry.values())).field:
        raise ValueError("context must be a finite tensor over the same field")
    if not isinstance(row, WitnessedOrder):
        raise ValueError("context lifting requires an ordinary restriction row")
    keys = tuple(dict.fromkeys(row.positive+row.negative))
    if not keys:
        raise ValueError("context export requires at least one registered row atom")
    if not isinstance(key_map, Mapping) or any(key not in key_map for key in keys):
        raise ValueError("provide a new tensor ID for each used row atom")
    names = tuple(key_map[key] for key in keys)
    if any(not isinstance(name, str) or not name or name in registry for name in names) or len(set(names)) != len(names):
        raise ValueError("lifted tensor IDs must be fresh unique nonempty strings")
    if any(key not in registry for key in keys):
        raise ValueError("ordinary row references an unregistered tensor")
    if max(len(row.positive), len(row.negative)) > limits.max_blocks:
        raise StateAssemblyLimit("block cap reached; construction incomplete")
    old_source, old_target = _shape(registry, row.positive), _shape(registry, row.negative)
    source_shape = tuple(a*b for a, b in zip(old_source, context_tensor.shape))
    target_shape = tuple(a*b for a, b in zip(old_target, context_tensor.shape))
    limits.tensor(source_shape)
    limits.tensor(target_shape)
    limits.maps(source_shape, target_shape)
    _check_row(row, registry, limits)
    f, maps = context_tensor.field, []
    for axis, mapping in enumerate(row.maps):
        width = context_tensor.shape[axis]
        rows = tuple(tuple(value if old_context == context else 0
                           for value in row_values for old_context in range(width))
                     for row_values in mapping.rows for context in range(width))
        maps.append(LinearMap(f, mapping.input_size*width, rows))
    items = tuple((key_map[key], registry[key].tensor_product(context_tensor)) for key in keys)
    lifted = WitnessedOrder(tuple(key_map[key] for key in row.positive),
                           tuple(key_map[key] for key in row.negative), tuple(maps),
                           f"tensor-context({row.label})")
    result = WitnessedOrderExport(lifted, items, key_pairs=tuple(zip(keys, names)))
    result.require_valid(limits)
    return result
