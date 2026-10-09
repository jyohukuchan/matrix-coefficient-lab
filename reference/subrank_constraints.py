"""Exact diagonal subrank restrictions for finite witnessed state rows.

An export witnesses n independent units <= T by three actual local maps.
It yields the additive-state lower bound state(T) >= n; n is a subrank
lower bound, not the tensor rank or an assertion of rank additivity.
Every export is checked at all tensor coefficients over the supplied field.
"""

from itertools import product

from .constraint_generators import WitnessedOrderExport
from .convolution import lagrange_basis
from .fields import FiniteField
from .maps import LinearMap, validate_local_maps
from .polynomials import Polynomial
from .sectors import convolution_tensor
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, matrix_multiplication_tensor
from .witnessed_constraints import (DEFAULT_LIMITS, StateAssemblyLimit,
                                    StateAssemblyLimits, WitnessedOrder)


def _natural(value, name):
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")


def _field(field):
    if not isinstance(field, FiniteField):
        raise ValueError("supply an explicitly represented finite field")


def _names(source_key, unit_key):
    if any(not isinstance(key, str) or not key for key in (source_key, unit_key)):
        raise ValueError("source and unit IDs must be nonempty strings")


def _preflight(source_shape, n, limits):
    """Bound dense tensors/maps and even bookkeeping for empty-axis inputs."""
    if not isinstance(limits, StateAssemblyLimits):
        raise ValueError("supply explicit state assembly limits")
    _natural(n, "diagonal size")
    if n > limits.max_blocks or limits.max_blocks < 1:
        raise StateAssemblyLimit("diagonal block cap reached; construction incomplete")
    target_shape = (n, n, n)
    limits.tensor(source_shape)
    limits.tensor(target_shape)
    limits.tensor((1, 1, 1))  # The inventory also contains the scalar unit.
    limits.maps(source_shape, target_shape)
    # LinearMap.column_supports also allocates one list per source coordinate,
    # including when dense tensor or map sizes vanish because n or an axis is 0.
    if sum(source_shape)+3*n > limits.max_constraint_entries:
        raise StateAssemblyLimit("axis bookkeeping cap reached; construction incomplete")


def _source_tensor(source):
    if isinstance(source, TensorScheme):
        return source.target
    if not isinstance(source, FiniteTensor):
        raise ValueError("source must be an exact finite tensor or tensor scheme")
    return source


def export_subrank_order(source, maps, source_key, *, unit_key="unit",
                         limits=DEFAULT_LIMITS):
    """Export a supplied diagonal lower restriction, after coefficient checking.

    LinearMap rows are [new diagonal coordinate][old source coordinate].
    A TensorScheme source is additionally required to be an exact scheme.
    Repeated negative unit IDs denote independent direct-sum occurrences,
    so they never collapse into the field scalar n, even in characteristic 2.
    """
    tensor = _source_tensor(source)
    _field(tensor.field)
    _names(source_key, unit_key)
    maps = validate_local_maps(tensor.field, tensor.shape, maps)
    sizes = tuple(mapping.output_size for mapping in maps)
    if len(set(sizes)) != 1:
        raise ValueError("subrank maps must have equal diagonal output sizes")
    n = sizes[0]
    _preflight(tensor.shape, n, limits)
    unit = FiniteTensor.unit(tensor.field)
    if source_key == unit_key and tensor != unit:
        raise ValueError("the unit ID must denote the actual scalar unit")
    if isinstance(source, TensorScheme):
        source.require_exact()
    items = ((source_key, tensor),) if source_key == unit_key else (
        (source_key, tensor), (unit_key, unit))
    row = WitnessedOrder((source_key,), (unit_key,)*n, maps, "diagonal-subrank-lower")
    exported = WitnessedOrderExport(row, items)
    exported.require_valid(limits)
    return exported


def export_selected_subrank_order(source, selections, source_key, *, unit_key="unit",
                                  limits=DEFAULT_LIMITS):
    """Select three coordinate lists, then check their complete diagonal tensor.

    This accepts FiniteTensor or TensorScheme sources. Selection alone does
    not establish a subrank bound: any unwanted off-diagonal coefficient is
    rejected by export_subrank_order.
    """
    tensor = _source_tensor(source)
    _field(tensor.field)
    _names(source_key, unit_key)
    selections = tuple(tuple(indices) for indices in selections)
    if len(selections) != 3 or len({len(indices) for indices in selections}) != 1:
        raise ValueError("supply three coordinate lists with equal diagonal sizes")
    _preflight(tensor.shape, len(selections[0]), limits)
    maps = tuple(LinearMap.selection(tensor.field, size, indices)
                 for size, indices in zip(tensor.shape, selections))
    return export_subrank_order(source, maps, source_key, unit_key=unit_key, limits=limits)


def export_convolution_subrank_order(field, a, b, n, source_key, *, unit_key="unit",
                                     nodes=None, limits=DEFAULT_LIMITS):
    """Witness n Unit <= C(a,b), for n <= min(a,b) and n available nodes.

    For V[u,i] = nodes[u]**i, input maps are (V^-1)^transpose,
    padded by zero columns for old input degrees >= n. The output map is
    Vout[w,k] = nodes[w]**k for every convolution output degree. Thus the
    transformed coefficient at (u,v,w) is L_u(nodes[w])*L_v(nodes[w]),
    exactly 1 for u=v=w and 0 otherwise. No sampled-input test is used.
    Canonical encodings enumerate extension elements, not integer casts.
    """
    _field(field)
    _names(source_key, unit_key)
    for value, name in ((a, "first input dimension"), (b, "second input dimension"),
                        (n, "diagonal size")):
        _natural(value, name)
    if n > min(a, b):
        raise ValueError("diagonal size exceeds a convolution input dimension")
    if n > field.order:
        raise ValueError("not enough field elements for the diagonal subrank nodes")
    shape = (a, b, max(0, a+b-1))
    _preflight(shape, n, limits)
    nodes = tuple(range(n)) if nodes is None else tuple(nodes)
    if len(nodes) != n:
        raise ValueError("supply exactly n distinct diagonal subrank nodes")
    # Lagrange coefficients are the columns of the inverse Vandermonde.
    basis = lagrange_basis(field, nodes)
    maps = tuple(LinearMap(field, size,
                          tuple(tuple(p.coefficient(i) if i < n else 0
                                      for i in range(size)) for p in basis))
                 for size in (a, b))
    output = LinearMap(field, shape[2],
                       tuple(tuple(field.pow(node, k) for k in range(shape[2]))
                             for node in nodes))
    exported = export_subrank_order(convolution_tensor(field, a, b), (*maps, output),
                                    source_key, unit_key=unit_key, limits=limits)
    return WitnessedOrderExport(exported.row, exported.inventory, nodes)


def export_matrix_subrank_order(field, d, source_key, *, n=None, unit_key="unit",
                                limits=DEFAULT_LIMITS):
    """Witness d Unit <= M_d by selecting the diagonal matrix coordinates.

    An optional n <= d selects only the first n diagonal positions. No
    distinct field elements are required: this works over every field.
    """
    _field(field)
    _names(source_key, unit_key)
    if type(d) is not int or d < 1:
        raise ValueError("matrix dimension must be a positive integer")
    n = d if n is None else n
    _natural(n, "diagonal size")
    if n > d:
        raise ValueError("diagonal size exceeds the matrix dimension")
    shape = (d*d,)*3
    _preflight(shape, n, limits)
    indices = tuple(i*d+i for i in range(n))
    tensor = matrix_multiplication_tensor(field, d, d, d)
    return export_selected_subrank_order(tensor, (indices,)*3, source_key,
                                         unit_key=unit_key, limits=limits)


def export_projective_convolution_subrank_order(field, a, b, n, source_key, *,
                                                unit_key="unit", nodes=None,
                                                limits=DEFAULT_LIMITS):
    """Witness n Unit <= C(a,b) using n-1 finite nodes and infinity.

    Requires n <= min(a,b,field.order+1). Finite input rows contain the
    ordinary Lagrange polynomials on the n-1 finite nodes, of degree <= n-2.
    The final input row is the monic vanishing polynomial of degree n-1.
    All input rows are padded by zeros at old degrees >= n. Finite output
    rows evaluate all source output coefficients. The infinity output row
    selects degree 2*(n-1), which need not be the highest source degree.

    Its selected coefficient is one only for the product of the two monic
    vanishing polynomials, and zero for all other row pairs. Infinity is a
    coordinate functional, never a field element or a None encoding. The
    export's nodes metadata contains only its finite nodes.
    """
    _field(field)
    _names(source_key, unit_key)
    for value, name in ((a, "first input dimension"), (b, "second input dimension"),
                        (n, "diagonal size")):
        _natural(value, name)
    if n > min(a, b):
        raise ValueError("diagonal size exceeds a convolution input dimension")
    if n > field.order+1:
        raise ValueError("not enough finite field elements plus infinity for the diagonal subrank nodes")
    if n == 0:
        return export_convolution_subrank_order(field, a, b, n, source_key,
                                                unit_key=unit_key, nodes=nodes, limits=limits)
    shape = (a, b, a+b-1)
    _preflight(shape, n, limits)
    nodes = tuple(range(n-1)) if nodes is None else tuple(nodes)
    if len(nodes) != n-1:
        raise ValueError("supply exactly n-1 distinct finite projective subrank nodes")
    finite_basis = lagrange_basis(field, nodes)
    infinity = Polynomial.constant(field, 1)
    for node in nodes:
        infinity = infinity*Polynomial(field, (field.neg(node), 1))
    basis = (*finite_basis, infinity)
    inputs = tuple(LinearMap(field, size,
                            tuple(tuple(p.coefficient(i) if i < n else 0
                                        for i in range(size)) for p in basis))
                   for size in (a, b))
    output_rows = tuple(tuple(field.pow(node, k) for k in range(shape[2])) for node in nodes)
    output_rows += (tuple(int(k == 2*(n-1)) for k in range(shape[2])),)
    output = LinearMap(field, shape[2], output_rows)
    exported = export_subrank_order(convolution_tensor(field, a, b), (*inputs, output),
                                    source_key, unit_key=unit_key, limits=limits)
    row = WitnessedOrder(exported.row.positive, exported.row.negative,
                         exported.row.maps, "projective-infinity-convolution-subrank-lower")
    return WitnessedOrderExport(row, exported.inventory, nodes)


def _lower_data(exported):
    if not isinstance(exported, WitnessedOrderExport):
        raise ValueError("product inputs must be witnessed ordinary subrank exports")
    row, registry = exported.row, exported.registry
    if len(row.positive) != 1 or row.positive[0] not in registry:
        raise ValueError("a subrank input must have exactly one source tensor")
    source = registry[row.positive[0]]
    if not isinstance(source, FiniteTensor):
        raise ValueError("a subrank input must register an exact finite source tensor")
    validate_local_maps(source.field, source.shape, row.maps)
    if any(key not in registry or registry[key] != FiniteTensor.unit(source.field)
           for key in row.negative):
        raise ValueError("a subrank input target must consist of independent scalar units")
    n = len(row.negative)
    if tuple(mapping.output_size for mapping in row.maps) != (n,)*3:
        raise ValueError("subrank input map sizes do not match the diagonal units")
    return source, n


def export_subrank_product_order(first, second, source_key, *, unit_key="unit",
                                 limits=DEFAULT_LIMITS):
    """Witness n*m Unit <= A tensor B from two checked lower restrictions.

    The tensor product of the diagonal targets is distributed explicitly:
    pair (u,v) labels the independent unit occurrence u*m+v on EVERY leg.
    Source coordinates use the same lexicographic pair convention. This
    constructs actual maps; it does not multiply values of an additive state.
    """
    _names(source_key, unit_key)
    a, n = _lower_data(first)
    b, m = _lower_data(second)
    if a.field != b.field:
        raise ValueError("subrank product inputs must use the same field representation")
    shape = tuple(x*y for x, y in zip(a.shape, b.shape))
    _preflight(shape, n*m, limits)
    # Validate originals only after the product allocation preflight.
    first.require_valid(limits)
    second.require_valid(limits)
    f, maps = a.field, []
    for axis in range(3):
        left, right = first.row.maps[axis], second.row.maps[axis]
        rows = tuple(tuple(f.mul(x, y) for x in left.rows[u] for y in right.rows[v])
                     for u, v in product(range(n), range(m)))
        maps.append(LinearMap(f, shape[axis], rows))
    return export_subrank_order(a.tensor_product(b), maps, source_key,
                                unit_key=unit_key, limits=limits)
