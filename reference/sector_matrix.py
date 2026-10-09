"""Checked matrix restrictions from proof sectors and independent columns.

Two branch words, (left,middle) and (right,middle), share the first matrix
coordinates and supply its two output columns. These are concrete coefficient
maps, not a low-rank matrix scheme or a finite catalytic witness.
"""

from .constraint_generators import WitnessedOrderExport
from .fields import FiniteField
from .maps import LinearMap
from .sectors import ThreeSectorConstruction
from .tensors import FiniteTensor, matrix_multiplication_tensor
from .witnessed_constraints import DEFAULT_LIMITS, StateAssemblyLimit, WitnessedOrder


def export_retained_two_copy_matrix_order(construction, source_key, target_key, *,
                                          limits=DEFAULT_LIMITS):
    """Export 2 R(a,h) -> M2 for a>=2,h>=1, using independent row copies.

    Each retained tensor supplies one matrix row: its middle output is one
    dot product; adding selected left/right outputs supplies the other.
    The second map duplicates the same B coordinates in both source copies.
    All copy occurrences and output additions are explicit local maps.
    """
    if not isinstance(construction, ThreeSectorConstruction) or construction.a < 2 or construction.h < 1:
        raise ValueError('supply a concrete sector with a>=2,h>=1')
    if any(not isinstance(key, str) or not key for key in (source_key, target_key)) or source_key == target_key:
        raise ValueError('source and target keys must be distinct nonempty strings')
    if limits.max_blocks < 2:
        raise StateAssemblyLimit('two-copy block cap reached')
    source_shape = tuple(2*size for size in construction.shape)
    limits.tensor(source_shape)
    limits.tensor((4, 4, 4))
    limits.maps(source_shape, (4, 4, 4))
    a, width, output = construction.shape
    middle, right = construction.middle_z_start, construction.right_start
    first = LinearMap.selection(construction.field, 2*a, (0, a-1, a, 2*a-1))
    y_indices = (middle, 0, construction.h, right)
    second = LinearMap(construction.field, 2*width,
        tuple(tuple(int(i % width == selected) for i in range(2*width)) for selected in y_indices))
    third_indices = ((middle,), (0, right+a-1), (output+middle,), (output, output+right+a-1))
    third = LinearMap(construction.field, 2*output,
        tuple(tuple(int(i in selected) for i in range(2*output)) for selected in third_indices))
    source = construction.retained
    target = matrix_multiplication_tensor(construction.field, 2, 2, 2)
    row = WitnessedOrder((source_key, source_key), (target_key,), (first, second, third),
                         'retained-two-copy-matrix-rows')
    exported = WitnessedOrderExport(row, ((source_key, source), (target_key, target)))
    exported.require_valid(limits)
    return exported


def export_matrix_column_sum_order(field, a, b, columns, source_key, target_key, *,
                                    limits=DEFAULT_LIMITS):
    """Export columns copies of M(a,b,1) -> M(a,b,columns) with actual maps.

    The first leg repeats each matrix entry in every independent source
    block. The other legs permute column-major source blocks to row-major
    matrix coordinates. Copies are integer occurrences, even if their count
    vanishes as a field scalar. No tensor-rank additivity is assumed.
    """
    if not isinstance(field, FiniteField):
        raise ValueError('supply an explicitly represented finite field')
    if any(type(n) is not int or n < 1 for n in (a, b, columns)):
        raise ValueError('matrix dimensions must be positive integers')
    if any(not isinstance(key, str) or not key for key in (source_key, target_key)) or source_key == target_key:
        raise ValueError('source and target keys must be distinct nonempty strings')
    if columns > limits.max_blocks:
        raise StateAssemblyLimit('column direct-sum block cap reached')
    source_shape = (columns*a*b, columns*b, columns*a)
    target_shape = (a*b, b*columns, a*columns)
    limits.tensor(source_shape)
    limits.tensor(target_shape)
    limits.maps(source_shape, target_shape)
    first = LinearMap(field, columns*a*b,
                      tuple(tuple(int(index % (a*b) == row) for index in range(columns*a*b))
                            for row in range(a*b)))
    second = LinearMap.selection(field, columns*b,
                                 tuple(k*b+j for j in range(b) for k in range(columns)))
    third = LinearMap.selection(field, columns*a,
                                tuple(k*a+i for i in range(a) for k in range(columns)))
    source = matrix_multiplication_tensor(field, a, b, 1)
    target = matrix_multiplication_tensor(field, a, b, columns)
    row = WitnessedOrder((source_key,)*columns, (target_key,), (first, second, third),
                         'matrix-independent-column-sum')
    exported = WitnessedOrderExport(row, ((source_key, source), (target_key, target)))
    exported.require_valid(limits)
    return exported


def export_retained_matrix_order(construction, source_key, target_key, *,
                                 limits=DEFAULT_LIMITS):
    """Export R(a,h)^2 -> M(a,a,2), for positive a,h over the given field.

    Each branch uses its local polynomial input coordinate zero. The middle
    branch reverses the second first-leg factor. The result retains both
    columns of an actual rectangular matrix product, without Fourier copies
    or interpolation. Tensor and map bounds are checked before allocation.
    """
    if not isinstance(construction, ThreeSectorConstruction):
        raise ValueError('supply a concrete three-sector construction')
    if not construction.a or not construction.h:
        raise ValueError('matrix extraction requires positive sector parameters')
    if any(not isinstance(key, str) or not key for key in (source_key, target_key)):
        raise ValueError('source and target keys must be nonempty strings')
    a, h = construction.a, construction.h
    source_shape = tuple(size*size for size in construction.shape)
    target_shape = (a*a, 2*a, 2*a)
    limits.tensor(source_shape)
    limits.tensor(target_shape)
    limits.maps(source_shape, target_shape)
    width, output = construction.shape[1:]
    right, middle = construction.right_start, construction.middle_z_start
    selections = (
        tuple(i*a+(a-1-j) for i in range(a) for j in range(a)),
        tuple((right*k)*width+h+j for j in range(a) for k in range(2)),
        tuple((i+right*k)*output+middle for i in range(a) for k in range(2)),
    )
    maps = tuple(LinearMap.selection(construction.field, size, selected)
                 for size, selected in zip(source_shape, selections))
    source = construction.retained.tensor_product(construction.retained)
    target = matrix_multiplication_tensor(construction.field, a, a, 2)
    if source_key == target_key and source != target:
        raise ValueError('different tensors require distinct inventory keys')
    inventory = ((source_key, source),) if source_key == target_key else (
        (source_key, source), (target_key, target))
    row = WitnessedOrder((source_key,), (target_key,), maps,
                         'retained-square-matrix-columns')
    exported = WitnessedOrderExport(row, inventory, exponent=2)
    exported.require_valid(limits)
    return exported


def export_star_matrix_order(field, source_key, target_key, *, limits=DEFAULT_LIMITS):
    """Export star^2 -> M2 with signed local maps, over any represented field.

    The star support is (0,0,1),(0,1,0),(1,0,2),(1,2,0). Its square sends
    A and (beta,u,v,W) to (beta*A,A*u,A^t*v,tr(A^t*W)). For B=(u,w), put
    beta=w0, v=(-w1,w0), W=w1*I. The trace correction to the transposed
    block recovers A*w. This uses the 2-by-2 adjugate identity and remains
    valid in characteristic two; signs are embedded field scalars.
    """
    if not isinstance(field, FiniteField):
        raise ValueError('supply an explicitly represented finite field')
    if any(not isinstance(key, str) or not key for key in (source_key, target_key)) or source_key == target_key:
        raise ValueError('source and target keys must be distinct nonempty strings')
    limits.tensor((4, 9, 9))
    limits.tensor((4, 4, 4))
    limits.maps((4, 9, 9), (4, 4, 4))
    support = {(0, 0, 1), (0, 1, 0), (1, 0, 2), (1, 2, 0)}
    star = FiniteTensor.from_function(field, (2, 3, 3),
                                    lambda i, j, k: int((i, j, k) in support))
    source = star.tensor_product(star)
    target = matrix_multiplication_tensor(field, 2, 2, 2)
    def sparse_map(rows):
        return LinearMap(field, 9, tuple(tuple(field.embed(row.get(i, 0)) for i in range(9))
                                        for row in rows))
    maps = (LinearMap.identity(field, 4),
            sparse_map(({1: 1}, {0: 1, 6: 1}, {2: 1}, {3: -1, 4: 1, 8: 1})),
            sparse_map(({3: 1}, {2: -1, 4: 1, 8: 1}, {6: 1}, {0: 1, 1: 1})))
    row = WitnessedOrder((source_key,), (target_key,), maps, 'star-square-matrix-trace-correction')
    exported = WitnessedOrderExport(row, ((source_key, source), (target_key, target)), exponent=2)
    exported.require_valid(limits)
    return exported
