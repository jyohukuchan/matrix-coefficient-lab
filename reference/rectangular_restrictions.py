"""Concrete sharing of a left operand across independent matrix blocks.

This finite restriction is used in catalyst lower-bound research. It gives
actual coefficient maps, not an implementation of an external rank theorem.
"""

from .catalyst_compiler import DEFAULT_LIMITS, CompiledRestriction
from .maps import LinearMap
from .tensors import FiniteTensor, matrix_multiplication_tensor, shared_first_tensor


def share_first_axis(tensor, copies, *, scalar_gain=0, limits=DEFAULT_LIMITS):
    """Merge the first leg of copies independent tensors, retaining scalar gains.

    Y/Z still have separate copy labels. The first tensor coordinates sum
    across source copies, while every scalar gain retains its own coordinates
    on all three legs. Useful for fixed-first-axis flattening amplification.
    """
    if not isinstance(tensor, FiniteTensor) or type(copies) is not int or copies < 1:
        raise ValueError("supply a finite tensor and a positive integer copy count")
    if type(scalar_gain) is not int or scalar_gain < 0:
        raise ValueError("scalar gain must be a nonnegative integer copy count")
    source_shape = tuple(copies*(size+scalar_gain) for size in tensor.shape)
    target_shape = (copies*scalar_gain+tensor.shape[0],
                    copies*scalar_gain+copies*tensor.shape[1],
                    copies*scalar_gain+copies*tensor.shape[2])
    limits.maps(source_shape, target_shape)
    field, unit = tensor.field, FiniteTensor.unit(tensor.field)
    block = FiniteTensor.direct_sum((unit,)*scalar_gain+(tensor,))
    source = FiniteTensor.direct_sum((block,)*copies)
    shared = shared_first_tensor((tensor,)*copies)
    target = FiniteTensor.direct_sum((unit,)*(copies*scalar_gain)+(shared,))
    maps = []
    for axis, width in enumerate(tensor.shape):
        stride = scalar_gain+width
        rows = []
        for h in range(copies):
            for gain in range(scalar_gain):
                selected = h*stride+gain
                rows.append(tuple(int(i == selected) for i in range(source_shape[axis])))
        if axis == 0:
            for out in range(width):
                selected = {h*stride+scalar_gain+out for h in range(copies)}
                rows.append(tuple(int(i in selected) for i in range(source_shape[axis])))
        else:
            for h in range(copies):
                for out in range(width):
                    selected = h*stride+scalar_gain+out
                    rows.append(tuple(int(i == selected) for i in range(source_shape[axis])))
        maps.append(LinearMap(field, source_shape[axis], rows))
    result = CompiledRestriction(source, target, tuple(maps))
    result.require_valid()
    return result


def share_matrix_left_operand(field, n, inner, m, copies, *, scalar_gain=0, limits=DEFAULT_LIMITS):
    """copies*(gain*unit ⊕ M(n,inner,m)) -> gain*copies*unit ⊕ M(n,inner,m*copies).

    Matrix X rows sum the same entry in every independent block. Y and Z
    selections regroup their column indices. Scalar blocks stay independent.
    The direct-sum labels in the source prevent unwanted cross products.
    """
    if any(type(value) is not int or value < 1 for value in (n, inner, m, copies)):
        raise ValueError("matrix dimensions and copy count must be positive integers")
    if type(scalar_gain) is not int or scalar_gain < 0:
        raise ValueError("scalar gain must be a nonnegative integer copy count")
    matrix_shape = (n*inner, inner*m, n*m)
    source_shape = tuple(copies*(size+scalar_gain) for size in matrix_shape)
    target_shape = (copies*scalar_gain+n*inner,
                    copies*scalar_gain+inner*m*copies,
                    copies*scalar_gain+n*m*copies)
    limits.maps(source_shape, target_shape)
    unit = FiniteTensor.unit(field)
    matrix = matrix_multiplication_tensor(field, n, inner, m)
    block = FiniteTensor.direct_sum((unit,)*scalar_gain+(matrix,))
    source = FiniteTensor.direct_sum((block,)*copies)
    rectangle = matrix_multiplication_tensor(field, n, inner, copies*m)
    target = FiniteTensor.direct_sum((unit,)*(copies*scalar_gain)+(rectangle,))
    rows = []
    stride = scalar_gain+matrix_shape[0]
    for block_index in range(copies):
        for gain in range(scalar_gain):
            selected = block_index*stride+gain
            rows.append(tuple(int(i == selected) for i in range(source_shape[0])))
    for x in range(matrix_shape[0]):
        selected = {block_index*stride+scalar_gain+x for block_index in range(copies)}
        rows.append(tuple(int(i in selected) for i in range(source_shape[0])))
    first = LinearMap(field, source_shape[0], rows)
    maps = [first]
    for axis, outer in ((1, inner), (2, n)):
        stride = scalar_gain+matrix_shape[axis]
        selected = [block_index*stride+gain for block_index in range(copies) for gain in range(scalar_gain)]
        selected.extend(block_index*stride+scalar_gain+row*m+column
                        for row in range(outer) for block_index in range(copies) for column in range(m))
        maps.append(LinearMap.selection(field, source_shape[axis], selected))
    result = CompiledRestriction(source, target, tuple(maps))
    result.require_valid()
    return result
