"""Recover the two non-shared branch labels using independent coordinate maps."""

from itertools import product

from .maps import LinearMap
from .tensors import FiniteTensor, _tensor_family, shared_first_tensor


def tag_shared_scheme(scheme, branches, y_labels, z_labels):
    """Insert branch tags without duplicating X, from BranchTagging.lean.

    Both Y and Z coordinates must determine the label of every supported
    branch coefficient. Check each branch before field addition, so overlaps
    or cancellation in small characteristic cannot hide invalid labels.
    """
    scheme.require_exact()
    branches = _tensor_family(branches, uniform=True)
    x, y, z = branches[0].shape
    count = len(branches)
    y_labels, z_labels = tuple(y_labels), tuple(z_labels)
    for labels, size in ((y_labels, y), (z_labels, z)):
        if len(labels) != size or any(type(h) is not int or not 0 <= h < count for h in labels):
            raise ValueError("coordinate labels must match the side axis and branch count")
    if scheme.field != branches[0].field or scheme.shape != (x, y, z):
        raise ValueError("source scheme must use the branches' field and ambient shape")
    for h, branch in enumerate(branches):
        for i, j, k in product(range(x), range(y), range(z)):
            if branch.coefficient(i, j, k) and (y_labels[j] != h or z_labels[k] != h):
                raise ValueError("both side coordinates must determine each supported branch label")
    total = FiniteTensor.from_function(scheme.field, (x, y, z),
                                       lambda i, j, k: scheme.field.sum(b.coefficient(i, j, k) for b in branches))
    if scheme.target != total:
        raise ValueError("source scheme must target the coefficientwise sum of branches")
    maps = (LinearMap.identity(scheme.field, x),
            LinearMap(scheme.field, y, tuple(tuple(int(old == j and y_labels[j] == h) for old in range(y))
                                             for h in range(count) for j in range(y))),
            LinearMap(scheme.field, z, tuple(tuple(int(old == k and z_labels[k] == h) for old in range(z))
                                             for h in range(count) for k in range(z))))
    result = scheme.restrict(*maps)
    if result.target != shared_first_tensor(branches):
        raise ValueError("tag insertion does not match the shared-input tensor")
    return result
