"""Remove independent scalar summands from an exact rank certificate.

Only scalar units have this elimination rule. It does NOT permit cancelling
a general catalyst from a tensor restriction or adding arbitrary tensor ranks.
"""

from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor


def eliminate_scalar_summands(scheme, core, copies=1):
    """Return a verified core decomposition with scheme.terms-copies terms.

    Require target == core direct-sum copies*unit on all three axes. For each
    final scalar output, solve for one rank-one term, then project the other
    two legs away from its scalar coordinate. This is ordinary exact linear
    elimination; all division is by a checked nonzero coefficient.
    """
    if not isinstance(scheme, TensorScheme) or not isinstance(core, FiniteTensor):
        raise ValueError("supply an exact tensor scheme and exact core tensor")
    if type(copies) is not int or copies < 0 or core.field != scheme.field:
        raise ValueError("scalar copies must be nonnegative and fields must match")
    expected_shape = tuple(size+copies for size in core.shape)
    if scheme.shape != expected_shape or copies > scheme.terms:
        raise ValueError("scheme cannot have this core and scalar-summand count")
    # Shape validation bounds copies before building any scalar block list.
    unit = FiniteTensor.unit(scheme.field)
    expected = FiniteTensor.direct_sum((core,)+(unit,)*copies)
    if scheme.target != expected:
        raise ValueError("target must be exactly the core plus independent scalar units")
    scheme.require_exact()
    result = scheme
    for remaining in range(copies, 0, -1):
        pivot = next((q for q, vector in enumerate(result.c) if vector[-1]), None)
        if pivot is None:
            raise ValueError("scalar output has no nonzero pivot coefficient")
        f = result.field
        inverse = f.inv(result.c[pivot][-1])
        aa, bb, cc = [], [], []
        for q, (a, b, c) in enumerate(zip(result.a, result.b, result.c)):
            if q == pivot:
                continue
            factor = f.mul(c[-1], inverse)
            aa.append(a[:-1])
            bb.append(b[:-1])
            cc.append(tuple(f.sub(value, f.mul(factor, old))
                            for value, old in zip(c[:-1], result.c[pivot][:-1])))
        target = FiniteTensor.direct_sum((core,)+(unit,)*(remaining-1))
        result = TensorScheme(target, aa, bb, cc)
        result.require_exact()
    return result
