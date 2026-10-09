"""Actual shared-dot quotients and small polynomial border-rank controls.

C_q has q+2 polynomial terms at leading degree three. That is a verified
degeneration, not an exact q+2-term rank scheme, a catalytic certificate, or
a 9/4 matrix generator. Base-field recovery explicitly pays for its nodes.
"""

from dataclasses import dataclass
from itertools import product
from math import prod

from .catalyst_compiler import CompilerLimits, CompilerResourceLimit
from .constraint_generators import WitnessedOrderExport
from .constructions import constant_weights, nonzero_nodes
from .fields import FiniteField
from .maps import LinearMap
from .mixed_dot_lengths import mixed_dot_length_tensor
from .polynomials import Polynomial
from .subrank_constraints import export_selected_subrank_order
from .tensor_degenerations import TensorDegeneration
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor
from .witnessed_constraints import (DEFAULT_LIMITS as DEFAULT_STATE_LIMITS,
                                    StateAssemblyLimit, StateAssemblyLimits, WitnessedOrder)


DEFAULT_COMPILER_LIMITS = CompilerLimits()


def _parameters(field, q):
    if not isinstance(field, FiniteField) or type(q) is not int or q < 1:
        raise ValueError("supply a represented field and positive integer shared-dot length")


def _state_preflight(field, q, source_shape, target_shape, limits):
    _parameters(field, q)
    if not isinstance(limits, StateAssemblyLimits):
        raise ValueError("supply explicit state assembly limits")
    limits.tensor(source_shape)
    limits.tensor(target_shape)
    limits.maps(source_shape, target_shape)
    if sum(source_shape)+sum(target_shape) > limits.max_constraint_entries:
        raise StateAssemblyLimit("shared-dot bookkeeping cap reached; construction incomplete")


def _names(first, second):
    if any(not isinstance(key, str) or not key for key in (first, second)) or first == second:
        raise ValueError("tensor IDs must be distinct nonempty strings")


def shared_dot_tensor(field, q, *, limits=DEFAULT_STATE_LIMITS):
    """C_q(0,i,i)=C_q(i,0,i)=C_q(i,i,0)=1 for i=1,...,q."""
    _parameters(field, q)
    shape = (q+1,)*3
    _state_preflight(field, q, shape, (0, 0, 0), limits)
    return _shared_dot_tensor_unchecked(field, q)


def _shared_dot_tensor_unchecked(field, q):
    shape = (q+1,)*3
    return FiniteTensor.from_function(field, shape,
        lambda i, j, k: int((i == 0 and j == k > 0) or (j == 0 and i == k > 0)
                           or (k == 0 and i == j > 0)))


def shared_dot_scheme(field, q, *, limits=DEFAULT_COMPILER_LIMITS):
    """An actual same-field rank upper of at most 2q+1; no optimality claim.

    Odd characteristic uses the average of (e0+ei)^3 and (e0-ei)^3,
    subtracting q*e0^3. If q casts to zero the last term is omitted. In
    characteristic two, a different cancellation identity avoids division.
    All returned coefficients are independently verified against C_q.
    """
    _parameters(field, q)
    if not isinstance(limits, CompilerLimits):
        raise ValueError("supply explicit compiler limits")
    terms = 2*q+1 if field.p == 2 else 2*q+int(field.embed(q) != 0)
    limits.scheme((q+1,)*3, terms)
    vector = lambda *indices: tuple(int(j in indices) for j in range(q+1))
    e0 = vector(0)
    aa, bb, cc = [], [], []
    if field.p == 2:
        p, w = vector(q), (0,)+(1,)*q
        for j in range(1, q):
            aa.extend((vector(q, j), vector(0, q, j)))
            bb.extend((vector(j), vector(0, j)))
            cc.extend((vector(0, j), vector(j)))
        aa.extend((p, e0, vector(0, q)))
        bb.extend((w, vector(0, q), e0))
        cc.extend((e0, p, w))
    else:
        half, minus_one = field.inv(field.embed(2)), field.neg(1)
        for i in range(1, q+1):
            for sign in (1, minus_one):
                v = tuple(1 if j == 0 else sign if j == i else 0 for j in range(q+1))
                aa.append(tuple(field.mul(half, x) for x in v))
                bb.append(v)
                cc.append(v)
        correction = field.neg(field.embed(q))
        if correction:
            aa.append((correction,)+(0,)*q)
            bb.append(e0)
            cc.append(e0)
    scheme = TensorScheme(_shared_dot_tensor_unchecked(field, q), aa, bb, cc)
    scheme.require_exact()
    return scheme


def export_shared_dot_quotient_order(field, q, source_key, target_key, *, limits=DEFAULT_STATE_LIMITS):
    """C_q <= S_q by actual local sums, checked at every coefficient."""
    _names(source_key, target_key)
    _parameters(field, q)
    source_shape, target_shape = (2*q+1,)*3, (q+1,)*3
    _state_preflight(field, q, source_shape, target_shape, limits)
    if limits.max_blocks < 3:
        raise StateAssemblyLimit("shared-dot source block cap reached; construction incomplete")
    size = 2*q+1
    selections = ((0,), (q,), (2*q,))
    pairs = (lambda i: (i, q+i), lambda i: (i-1, q+i), lambda i: (i-1, q+i-1))
    maps = tuple(LinearMap(field, size,
        tuple(tuple(int(column in coordinates) for column in range(size))
              for coordinates in (selections[axis],)+tuple(pairs[axis](i) for i in range(1, q+1))))
        for axis in range(3))
    source, target = mixed_dot_length_tensor(field, q, limits=limits), shared_dot_tensor(field, q, limits=limits)
    row = WitnessedOrder((source_key,), (target_key,), maps, "shared-dot-quotient")
    exported = WitnessedOrderExport(row, ((source_key, source), (target_key, target)))
    exported.require_valid(limits)
    return exported


def _polynomial_preflight(shape, terms, degree, limits):
    if not isinstance(limits, CompilerLimits):
        raise ValueError("supply explicit compiler limits")
    limits.tensor(shape)
    # Includes a conservative bound for local polynomial coefficient storage
    # and the complete formal tensor-coordinate check before construction.
    count = (degree+1)*(terms*sum(shape)+prod(shape))
    if terms > limits.max_terms or count > limits.max_scheme_scalars:
        raise CompilerResourceLimit("polynomial coefficient cap reached; construction incomplete")


def require_shared_dot_identity(degeneration, *, limits=DEFAULT_COMPILER_LIMITS):
    """Verify every formal coefficient, including the degree-four/six noise."""
    if not isinstance(degeneration, TensorDegeneration):
        raise ValueError("supply an actual polynomial tensor degeneration")
    field = degeneration.field
    if len(set(degeneration.shape)) != 1 or degeneration.shape[0] < 2 or degeneration.leading != 3:
        raise ValueError("shared-dot degeneration must have equal axes and leading degree three")
    _polynomial_preflight(degeneration.shape, degeneration.terms, degeneration.degree_bound, limits)
    q = degeneration.shape[0]-1
    for i, j, k in product(range(q+1), repeat=3):
        target_coefficient = int((i == 0 and j == k > 0) or (j == 0 and i == k > 0)
                                 or (k == 0 and i == j > 0))
        if degeneration.target.coefficient(i, j, k) != target_coefficient:
            raise ValueError("degeneration must target the actual shared-dot tensor")
        one_zero = int((i == 0)+(j == 0)+(k == 0) == 1)
        fourth = field.sub(int(i == j == k > 0), one_zero)
        sixth = field.neg(int(i > 0 and j > 0 and k > 0))
        expected = Polynomial(field, (0, 0, 0, target_coefficient, fourth, 0, sixth))
        if degeneration.tensor_polynomial(i, j, k) != expected:
            raise ValueError(f"shared-dot full formal identity fails at {(i, j, k)}")
    degeneration.require_valid()


def shared_dot_degeneration(field, q, *, limits=DEFAULT_COMPILER_LIMITS):
    """The actual q+2 terms of t sum_i(e0+t ei)^3 -(e0+t²W)^3 +(1-qt)e0³."""
    _parameters(field, q)
    shape, terms = (q+1,)*3, q+2
    _polynomial_preflight(shape, terms, 6, limits)
    zero, one, t = Polynomial(field), Polynomial.constant(field, 1), Polynomial.monomial(field, 1)
    aa, bb, cc = [], [], []
    for i in range(1, q+1):
        vector = tuple(one if j == 0 else t if j == i else zero for j in range(q+1))
        aa.append(tuple(t*p for p in vector))
        bb.append(vector)
        cc.append(vector)
    w = (one,)+(Polynomial.monomial(field, 2),)*q
    aa.append(tuple(-p for p in w))
    bb.append(w)
    cc.append(w)
    aa.append((Polynomial(field, (1, field.neg(field.embed(q)))),)+(zero,)*q)
    bb.append((one,)+(zero,)*q)
    cc.append((one,)+(zero,)*q)
    target = _shared_dot_tensor_unchecked(field, q)
    degeneration = TensorDegeneration(target, 3, aa, bb, cc)
    require_shared_dot_identity(degeneration, limits=limits)
    return degeneration


def shared_dot_w_degeneration(field, *, limits=DEFAULT_COMPILER_LIMITS):
    """The two-term (e1+t e0)^3-e1^3 family has leading t*C1 in every field."""
    _parameters(field, 1)
    _polynomial_preflight((2,)*3, 2, 3, limits)
    zero, one = Polynomial(field), Polynomial.constant(field, 1)
    vector = (Polynomial.monomial(field, 1), one)
    e1 = (zero, one)
    degeneration = TensorDegeneration(_shared_dot_tensor_unchecked(field, 1), 1,
        (vector, (zero, -one)), (vector, e1), (vector, e1))
    degeneration.require_valid()
    return degeneration


@dataclass(frozen=True)
class SharedDotPowerRecovery:
    exponent: int
    degeneration: TensorDegeneration
    scheme: TensorScheme
    nodes: tuple[int, ...]


def recover_shared_dot_power(field, q, exponent=1, *, nodes=None, limits=DEFAULT_COMPILER_LIMITS):
    """Recover a small powered exact scheme over this same field using 3e+1 nodes."""
    _parameters(field, q)
    if type(exponent) is not int or exponent < 0:
        raise ValueError("the exponent must be a nonnegative integer")
    if not isinstance(limits, CompilerLimits):
        raise ValueError("supply explicit compiler limits")
    # (q+1)^(3e) >= 2^(3e), so this guard precedes powers and long loops.
    if 3*exponent >= limits.max_tensor_coefficients.bit_length():
        raise CompilerResourceLimit("shared-dot power cap reached; construction incomplete")
    size, terms, count = (q+1)**exponent, (q+2)**exponent, 3*exponent+1
    shape = (size,)*3
    _polynomial_preflight(shape, terms, 6*exponent, limits)
    selected = nonzero_nodes(field, count) if nodes is None else tuple(nodes)
    if len(selected) < count:
        raise ValueError("shared-dot recovery requires at least 3*exponent+1 distinct nonzero nodes")
    limits.scheme(shape, len(selected)*terms)
    constant_weights(field, selected)  # validate actual nodes before any polynomial powers
    if exponent == 0:
        unit_scheme = TensorScheme.from_tensor(FiniteTensor.unit(field))
        degeneration = TensorDegeneration.from_scheme(unit_scheme)
    else:
        _polynomial_preflight((q+1,)*3, q+2, 6, limits)
        degeneration = shared_dot_degeneration(field, q, limits=limits).tensor_power(exponent)
    degeneration.require_valid()
    scheme = degeneration.recover_scheme(selected)
    scheme.require_exact()
    return SharedDotPowerRecovery(exponent, degeneration, scheme, selected)


def export_shared_dot_subrank_order(field, q, source_key, *, unit_key="unit", limits=DEFAULT_STATE_LIMITS):
    """Three actual independent units <= C_q for q>=3; no subrank optimality."""
    _parameters(field, q)
    _names(source_key, unit_key)
    if q < 3:
        raise ValueError("these three-unit selectors require q>=3")
    _state_preflight(field, q, (q+1,)*3, (3,)*3, limits)
    if limits.max_blocks < 3:
        raise StateAssemblyLimit("shared-dot unit block cap reached; construction incomplete")
    return export_selected_subrank_order(shared_dot_tensor(field, q, limits=limits),
                                         ((0, 2, 3), (1, 0, 3), (1, 2, 0)), source_key,
                                         unit_key=unit_key, limits=limits)
