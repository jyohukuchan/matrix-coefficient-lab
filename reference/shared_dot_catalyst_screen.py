"""Finite coefficient controls for independently audited shared-dot screens.

The all-field rational-function argument and slice-space geometry are
in research/shared-dot-catalysts.md. These controls check formal polynomials
and actual matrices, not the unbounded argument or a Lean theorem. Screens
never become ordinary normalized-state constraints.
"""

from dataclasses import dataclass

from .catalyst_compiler import CompilerLimits, CompilerResourceLimit
from .catalyst_search_experiment import field_matrix_rank
from .fields import FiniteField
from .maps import LinearMap
from .shared_dot_orders import shared_dot_degeneration, shared_dot_tensor, shared_dot_w_degeneration
from .shared_dot_koszul import (SharedDotKoszulLimit, SharedDotKoszulLimits,
                               verify_shared_dot_three_integer_minor, verify_shared_dot_two_integer_transform)
from .tensors import FiniteTensor, matrix_multiplication_tensor
from .witnessed_constraints import DEFAULT_LIMITS, StateAssemblyLimit, StateAssemblyLimits


def _parameters(field, q, gain, repetitions, limits):
    if not isinstance(field, FiniteField) or any(type(v) is not int or v < 1 for v in (q, gain, repetitions)):
        raise ValueError("supply a represented field and positive integer length, gain, repetitions")
    if not isinstance(limits, StateAssemblyLimits):
        raise ValueError("supply explicit state assembly limits")
    n = repetitions*(gain+4*(q+1))
    limits.tensor((n,)*3)
    limits.maps((n,)*3, (n,)*3)
    if 3*n > limits.max_constraint_entries or max(gain+1, repetitions) > limits.max_blocks:
        raise StateAssemblyLimit("shared-dot commutator bookkeeping or block cap reached")


def _slice(tensor, vector):
    field = tensor.field
    support = tuple((i, value) for i, value in enumerate(vector) if value)
    _, y, z = tensor.shape
    return LinearMap(field, y, tuple(tuple(field.sum(field.mul(value, tensor.coefficient(i, j, k))
        for i, value in support) for j in range(y)) for k in range(z)))


def _inverse(mapping):
    field, n = mapping.field, mapping.input_size
    if mapping.output_size != n:
        raise ValueError("pivot slice must be square")
    rows = [list(row)+[int(i == j) for j in range(n)] for i, row in enumerate(mapping.rows)]
    for column in range(n):
        pivot = next((i for i in range(column, n) if rows[i][column]), None)
        if pivot is None:
            raise ValueError("pivot slice is singular")
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = field.inv(rows[column][column])
        rows[column] = [field.mul(scale, value) for value in rows[column]]
        for i in range(n):
            if i != column and (scale := rows[i][column]):
                rows[i] = [field.sub(a, field.mul(scale, b)) for a, b in zip(rows[i], rows[column])]
    return LinearMap(field, n, tuple(tuple(row[n:]) for row in rows))


@dataclass(frozen=True)
class SharedDotCommutatorControl:
    auxiliary: FiniteTensor
    tensor: FiniteTensor
    length: int
    gain: int
    repetitions: int
    vectors: tuple[tuple[int, ...], ...]
    slices: tuple[LinearMap, ...]
    normalized: tuple[LinearMap, ...]
    commutator: LinearMap
    identity_rank: int
    commutator_rank: int

    @property
    def rank_lower_control(self):
        return self.identity_rank+(self.commutator_rank+1)//2

    def require_valid(self, limits=DEFAULT_LIMITS):
        field, q, gain, repetitions = self.tensor.field, self.length, self.gain, self.repetitions
        _parameters(field, q, gain, repetitions, limits)
        auxiliary = shared_dot_tensor(field, q, limits=limits)
        core = matrix_multiplication_tensor(field, 2, 2, 2).tensor_product(auxiliary)
        block = FiniteTensor.direct_sum((FiniteTensor.unit(field),)*gain+(core,))
        expected = FiniteTensor.direct_sum((block,)*repetitions)
        if self.auxiliary != auxiliary or self.tensor != expected:
            raise ValueError("shared-dot control coefficient tensor was corrupted")
        if len(self.vectors) != 3 or any(len(v) != expected.shape[0] for v in self.vectors):
            raise ValueError("shared-dot slice control vector has the wrong size")
        for vector in self.vectors:
            for value in vector:
                field.check(value)
        slices = tuple(_slice(expected, vector) for vector in self.vectors)
        n = expected.shape[1]
        inverse = _inverse(slices[0])
        normalized = tuple(inverse.compose(mapping) for mapping in slices)
        if normalized[0] != LinearMap.identity(field, n):
            raise ValueError("shared-dot pivot slice failed actual inversion")
        bc, cb = normalized[1].compose(normalized[2]), normalized[2].compose(normalized[1])
        commutator = LinearMap(field, n, tuple(tuple(field.sub(a, b) for a, b in zip(left, right))
            for left, right in zip(bc.rows, cb.rows)))
        identity_rank = field_matrix_rank(field, slices[0].rows, n)
        commutator_rank = field_matrix_rank(field, commutator.rows, n)
        if (self.slices != slices or self.normalized != normalized or self.commutator != commutator
                or self.identity_rank != identity_rank or self.commutator_rank != commutator_rank
                or identity_rank != repetitions*(gain+4*(q+1))
                or commutator_rank != repetitions*4*(q+1)):
            raise ValueError("shared-dot actual commutator ranks disagree with the asserted controls")


def shared_dot_commutator_control(field, q, *, gain=1, repetitions=1, limits=DEFAULT_LIMITS):
    """Actual invertible slice and E12/E21 commutator, including characteristic two."""
    _parameters(field, q, gain, repetitions, limits)
    auxiliary = shared_dot_tensor(field, q, limits=limits)
    core = matrix_multiplication_tensor(field, 2, 2, 2).tensor_product(auxiliary)
    block = FiniteTensor.direct_sum((FiniteTensor.unit(field),)*gain+(core,))
    tensor = FiniteTensor.direct_sum((block,)*repetitions)
    n, size = tensor.shape[0], q+1
    vectors = [[0]*n for _ in range(3)]
    for copy in range(repetitions):
        start = copy*(gain+4*size)
        for i in range(gain):
            vectors[0][start+i] = 1
        for matrix_index in (0, 3):
            vectors[0][start+gain+matrix_index*size] = 1
            vectors[0][start+gain+matrix_index*size+1] = 1
        for target, matrix_index in ((1, 1), (2, 2)):
            vectors[target][start+gain+matrix_index*size] = 1
            vectors[target][start+gain+matrix_index*size+1] = 1
    vectors = tuple(tuple(vector) for vector in vectors)
    slices = tuple(_slice(tensor, vector) for vector in vectors)
    inverse = _inverse(slices[0])
    normalized = tuple(inverse.compose(mapping) for mapping in slices)
    bc, cb = normalized[1].compose(normalized[2]), normalized[2].compose(normalized[1])
    commutator = LinearMap(field, n, tuple(tuple(field.sub(a, b) for a, b in zip(left, right))
        for left, right in zip(bc.rows, cb.rows)))
    control = SharedDotCommutatorControl(auxiliary, tensor, q, gain, repetitions, vectors,
        slices, normalized, commutator, field_matrix_rank(field, slices[0].rows, n),
        field_matrix_rank(field, commutator.rows, n))
    control.require_valid(limits)
    return control


def screen_shared_dot_catalyst(catalyst, *, auxiliary, m=1, limits=DEFAULT_LIMITS,
                              compiler_limits=CompilerLimits()):
    """Necessary analytic screen for any finite D and literal shared C_q.

    Integral Koszul certificates exclude q=2,3 for every finite D. Rational
    commutator transfer excludes q>=4 and the two-term q=1 degeneration.
    Passing either screen does not establish existence.
    """
    if not isinstance(catalyst, FiniteTensor) or not isinstance(auxiliary, FiniteTensor) or catalyst.field != auxiliary.field:
        raise ValueError("supply actual catalyst and auxiliary tensors over the same represented field")
    if type(m) is not int or m < 1:
        raise ValueError("gain must be a positive integer")
    if not isinstance(limits, StateAssemblyLimits) or not isinstance(compiler_limits, CompilerLimits):
        raise ValueError("supply explicit state/compiler limits")
    report = {"status": "unclassified", "reasons": [], "gain": m,
        "scope": "D + m Unit + M2 tensor C_q <= D + 5 C_q",
        "verification": "analytic arguments in research/shared-dot-catalysts.md; exact finite controls; not Lean",
        "finite_sampling_proves_all_fields": False}
    if len(set(auxiliary.shape)) != 1 or auxiliary.shape[0] < 2:
        report['reasons'].append("auxiliary is outside literal shared-dot coordinates")
        return report
    q, field = auxiliary.shape[0]-1, auxiliary.field
    try:
        limits.tensor(auxiliary.shape)
        if auxiliary != shared_dot_tensor(field, q, limits=limits):
            report['reasons'].append("auxiliary coefficients are not the shared-dot tensor")
            return report
        control = shared_dot_commutator_control(field, q, gain=m, limits=limits)
        degeneration = (shared_dot_w_degeneration(field, limits=compiler_limits) if q == 1 else
                        shared_dot_degeneration(field, q, limits=compiler_limits))
        koszul_limits = SharedDotKoszulLimits(max_matrix_entries=limits.max_map_entries,
                                             max_operations=limits.max_constraint_entries)
        universal = (verify_shared_dot_three_integer_minor(limits=koszul_limits) if q == 3 else
                     verify_shared_dot_two_integer_transform(limits=koszul_limits) if q == 2 else None)
    except (StateAssemblyLimit, CompilerResourceLimit, SharedDotKoszulLimit) as error:
        report.update(status="resource_cap", reasons=[str(error)])
        return report
    report.update(length=q, necessary_gain_boundary=0,
        coefficient_controls={"invertible_slice_rank": control.identity_rank,
            "commutator_rank": control.commutator_rank,
            "rank_lower_control": control.rank_lower_control,
            "source_border_terms_per_five_copies": 5*degeneration.terms,
            "formal_leading_degree": degeneration.leading,
            "sampling_scope": "finite coefficient controls; rational-function iteration is analytic"})
    if q in (2, 3):
        source_slope = 5*degeneration.terms
        report.update(status="analytically_excluded", exclusion="all_field_integral_koszul",
            reasons=["the integral unimodular Koszul certificate and additive iteration exclude every finite catalyst"],
            arbitrary_finite_catalyst=True,
            integral_minor={"size": universal.minor_size, "determinant": universal.minor_determinant,
                "rank_one_factor": universal.rank_one_factor, "target_border_lower": universal.border_rank_lower,
                "iteration_inequality": f"{universal.minor_size}*n <= 6*(C_D+{source_slope}*n), impossible for unbounded n"})
    elif q == 1 or q+m > 4:
        report.update(status="analytically_excluded", exclusion="all_field_border_commutator",
            reasons=["the target commutator slope exceeds the supplied source polynomial slope"],
            arbitrary_finite_catalyst=True)
    elif not any(catalyst.coefficients) and q >= 3:
        report.update(status="analytically_excluded", exclusion="zero_catalyst_slice_space",
            reasons=["the audited all-characteristic shared-dot slice-space argument excludes D=0"])
    else:
        report['reasons'].append("these necessary screens do not decide this catalyst")
    return report
