"""Necessary screen for one dot_2 in each of the three singleton orientations.

The analytic exclusions use the slice-space arguments in
research/variable-auxiliaries.md and the rank/substitution argument in
research/mixed-dot-rank-catalysts.md, not a finite-field search or a Lean theorem.
This module checks its coefficient formulas and finite rank-profile arithmetic.
Unsupported auxiliaries remain unclassified; the literal auxiliary is excluded
for every finite catalyst. Rank bounds here are
catalyst screens, never ordinary normalized-state or subrank constraints.
"""

from dataclasses import dataclass
from itertools import product
from math import prod

from .catalyst_search_experiment import field_matrix_rank
from .fields import FiniteField
from .schemes import strassen_scheme
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, matrix_multiplication_tensor


class MixedDotScreenLimit(ValueError):
    """A declared resource cap prevents completing the controls."""


@dataclass(frozen=True)
class MixedDotScreenBudget:
    max_tensor_entries: int = 100_000
    max_matrix_entries: int = 100_000
    max_profile_cases: int = 256
    max_slice_cases: int = 4096

    def __post_init__(self):
        if any(type(value) is not int or value < 0 for value in
               (self.max_tensor_entries, self.max_matrix_entries,
                self.max_profile_cases, self.max_slice_cases)):
            raise ValueError("screen caps must be nonnegative integers")


DEFAULT_BUDGET = MixedDotScreenBudget()
KINDS = ("zero", "unit", "matrix")


def _budget(budget):
    if not isinstance(budget, MixedDotScreenBudget):
        raise ValueError("supply an explicit mixed-dot screen budget")


def _tensor_cap(shape, budget):
    if prod(shape) > budget.max_tensor_entries or sum(shape) > budget.max_tensor_entries:
        raise MixedDotScreenLimit("dense tensor or coordinate cap reached; controls incomplete")


def _kind(kind):
    if kind not in KINDS:
        raise ValueError("catalyst kind must be zero, unit or matrix")


def mixed_dot_tensor(field, *, budget=DEFAULT_BUDGET):
    """Actual (5,5,5) direct sum, in singleton-leg order first/second/third."""
    if not isinstance(field, FiniteField):
        raise ValueError("supply an explicitly represented finite field")
    _budget(budget)
    _tensor_cap((5, 5, 5), budget)
    dot = FiniteTensor.from_function(field, (1, 2, 2), lambda i, j, k: int(j == k))
    return FiniteTensor.direct_sum((dot, dot.permute_axes((1, 0, 2)), dot.permute_axes((1, 2, 0))))


def exact_slice_rank(tensor, vector, *, budget=DEFAULT_BUDGET):
    """Actual first-leg contraction and exact rank, not a generic-rank guess."""
    if not isinstance(tensor, FiniteTensor):
        raise ValueError("slice source must be an exact finite tensor")
    _budget(budget)
    _tensor_cap(tensor.shape, budget)
    if budget.max_slice_cases < 1:
        raise MixedDotScreenLimit("slice-control cap reached; controls incomplete")
    x, y, z = tensor.shape
    if y*z > budget.max_matrix_entries:
        raise MixedDotScreenLimit("slice matrix cap reached; controls incomplete")
    vector = tuple(vector)
    if len(vector) != x:
        raise ValueError("slice vector has the wrong first coordinate size")
    for value in vector:
        tensor.field.check(value)
    support = tuple((i, value) for i, value in enumerate(vector) if value)
    f = tensor.field
    rows = tuple(tuple(f.sum(f.mul(value, coefficient) for i, value in support
                             if (coefficient := tensor.coefficient(i, j, k)))
                       for k in range(z)) for j in range(y))
    return field_matrix_rank(f, rows, z)


def mixed_dot_core_input(field, A, B, C):
    """Encode A(2x2), B(4x2), C(2x4) in M2 tensor S's first coordinates.

    The first coordinate is matrix_index*5 + dot_coordinate. B's row is
    (matrix_row,dot_coordinate), while C's column is (matrix_inner,dot_coordinate).
    Their actual slice rank is 4 rank(A)+2 rank(B)+2 rank(C).
    """
    if not isinstance(field, FiniteField):
        raise ValueError("supply an explicitly represented finite field")
    matrices = tuple(tuple(tuple(row) for row in matrix) for matrix in (A, B, C))
    for matrix, shape in zip(matrices, ((2, 2), (4, 2), (2, 4))):
        if len(matrix) != shape[0] or any(len(row) != shape[1] for row in matrix):
            raise ValueError("core matrices must have shapes (2,2), (4,2), (2,4)")
        for row in matrix:
            for value in row:
                field.check(value)
    A, B, C = matrices
    vector = [0]*20
    for i, j in product(range(2), repeat=2):
        coordinate = (2*i+j)*5
        vector[coordinate] = A[i][j]
        for r in range(2):
            vector[coordinate+1+r] = B[2*i+r][j]
            vector[coordinate+3+r] = C[i][2*j+r]
    return tuple(vector)


def weighted_linear_dimension_maximum(kind, *, rank_bound=None, budget=DEFAULT_BUDGET):
    """Enumerate analytic generic-rank profiles, not field-valued samples.

    Maximal linear dimensions for projection ranks 0/1/2 are f=(0,2,4)
    for 2x2 and t=(0,4,8) for 4x2 or 2x4 matrices. These bounds use the
    all-field common-factor classification from the note. The returned
    maximizers describe this weighted upper estimate, not sampled spaces.
    """
    _kind(kind)
    _budget(budget)
    default_bound = {"zero": 10, "unit": 11, "matrix": 14}[kind]
    rank_bound = default_bound if rank_bound is None else rank_bound
    if type(rank_bound) is not int or rank_bound < 0:
        raise ValueError("slice rank bound must be a nonnegative integer")
    domains = (range(3), range(3), range(3), range(2))
    if kind == "unit":
        domains += (range(2),)
    elif kind == "matrix":
        domains += (range(3),)
    count = prod(len(domain) for domain in domains)
    if count > budget.max_profile_cases:
        raise MixedDotScreenLimit("rank-profile cap reached; controls incomplete")
    f, t, maximum, maximizers = (0, 2, 4), (0, 4, 8), -1, []
    for ranks in product(*domains):
        a, b, c, g = ranks[:4]
        rank = 4*a+2*b+2*c+g
        dimension = f[a]+t[b]+t[c]+g
        if kind == "unit":
            rank += ranks[4]
            dimension += ranks[4]
        elif kind == "matrix":
            rank += 2*ranks[4]
            dimension += f[ranks[4]]
        if rank <= rank_bound:
            if dimension > maximum:
                maximum, maximizers = dimension, [ranks]
            elif dimension == maximum:
                maximizers.append(ranks)
    return {"rank_bound": rank_bound, "maximum_dimension": maximum,
            "maximizing_profiles": tuple(maximizers), "checked_profiles": count,
            "profile_order": ("A", "B", "C", "gain") +
                             (() if kind == "zero" else ("D_scalar" if kind == "unit" else "D_matrix",))}


def weighted_diagonal_dimension_maximum(n, *, budget=DEFAULT_BUDGET):
    """The same upper estimate for n diagonal catalyst scalars and one gain.

    For each of 27 matrix-rank profiles, choose the largest admissible
    scalar projection count analytically. This avoids enumerating 2**(n+1)
    scalar profiles; it is a combinatorial bound, not a finite-field sample.
    """
    if type(n) is not int or n < 0:
        raise ValueError("diagonal size must be a nonnegative integer")
    _budget(budget)
    if budget.max_profile_cases < 27:
        raise MixedDotScreenLimit("rank-profile cap reached; controls incomplete")
    f, t, maximum, maximizers = (0, 2, 4), (0, 4, 8), -1, []
    for a, b, c in product(range(3), repeat=3):
        scalars = min(n+1, n+10-4*a-2*b-2*c)
        if scalars < 0:
            continue
        dimension = f[a]+t[b]+t[c]+scalars
        if dimension > maximum:
            maximum, maximizers = dimension, [(a, b, c, scalars)]
        elif dimension == maximum:
            maximizers.append((a, b, c, scalars))
    return {"rank_bound": n+10, "maximum_dimension": maximum,
            "maximizing_profiles": tuple(maximizers), "checked_profiles": 27,
            "profile_order": ("A", "B", "C", "active_scalar_projections")}


@dataclass(frozen=True)
class MixedDotCoefficientControl:
    auxiliary: FiniteTensor
    catalyst: FiniteTensor
    source: FiniteTensor
    target: FiniteTensor
    auxiliary_scheme: TensorScheme
    target_scheme: TensorScheme
    source_vectors: tuple[tuple[int, ...], ...]
    source_slice_ranks: tuple[int, ...]
    target_vector: tuple[int, ...]
    target_slice_rank: int
    unit_exception_coordinates: tuple[int, ...] = ()
    unit_exception_vector: tuple[int, ...] = ()
    unit_exception_rank: int | None = None


def mixed_dot_coefficient_control(field, kind, *, diagonal_size=None, budget=DEFAULT_BUDGET):
    """Build exact one-gain source/target controls for a supported small D."""
    if kind == "diagonal":
        if type(diagonal_size) is not int or diagonal_size < 2:
            raise ValueError("diagonal controls require an integer size at least two")
    else:
        _kind(kind)
        if diagonal_size is not None:
            raise ValueError("diagonal_size is only supported for diagonal catalyst controls")
    _budget(budget)
    if not isinstance(field, FiniteField):
        raise ValueError("supply an explicitly represented finite field")
    size = diagonal_size if kind == "diagonal" else {"zero": 0, "unit": 1, "matrix": 4}[kind]
    for shape in ((5,)*3, (4,)*3, (20,)*3, (size+25,)*3, (size+21,)*3):
        _tensor_cap(shape, budget)
    if (size+25)**2 > budget.max_matrix_entries:
        raise MixedDotScreenLimit("control matrix cap reached; controls incomplete")
    target_terms = 43+(size if kind == "diagonal" else {"zero": 0, "unit": 1, "matrix": 7}[kind])
    if 3*target_terms*(size+21) > budget.max_matrix_entries:
        raise MixedDotScreenLimit("scheme coefficient matrix cap reached; controls incomplete")
    if budget.max_slice_cases < (5 if kind == "unit" else 4):
        raise MixedDotScreenLimit("slice-control cap reached; controls incomplete")
    S, unit = mixed_dot_tensor(field, budget=budget), FiniteTensor.unit(field)
    matrix_scheme = strassen_scheme(field).as_tensor_scheme()
    matrix = matrix_scheme.target
    scalar_scheme = TensorScheme.from_tensor(unit)
    s_scheme = TensorScheme.from_tensor(S)
    s_scheme.require_exact()
    core = matrix.tensor_product(S)
    core_scheme = matrix_scheme.tensor_product(s_scheme)
    if kind == "zero":
        D, d_vector, d_schemes = FiniteTensor.zero(field, (0, 0, 0)), (), ()
    elif kind == "unit":
        D, d_vector, d_schemes = unit, (1,), (scalar_scheme,)
    elif kind == "matrix":
        D, d_vector, d_schemes = matrix, (1, 0, 0, 1), (matrix_scheme,)
    else:
        D = FiniteTensor.from_function(field, (size,)*3, lambda i, j, k: int(i == j == k))
        d_vector, d_schemes = (1,)*size, (TensorScheme.from_tensor(D),)
    source = FiniteTensor.direct_sum(((D,) if size else ())+(S,)*5)
    target = FiniteTensor.direct_sum(((D,) if size else ())+(unit, core))
    target_scheme = TensorScheme.direct_sum(d_schemes+(scalar_scheme, core_scheme))
    target_scheme.require_exact()
    if s_scheme.terms != 6 or target_scheme.target != target or target_scheme.terms != target_terms:
        raise ValueError("supplied mixed-dot rank upper controls have the wrong target or term count")
    source_vectors = tuple(d_vector+tuple(value for copy in range(5)
                                         for value in ((int(copy >= killed), 1, 1, 1, 1)))
                           for killed in (0, 2, 5))
    source_ranks = tuple(exact_slice_rank(source, vector, budget=budget) for vector in source_vectors)
    A = ((1, 0), (0, 1))
    B = ((1, 0), (0, 1), (0, 0), (0, 0))
    C = ((1, 0, 0, 0), (0, 1, 0, 0))
    core_vector = mixed_dot_core_input(field, A, B, C)
    target_vector = d_vector+(1,)+core_vector
    target_rank = exact_slice_rank(target, target_vector, budget=budget)
    extra = size if kind == "diagonal" else {"zero": 0, "unit": 1, "matrix": 4}[kind]
    expected_target = 17+extra
    if source_ranks != (20+extra, 16+extra, 10+extra) or target_rank != expected_target:
        raise ValueError("actual coefficient slices disagree with the mixed-dot rank formulas")
    coordinates, exception_vector, exception_rank = (), (), None
    if kind == "unit":
        # Pure common-image A plane plus all eight B and eight C coordinates;
        # both scalar coordinates d,g are zero. This is an actual dimension-18
        # first space with a witnessed slice of rank 12, exceeding eleven.
        coordinates = tuple(2+i for i in range(20) if i not in (10, 15))
        exception_vector = (0, 0)+mixed_dot_core_input(field, ((1, 0), (0, 0)), B, C)
        exception_rank = exact_slice_rank(target, exception_vector, budget=budget)
        if exception_rank != 12 or any(exception_vector[i] for i in range(22) if i not in coordinates):
            raise ValueError("unit exceptional space control failed")
    return MixedDotCoefficientControl(S, D, source, target, s_scheme, target_scheme,
                                      source_vectors, source_ranks, target_vector, target_rank,
                                      coordinates, exception_vector, exception_rank)


def _canonical_diagonal_size(tensor):
    n = tensor.shape[0]
    if tensor.shape != (n,)*3:
        return None
    if n == 0:
        return 0
    for index, value in enumerate(tensor.coefficients):
        ij, k = divmod(index, n)
        i, j = divmod(ij, n)
        if value != int(i == j == k):
            return None
    return n


def _canonical_same_orientation_dots(tensor):
    """Recognize the literal direct sum; no basis-change inference."""
    for axis in range(3):
        n = tensor.shape[axis]
        others = tuple(i for i in range(3) if i != axis)
        if not n or any(tensor.shape[i] != 2*n for i in others):
            continue
        y, z = tensor.shape[1:]
        for index, value in enumerate(tensor.coefficients):
            ij, k = divmod(index, z)
            i, j = divmod(ij, y)
            coordinates = (i, j, k)
            a, b = (coordinates[leg] for leg in others)
            if value != int(a == b and a//2 == coordinates[axis]):
                break
        else:
            return n, axis
    return None


def _canonical_mixed_orientation_dots(tensor):
    total, remainder = divmod(sum(tensor.shape), 5)
    if remainder or total == 0:
        return None
    counts = tuple(2*total-size for size in tensor.shape)
    if any(n < 0 for n in counts) or sum(counts) != total:
        return None
    support, offsets = set(), [0, 0, 0]
    y, z = tensor.shape[1:]
    for axis, count in enumerate(counts):
        for _ in range(count):
            for coordinate in range(2):
                triple = tuple(offsets[leg]+(0 if leg == axis else coordinate) for leg in range(3))
                support.add((triple[0]*y+triple[1])*z+triple[2])
            offsets = [offsets[leg]+(1 if leg == axis else 2) for leg in range(3)]
    if any(value != int(index in support) for index, value in enumerate(tensor.coefficients)):
        return None
    return counts


def _dot_catalyst_controls(tensor, axis, budget):
    from .mixed_dot_lengths import mixed_dot_length_coefficient_control
    from .witnessed_constraints import StateAssemblyLimit, StateAssemblyLimits
    if budget.max_slice_cases < 4:
        raise MixedDotScreenLimit("four slice controls exceed the declared cap")
    axes = (axis,)+tuple(i for i in range(3) if i != axis)
    limits = StateAssemblyLimits(max_tensor_entries=budget.max_tensor_entries,
                                max_map_entries=budget.max_matrix_entries,
                                max_constraint_entries=budget.max_tensor_entries)
    try:
        control = mixed_dot_length_coefficient_control(
            tensor.field, 2, catalyst=tensor.permute_axes(axes),
            catalyst_input=(1,)*tensor.shape[axis], limits=limits)
        control.require_valid(limits=limits)
    except StateAssemblyLimit as error:
        raise MixedDotScreenLimit(str(error)) from error
    return {"source_slice_ranks_full_two_zero_all_zero": control.source_slice_ranks,
            "target_full_rank_slice": control.target_slice_rank,
            "sampling_scope": "explicit coefficient witnesses only; generic geometry is analytic"}


def _all_catalyst_rank_control(tensor, gain, budget):
    from .mixed_dot_rank import mixed_dot_rank_control
    from .witnessed_constraints import StateAssemblyLimit, StateAssemblyLimits
    if budget.max_profile_cases < 1 or budget.max_slice_cases < 1:
        raise MixedDotScreenLimit("rank-control bookkeeping cap reached")
    limits = StateAssemblyLimits(max_tensor_entries=budget.max_tensor_entries,
                                max_map_entries=budget.max_matrix_entries,
                                max_constraint_entries=budget.max_tensor_entries)
    try:
        control = mixed_dot_rank_control(tensor.field, limits=limits)
    except StateAssemblyLimit as error:
        raise MixedDotScreenLimit(str(error)) from error
    upper = sum(value != 0 for value in tensor.coefficients)
    copies = upper//gain+1
    return {"core_first_flattening_rank": control.first_flattening_rank,
            "rectangle_shape": control.rectangle.shape,
            "analytic_core_rank_lower_per_copy": control.analytic_core_rank_lower,
            "catalyst_coordinate_rank_upper": upper,
            "contradicting_repetitions": copies,
            "target_rank_lower": (30+gain)*copies,
            "source_rank_upper": upper+30*copies,
            "iteration_materialized": False,
            "scope": "replicated rectangular kernel-map proof and substitution; finite controls are not tensor-rank samples"}


def screen_mixed_dot_catalyst(tensor, *, auxiliary=None, m=1, budget=DEFAULT_BUDGET):
    """Analytic screen for canonical S and every finite catalyst.

    analytically_excluded cites the note's analytic proof. Coefficient
    controls validate its formulas; their finite samples do not prove that
    proof. Unsupported S remains unclassified and supplies no existence
    evidence. Positive gains reduce to one by an output restriction.
    """
    if not isinstance(tensor, FiniteTensor) or not isinstance(tensor.field, FiniteField):
        raise ValueError("supply an exact finite catalyst tensor")
    _budget(budget)
    if type(m) is not int or m < 1:
        raise ValueError("gain must be a positive integer")
    if auxiliary is not None and (not isinstance(auxiliary, FiniteTensor) or auxiliary.field != tensor.field):
        raise ValueError("auxiliary must be an exact tensor over the catalyst field")
    report = {"status": "unclassified", "reasons": [],
              "scope": "D + m Unit + M2 tensor S <= D + 5 S; m>=1; exact three-oriented dot_2 S",
              "field": {"p": tensor.field.p, "modulus": tensor.field.modulus},
              "verification": "analytic arguments in research/variable-auxiliaries.md and research/mixed-dot-rank-catalysts.md; exact coefficient controls; not Lean",
              "finite_sampling_proves_all_fields": False, "gain": m}
    try:
        _tensor_cap(tensor.shape, budget)
        if auxiliary is not None:
            _tensor_cap(auxiliary.shape, budget)
        expected_s = mixed_dot_tensor(tensor.field, budget=budget)
        if auxiliary is not None and auxiliary != expected_s:
            report["reasons"].append("auxiliary is outside the exact supported coefficient tensor")
            return report
        report.update(arbitrary_finite_catalyst=True,
                      all_catalyst_rank_control=_all_catalyst_rank_control(tensor, m, budget))
        dots = _canonical_same_orientation_dots(tensor)
        if dots is not None:
            n, axis = dots
            controls = _dot_catalyst_controls(tensor, axis, budget)
            report.update(status="analytically_excluded", catalyst_kind="same_orientation_dots",
                          dot_copies=n, singleton_leg=axis, coefficient_controls=controls,
                          reasons=["matrix direct-sum rank additivity and the core rank lower bound exclude this literal catalyst"],
                          rank_obstruction={"catalyst_rank": 2*n,
                              "source_rank": 2*n+30, "target_rank_lower": 2*n+30+m,
                              "scope": "analytic all-field rank bound; not an ordinary state row"})
            return report
        counts = _canonical_mixed_orientation_dots(tensor)
        if counts is not None:
            axis = counts.index(min(counts))
            controls = _dot_catalyst_controls(tensor, axis, budget)
            rank_d = 2*sum(counts)
            report.update(status="analytically_excluded", catalyst_kind="mixed_orientation_dots",
                dot_counts=counts, coefficient_controls=controls,
                reasons=["matrix direct-sum rank additivity and the core rank lower bound exclude all orientation counts"],
                chosen_first_orientation=axis,
                rank_obstruction={"catalyst_rank": rank_d, "source_rank": rank_d+30,
                    "target_rank_lower": rank_d+30+m,
                    "scope": "analytic all-field rank bound; not an ordinary state row"})
            return report
        diagonal_size = _canonical_diagonal_size(tensor)
        if diagonal_size == 0:
            kind = "zero"
        elif diagonal_size == 1:
            kind = "unit"
        elif diagonal_size is not None:
            kind = "diagonal"
        elif tensor.shape == (4,)*3 and tensor == matrix_multiplication_tensor(tensor.field, 2, 2, 2):
            kind = "matrix"
        else:
            report.update(status="analytically_excluded", catalyst_kind="arbitrary_finite",
                reasons=["replicated rectangular rank and substitution exclude this auxiliary for every finite catalyst"])
            return report
        dimensions = (weighted_diagonal_dimension_maximum(diagonal_size, budget=budget)
                      if kind == "diagonal" else weighted_linear_dimension_maximum(kind, budget=budget))
        control = mixed_dot_coefficient_control(tensor.field, kind,
                    diagonal_size=diagonal_size if kind == "diagonal" else None, budget=budget)
    except MixedDotScreenLimit as error:
        report.update(status="resource_cap", reasons=[str(error)])
        return report
    span_bound = 5 if kind == "diagonal" else {"zero": 3, "unit": 4, "matrix": 3}[kind]
    report.update(status="analytically_excluded", catalyst_kind=kind,
                  reasons=["mixed-dot slice-space argument excludes this exact supported catalyst"],
                  target_first_dimension=control.target.shape[0],
                  scalar_form_span_bound=span_bound,
                  common_kernel_dimension_lower=control.target.shape[0]-span_bound,
                  dimension_control=dimensions,
                  coefficient_controls={"source_slice_ranks_full_two_zero_all_zero": control.source_slice_ranks,
                                        "target_full_rank_slice": control.target_slice_rank,
                                        "auxiliary_upper_terms": control.auxiliary_scheme.terms,
                                        "target_upper_terms": control.target_scheme.terms,
                                        "sampling_scope": "explicit coefficient witnesses only; generic geometry is analytic"})
    if kind == "unit":
        report["unit_exception"] = {"dimension": len(control.unit_exception_coordinates),
                                    "witnessed_slice_rank": control.unit_exception_rank,
                                    "source_slice_bound": 11}
    if kind == "diagonal":
        report["diagonal_size"] = diagonal_size
        report["diagonal_argument"] = {
            "no_pure_A_case": "core generic rank sixteen plus the surviving scalar axes exceeds n+10",
            "pure_A_case": "the remaining n-2 scalar axes cannot contain n distinct source-D scalar forms",
            "scope": "analytic all-n argument; coefficient controls do not establish the geometry"}
    return report
