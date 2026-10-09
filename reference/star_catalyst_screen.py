"""Exact computational consequences of audited singular-star obstructions.

Scope: D + m Unit + M2 tensor S_star <= D + 5 S_star, m>0.
The mathematical implications are documented in research/star-*.md and
research/variable-auxiliaries.md (the stronger three-plane refinement), and have
not been formalized in Lean. Passing this screen supplies no catalyst maps.
Finite-point slice ranks are explicitly not generic slice ranks. The separate
universal screen uses the integral certificate in research/star-koszul-catalysts.md.
"""

from dataclasses import dataclass
from itertools import product
from math import prod

from .catalyst_search_experiment import field_matrix_rank
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor


@dataclass(frozen=True)
class StarScreenBudget:
    max_tensor_entries: int = 100_000
    max_slices: int = 4096

    def __post_init__(self):
        if any(type(v) is not int or v < 0 for v in (self.max_tensor_entries, self.max_slices)):
            raise ValueError('screen caps must be nonnegative integers')


def _flatten_rows(tensor, axis):
    others = tuple(i for i in range(3) if i != axis)
    rows = []
    for i in range(tensor.shape[axis]):
        row = []
        for j, k in product(range(tensor.shape[others[0]]), range(tensor.shape[others[1]])):
            coordinate = [0, 0, 0]
            coordinate[axis], coordinate[others[0]], coordinate[others[1]] = i, j, k
            row.append(tensor.coefficient(*coordinate))
        rows.append(tuple(row))
    return tuple(rows), tensor.shape[others[0]]*tensor.shape[others[1]]


def screen_star_catalyst_geometry(tensor, *, supplied_scheme=None, budget=StarScreenBudget()):
    """Earlier necessary geometry controls, retained independently for diagnostics."""
    if not isinstance(tensor, FiniteTensor) or not isinstance(budget, StarScreenBudget):
        raise ValueError('supply an exact finite tensor and screen budget')
    report = {'status': 'not_excluded', 'reasons': [],
              'scope': 'fixed S_star,d=2,k=5,positive gain; restriction coefficients over this represented field',
              'field': {'p': tensor.field.p, 'modulus': tensor.field.modulus},
              'verification': 'exact arithmetic under independently audited, unformalized mathematical lemmas'}
    if prod(tensor.shape) > budget.max_tensor_entries or sum(tensor.shape) > budget.max_tensor_entries:
        report.update(status='resource_cap', reasons=['input tensor exceeds screen budget'])
        return report
    if not all(tensor.shape):
        report.update(status='excluded', reasons=['zero-dimensional catalyst'], flattening_ranks=(0, 0, 0))
        return report
    rows = [_flatten_rows(tensor, axis) for axis in range(3)]
    ranks = tuple(field_matrix_rank(tensor.field, matrix, width) for matrix, width in rows)
    h, b, c = ranks
    report['flattening_ranks'] = ranks
    reasons = report['reasons']
    if h < 6:
        reasons.append('slice-space refinement requires first flattening rank at least six')
    if min(b, c) < 6:
        reasons.append('three-plane refinement requires other flattening ranks at least six')
    if min(b, c) <= 7 and max(b, c) < 10:
        reasons.append('generic slice rank at most seven requires another flattening rank at least ten')
    q = tensor.field.order
    # Count rank-one first-factor forms on the actual coefficient-field
    # three-plane. This is a tensor-rank bound, never an LP state lower row.
    required_rank = max(6,
                        (4*(q*q+q+1)+q*q-1)//(q*q))
    report['required_same_field_tensor_rank'] = required_rank
    if supplied_scheme is not None:
        if not isinstance(supplied_scheme, TensorScheme) or supplied_scheme.target != tensor:
            raise ValueError('supplied rank upper certificate must have this exact target')
        supplied_scheme.require_exact()
        report['supplied_upper_terms'] = supplied_scheme.terms
        if supplied_scheme.terms < required_rank:
            reasons.append('verified same-field rank upper is below the necessary rank')
    if reasons:
        report['status'] = 'excluded'
        return report

    # Use a basis of the actual slice image, rather than counting nonzero
    # vectors in an input space with a flattening kernel.
    basis = []
    for row in rows[0][0]:
        if field_matrix_rank(tensor.field, basis+[row], rows[0][1]) > len(basis):
            basis.append(row)
    tested, maximum = 0, 0
    histogram = {}
    total = tensor.field.order**h-1
    for vector in product(range(tensor.field.order), repeat=h):
        if not any(vector):
            continue
        if tested >= budget.max_slices:
            break
        coefficients = tuple(tensor.field.sum(tensor.field.mul(v, row[index]) for v, row in zip(vector, basis))
                             for index in range(rows[0][1]))
        matrix = tuple(coefficients[j*tensor.shape[2]:(j+1)*tensor.shape[2]] for j in range(tensor.shape[1]))
        rank = field_matrix_rank(tensor.field, matrix, tensor.shape[2])
        tested += 1
        maximum = max(maximum, rank)
        histogram[rank] = histogram.get(rank, 0)+1
    report['base_field_slices'] = {'checked': tested, 'total_nonzero': total,
                                  'complete': tested == total, 'histogram': histogram,
                                  'maximum_observed': maximum, 'generic_rank': 'not_computed'}
    if tested == total and maximum < 4:
        reasons.append('no base-field slice has the necessary matrix rank four')
        report['status'] = 'excluded'
    return report


def screen_star_catalyst(tensor, *, supplied_scheme=None, budget=StarScreenBudget(),
                        m=1, koszul_limits=None):
    """All-field analytic exclusion for every finite catalyst of the literal star.

    Integer components and the complete source polynomial are replayed. The
    unbounded rational-function iteration is a mathematical proof, not a
    finite-field enumeration or a new Lean theorem.
    """
    from .star_degeneration import star_polynomial_degeneration
    from .star_koszul import StarKoszulLimit, StarKoszulLimits, verify_star_integer_components
    if not isinstance(tensor, FiniteTensor) or not isinstance(budget, StarScreenBudget):
        raise ValueError('supply an exact finite tensor and screen budget')
    if type(m) is not int or m < 1:
        raise ValueError('gain must be a positive integer')
    limits = StarKoszulLimits() if koszul_limits is None else koszul_limits
    if not isinstance(limits, StarKoszulLimits):
        raise ValueError('supply explicit star Koszul limits')
    report = {'status': 'not_excluded', 'reasons': [], 'gain': m,
        'scope': 'fixed literal S_star,d=2,k=5,positive gain; every finite catalyst and every field',
        'field': {'p': tensor.field.p, 'modulus': tensor.field.modulus},
        'verification': 'integral certificate and analytic iteration in research/star-koszul-catalysts.md; not Lean',
        'finite_sampling_proves_all_fields': False}
    if prod(tensor.shape) > budget.max_tensor_entries or sum(tensor.shape) > budget.max_tensor_entries:
        report.update(status='resource_cap', reasons=['input tensor exceeds screen budget'])
        return report
    if supplied_scheme is not None:
        if not isinstance(supplied_scheme, TensorScheme) or supplied_scheme.target != tensor:
            raise ValueError('supplied rank upper certificate must have this exact target')
        supplied_scheme.require_exact()
        report['supplied_upper_terms'] = supplied_scheme.terms
    try:
        control = verify_star_integer_components(limits=limits)
    except StarKoszulLimit as error:
        report.update(status='resource_cap', reasons=[str(error)])
        return report
    polynomial = star_polynomial_degeneration(tensor.field)
    polynomial.require_valid()
    upper = sum(value != 0 for value in tensor.coefficients)
    gap = control.rank_lower-15*control.rank_one_factor
    copies = control.rank_one_factor*upper//gap+1
    report.update(status='excluded', arbitrary_finite_catalyst=True,
        reasons=['integral Koszul components and rational-function iteration exclude every finite catalyst'],
        coefficient_controls={'formal_source_terms': polynomial.terms,
            'formal_leading_degree': polynomial.leading,
            'integral_matrix_rank_lower': control.rank_lower,
            'rank_one_factor': control.rank_one_factor,
            'target_border_rank_lower': control.border_rank_lower},
        iteration={'catalyst_coordinate_rank_upper': upper, 'contradicting_repetitions': copies,
            'target_matrix_rank_lower': control.rank_lower*copies,
            'source_matrix_rank_upper': control.rank_one_factor*(upper+15*copies),
            'materialized': False})
    return report
