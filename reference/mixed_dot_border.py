"""Exact finite controls for the universal mixed-dot border-rank-six proof.

The adjugate polynomial identity for rank at most five is derived in
research/mixed-dot-border-rank.md. Finite tests verify its coefficient
witness; they do not establish the universal closure argument.
"""

from itertools import permutations

from .fields import FiniteField
from .maps import LinearMap
from .mixed_dot_catalyst_screen import MixedDotScreenBudget, MixedDotScreenLimit, mixed_dot_tensor
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor


def adjugate_slice_relation(tensor, vectors, *, budget=MixedDotScreenBudget()):
    """Return A_x adj(A_y) A_z - A_z adj(A_y) A_x, with second-leg rows."""
    if not isinstance(tensor, FiniteTensor) or tensor.shape[1:] != (5, 5):
        raise ValueError("supply an actual tensor with second and third dimensions five")
    if not isinstance(budget, MixedDotScreenBudget):
        raise ValueError("supply an explicit mixed-dot budget")
    if len(tensor.coefficients) > budget.max_tensor_entries or sum(tensor.shape) > budget.max_tensor_entries:
        raise MixedDotScreenLimit("adjugate tensor cap reached")
    if budget.max_matrix_entries < 25 or budget.max_slice_cases < 3:
        raise MixedDotScreenLimit("adjugate matrix or slice cap reached")
    vectors = tuple(tuple(vector) for vector in vectors)
    if len(vectors) != 3 or any(len(vector) != tensor.shape[0] for vector in vectors):
        raise ValueError("supply three first-slice vectors of the actual first dimension")
    field = tensor.field
    for vector in vectors:
        for value in vector:
            field.check(value)
    slices = tuple(LinearMap(field, 5, tuple(tuple(field.sum(field.mul(value, tensor.coefficient(i, j, k))
        for i, value in enumerate(vector) if value) for k in range(5)) for j in range(5))) for vector in vectors)

    def determinant4(rows):
        values = []
        for order in permutations(range(4)):
            term = field.embed(-1 if sum(order[i] > order[j] for i in range(4) for j in range(i+1, 4)) % 2 else 1)
            for i, j in enumerate(order):
                term = field.mul(term, rows[i][j])
            values.append(term)
        return field.sum(values)

    middle = slices[1].rows
    adjugate = LinearMap(field, 5, tuple(tuple(field.mul(field.embed((-1)**(i+j)),
        determinant4(tuple(tuple(middle[r][c] for c in range(5) if c != i)
                           for r in range(5) if r != j))) for j in range(5)) for i in range(5)))
    first = slices[0].compose(adjugate).compose(slices[2])
    second = slices[2].compose(adjugate).compose(slices[0])
    bracket = LinearMap(field, 5, tuple(tuple(field.sub(a, b) for a, b in zip(left, right))
        for left, right in zip(first.rows, second.rows)))
    return slices, adjugate, bracket


def mixed_dot_border_control(field, *, budget=MixedDotScreenBudget()):
    """Actual obstruction entry -1 and supplied six-term coefficient upper."""
    if not isinstance(field, FiniteField):
        raise ValueError("supply a represented finite field")
    tensor = mixed_dot_tensor(field, budget=budget)
    vectors = ((0, 1, 0, 0, 0), (1, 0, 1, 0, 1), (0, 0, 0, 1, 0))
    slices, adjugate, bracket = adjugate_slice_relation(tensor, vectors, budget=budget)
    expected_adjugate = tuple(tuple(field.neg(1) if (i, j) == (2, 3) else 0 for j in range(5)) for i in range(5))
    expected_bracket = tuple(tuple(field.neg(1) if (i, j) == (2, 4) else 0 for j in range(5)) for i in range(5))
    if adjugate.rows != expected_adjugate or bracket.rows != expected_bracket:
        raise ValueError("actual mixed-dot adjugate coefficient witness differs from the universal proof")
    scheme = TensorScheme.from_tensor(tensor)
    scheme.require_exact()
    if scheme.terms != 6:
        raise ValueError("actual mixed-dot upper control has the wrong term count")
    return {"status": "checked", "exact_rank_upper_terms": 6, "border_rank": 6,
        "obstruction_coordinate": (2, 4), "obstruction_coefficient": bracket.rows[2][4],
        "slices": slices, "adjugate": adjugate, "bracket": bracket,
        "verification": "exact finite coefficients plus analytic adjugate closure proof; not Lean"}
