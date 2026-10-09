"""Exact controls for the complete-convolution output/commutator obstruction.

Run: python3 -m reference.research.convolution_commutator_experiment

The independently audited argument is elementary, but not a Lean theorem.
For an output-concise tensor in X tensor Y tensor W with dim(W)=Z,
suppose one slice identifies Y with V <= W, dim(Y)=dim(V)=N, and two
other slices have image in V with normalized commutator rank c. Then
tensor rank >= Z + ceil(c/2): in any rank-one decomposition, select Z-N
output factors independent modulo V and quotient by their span K. This
removes Z-N terms, while preserving the three slices under W/K ~= V.
The remaining square tensor has rank >= N+ceil(c/2).

The square bound follows by writing the normalized slices as U D V and
U E V, where UV=I and D,E are diagonal. Their commutator is a sum of
two products through I-VU, of rank r-N. If some decomposition coefficients
at the chosen identity slice vanish, perturb that slice over an infinite
extension; its inverse is regular at the chosen slice, so the commutator
minor bound specializes there. This addresses zero weights and finite fields.

For q copies of m Unit + M2 tensor C(a,b), these controls check Z=q(m+4L),
N=q(m+4b), c=4qb, L=a+b-1. Swapping input legs replaces b by a. Thus
rank >= q(m+4L+2max(a,b)). A hypothetical finite catalyst for k=5 can be
reused q times and then projected away, giving the upper bound
rank(D)+5qL over an infinite closure. Since m+2max(a,b)-L>0, this is
impossible for unbounded q. No rank direct-sum additivity, rank cancellation,
or external rectangular-rank theorem is used. This program checks the finite
linear algebra and supplied rank upper controls, not the unbounded argument.
"""

from dataclasses import dataclass
from itertools import product

from reference.fields import FiniteField
from reference.maps import LinearMap
from reference.sectors import convolution_tensor
from reference.tensor_schemes import TensorScheme
from reference.tensors import FiniteTensor, matrix_multiplication_tensor
from reference.witnessed_constraints import (DEFAULT_LIMITS, StateAssemblyLimit,
                                            StateAssemblyLimits)


def _rank(field, rows, columns):
    """Exact sparse Gaussian elimination, including extension encodings."""
    pivots = {}
    for row in rows:
        if len(row) != columns:
            raise ValueError("rank matrix rows have inconsistent lengths")
        values = {i: field.check(value) for i, value in enumerate(row) if value}
        while values:
            column = min(values)
            if column not in pivots:
                inverse = field.inv(values[column])
                pivots[column] = {i: field.mul(value, inverse) for i, value in values.items()}
                break
            scalar = values[column]
            for i, value in pivots[column].items():
                updated = field.sub(values.get(i, 0), field.mul(scalar, value))
                if updated:
                    values[i] = updated
                else:
                    values.pop(i, None)
    return len(pivots)


def output_flattening_rank(tensor, *, limits=DEFAULT_LIMITS):
    """Compute the actual output flattening rank from all coefficients."""
    if not isinstance(tensor, FiniteTensor) or not isinstance(limits, StateAssemblyLimits):
        raise ValueError("supply an exact finite tensor and state assembly limits")
    limits.tensor(tensor.shape)
    x, y, z = tensor.shape
    rows = (tuple(tensor.coefficient(i, j, k) for i in range(x) for j in range(y))
            for k in range(z))
    return _rank(tensor.field, rows, x*y)


def slice_map(tensor, vector, *, limits=DEFAULT_LIMITS):
    """Contraction of first leg, as [output coordinate][input column]."""
    if not isinstance(tensor, FiniteTensor) or not isinstance(limits, StateAssemblyLimits):
        raise ValueError("supply an exact finite tensor and state assembly limits")
    limits.tensor(tensor.shape)
    limits.maps((tensor.shape[1],)*3, (tensor.shape[2],)*3)
    if sum(tensor.shape) > limits.max_constraint_entries:
        raise StateAssemblyLimit("slice coordinate cap reached; construction incomplete")
    vector = tuple(vector)
    if len(vector) != tensor.shape[0]:
        raise ValueError("slice input vector has the wrong length")
    for value in vector:
        tensor.field.check(value)
    f, (_, y, z) = tensor.field, tensor.shape
    support = tuple((i, value) for i, value in enumerate(vector) if value)
    rows = tuple(tuple(f.sum(f.mul(value, coefficient)
                             for i, value in support
                             if (coefficient := tensor.coefficient(i, j, k)))
                       for j in range(y)) for k in range(z))
    return LinearMap(f, y, rows)


def _preflight(a, b, m, q, limits):
    if not isinstance(limits, StateAssemblyLimits):
        raise ValueError("supply explicit state assembly limits")
    if q*(m+1) > limits.max_blocks:
        raise StateAssemblyLimit("commutator block cap reached; construction incomplete")
    L = a+b-1
    shape = tuple(q*(m+4*size) for size in (a, b, L))
    for candidate in ((1, 1, 1), (4, 4, 4), (a, b, L), (4*a, 4*b, 4*L),
                      (m+4*a, m+4*b, m+4*L), shape):
        limits.tensor(candidate)
    n, z = max(shape[:2]), shape[2]
    limits.maps((n,)*3, (z,)*3)  # Three full slices and quotient-sized matrices.
    limits.maps((n,)*3, (n,)*3)  # Normalized slices and commutator.
    # A naive supplied scheme has q*(m+8ab) terms. This check precedes
    # its coefficient-family allocation and verification.
    if q*(m+8*a*b)*sum(shape) > limits.max_constraint_entries:
        raise StateAssemblyLimit("supplied scheme coefficient cap reached; construction incomplete")
    return shape


def mixed_output_quotient(field, selected, output_size, *, mixed=False, limits=DEFAULT_LIMITS):
    """Explicit W/K -> V coordinates, allowing a genuinely mixed complement.

    For every discarded coordinate d_j, K contains
    e_d_j + sum_i H[i,j] e_selected_i; Q(e_d_j)=-H[:,j].
    Q is the identity on the selected V coordinates. Canonical field
    encodings in H include nonbase elements when working over F4.
    """
    if not isinstance(field, FiniteField) or not isinstance(limits, StateAssemblyLimits):
        raise ValueError("supply an exact finite field and state assembly limits")
    if type(output_size) is not int or output_size < 0 or type(mixed) is not bool:
        raise ValueError("output size must be nonnegative and mixed must be boolean")
    if output_size*output_size > limits.max_map_entries:
        raise StateAssemblyLimit("quotient and complement map cap reached; construction incomplete")
    if output_size > limits.max_constraint_entries:
        raise StateAssemblyLimit("quotient coordinate cap reached; construction incomplete")
    selected = tuple(selected)
    if len(set(selected)) != len(selected) or any(
            type(i) is not int or not 0 <= i < output_size for i in selected):
        raise ValueError("selected output coordinates must be distinct and in range")
    chosen = set(selected)
    discarded = tuple(i for i in range(output_size) if i not in chosen)
    h = tuple(tuple((i+2*j+1) % field.order if mixed else 0
                    for j in range(len(discarded))) for i in range(len(selected)))
    rows = []
    for i, position in enumerate(selected):
        row = [0]*output_size
        row[position] = 1
        for j, old in enumerate(discarded):
            row[old] = field.neg(h[i][j])
        rows.append(tuple(row))
    complement = []
    for j, position in enumerate(discarded):
        column = [0]*output_size
        column[position] = 1
        for i, old in enumerate(selected):
            column[old] = h[i][j]
        complement.append(tuple(column))
    return LinearMap(field, output_size, tuple(rows)), tuple(complement)


@dataclass(frozen=True)
class ConvolutionCommutatorControl:
    tensor: FiniteTensor
    a: int
    b: int
    m: int
    q: int
    swapped_inputs: bool
    slice_vectors: tuple[tuple[int, ...], ...]
    slices: tuple[LinearMap, LinearMap, LinearMap]
    image_coordinates: tuple[int, ...]
    quotient: LinearMap
    complement: tuple[tuple[int, ...], ...]
    normalized_slices: tuple[LinearMap, LinearMap, LinearMap]
    commutator: LinearMap
    output_rank: int
    identity_rank: int
    commutator_rank: int
    supplied_scheme: TensorScheme

    @property
    def lower_bound(self):
        return self.output_rank+(self.commutator_rank+1)//2

    @property
    def comparison_gap(self):
        """Per-q rank gap against 5 copies of C; not a finite-field upper."""
        return self.lower_bound-5*self.q*(self.a+self.b-1)


def convolution_commutator_control(field, a, b, *, m=1, q=1,
                                    swapped_inputs=False, mixed_complement=False,
                                    limits=DEFAULT_LIMITS):
    """Construct and check actual slices and a supplied exact naive scheme."""
    if not isinstance(field, FiniteField):
        raise ValueError("supply an explicitly represented finite field")
    if any(type(value) is not int or value < 1 for value in (a, b, q)):
        raise ValueError("convolution dimensions and repetition count must be positive integers")
    if type(m) is not int or m < 0:
        raise ValueError("unit count must be a nonnegative integer")
    if type(swapped_inputs) is not bool or type(mixed_complement) is not bool:
        raise ValueError("input swap and mixed complement flags must be booleans")
    _preflight(a, b, m, q, limits)
    unit = FiniteTensor.unit(field)
    convolution = convolution_tensor(field, a, b)
    core = matrix_multiplication_tensor(field, 2, 2, 2).tensor_product(convolution)
    block = FiniteTensor.direct_sum((unit,)*m+(core,))
    tensor = FiniteTensor.direct_sum((block,)*q)
    if swapped_inputs:
        tensor = tensor.permute_axes((1, 0, 2))
    first_size, column_degree = (b, a) if swapped_inputs else (a, b)
    L = a+b-1
    x, n, z = tensor.shape
    vectors = [[0]*x for _ in range(3)]
    selected = []
    for copy in range(q):
        start = copy*(m+4*first_size)
        output_start = copy*(m+4*L)
        for scalar in range(m):
            vectors[0][start+scalar] = 1
            selected.append(output_start+scalar)
        vectors[0][start+m] = vectors[0][start+m+3*first_size] = 1
        vectors[1][start+m+first_size] = 1  # E12 at polynomial degree 0.
        vectors[2][start+m+2*first_size] = 1  # E21 at polynomial degree 0.
        selected.extend(output_start+m+matrix*L+degree
                        for matrix in range(4) for degree in range(column_degree))
    slices = tuple(slice_map(tensor, vector, limits=limits) for vector in vectors)
    selected = tuple(selected)
    if len(selected) != n:
        raise ValueError("identity image coordinate count disagrees with input column size")
    selected_set = set(selected)
    position = {old: i for i, old in enumerate(selected)}
    embedding = LinearMap(field, n,
                          tuple(tuple(int(position.get(k) == i) for i in range(n))
                                for k in range(z)))
    if slices[0] != embedding:
        raise ValueError("the actual identity slice is not the claimed column embedding")
    if any(any(row) for mapping in slices for k, row in enumerate(mapping.rows)
           if k not in selected_set):
        raise ValueError("a selected slice has image outside the identity image V")
    quotient, complement = mixed_output_quotient(field, selected, z, mixed=mixed_complement, limits=limits)
    if quotient.compose(embedding) != LinearMap.identity(field, n):
        raise ValueError("output quotient is not the identity on V")
    if any(any(quotient.apply(column)) for column in complement):
        raise ValueError("output quotient does not annihilate its claimed complement")
    # A basis of K is visibly independent modulo V: its discarded coordinate
    # matrix is the identity. K therefore has dimension Z-N and K intersect V=0.
    if len(complement) != z-n:
        raise ValueError("output complement has the wrong dimension")
    normalized = tuple(quotient.compose(mapping) for mapping in slices)
    if normalized[0] != LinearMap.identity(field, n):
        raise ValueError("normalized identity slice is not I_N")
    bc, cb = normalized[1].compose(normalized[2]), normalized[2].compose(normalized[1])
    commutator = LinearMap(field, n,
                           tuple(tuple(field.sub(left, right) for left, right in zip(brow, crow))
                                 for brow, crow in zip(bc.rows, cb.rows)))
    output_rank = output_flattening_rank(tensor, limits=limits)
    identity_rank = _rank(field, slices[0].rows, n)
    commutator_rank = _rank(field, commutator.rows, n)
    if (output_rank, identity_rank, commutator_rank) != (z, n, 4*q*column_degree):
        raise ValueError("actual flattening or commutator rank differs from the audited formula")
    scheme = TensorScheme.from_tensor(tensor)
    scheme.require_exact()
    lower = z+(commutator_rank+1)//2
    if scheme.terms != q*(m+8*a*b) or lower > scheme.terms:
        raise ValueError("audited lower bound exceeds the supplied exact rank upper control")
    return ConvolutionCommutatorControl(tensor, a, b, m, q, swapped_inputs,
                                        tuple(tuple(vector) for vector in vectors), slices,
                                        selected, quotient, complement, normalized,
                                        commutator, output_rank, identity_rank,
                                        commutator_rank, scheme)


def main():
    fields = (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)), FiniteField(5))
    cases = 0
    for field, a, b, q, swap in product(fields, range(2, 5), range(2, 5), (1, 2), (False, True)):
        control = convolution_commutator_control(field, a, b, q=q, swapped_inputs=swap)
        column_degree = a if swap else b
        assert control.lower_bound == q*(1+4*(a+b-1)+2*column_degree)
        if column_degree == max(a, b):
            assert control.comparison_gap > 0
        cases += 1
    for field in fields:
        convolution_commutator_control(field, 2, 4, q=2, mixed_complement=True)
    print(f"Verified {cases} complete-convolution slice/flattening/commutator controls and 4 mixed complements.")
    print("Fields F2,F3,F4,F5; a,b=2..4; q=1,2; both input orientations; exact supplied upper schemes.")
    print("The audited rank bound is q*(m+4*(a+b-1)+2*max(a,b)); no catalyst witness or Lean theorem is produced.")


if __name__ == "__main__":
    main()
