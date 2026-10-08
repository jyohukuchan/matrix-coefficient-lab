"""Bounded exact F2/F4 restrictions and optional unseeded SAT experiments.

This standalone research experiment is not a near-9/4 coefficient generator.
The default demonstration uses only the standard library. Optional z3 is loaded
only inside SAT search functions; no known coefficient scheme is planted.
"""

from dataclasses import dataclass
from collections import Counter
from itertools import combinations, product
from math import comb
from random import Random
import argparse
import json
import sys

from .fields import FiniteField
from .maps import LinearMap
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, matrix_multiplication_tensor


@dataclass(frozen=True)
class SearchResult:
    status: str
    reason: str
    checked: int
    total_candidates: int | None
    maps: tuple[LinearMap, LinearMap, LinearMap] | None = None
    scheme: TensorScheme | None = None


def gf2_rank(rows, columns):
    """Rank of bit-packed rows; ignore bits outside the declared matrix."""
    basis = {}
    mask = (1 << columns)-1
    for row in rows:
        row &= mask
        while row:
            pivot = row.bit_length()-1
            if pivot not in basis:
                basis[pivot] = row
                break
            row ^= basis[pivot]
    return len(basis)


def solve_many_gf2(rows, right_sides, variables, outputs):
    """Solve all RHS columns together; return one solution with free values zero.

    Each row and RHS row is bit-packed. The returned tuple contains one bit
    mask of values per RHS column, or None if any RHS is inconsistent.
    """
    rows, right_sides = tuple(rows), tuple(right_sides)
    if len(rows) != len(right_sides):
        raise ValueError("matrix and right-hand sides must have equal row counts")
    if any(row < 0 or row >> variables for row in rows):
        raise ValueError("matrix row does not fit declared variable count")
    if any(row < 0 or row >> outputs for row in right_sides):
        raise ValueError("right-hand side row does not fit declared output count")
    augmented = [row | (rhs << variables) for row, rhs in zip(rows, right_sides)]
    position = 0
    pivots = []
    for variable in range(variables):
        selected = next((i for i in range(position, len(augmented))
                         if augmented[i] >> variable & 1), None)
        if selected is None:
            continue
        augmented[position], augmented[selected] = augmented[selected], augmented[position]
        for i in range(len(augmented)):
            if i != position and augmented[i] >> variable & 1:
                augmented[i] ^= augmented[position]
        pivots.append(variable)
        position += 1
    mask = (1 << variables)-1
    if any((row & mask) == 0 and row >> variables for row in augmented):
        return None
    solutions = [0]*outputs
    for i, variable in enumerate(pivots):
        for output in range(outputs):
            solutions[output] |= ((augmented[i] >> (variables+output)) & 1) << variable
    return tuple(solutions)


def solve_many_field(field, rows, right_sides, variables, outputs):
    """Exact finite-field Gaussian elimination, with all RHS columns together."""
    if any(type(value) is not int or value < 0 for value in (variables, outputs)):
        raise ValueError("Gaussian variable and output counts must be nonnegative integers")
    rows, right_sides = tuple(tuple(row) for row in rows), tuple(tuple(row) for row in right_sides)
    if len(rows) != len(right_sides) or any(len(row) != variables for row in rows) or any(len(row) != outputs for row in right_sides):
        raise ValueError("matrix/RHS dimensions do not match the declared shape")
    augmented = [list(row+rhs) for row, rhs in zip(rows, right_sides)]
    for row in augmented:
        for value in row:
            field.check(value)
    position = 0
    pivots = []
    for variable in range(variables):
        selected = next((i for i in range(position, len(augmented)) if augmented[i][variable]), None)
        if selected is None:
            continue
        augmented[position], augmented[selected] = augmented[selected], augmented[position]
        inverse = field.inv(augmented[position][variable])
        augmented[position] = [field.mul(value, inverse) for value in augmented[position]]
        for i in range(len(augmented)):
            if i != position and augmented[i][variable]:
                scalar = augmented[i][variable]
                augmented[i] = [field.sub(value, field.mul(scalar, pivot))
                                for value, pivot in zip(augmented[i], augmented[position])]
        pivots.append(variable)
        position += 1
    if any(not any(row[:variables]) and any(row[variables:]) for row in augmented):
        return None
    solution = [[0]*variables for _ in range(outputs)]
    for i, variable in enumerate(pivots):
        for output in range(outputs):
            solution[output][variable] = augmented[i][variables+output]
    return tuple(tuple(row) for row in solution)


def field_matrix_rank(field, rows, columns):
    """Exact rank, including extension-field coefficients in canonical encoding."""
    if type(columns) is not int or columns < 0:
        raise ValueError("matrix column count must be a nonnegative integer")
    rows = [list(row) for row in rows]
    if any(len(row) != columns for row in rows):
        raise ValueError("rank matrix rows have the wrong length")
    for row in rows:
        for value in row:
            field.check(value)
    position = 0
    for column in range(columns):
        selected = next((i for i in range(position, len(rows)) if rows[i][column]), None)
        if selected is None:
            continue
        rows[position], rows[selected] = rows[selected], rows[position]
        inverse = field.inv(rows[position][column])
        rows[position] = [field.mul(value, inverse) for value in rows[position]]
        for i in range(position+1, len(rows)):
            if rows[i][column]:
                scalar = rows[i][column]
                rows[i] = [field.sub(value, field.mul(scalar, pivot))
                           for value, pivot in zip(rows[i], rows[position])]
        position += 1
    return position


def flattening_ranks(tensor):
    if tensor.field not in (FiniteField(2), FiniteField(2, (1, 1, 1))):
        raise ValueError("this experiment supports explicitly represented F2 or F4")
    binary = tensor.field == FiniteField(2)
    result = []
    for axis in range(3):
        sizes = tensor.shape
        others = tuple(i for i in range(3) if i != axis)
        width = sizes[others[0]]*sizes[others[1]]
        rows = []
        for i in range(sizes[axis]):
            row = 0
            field_row = []
            for j, k in product(range(sizes[others[0]]), range(sizes[others[1]])):
                coordinate = [0, 0, 0]
                coordinate[axis] = i
                coordinate[others[0]], coordinate[others[1]] = j, k
                value = tensor.coefficient(*coordinate)
                if binary:
                    row |= value << (j*sizes[others[1]]+k)
                field_row.append(value)
            rows.append(row if binary else field_row)
        result.append(gf2_rank(rows, width) if binary
                      else field_matrix_rank(tensor.field, rows, width))
    return tuple(result)


def _map_from_bits(field, inputs, outputs, bits):
    return LinearMap(field, inputs,
                     tuple(tuple((bits >> (x*inputs+i)) & 1 for i in range(inputs))
                           for x in range(outputs)))


def solve_third_leg(source, target, left, right):
    """Given two actual local maps, solve for the third and verify all coefficients."""
    if source.field not in (FiniteField(2), FiniteField(2, (1, 1, 1))) or target.field != source.field:
        raise ValueError("source and target must use the same explicitly represented F2 or F4")
    if (left.input_size, left.output_size) != (source.shape[0], target.shape[0]):
        raise ValueError("left map has incorrect dimensions")
    if (right.input_size, right.output_size) != (source.shape[1], target.shape[1]):
        raise ValueError("right map has incorrect dimensions")
    if left.field != source.field or right.field != source.field:
        raise ValueError("local maps have the wrong field")
    source_z, target_z = source.shape[2], target.shape[2]
    if source.field != FiniteField(2):
        field = source.field
        rows, rhs = [], []
        for x, y in product(range(target.shape[0]), range(target.shape[1])):
            rows.append(tuple(field.sum(field.mul(a, field.mul(b, source.coefficient(i, j, z)))
                                        for (i, a), (j, b) in product(left.row_supports[x], right.row_supports[y]))
                              for z in range(source_z)))
            rhs.append(tuple(target.coefficient(x, y, z) for z in range(target_z)))
        solution = solve_many_field(field, rows, rhs, source_z, target_z)
        if solution is None:
            return None
        third = LinearMap(field, source_z, solution)
        if source.restrict(left, right, third) != target:
            raise AssertionError("extension-field solver fails the full restriction coefficient check")
        return third
    rows, rhs = [], []
    for x, y in product(range(target.shape[0]), range(target.shape[1])):
        row = 0
        for z in range(source_z):
            value = 0
            for (i, a), (j, b) in product(left.row_supports[x], right.row_supports[y]):
                value ^= a*b*source.coefficient(i, j, z)
            row |= value << z
        rows.append(row)
        rhs.append(sum(target.coefficient(x, y, z) << z for z in range(target_z)))
    solution = solve_many_gf2(rows, rhs, source_z, target_z)
    if solution is None:
        return None
    third = LinearMap(source.field, source_z,
                      tuple(tuple(bits >> z & 1 for z in range(source_z)) for bits in solution))
    if source.restrict(left, right, third) != target:
        raise AssertionError("linear solver failed the independent coefficient restriction check")
    return third


def bounded_restriction_search(source, target, *, max_pairs=10000):
    """Enumerate ALL first-two-leg maps in bit order, with an explicit pair cap.

    'exhausted' is nonexistence only for these exact fixed source/target tensors.
    'resource_cap' never implies nonexistence. 'filtered' is an exact necessary
    flattening obstruction. This search does not range over catalyst tensors.
    """
    if type(max_pairs) is not int or max_pairs < 0:
        raise ValueError("pair cap must be a nonnegative integer")
    if source.field != FiniteField(2) or target.field != source.field:
        raise ValueError("source and target must be over F2")
    left_bits = source.shape[0]*target.shape[0]
    right_bits = source.shape[1]*target.shape[1]
    total = 1 << (left_bits+right_bits)
    source_ranks, target_ranks = flattening_ranks(source), flattening_ranks(target)
    if any(t > s for s, t in zip(source_ranks, target_ranks)):
        return SearchResult("filtered", "target flattening rank exceeds source rank", 0, total)
    checked = 0
    for bits_left in range(1 << left_bits):
        left = _map_from_bits(source.field, source.shape[0], target.shape[0], bits_left)
        for bits_right in range(1 << right_bits):
            if checked == max_pairs:
                return SearchResult("resource_cap", "pair cap reached; search is incomplete", checked, total)
            right = _map_from_bits(source.field, source.shape[1], target.shape[1], bits_right)
            checked += 1
            # Necessary map-rank conditions, independent of chosen C.
            left_rows = [sum(value << i for i, value in enumerate(row)) for row in left.rows]
            right_rows = [sum(value << i for i, value in enumerate(row)) for row in right.rows]
            if (gf2_rank(left_rows, source.shape[0]) < target_ranks[0]
                    or gf2_rank(right_rows, source.shape[1]) < target_ranks[1]):
                continue
            third = solve_third_leg(source, target, left, right)
            if third is not None:
                return SearchResult("found", "actual local maps independently coefficient-verified",
                                    checked, total, (left, right, third))
    return SearchResult("exhausted", "all map pairs exhausted for this fixed source and target", checked, total)


def unseeded_square_sat(n=2, rank=7, *, timeout_ms=30000, solver_path=None):
    """Solve formal coefficient equations without importing or planting a scheme.

    Sound symmetry: normalize nonzero F2 input pairs (already scalar-normalized),
    merge repeated input pairs, sort distinct active pairs, and pad with zero
    terms at the end. Thus UNSAT concerns rank AT MOST the supplied rank.
    Optional z3 dependency stays outside the ordinary reference implementation.
    """
    if type(n) is not int or n < 1 or type(rank) is not int or rank < 1:
        raise ValueError("matrix size and rank must be positive integers")
    if type(timeout_ms) is not int or timeout_ms < 1:
        raise ValueError("solver timeout must be a positive integer")
    if solver_path is not None:
        sys.path.insert(0, solver_path)
    try:
        import z3
    except ImportError:
        return SearchResult("unavailable", "optional z3-solver is not installed", 0, None)
    width = n*n
    target = matrix_multiplication_tensor(FiniteField(2), n, n, n)
    aa = [z3.BitVec(f"a_{t}", width) for t in range(rank)]
    bb = [z3.BitVec(f"b_{t}", width) for t in range(rank)]
    cc = [z3.BitVec(f"c_{t}", width) for t in range(rank)]
    solver = z3.Solver()
    solver.set(timeout=timeout_ms, random_seed=0)
    active = [a != 0 for a in aa]
    for t in range(rank):
        solver.add(active[t] == (bb[t] != 0), active[t] == (cc[t] != 0))
        if t+1 < rank:
            solver.add(z3.Implies(active[t+1], active[t]))
            solver.add(z3.Implies(active[t+1],
                                 z3.ULT(z3.Concat(aa[t], bb[t]), z3.Concat(aa[t+1], bb[t+1]))))
    bit = lambda vector, i: z3.Extract(i, i, vector) == 1
    for x, y, z in product(range(width), repeat=3):
        values = [z3.And(bit(a, x), bit(b, y), bit(c, z)) for a, b, c in zip(aa, bb, cc)]
        parity = values[0]
        for value in values[1:]:
            parity = z3.Xor(parity, value)
        solver.add(parity == bool(target.coefficient(x, y, z)))
    outcome = solver.check()
    if outcome == z3.unsat:
        return SearchResult("exhausted", "complete normalized coefficient encoding is UNSAT", 0, None)
    if outcome != z3.sat:
        return SearchResult("resource_cap", f"solver incomplete: {solver.reason_unknown()}", 0, None)
    model = solver.model()
    family = lambda vectors: tuple(tuple((model.eval(v).as_long() >> i) & 1 for i in range(width))
                                  for v in vectors)
    a, b, c = family(aa), family(bb), family(cc)
    scheme = TensorScheme(target, a, b, c)
    scheme.require_exact()
    # Independently eliminate C after SAT chose A,B; do not rely on SAT's C.
    diagonal = diagonal_tensor(target.field, rank)
    left = LinearMap(target.field, rank, tuple(tuple(a[t][i] for t in range(rank)) for i in range(width)))
    right = LinearMap(target.field, rank, tuple(tuple(b[t][i] for t in range(rank)) for i in range(width)))
    third = solve_third_leg(diagonal, target, left, right)
    if third is None:
        raise AssertionError("SAT factors do not admit an independently solved output family")
    recovered_c = tuple(tuple(third.rows[z][t] for z in range(width)) for t in range(rank))
    recovered = TensorScheme(target, a, b, recovered_c)
    recovered.require_exact()
    return SearchResult("found", "unseeded SAT certificate, independently recovered C and all coefficients verified",
                        1, None, (left, right, third), recovered)


def bounded_restriction_sat(source, target, *, timeout_ms=60000, max_memory_mb=512,
                            max_products=50000, max_map_bits=4096,
                            max_boolean_products=250000,
                            identical_blocks=None, solver_path=None):
    """Optional exact full-map SAT with explicit construction/time/memory caps.

    UNSAT excludes only this fixed source and target over its supplied F2/F4.
    F4 uses two canonical power-basis bits with alpha^2=alpha+1. If identical_blocks
    is supplied, verify the source is exactly that many identical direct-sum
    blocks, then sort their first-leg map blocks. This symmetry is sound.
    """
    if source.field not in (FiniteField(2), FiniteField(2, (1, 1, 1))) or target.field != source.field:
        raise ValueError("source and target must use the same explicitly represented F2 or F4")
    field = source.field
    degree = field.degree
    if any(type(value) is not int or value < 1 for value in
           (timeout_ms, max_memory_mb, max_products, max_map_bits, max_boolean_products)):
        raise ValueError("solver resource caps must be positive integers")
    map_bits = degree*sum(a*b for a, b in zip(source.shape, target.shape))
    if map_bits > max_map_bits:
        return SearchResult("resource_cap", f"encoding needs {map_bits} map bits; construction cap exceeded", 0, None)
    if not any(target.coefficients):
        maps = tuple(LinearMap(source.field, inputs, ((0,)*inputs,)*outputs)
                     for inputs, outputs in zip(source.shape, target.shape))
        if source.restrict(*maps) != target:
            raise AssertionError("zero maps fail the exact zero target check")
        return SearchResult("found", "exact zero target restricted by explicit zero maps", 0, None, maps)
    if any(t > s for s, t in zip(flattening_ranks(source), flattening_ranks(target))):
        return SearchResult("filtered", "target flattening rank exceeds source rank", 0, None)
    support = tuple((x, y, z, source.coefficient(x, y, z)) for x, y, z in product(*(range(d) for d in source.shape))
                    if source.coefficient(x, y, z))
    products_count = len(support)*target.shape[0]*target.shape[1]*target.shape[2]
    if products_count > max_products:
        return SearchResult("resource_cap", f"encoding would require {products_count} cubic products; construction cap exceeded", 0, None)
    cubic_terms = {}
    for scalar in {value for x, y, z, value in support}:
        cubic_terms[scalar] = tuple((a, b, c, field.mul(scalar, field.mul(field.mul(1 << a, 1 << b), 1 << c)))
                                    for a, b, c in product(range(degree), repeat=3))
    boolean_products = len(target.coefficients)*sum(
        sum(value.bit_count() for a, b, c, value in cubic_terms[scalar]) for x, y, z, scalar in support)
    if boolean_products > max_boolean_products:
        return SearchResult("resource_cap", f"encoding needs {boolean_products} Boolean cubic products; construction cap exceeded", 0, None)
    block_width = None
    if identical_blocks is not None:
        if type(identical_blocks) is not int or identical_blocks < 1 or any(d % identical_blocks for d in source.shape):
            raise ValueError("identical block count must divide every source dimension")
        shape = tuple(d//identical_blocks for d in source.shape)
        block = FiniteTensor.from_function(source.field, shape, source.coefficient)
        if FiniteTensor.direct_sum((block,)*identical_blocks) != source:
            raise ValueError("source is not the supplied number of identical direct-sum blocks")
        block_width = shape[0]
    if solver_path is not None:
        sys.path.insert(0, solver_path)
    try:
        import z3
    except ImportError:
        return SearchResult("unavailable", "optional z3-solver is not installed", 0, None)
    solver = z3.Solver()
    solver.set(timeout=timeout_ms, max_memory=max_memory_mb, random_seed=0)
    vectors = tuple(tuple(z3.BitVec(f"map_{axis}_{row}", degree*inputs) for row in range(outputs))
                    for axis, (inputs, outputs) in enumerate(zip(source.shape, target.shape)))
    bits = tuple(tuple(tuple(tuple(z3.Extract(degree*i+component, degree*i+component, row) == 1
                                  for component in range(degree)) for i in range(inputs)) for row in family)
                 for inputs, family in zip(source.shape, vectors))
    for axis, family in enumerate(vectors):
        for row, vector in enumerate(family):
            if any(target.coefficient(*index) for index in product(*(range(d) for d in target.shape))
                   if index[axis] == row):
                solver.add(vector != 0)
    if block_width:
        keys = []
        for block in range(identical_blocks):
            components = [z3.Extract(degree*(block*block_width+i)+degree-1,
                                     degree*(block*block_width+i), row)
                          for row in vectors[0] for i in range(block_width)]
            keys.append(components[0] if len(components) == 1 else z3.Concat(*components))
        solver.add(*(z3.ULE(first, second) for first, second in zip(keys, keys[1:])))
    for x, y, z in product(*(range(d) for d in target.shape)):
        terms = []
        for i, j, k, scalar in support:
            terms.extend((z3.And(bits[0][x][i][a], bits[1][y][j][b], bits[2][z][k][c]), value)
                         for a, b, c, value in cubic_terms[scalar])
        for component in range(degree):
            parity = z3.BoolVal(False)
            for expression, scalar in terms:
                if scalar >> component & 1:
                    parity = z3.Xor(parity, expression)
            solver.add(parity == bool(target.coefficient(x, y, z) >> component & 1))
    outcome = solver.check()
    if outcome == z3.unsat:
        return SearchResult("exhausted", "full coefficient restriction encoding UNSAT for these exact tensors", 1, None)
    if outcome != z3.sat:
        return SearchResult("resource_cap", f"solver incomplete: {solver.reason_unknown()}", 1, None)
    model = solver.model()
    maps = tuple(LinearMap(source.field, inputs,
                           tuple(tuple(model.eval(row, model_completion=True).as_long() >> (degree*i) & (field.order-1)
                                       for i in range(inputs)) for row in family))
                 for inputs, family in zip(source.shape, vectors))
    if source.restrict(*maps) != target:
        raise AssertionError("SAT model fails the independent full restriction coefficient check")
    third = solve_third_leg(source, target, maps[0], maps[1])
    if third is None:
        raise AssertionError("SAT input maps fail independent third-leg elimination")
    return SearchResult("found", "actual SAT restriction independently verified; third leg re-solved", 1, None,
                        (maps[0], maps[1], third))


def unseeded_convolution_sat(field, a, b, rank, *, timeout_ms=10000, solver_path=None):
    """Find formal F2/F4 convolution coefficients from an unseeded scalar source.

    No projective-interpolation factors or other algorithm are supplied to SAT.
    A cap or UNSAT result concerns only the specified field and rank threshold.
    """
    from .sectors import convolution_tensor

    if type(rank) is not int or rank < 0:
        raise ValueError("rank threshold must be a nonnegative integer")
    target = convolution_tensor(field, a, b)
    result = bounded_restriction_sat(diagonal_tensor(field, rank), target,
                                    timeout_ms=timeout_ms, solver_path=solver_path)
    if result.status != "found":
        return result
    left, right, third = result.maps
    scheme = TensorScheme(target,
                          tuple(tuple(left.rows[x][t] for x in range(a)) for t in range(rank)),
                          tuple(tuple(right.rows[y][t] for y in range(b)) for t in range(rank)),
                          tuple(tuple(third.rows[z][t] for z in range(target.shape[2]))
                                for t in range(rank)))
    scheme.require_exact()
    return SearchResult("found", "unseeded convolution factors with independently recovered C and all coefficients verified",
                        result.checked, result.total_candidates, result.maps, scheme)


def bounded_s233_catalyst_experiment(*, pencil_kind="regular", timeout_ms=60000, solver_path=None):
    """One fixed D=0, d=2, k=5 candidate; no claim about arbitrary catalysts."""
    from .catalysts import CatalyticCertificate

    f = FiniteField(2)
    if pencil_kind == "regular":
        S = FiniteTensor.from_function(f, (2, 3, 3), lambda x, y, z: int(z == y+x))
        description = "D=0,S=(I3,J3),d=2,k=5,m=1"
    elif pencil_kind == "singular":
        support = ((0, 0, 0), (0, 1, 2), (1, 0, 1), (1, 2, 2))
        S = FiniteTensor.from_function(f, (2, 3, 3), lambda x, y, z: int((x, y, z) in support))
        description = "D=0,S=[[x,y,0],[0,0,x],[0,0,y]],d=2,k=5,m=1"
    else:
        raise ValueError("pencil kind must be regular or singular")
    scheme_search = bounded_tensor_rank_search(S, max_rank=4, max_subsets=7315)
    if scheme_search.status != "found" or scheme_search.scheme.terms != 4:
        raise AssertionError("S233 does not have the expected certified four-term upper bound")
    lower = output_space_rank_one_filter(S)
    if lower["tensor_rank_lower_bound"] != 4:
        raise AssertionError("S233 output-space filter did not certify the rank-four lower bound")
    D = FiniteTensor.zero(f, (0, 0, 0))
    source = FiniteTensor.direct_sum((S,)*5)
    target = FiniteTensor.direct_sum((FiniteTensor.unit(f), matrix_multiplication_tensor(f, 2, 2, 2).tensor_product(S)))
    result = bounded_restriction_sat(source, target, timeout_ms=timeout_ms, max_memory_mb=512,
                                    max_products=50000, identical_blocks=5, solver_path=solver_path)
    verified = False
    if result.status == "found":
        certificate = CatalyticCertificate(2, 5, 1, D, S, result.maps)
        certificate.require_valid()
        verified = True
    return {"field": "F2", "fixed_candidate": description,
            "source_shape": source.shape, "target_shape": target.shape,
            "map_coefficient_bits": sum(a*b for a, b in zip(source.shape, target.shape)),
            "coefficient_equations": len(target.coefficients),
            "cubic_product_occurrences": sum(bool(x) for x in source.coefficients)*len(target.coefficients),
            "timeout_ms": timeout_ms, "solver_memory_cap_mb": 512,
            "actual_catalytic_certificate_verified": verified,
            "S_coefficients": S.coefficients, "S_exact_scheme": _summary(scheme_search),
            "S_rank_lower_filter": lower, "result": _summary(result),
            "limitations": "Exhaustion excludes only this fixed D=0 tensor over F2; a cap or unknown gives no nonexistence conclusion."}


def bounded_singular_rectangular_embedding(*, field=None, timeout_ms=60000, solver_path=None):
    """Historical bounded search of a now analytically excluded restriction.

    The common-image obstruction in research/pencil-frontier.md excludes
    this restriction over every field. Keep the old experiment reproducible,
    but do not spend another dense search budget on it. A claimed SAT model
    still receives the complete independent coefficient check.
    """
    f = FiniteField(2) if field is None else field
    if f not in (FiniteField(2), FiniteField(2, (1, 1, 1))):
        raise ValueError("embedding experiment supports explicit F2 or F4")
    support = ((0, 0, 0), (0, 1, 2), (1, 0, 1), (1, 2, 2))
    S = FiniteTensor.from_function(f, (2, 3, 3), lambda x, y, z: int((x, y, z) in support))
    source = matrix_multiplication_tensor(f, 2, 2, 2).tensor_product(S)
    target = matrix_multiplication_tensor(f, 3, 2, 4)
    result = bounded_restriction_sat(source, target, timeout_ms=timeout_ms, max_memory_mb=512,
                                    max_products=50000, solver_path=solver_path)
    augmented_maps = None
    if result.status == "found":
        augmented_maps = tuple(LinearMap(f, mapping.input_size+1,
                                         ((1,)+(0,)*mapping.input_size,)
                                         +tuple((0,)+row for row in mapping.rows)) for mapping in result.maps)
        full_source = FiniteTensor.direct_sum((FiniteTensor.unit(f), source))
        full_target = FiniteTensor.direct_sum((FiniteTensor.unit(f), target))
        if full_source.restrict(*augmented_maps) != full_target:
            raise AssertionError("scalar-identity augmentation does not preserve the full restriction")
    field_products = sum(bool(value) for value in source.coefficients)*len(target.coefficients)
    return {"field": f"F{f.order}", "field_modulus": f.modulus,
            "fixed_embedding": "M2 tensor singular S233 -> M(3,2,4)",
            "source_shape": source.shape, "target_shape": target.shape,
            "map_field_variables": sum(a*b for a, b in zip(source.shape, target.shape)),
            "map_coefficient_bits": f.degree*sum(a*b for a, b in zip(source.shape, target.shape)),
            "coefficient_equations": len(target.coefficients),
            "bit_coefficient_equations": f.degree*len(target.coefficients),
            "cubic_product_occurrences": field_products,
            "Boolean_cubic_product_occurrences": field_products*(1 if f.degree == 1 else 11),
            "timeout_ms": timeout_ms, "solver_memory_cap_mb": 512,
            "actual_embedding_verified": result.status == "found",
            "scalar_gain_identity_verified": augmented_maps is not None,
            "scalar_augmented_maps": None if augmented_maps is None else [mapping.rows for mapping in augmented_maps],
            "S_coefficients": S.coefficients, "result": _summary(result),
            "limitations": "A cap gives no embedding/nonexistence conclusion. A verified embedding is a lower-bound witness, not a positive catalyst."}


def diagonal_tensor(field, size):
    return FiniteTensor.from_function(field, (size, size, size), lambda x, y, z: int(x == y == z))


def output_space_rank_one_filter(target):
    """An exact finite F2 obstruction to attaining the output flattening rank.

    A decomposition with exactly R=dim(output slice space) terms must have
    every input outer product in that space. If its complete rank-one inventory
    does not span the space, rank must exceed R. This is not a general solver.
    """
    if target.field != FiniteField(2):
        raise ValueError("rank-one inventory supports F2 only")
    a, b, c = target.shape
    slices = [sum(target.coefficient(x, y, z) << (x*b+y)
                  for x, y in product(range(a), range(b))) for z in range(c)]
    rank = gf2_rank(slices, a*b)
    pairs, columns = [], []
    for left, right in product(range(1, 1 << a), range(1, 1 << b)):
        column = sum((left >> x & 1)*(right >> y & 1) << (x*b+y)
                     for x, y in product(range(a), range(b)))
        if gf2_rank(slices+[column], a*b) == rank:
            pairs.append((left, right))
            columns.append(column)
    span_rank = gf2_rank(columns, a*b)
    return {"output_flattening_rank": rank, "rank_one_pairs_in_output_space": pairs,
            "rank_one_span_rank": span_rank,
            "tensor_rank_lower_bound": rank+int(span_rank < rank)}


def bounded_tensor_rank_search(target, *, max_rank, max_subsets=10000, max_input_pairs=10000):
    """Unseeded normalized-input enumeration for an arbitrary tiny F2 tensor."""
    if target.field != FiniteField(2):
        raise ValueError("bounded coefficient search supports F2 only")
    if any(type(value) is not int or value < 0 for value in (max_rank, max_subsets, max_input_pairs)):
        raise ValueError("rank and search caps must be nonnegative integers")
    if not any(target.coefficients):
        scheme = TensorScheme(target, (), (), ())
        scheme.require_exact()
        return SearchResult("found", "exact zero target with zero terms", 0, 0, scheme=scheme)
    a, b, c = target.shape
    pair_count = ((1 << a)-1)*((1 << b)-1)
    if pair_count > max_input_pairs:
        return SearchResult("resource_cap", "input-factor inventory exceeds construction cap", 0, None)
    vectors = lambda n: tuple(tuple(mask >> i & 1 for i in range(n)) for mask in range(1, 1 << n))
    pairs = tuple(product(vectors(a), vectors(b)))
    lower = max(flattening_ranks(target))
    if lower > max_rank:
        return SearchResult("filtered", "flattening rank exceeds maximum permitted rank", 0, 0)
    total = sum(comb(len(pairs), rank) for rank in range(lower, min(max_rank, len(pairs))+1))
    rhs = tuple(sum(target.coefficient(x, y, z) << z for z in range(c))
                for x, y in product(range(a), range(b)))
    checked = 0
    for rank in range(lower, min(max_rank, len(pairs))+1):
        for subset in combinations(pairs, rank):
            if checked == max_subsets:
                return SearchResult("resource_cap", "coefficient subset cap reached; search incomplete", checked, total)
            checked += 1
            rows = tuple(sum(left[x]*right[y] << term for term, (left, right) in enumerate(subset))
                         for x, y in product(range(a), range(b)))
            output = solve_many_gf2(rows, rhs, rank, c)
            if output is None:
                continue
            scheme = TensorScheme(target, tuple(left for left, right in subset),
                                  tuple(right for left, right in subset),
                                  tuple(tuple(output[z] >> term & 1 for z in range(c))
                                        for term in range(rank)))
            scheme.require_exact()
            return SearchResult("found", "unseeded coefficient family independently verified", checked, total,
                                scheme=scheme)
    return SearchResult("exhausted", "all normalized subsets through maximum rank exhausted", checked, total)


def koszul_rank_bound(target, degree, *, max_bit_entries=100000000):
    """Exact characteristic-two exterior flattening; report an explicit cap.

    T: Lambda^p(A) tensor B* -> Lambda^(p+1)(A) tensor C.
    Matrix rank divided by choose(dimA-1,p) bounds ordinary tensor rank.
    Bit-packed elimination operates on true coefficients, not sample inputs.
    """
    if target.field != FiniteField(2):
        raise ValueError("this exact exterior-flattening experiment supports F2 only")
    a, b, c = target.shape
    if type(degree) is not int or not 0 <= degree < a:
        raise ValueError("Koszul degree must lie between zero and dimA-1")
    matrix_rows, matrix_columns = comb(a, degree+1)*c, comb(a, degree)*b
    if matrix_rows*matrix_columns > max_bit_entries:
        return {"status": "resource_cap", "matrix_shape": (matrix_rows, matrix_columns)}
    inputs = tuple(combinations(range(a), degree))
    outputs = tuple(combinations(range(a), degree+1))
    positions = {subset: i for i, subset in enumerate(outputs)}
    rows = [0]*matrix_rows
    support = tuple((x, y, z) for x, y, z in product(range(a), range(b), range(c))
                    if target.coefficient(x, y, z))
    for column, subset in enumerate(inputs):
        for x, y, z in support:
            if x not in subset:
                row = positions[tuple(sorted(subset+(x,)))]*c+z
                rows[row] ^= 1 << (column*b+y)
    rank = gf2_rank(rows, matrix_columns)
    factor = comb(a-1, degree)
    return {"status": "checked", "matrix_shape": (matrix_rows, matrix_columns),
            "matrix_rank": rank, "rank_one_factor": factor,
            "tensor_rank_lower_bound": (rank+factor-1)//factor}


def convolution_family_checks():
    """Actual C22 restriction and rank checks, plus still-inconclusive filters."""
    from .projective_convolution import projective_convolution_scheme
    from .sectors import convolution_tensor

    f = FiniteField(2)
    C22 = convolution_tensor(f, 2, 2)
    C23 = convolution_tensor(f, 2, 3)
    C33 = convolution_tensor(f, 3, 3)
    C22_scheme = TensorScheme(C22, ((1, 0), (1, 1), (0, 1)),
                              ((1, 0), (1, 1), (0, 1)), ((1, 1, 0), (0, 1, 0), (0, 1, 1)))
    C22_scheme.require_exact()
    dual = FiniteTensor.from_function(f, (2, 2, 2), lambda i, j, k: int(i+j == k))
    maps = (LinearMap.identity(f, 2), LinearMap.identity(f, 2), LinearMap.selection(f, 3, (0, 1)))
    if C22.restrict(*maps) != dual:
        raise AssertionError("C22 does not restrict to actual dual-number multiplication")
    if C22_scheme.restrict(*maps).target != dual:
        raise AssertionError("scheme restriction does not preserve the dual-number tensor")
    rank4 = bounded_tensor_rank_search(C23, max_rank=4, max_subsets=5985)
    rank5 = bounded_tensor_rank_search(C23, max_rank=5, max_subsets=26334)
    if rank4.status != "exhausted" or rank5.status != "found" or rank5.scheme.terms != 5:
        raise AssertionError("unexpected complete C23 rank-search result")
    # Six-product diagonal/cross identity supplies a certified C33 upper bound.
    aa = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1))
    cc = ((1, 1, 1, 0, 0), (0, 1, 1, 1, 0), (0, 0, 1, 1, 1),
          (0, 1, 0, 0, 0), (0, 0, 1, 0, 0), (0, 0, 0, 1, 0))
    C33_scheme = TensorScheme(C33, aa, aa, cc)
    C33_scheme.require_exact()
    filter33 = output_space_rank_one_filter(C33)
    if filter33["tensor_rank_lower_bound"] != 6:
        raise AssertionError("expected the exact C33 rank-six certificate")
    # Reusable exact projective interpolation, finite nodes plus infinity.
    extension = FiniteField(2, (1, 1, 1))
    C23_extension = projective_convolution_scheme(extension, 2, 3)
    C33_extension = projective_convolution_scheme(extension, 3, 3)
    M = matrix_multiplication_tensor(f, 2, 2, 2)
    koszul = {}
    for name, tensor in (("C22", C22), ("C23", C23), ("C33", C33)):
        target = FiniteTensor.direct_sum((FiniteTensor.unit(f), M.tensor_product(tensor))).permute_axes((0, 2, 1))
        koszul[name] = [koszul_rank_bound(target, degree) for degree in (1, 2, 3)]
    rectangle = FiniteTensor.direct_sum((FiniteTensor.unit(f), matrix_multiplication_tensor(f, 4, 2, 4)))
    koszul["Unit_plus_M424"] = [koszul_rank_bound(rectangle.permute_axes((0, 2, 1)), degree)
                                 for degree in (1, 2, 3, 4)]
    rng = Random(0)
    projections = []
    for trial in range(12):
        while True:
            masks = tuple(rng.randrange(1 << 9) for _ in range(8))
            if gf2_rank(masks, 9) == 8:
                break
        mapping = LinearMap(f, 9, tuple(tuple(mask >> i & 1 for i in range(9)) for mask in masks))
        projected = rectangle.permute_axes((0, 2, 1)).restrict(mapping, LinearMap.identity(f, 17), LinearMap.identity(f, 9))
        projections.append({"trial": trial, "row_masks": masks, "bound": koszul_rank_bound(projected, 1)})
    return {"C22_exact_rank": 3, "C22_restriction_to_dual_numbers_verified": True,
            "C22_dual_restriction_maps": [mapping.rows for mapping in maps],
            "C23_rank4_complete_search": _summary(rank4),
            "C23_rank5_unseeded_search": _summary(rank5),
            "C23_output_space_filter": output_space_rank_one_filter(C23),
            "C33_exact_rank": 6, "C33_output_space_filter": filter33,
            "F4_C23_certified_terms": C23_extension.terms,
            "F4_C33_projective_certified_terms": C33_extension.terms,
            "koszul_bounds": koszul,
            "M424_projection_seed": 0, "M424_first_leg_projections": projections,
            "limitations": "These ranks and maps do not certify a useful positive catalyst; C22 exclusion invokes separate classical algebra lower bounds."}


def cube_inventory():
    """Complete 256-tensor F2 inventory, with independently certified exact ranks.

    At each rank, every normalized input-pair subset is tried. A found scheme
    is independently verified, and every smaller rank was fully exhausted.
    The eight GL2(F2)^3 orbits are enumerated, not inferred from a classification
    theorem. These are finite cube results, not catalytic restrictions.
    """
    f = FiniteField(2)
    vectors = ((1, 0), (0, 1), (1, 1))
    candidates = tuple(product(vectors, repeat=2))
    subsets = {rank: tuple(combinations(candidates, rank)) for rank in range(4)}
    tensors, metadata = {}, {}
    counts, concise, pencils = Counter(), Counter(), Counter()
    for code in range(256):
        target = FiniteTensor(f, (2, 2, 2), tuple(code >> i & 1 for i in range(8)))
        tensors[code] = target
        rhs = tuple(sum(target.coefficient(x, y, z) << z for z in range(2))
                    for x, y in product(range(2), repeat=2))
        found = None
        for rank in range(4):
            for pairs in subsets[rank]:
                rows = tuple(sum(a[x]*b[y] << term for term, (a, b) in enumerate(pairs))
                             for x, y in product(range(2), repeat=2))
                output = solve_many_gf2(rows, rhs, rank, 2)
                if output is None:
                    continue
                found = TensorScheme(target, tuple(a for a, b in pairs), tuple(b for a, b in pairs),
                                     tuple(tuple(output[z] >> term & 1 for z in range(2))
                                           for term in range(rank)))
                found.require_exact()
                break
            if found is not None:
                break
        if found is None:
            raise AssertionError("inventory did not find a scheme through rank three")
        for left, right in product(product(range(2), repeat=2), repeat=2):
            if found.apply(left, right) != target.contract(left, right):
                raise AssertionError("inventory certificate fails an input pair")
        ranks = flattening_ranks(target)
        counts[found.terms] += 1
        invertible = 0
        for mask in (1, 2, 3):
            matrix = tuple(tuple(sum((mask >> z & 1)*target.coefficient(x, y, z)
                                     for z in range(2)) % 2 for y in range(2)) for x in range(2))
            invertible += matrix[0][0]*matrix[1][1] ^ matrix[0][1]*matrix[1][0]
        if ranks == (2, 2, 2):
            concise[found.terms] += 1
            pencils[(found.terms, invertible)] += 1
        metadata[code] = (found.terms, ranks, invertible)
    gl = tuple(_map_from_bits(f, 2, 2, bits) for bits in range(16)
               if gf2_rank((bits & 3, bits >> 2), 2) == 2)
    unseen = set(range(256))
    orbits = []
    while unseen:
        seed = min(unseen)
        orbit = set()
        for maps in product(gl, repeat=3):
            transformed = tensors[seed].restrict(*maps)
            code = sum(value << i for i, value in enumerate(transformed.coefficients))
            orbit.add(code)
        if any(metadata[code] != metadata[seed] for code in orbit):
            raise AssertionError("orbit enumeration changes exact rank, flattenings, or pencil count")
        if not orbit <= unseen:
            raise AssertionError("invertible-map orbits unexpectedly overlap")
        unseen -= orbit
        rank, ranks, invertible = metadata[seed]
        orbits.append({"canonical_code": seed, "members": len(orbit), "exact_rank": rank,
                       "flattening_ranks": ranks, "invertible_output_slices": invertible})
    return {"status": "exhausted", "tensors": 256,
            "all_found_coefficients_independently_verified": True,
            "all_16_input_pairs_verified_per_tensor": True,
            "exact_rank_counts": dict(sorted(counts.items())),
            "concise_exact_rank_counts": dict(sorted(concise.items())),
            "concise_pencil_counts": [{"rank": rank, "invertible_slices": inv, "count": count}
                                       for (rank, inv), count in sorted(pencils.items())],
            "invertible_local_maps_per_leg": len(gl), "orbits": orbits}


def _summary(result):
    output = {"status": result.status, "reason": result.reason,
              "checked": result.checked, "total_candidates": result.total_candidates}
    if result.maps is not None:
        output["maps"] = [mapping.rows for mapping in result.maps]
    if result.scheme is not None:
        output["scheme"] = {"field_characteristic": result.scheme.field.p,
                            "field_modulus": result.scheme.field.modulus,
                            "shape": result.scheme.shape, "terms": result.scheme.terms,
                            "a": result.scheme.a, "b": result.scheme.b, "c": result.scheme.c,
                            "all_coefficients_verified": result.scheme.verify()}
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("demo", "sat", "convolution-sat", "inventory", "families", "catalyst", "embedding"), default="demo")
    parser.add_argument("--size", type=int, default=2)
    parser.add_argument("--rank", type=int, default=7)
    parser.add_argument("--a", type=int, default=2, help="first convolution input dimension")
    parser.add_argument("--b", type=int, default=3, help="second convolution input dimension")
    parser.add_argument("--field", choices=("F2", "F4"), default="F2")
    parser.add_argument("--pencil-kind", choices=("regular", "singular"), default="regular")
    parser.add_argument("--timeout-ms", type=int, default=30000)
    parser.add_argument("--solver-path")
    parser.add_argument("--output", help="optional JSON certificate/report path")
    args = parser.parse_args()
    if args.mode == "sat":
        results = {f"unseeded_{args.size}x{args.size}_rank_at_most_{args.rank}": _summary(unseeded_square_sat(
            args.size, args.rank, timeout_ms=args.timeout_ms, solver_path=args.solver_path))}
    elif args.mode == "convolution-sat":
        field = FiniteField(2) if args.field == "F2" else FiniteField(2, (1, 1, 1))
        results = {f"unseeded_C{args.a}{args.b}_rank_at_most_{args.rank}": _summary(unseeded_convolution_sat(
            field, args.a, args.b, args.rank, timeout_ms=args.timeout_ms, solver_path=args.solver_path))}
    elif args.mode == "inventory":
        results = cube_inventory()
    elif args.mode == "families":
        results = convolution_family_checks()
    elif args.mode == "catalyst":
        results = bounded_s233_catalyst_experiment(pencil_kind=args.pencil_kind,
                                                  timeout_ms=args.timeout_ms, solver_path=args.solver_path)
    elif args.mode == "embedding":
        field = FiniteField(2) if args.field == "F2" else FiniteField(2, (1, 1, 1))
        results = bounded_singular_rectangular_embedding(field=field, timeout_ms=args.timeout_ms, solver_path=args.solver_path)
    else:
        f = FiniteField(2)
        # Genuine unseeded tiny local restriction. Every pair is enumerated.
        source = diagonal_tensor(f, 4)
        target = matrix_multiplication_tensor(f, 2, 1, 2)
        found = bounded_restriction_search(source, target, max_pairs=65536)
        assert found.status == "found"
        # Fixed 2x2x2 target with three active support triples. Exhaustion
        # establishes that diagonal rank2 cannot produce it, despite all
        # three individual flattening ranks being two.
        impossible = FiniteTensor.from_function(f, (2, 2, 2),
                                                lambda x, y, z: int((x, y, z) in ((0, 0, 1), (0, 1, 0), (1, 0, 0))))
        exhausted = bounded_restriction_search(diagonal_tensor(f, 2), impossible, max_pairs=256)
        assert exhausted.status == "exhausted"
        capped = bounded_restriction_search(source, target, max_pairs=1)
        assert capped.status == "resource_cap"
        filtered = bounded_restriction_search(diagonal_tensor(f, 1), target, max_pairs=100)
        assert filtered.status == "filtered"
        results = {"found": _summary(found), "exhausted": _summary(exhausted),
                   "resource_cap": _summary(capped), "flattening_filter": _summary(filtered)}
    text = json.dumps(results, indent=2)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as report:
            report.write(text+"\n")
    print(text)


if __name__ == "__main__":
    main()
