"""Small exact matrix checks for the catalyst-obstructions research note.

Run with python3 -m reference.research.catalyst_obstructions_experiment.
This checks flattenings and commutators, not catalyst existence/nonexistence
by exhaustive search and not a tensor-rank decomposition.
"""

from itertools import combinations, product
from math import comb

from reference.fields import FiniteField
from reference.maps import LinearMap
from reference.tensor_schemes import TensorScheme
from reference.tensors import FiniteTensor, matrix_multiplication_tensor


def matrix_rank(field, matrix):
    rows = [list(row) for row in matrix]
    if not rows:
        return 0
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("ragged matrix")
    rank = 0
    for column in range(width):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = field.inv(rows[rank][column])
        rows[rank] = [field.mul(inverse, value) for value in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][column]:
                scale = rows[i][column]
                rows[i] = [field.sub(a, field.mul(scale, b))
                           for a, b in zip(rows[i], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def flattening_ranks(tensor):
    ranks = []
    for axis in range(3):
        other = [i for i in range(3) if i != axis]
        matrix = []
        for row in range(tensor.shape[axis]):
            entries = []
            for pair in product(*(range(tensor.shape[i]) for i in other)):
                coordinate = [0, 0, 0]
                coordinate[axis] = row
                for i, value in zip(other, pair):
                    coordinate[i] = value
                entries.append(tensor.coefficient(*coordinate))
            matrix.append(entries)
        ranks.append(matrix_rank(tensor.field, matrix))
    return tuple(ranks)


def product_matrix(field, left, right):
    return [[field.sum(field.mul(a, right[h][j]) for h, a in enumerate(row))
             for j in range(len(right[0]))] for row in left]


def slice_matrix(tensor, coefficients):
    return [[tensor.field.sum(tensor.field.mul(coefficients[x], tensor.coefficient(x, y, z))
                              for x in range(tensor.shape[0]))
             for z in range(tensor.shape[2])] for y in range(tensor.shape[1])]


def dot_tensor(field, rank):
    return FiniteTensor.from_function(field, (1, rank, rank), lambda x, y, z: int(y == z))


def cube_orbit_partition():
    """Exhaust all F2 cubes under all six independent GL2 changes per leg."""
    field = FiniteField(2)
    gl = tuple(m for m in product((0, 1), repeat=4)
               if (m[0]*m[3]) ^ (m[1]*m[2]))
    assert len(gl) == 6
    maps = tuple(LinearMap(field, 2, (m[:2], m[2:])) for m in gl)

    def tensor(mask):
        return FiniteTensor(field, (2, 2, 2), tuple((mask >> i) & 1 for i in range(8)))

    def mask(T):
        return sum(value << i for i, value in enumerate(T.coefficients))

    concise = {x for x in range(256) if flattening_ranks(tensor(x)) == (2, 2, 2)}
    expected = {"split": (129, 108), "dual_numbers": (22, 54), "F4": (233, 12)}
    orbits = {}
    covered = set()
    for name, (seed, size) in expected.items():
        orbit = {mask(tensor(seed).restrict(*triple)) for triple in product(maps, repeat=3)}
        assert len(orbit) == size and orbit <= concise and not orbit & covered
        covered |= orbit
        orbits[name] = orbit
    assert covered == concise and len(concise) == 174

    # Independent explicit rank upper certificates for the algebra representatives.
    for name, u, v in (("dual_numbers", 0, 0), ("F4", 1, 1)):
        table = FiniteTensor(field, (2, 2, 2), (1, 0, 0, 1, 0, 1, u, v))
        scheme = TensorScheme(table, ((1, 0), (0, 1), (1, 1)),
                              ((1, 0), (0, 1), (1, 1)),
                              ((1, 1), (u, field.sub(v, 1)), (0, 1)))
        assert scheme.verify() and scheme.terms == 3
        if name == "dual_numbers":
            swap = LinearMap(field, 2, ((0, 1), (1, 0)))
            identity = LinearMap.identity(field, 2)
            assert table.restrict(identity, identity, swap) == tensor(22)
        else:
            extension = FiniteField(2, (1, 1, 1))
            assert table == FiniteTensor.from_function(
                field, (2, 2, 2), lambda x, y, z: extension.coordinates(
                    extension.mul(1 << x, 1 << y))[z])
            assert table == tensor(233)
    assert TensorScheme.from_tensor(tensor(129)).terms == 2
    return {name: len(orbit) for name, orbit in orbits.items()}


def auxiliary_223_partition():
    """Exhaust the two GL2×GL2×GL3 orbits of concise F2 (2,2,3) tensors."""
    field = FiniteField(2)

    def invertible_maps(size):
        result = []
        for entries in product((0, 1), repeat=size*size):
            rows = tuple(entries[i*size:(i+1)*size] for i in range(size))
            if matrix_rank(field, rows) == size:
                result.append(LinearMap(field, size, rows))
        return tuple(result)

    gl2, gl3 = invertible_maps(2), invertible_maps(3)
    assert len(gl2) == 6 and len(gl3) == 168

    def tensor(mask):
        return FiniteTensor(field, (2, 2, 3), tuple((mask >> i) & 1 for i in range(12)))

    concise = {x for x in range(4096) if flattening_ranks(tensor(x)) == (2, 2, 3)}
    covered = set()
    counts = {}
    for name, seed, size in (("rank_one_kernel", 273, 1512), ("rank_two_kernel", 2193, 1008)):
        T = tensor(seed)
        orbit = {sum(value << i for i, value in enumerate(T.restrict(*maps).coefficients))
                 for maps in product(gl2, gl2, gl3)}
        assert len(orbit) == size and not covered & orbit and orbit <= concise
        covered |= orbit
        counts[name] = size
        if name == "rank_one_kernel":
            scheme = TensorScheme.from_tensor(T)
            output = LinearMap(field, 3, ((1, 0, 0), (0, 1, 1)))
        else:
            scheme = TensorScheme(T, ((1, 0), (0, 1), (1, 1)),
                                  ((1, 0), (0, 1), (1, 1)),
                                  ((1, 1, 0), (0, 1, 1), (0, 1, 0)))
            output = LinearMap(field, 3, ((1, 0, 0), (0, 1, 0)))
        assert scheme.verify() and scheme.terms == 3
        identity = LinearMap.identity(field, 2)
        dual_numbers = FiniteTensor(field, (2, 2, 2), (1, 0, 0, 1, 0, 1, 0, 0))
        assert T.restrict(identity, identity, output) == dual_numbers
    assert covered == concise and len(concise) == 2520
    return counts


def gf2_koszul_rank(tensor, p):
    """Exact GF2 Koszul matrix, B tensor wedge^p(A) -> C tensor wedge^(p+1)(A)."""
    if tensor.field != FiniteField(2):
        raise ValueError("this bit-matrix experiment is only over F2")
    a, b, c = tensor.shape
    inputs = tuple(combinations(range(a), p))
    outputs = tuple(combinations(range(a), p+1))
    out_index = {subset: i for i, subset in enumerate(outputs)}
    rows = [0]*(len(outputs)*c)
    supports = tuple(tuple((j, k) for j in range(b) for k in range(c)
                           if tensor.coefficient(i, j, k)) for i in range(a))
    for h, subset in enumerate(inputs):
        for i in range(a):
            if i in subset:
                continue
            output = out_index[tuple(sorted((*subset, i)))]
            for j, k in supports[i]:
                rows[output*c+k] ^= 1 << (h*b+j)
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length()-1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots), len(rows), len(inputs)*b, comb(a-1, p)


def verify_rectangular_shared_merge():
    """Full coefficient check: q independent M424 copies restrict to M42(4q)."""
    for field, q in product((FiniteField(2), FiniteField(3)), (1, 2)):
        factor = matrix_multiplication_tensor(field, 4, 2, 4)
        copies = FiniteTensor.direct_sum((factor,)*q)
        first = LinearMap(field, 8*q,
                          tuple(tuple(int(i % 8 == x) for i in range(8*q)) for x in range(8)))
        second_indices = tuple(block*8+j*4+h
                               for j in range(2) for block in range(q) for h in range(4))
        third_indices = tuple(block*16+i*4+h
                              for i in range(4) for block in range(q) for h in range(4))
        second = LinearMap.selection(field, 8*q, second_indices)
        third = LinearMap.selection(field, 16*q, third_indices)
        expected = matrix_multiplication_tensor(field, 4, 2, 4*q)
        assert copies.restrict(first, second, third) == expected


def main():
    fields = (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)), FiniteField(5))
    cases = 0
    for field, r, copies in product(fields, (1, 2), (1, 2)):
        unit = FiniteTensor.unit(field)
        S = dot_tensor(field, r)
        T = matrix_multiplication_tensor(field, 2, 2, 2).tensor_product(S)
        assert flattening_ranks(S) == (1, r, r)
        assert flattening_ranks(T) == (4, 4*r, 4*r)
        source = FiniteTensor.direct_sum((unit,)+(S,)*5)
        target = FiniteTensor.direct_sum((unit, unit, T))
        assert flattening_ranks(source) == (6, 1+5*r, 1+5*r)
        assert flattening_ranks(target) == (6, 2+4*r, 2+4*r)

        repeated = FiniteTensor.direct_sum((unit, T)*copies)
        identity = [0]*repeated.shape[0]
        off12 = identity.copy()
        off21 = identity.copy()
        for q in range(copies):
            start = 5*q
            identity[start] = identity[start+1] = identity[start+4] = 1
            off12[start+2] = off21[start+3] = 1
        A, B, C = (slice_matrix(repeated, v) for v in (identity, off12, off21))
        size = copies*(1+4*r)
        assert A == [[int(i == j) for j in range(size)] for i in range(size)]
        BC, CB = product_matrix(field, B, C), product_matrix(field, C, B)
        commutator = [[field.sub(a, b) for a, b in zip(x, y)] for x, y in zip(BC, CB)]
        commutator_rank = matrix_rank(field, commutator)
        assert commutator_rank == 4*copies*r
        lower_bound = size + (commutator_rank+1)//2
        assert lower_bound == copies*(1+6*r)
        assert lower_bound > 5*copies*r
        cases += 1
    print(f"Verified {cases} exact flattening/commutator cases over F2, F3, F4, F5.")
    print("For S=dot_r, repeated target rank lower bound is q*(1+6*r),")
    print("while a d=2,k=5 catalyst would give upper bound rank(D)+5*q*r.")
    print("These checks produce no catalyst coefficient witness.")
    sizes = cube_orbit_partition()
    print(f"Exhausted all 256 F2 cubes: concise orbit sizes {sizes}, total 174.")
    print("Algebra representative rank-upper certificates are verified; the")
    print("Alder-Strassen lower bound used to exclude them is a cited theorem.")
    sizes = auxiliary_223_partition()
    print(f"Exhausted all 4096 F2 (2,2,3) tensors: concise orbit sizes {sizes}, total 2520.")
    print("Both representatives have verified rank3 schemes and exact maps to dual numbers.")
    field = FiniteField(2)
    rectangular = matrix_multiplication_tensor(field, 4, 2, 4)
    target = FiniteTensor.direct_sum((FiniteTensor.unit(field), rectangular))
    rank, rows, columns, denominator = gf2_koszul_rank(target, 5)
    assert (rank, rows, columns, denominator) == (960, 1428, 1134, 56)
    print(f"Unit plus M(4,2,4): Koszul p5 rank{rank}/{denominator},")
    print("giving tensor-rank lower18; this Koszul check alone does not exclude rank20.")
    verify_rectangular_shared_merge()
    print("Verified shared-first merge of independent M424 copies into M42(4q),")
    print("at every coefficient over F2/F3 for q1/q2; lower24q uses the cited theorem.")


if __name__ == "__main__":
    main()
