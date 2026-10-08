"""Exact finite evidence for the D=0 singular-star obstruction.

Run with python3 -m reference.research.star_slice_obstruction_experiment.
The general proof is in star-zero-catalyst.md; these checks do not search maps
and do not exclude a nonzero arbitrary catalyst D.
"""

from itertools import product

from reference.fields import FiniteField
from .catalyst_obstructions_experiment import matrix_rank


def full_slice_count(q):
    return q**8-2*q**5-q**4+2*q**3


def bad_slice_count(q):
    return 2*q**5+q**4-2*q**3


def hyperplane_full_source_count(q, active_blocks):
    a = q*q-1
    return a**(5-active_blocks)*(a**active_blocks+(q-1)*(-1)**active_blocks)//q


def main():
    fields = (FiniteField(2), FiniteField(3), FiniteField(2, (1, 1, 1)))
    for field in fields:
        q, full = field.order, 0
        for entries in product(range(q), repeat=8):
            A = (entries[:2], entries[2:4])
            B = (entries[4:6], entries[6:8])
            stacked = (*A, *B)
            horizontal = tuple(a+b for a, b in zip(A, B))
            if matrix_rank(field, stacked) == matrix_rank(field, horizontal) == 2:
                full += 1
        if full != full_slice_count(q) or q**8-full != bad_slice_count(q):
            raise RuntimeError("independent slice enumeration disagrees with the formula")
        if q >= 3 and not bad_slice_count(q) < q**6:
            raise RuntimeError("the finite-extension fiber inequality fails")
        print(f"F_{q}: all {q**8} matrix pairs checked; full rank-eight slices {full},"
              f" bad pairs {q**8-full}; affine fiber budget {q**6}.")

    f = fields[0]
    hist = [0]*5
    nonzero_pairs = ((0, 1), (1, 0), (1, 1))
    for pairs in product(nonzero_pairs, repeat=5):
        value = 0
        for t, pair in enumerate(pairs):
            value = f.add(value, pair[0])
            if value == 0:
                hist[t] += 1
    if hist != [hyperplane_full_source_count(2, t) for t in range(1, 6)]:
        raise RuntimeError("binary hyperplane count disagrees with enumeration")
    if any(value >= full_slice_count(2) for value in hist):
        raise RuntimeError("binary rank-nine counting contradiction is absent")
    print(f"F_2 hyperplane source rank-ten counts: {hist}; target rank-nine count192.")

    # Check the two rank-seven subspaces for EVERY binary hyperplane normal.
    for normal in product(range(2), repeat=10):
        if not any(normal):
            continue
        active = next(i for i in range(5) if any(normal[2*i:2*i+2]))
        i, j = [h for h in range(5) if h != active][:2]
        coordinate = lambda h: tuple(int(k == h) for k in range(10))
        left = (normal, coordinate(2*i), coordinate(2*i+1))
        right = (normal, coordinate(2*j), coordinate(2*j+1))
        joined = left+right[1:]
        if (10-matrix_rank(f, left), 10-matrix_rank(f, right),
            10-matrix_rank(f, joined)) != (7, 7, 5):
            raise RuntimeError("the hyperplane spanning dimensions fail")
    print("All 1,023 binary hyperplanes have the required 7+7−5=9 spanning pair.")
    print("No catalyst maps were searched. This checks the D=0 proof ingredients only.")


if __name__ == "__main__":
    main()
