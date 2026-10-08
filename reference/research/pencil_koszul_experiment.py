"""Exact, bounded GF2 checks on four concise (2,3,3) pencil auxiliaries.

This produces verified auxiliary decompositions and projected Koszul ranks,
not a useful catalyst. Shared-first repetition checks the matrix-rank
amplification without assuming tensor-rank direct-sum additivity.
"""

from collections import Counter
from fractions import Fraction
from random import Random

from ..catalyst_search_experiment import (bounded_tensor_rank_search, gf2_rank,
                                          koszul_rank_bound)
from ..fields import FiniteField
from ..maps import LinearMap
from ..rectangular_restrictions import share_first_axis
from ..tensors import FiniteTensor, matrix_multiplication_tensor


FORMS = {
    "Jordan3": (((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                ((0, 1, 0), (0, 0, 1), (0, 0, 0))),
    "Jordan2_plus_0": (((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                       ((0, 1, 0), (0, 0, 0), (0, 0, 0))),
    "Jordan2_plus_1": (((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                       ((0, 1, 0), (0, 0, 0), (0, 0, 1))),
    "L1_plus_L1_transpose": (((1, 0, 0), (0, 0, 1), (0, 0, 0)),
                              ((0, 1, 0), (0, 0, 0), (0, 0, 1))),
}


def first_projection(tensor, masks):
    width = tensor.shape[0]
    if any(type(mask) is not int or not 0 < mask < (1 << width) for mask in masks):
        raise ValueError("projection masks must fit the input axis")
    if gf2_rank(masks, width) != len(masks):
        raise ValueError("projection must have full row rank")
    mapping = LinearMap(tensor.field, width,
                        tuple(tuple(mask >> i & 1 for i in range(width)) for mask in masks))
    return tensor.restrict(mapping, LinearMap.identity(tensor.field, tensor.shape[1]),
                           LinearMap.identity(tensor.field, tensor.shape[2]))


def hyperplane_quotients(width):
    """One representative for every one-dimensional kernel over GF2."""
    for kernel in range(1, 1 << width):
        pivot = (kernel & -kernel).bit_length()-1
        yield tuple((1 << j) | ((kernel >> j & 1) << pivot)
                    for j in range(width) if j != pivot)


def random_quotients(width, output, trials, seed=0):
    random = Random(seed)
    for _ in range(trials):
        while True:
            masks = tuple(random.randrange(1, 1 << width) for _ in range(output))
            if gf2_rank(masks, width) == output:
                yield masks
                break


def projection_checks(tensor, projections, degree):
    histogram, best = Counter(), None
    for masks in projections:
        projected = first_projection(tensor, masks)
        result = koszul_rank_bound(projected, degree, max_bit_entries=2_000_000)
        assert result["status"] == "checked"
        histogram[result["matrix_rank"]] += 1
        if best is None or result["matrix_rank"] > best["rank"]:
            best = {"rank": result["matrix_rank"], "factor": result["rank_one_factor"],
                    "masks": masks, "tensor": projected}
    return histogram, best


def verify_shared_amplification(tensor, degree, expected_rank):
    """Keep two units, merge first tensor coordinates, then isolate the core."""
    merged = share_first_axis(tensor, 2, scalar_gain=1)
    merged.require_valid()
    for scalar in range(2):
        assert merged.target.coefficient(scalar, scalar, scalar) == 1
    core = merged.target.restrict(*(LinearMap.selection(tensor.field, size, range(2, size))
                                   for size in merged.target.shape))
    result = koszul_rank_bound(core, degree, max_bit_entries=10_000_000)
    assert result["status"] == "checked" and result["matrix_rank"] == 2*expected_rank
    single = koszul_rank_bound(tensor, degree, max_bit_entries=2_000_000)
    assert result["rank_one_factor"] == single["rank_one_factor"]
    return result


def main():
    field = FiniteField(2)
    matrix = matrix_multiplication_tensor(field, 2, 2, 2)
    for name, slices in FORMS.items():
        auxiliary = FiniteTensor.from_function(field, (2, 3, 3),
                                               lambda x, y, z: slices[x][y][z])
        search = bounded_tensor_rank_search(auxiliary, max_rank=4, max_subsets=10_000)
        assert search.status == "found" and search.scheme.terms == 4
        search.scheme.require_exact()
        tensor = matrix.tensor_product(auxiliary)
        unprojected = [koszul_rank_bound(tensor, p, max_bit_entries=10_000_000)
                       for p in range(1, 7)]
        assert all(result["status"] == "checked" for result in unprojected)
        seven_hist, seven = projection_checks(tensor, hyperplane_quotients(8), 3)
        five_hist, five = projection_checks(tensor, random_quotients(8, 5, 256), 2)
        assert sum(seven_hist.values()) == 255 and sum(five_hist.values()) == 256
        verify_shared_amplification(seven["tensor"], 3, seven["rank"])
        verify_shared_amplification(five["tensor"], 2, five["rank"])
        ratios = [Fraction(result["matrix_rank"], result["rank_one_factor"])
                  for result in unprojected]
        ratios.extend((Fraction(seven["rank"], seven["factor"]),
                       Fraction(five["rank"], five["factor"])))
        print(name)
        print("  verified auxiliary rank4 scheme after", search.checked, "subsets")
        print("  unprojected ratios:", ", ".join(str(ratio) for ratio in ratios[:6]))
        print("  all255 codim1 quotient rank histogram:", dict(sorted(seven_hist.items())))
        print("  seed0 256 dimension5 quotient histogram:", dict(sorted(five_hist.items())))
        print("  best projected ratios:", ratios[-2], ratios[-1])
        print("  first dimension5 witness masks:", five["masks"])
        print("  q2 coefficient merge and exact Koszul rank doubling verified")
        assert max(ratios) <= 19
    print("No checked ratio exceeds19; these Koszul bounds do not exclude a rank20q budget.")
    print("Regular-pencil exclusion uses the separate commuting-source argument, not these ranks.")


if __name__ == "__main__":
    main()
