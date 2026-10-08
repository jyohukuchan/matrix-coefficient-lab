"""Run with: python3 -m reference.demo (from the repository root)."""

from .constructions import (fourier_filter, nonzero_nodes, primitive_root,
                            recover_leading_coefficient, separation_period)
from .fields import FiniteField
from .schemes import descend, naive_multiply, strassen_scheme


def main():
    print("Exact finite-field reference demo; this does not generate a 9/4 scheme.")
    field = FiniteField(2, (1, 1, 1))  # F_4 = F_2[X]/(X^2+X+1)
    # 2 encodes X, so rescaling introduces genuine extension coefficients.
    seed = strassen_scheme(field).rescale_terms(2, 1)
    if not seed.verify():
        raise RuntimeError("seed coefficient verification failed")
    print(f"F_4 test seed: {seed.terms} terms, all tensor coefficients verified.")

    for exponent in (0, 1, 2):
        powered = seed.tensor_power(exponent)
        projected = descend(powered)
        if not projected.verify():
            raise RuntimeError("descended coefficient verification failed")
        print(f"Power {exponent}: size {projected.n}, {powered.terms} extension terms"
              f" -> {projected.terms} base-field terms (fixed overhead 4).")

    # Deliberately compare the two orders: the wrong order repeats overhead.
    projected = descend(seed.tensor_power(2))
    repeated = descend(seed).tensor_power(2)
    print(f"At power 2: descend after powering = {projected.terms} terms;"
          f" power after descending = {repeated.terms} terms.")
    left = [[1, 0, 1, 1], [0, 1, 1, 0], [1, 1, 0, 1], [1, 0, 0, 1]]
    right = [[0, 1, 1, 0], [1, 1, 0, 1], [1, 0, 1, 1], [0, 1, 0, 1]]
    output = projected.multiply(left, right)
    recursive = descend(seed).multiply_recursive(left, right)
    expected = naive_multiply(projected.field, left, right)
    if output != expected or recursive != expected:
        raise RuntimeError("matrix multiplication check failed")
    print(f"F_2 4x4 product, exact scheme and recursive block evaluation agree: {output}")

    coefficients = (0, 0, 2, 3, 1)  # X^2 * (alpha + (alpha+1)*X + X^2)
    nodes = nonzero_nodes(field, 3)
    recovered = recover_leading_coefficient(field, coefficients, 2, nodes)
    if recovered != 2:
        raise RuntimeError("interpolation check failed")
    print(f"F_4 interpolation at nonzero nodes {nodes}: recovered encoded coefficient {recovered}.")

    for roots_field in (FiniteField(2, (1, 1, 0, 0, 1)), FiniteField(5, (2, 0, 1))):
        period = separation_period(roots_field, 1)
        root = primitive_root(roots_field, period)
        exponents = range(-period, period + 1)
        values = [fourier_filter(roots_field, period, root, e) for e in exponents]
        if values != [int(e % period == 0) for e in exponents]:
            raise RuntimeError("Fourier check failed")
        print(f"F_{roots_field.order} Fourier: period {period}, encoded root {root};"
              f" filter selects precisely multiples of {period} (negative exponents included).")


if __name__ == "__main__":
    main()
