"""Actual finite dual -> catalyst -> coefficients: python3 -m reference.witnessed_demo."""

import json
from pathlib import Path

from .catalyst_compiler import compile_absorption, compile_powered_catalyst
from .extraction import simplify_scheme
from .fields import FiniteField
from .schemes import naive_multiply, naive_scheme
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, matrix_multiplication_tensor
from .witnessed_constraints import (DetectorConstraint, WitnessedStateSystem,
                                    assemble_catalyst, clear_rational_dual,
                                    find_nonnegative_dual, rank_upper_constraint)


def main():
    # This is the output of the unseeded experiment, now a supplied checked
    # input certificate. Replaying the demo does not run or require a solver.
    artifact = Path(__file__).with_name("research") / "unseeded-f2-rank7.json"
    data = json.loads(artifact.read_text())
    field = FiniteField(2)
    matrix = TensorScheme(matrix_multiplication_tensor(field, 2, 2, 2), data["a"], data["b"], data["c"])
    matrix.require_exact()
    system = WitnessedStateSystem(2, 8, (("unit", FiniteTensor.unit(field)), ("matrix", matrix.target)),
                                 (rank_upper_constraint("unit", "matrix", matrix),
                                  DetectorConstraint("unit", "matrix")))
    dual = find_nonnegative_dual(system)
    if dual.status != "found":
        raise RuntimeError("finite toy dual was not found")
    integer = clear_rational_dual(system, dual.coefficients)
    certificate = assemble_catalyst(integer)
    D = TensorScheme.from_tensor(certificate.D)
    S = TensorScheme.from_tensor(certificate.S)
    powered = compile_powered_catalyst(certificate, D)
    seed = naive_scheme(field, 1, 1, 1)
    print("Witnessed finite dual -> actual catalytic maps -> actual coefficient arrays.")
    print("Known seven-term input certificate from an unseeded search; no 9/4 catalyst.")
    print(f"  Dual weights {integer.weights}, integer gain {integer.m}; d=2,k=8,D has {D.terms} terms.")
    for _ in range(2):
        step = compile_absorption(powered, seed, S)
        compact = simplify_scheme(step.matrix_scheme.as_tensor_scheme()).to_matrix_scheme(
            step.matrix_scheme.n, step.matrix_scheme.n, step.matrix_scheme.n)
        n = compact.n
        left = [[(i+j) % 2 for j in range(n)] for i in range(n)]
        right = [[int(i == j) for j in range(n)] for i in range(n)]
        if compact.multiply(left, right) != naive_multiply(field, left, right):
            raise RuntimeError("compiled dual disagrees with ordinary matrix multiplication")
        print(f"  Size {n}: {step.matrix_scheme.terms} generated terms -> {compact.terms} after exact merging;"
              " every tensor coefficient verified.")
        seed = compact


if __name__ == "__main__":
    main()
