"""Run actual catalytic coefficient compilation: python3 -m reference.catalyst_demo."""

from .catalyst_compiler import (CompilerResourceLimit, compile_absorption, compile_gain_absorption,
                               compile_gain_powered_catalyst,
                               compile_powered_catalyst, coordinate_catalyst)
from .schemes import naive_multiply, naive_scheme
from .tensor_schemes import TensorScheme


def main():
    certificate = coordinate_catalyst()
    D, S = TensorScheme.from_tensor(certificate.D), TensorScheme.from_tensor(certificate.S)
    powered = compile_powered_catalyst(certificate, D)
    seed = naive_scheme(certificate.field, 1, 1, 1)
    print("Actual catalytic maps and coefficient arrays; d=2,k=9 does not supply a 9/4 gap.")
    for _ in range(2):
        step = compile_absorption(powered, seed, S)
        size = step.matrix_scheme.n
        left = [[(i+j) % 2 for j in range(size)] for i in range(size)]
        right = [[(i+1) % 2 for j in range(size)] for i in range(size)]
        if step.matrix_scheme.multiply(left, right) != naive_multiply(certificate.field, left, right):
            raise RuntimeError("compiled catalyst scheme fails direct matrix multiplication")
        print(f"  {seed.n}x{seed.n} seed with {seed.terms} terms -> {size}x{size} with"
              f" {step.matrix_scheme.terms} terms; all maps and coefficients verified.")
        seed = step.matrix_scheme
    powered_gain = compile_gain_powered_catalyst(certificate, D, 2)
    step = compile_gain_absorption(powered_gain, naive_scheme(certificate.field, 1, 1, 1),
                                   S, powered_gain.catalyst_scheme)
    print(f"  Powered positive gains b=2: {powered_gain.gain} independent units ->"
          f" exact M_{powered_gain.matrix_size} with {step.matrix_scheme.terms} terms.")
    power_two = compile_powered_catalyst(certificate, D, 2)
    step = compile_absorption(power_two, naive_scheme(certificate.field, 1, 1, 1), S)
    print(f"  Powered catalyst b=2: {power_two.copies} S copies -> M_{power_two.matrix_size};"
          f" explicit composed maps yield {step.matrix_scheme.terms} terms.")
    try:
        compile_absorption(powered, seed, S)
    except CompilerResourceLimit as error:
        print(f"  Further dense construction stopped: {error}.")
    seed = naive_scheme(certificate.field, 1, 1, 1)
    print("  Retain positive scalar gains and remove only the independent scalar summands:")
    for _ in range(2):
        step = compile_gain_absorption(certificate, seed, S, D)
        print(f"    Size {step.matrix_scheme.n}: {step.raw_gain_scheme.terms} raw terms,"
              f" {step.gain_removed} gains removed -> {step.matrix_scheme.terms} exact matrix terms.")
        seed = step.matrix_scheme


if __name__ == "__main__":
    main()
