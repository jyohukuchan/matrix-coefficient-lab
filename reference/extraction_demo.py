"""Run the proof-derived boundary matrix path: python3 -m reference.extraction_demo."""

from .extraction import proof_matrix_pipeline
from .schemes import naive_multiply


def main():
    construction = proof_matrix_pipeline()
    scheme = construction.matrix_scheme
    left, right = [[1, 2], [7, 13]], [[9, 4], [3, 11]]
    result = scheme.multiply(left, right)
    if result != naive_multiply(scheme.field, left, right):
        raise RuntimeError("generated matrix scheme disagrees with direct multiplication")
    print("Proof-derived finite boundary path; no 9/4 or speed claim.")
    print("  " + " -> ".join(f"{name}: {terms} terms" for name, terms in construction.stage_terms))
    print(f"  Fourier period {construction.period}, encoded root {construction.root}.")
    print("  Branch counts before / after proportional merging / singleton factorization:")
    for branch, raw, merged, factored in zip(("left", "middle", "right"), construction.branch_schemes,
                                            construction.merged_branches, construction.factored_branches):
        print(f"    {branch}: {raw.terms} / {merged.terms} / {factored.terms}")
    print(f"  Exact 2x2 matrix multiplication over F_16: {result}; all 64 tensor coefficients verified.")
    for exponent in (1, 2):
        base = construction.descend_power(exponent)
        base.require_exact()
        print(f"  Matrix power {exponent}, then one fixed F_16/F_2 descent:"
              f" size {base.n}, {base.terms} terms (fixed overhead 16).")


if __name__ == "__main__":
    main()
