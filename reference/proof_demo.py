"""Run the finite proof constructions: python3 -m reference.proof_demo."""

from .convolution import convolution_scheme
from .fields import FiniteField
from .schemes import descend
from .sectors import ThreeSectorConstruction
from .separation import FiniteSeparation


def _check_bilinear(scheme):
    scheme.require_exact()
    left = tuple((i+1) % scheme.field.order for i in range(scheme.shape[0]))
    right = tuple((5*j+1) % scheme.field.order for j in range(scheme.shape[1]))
    if scheme.apply(left, right) != scheme.target.contract(left, right):
        raise RuntimeError("proof construction does not match direct bilinear contraction")


def main():
    print("Finite proof constructions; no speed or 9/4 matrix scheme claim.")
    field = FiniteField(2, (1, 1, 0, 1))  # F_8
    c = ThreeSectorConstruction(field, 2, 2)
    source = convolution_scheme(field, c.a, c.source_width)
    exact = c.recover_power(2, source_scheme=source)
    base = descend(exact)
    _check_bilinear(source)
    _check_bilinear(exact)
    _check_bilinear(base)
    print(f"C(2,7) evaluation/interpolation: {source.terms} source terms (coordinate decomposition: 14).")
    print(f"  Proof regions -> power 2 -> 3 nonzero nodes: {exact.terms} F_8 terms"
          f" -> {base.terms} F_2 terms, shape {base.shape}; all coefficients verified.")

    field = FiniteField(2, (1, 1, 0, 0, 1))  # F_16 supports order 15 and five nodes
    c = ThreeSectorConstruction(field, 2, 1)
    source = convolution_scheme(field, c.a, c.source_width)
    retained = c.recover_power(1, source_scheme=source)
    tagged = c.tag_scheme(retained)
    separation = FiniteSeparation(c.branches)
    if tagged.target != separation.source:
        raise RuntimeError("tagged regions do not match the separation source")
    exact = separation.recover_scheme(source_scheme=tagged)
    direct = separation.as_direct_sum(exact)
    base = descend(direct)
    for scheme in (tagged, exact, direct, base):
        _check_bilinear(scheme)
    print(f"C(2,4) -> retained regions: {source.terms} -> {retained.terms} terms;"
          f" independent side tags preserve {tagged.shape[0]} shared first-input coordinates.")
    print(f"  Fourier period {separation.period}, encoded primitive root {separation.root}:"
          f" {separation.period} complete source copies ->"
          f" {separation.period*tagged.terms} projected terms.")
    print(f"  Integer square weights, polynomial shift {separation.leading},"
          f" normalized degree <= {separation.normalized_degree_bound}:"
          f" 5 nonzero nodes -> {exact.terms} exact F_16 terms.")
    print(f"  Full direct sum of 3 branches tensor dot_3: shape {direct.shape},"
          f" {base.terms} F_2 terms; all coefficients and direct contraction agree.")


if __name__ == "__main__":
    main()
