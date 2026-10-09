"""Node-free exact extraction from a supplied formal tensor degeneration.

Each polynomial rank-one term contributes coefficient vectors at degree
triples i+j+k=leading. No field evaluations, interpolation nodes or extension
are involved. The resulting length is an upper bound, not an optimal-rank
or performance claim.
"""

from itertools import product
from math import prod

from .catalyst_compiler import CompilerLimits, CompilerResourceLimit
from .tensor_degenerations import TensorDegeneration
from .tensor_schemes import TensorScheme


DEFAULT_LIMITS = CompilerLimits()


def _degree_triples(supports, leading):
    """Inspect the smallest two supports, with constant-time membership in the third."""
    order = sorted(range(3), key=lambda axis: len(supports[axis]))
    first, second, third = order
    for i, j in product(supports[first], supports[second]):
        k = leading-i-j
        if k in supports[third]:
            values = [0, 0, 0]
            values[first], values[second], values[third] = i, j, k
            yield tuple(values)


def _preflight(degeneration, limits):
    """Bound validation, degree matching and output before formal multiplication.

    max_scheme_scalars also conservatively bounds stored input polynomial
    coefficients, the formal coordinate accumulator, degree-pair inspections,
    formal coefficient operations, and coordinate bookkeeping. These are counts, not elapsed-time or peak
    Python-memory estimates. Both supplied and extracted term counts are capped.
    """
    shape, leading = degeneration.shape, degeneration.leading
    limits.tensor(shape)
    if sum(shape) > limits.max_scheme_scalars:
        raise CompilerResourceLimit("coordinate bookkeeping cap reached; extraction incomplete")
    if leading+1 > limits.max_scheme_scalars:
        raise CompilerResourceLimit("leading-degree verification cap reached; extraction incomplete")
    limits.scheme(shape, degeneration.terms)
    families = degeneration.a, degeneration.b, degeneration.c
    input_slots = sum(len(polynomial.coefficients) for family in families for row in family for polynomial in row)
    if input_slots > limits.max_scheme_scalars:
        raise CompilerResourceLimit("input polynomial coefficient cap reached; extraction incomplete")
    degree_bound = degeneration.degree_bound
    formal_slots = (degree_bound+1)*prod(shape)
    if input_slots+formal_slots > limits.max_scheme_scalars:
        raise CompilerResourceLimit("formal polynomial validation cap reached; extraction incomplete")
    all_supports, inspected_pairs, validation_work, output_terms = [], 0, 0, 0
    for rows in zip(*families):
        counts = tuple(tuple(sum(bool(value) for value in polynomial.coefficients) for polynomial in row)
                       for row in rows)
        nonzero = tuple(sum(count) for count in counts)
        active = tuple(sum(bool(n) for n in count) for count in counts)
        slots = tuple(sum(len(polynomial.coefficients) for polynomial in row) for row in rows)
        # TensorDegeneration multiplies each active coordinate's first two
        # polynomials and then the third, before accumulating coefficients.
        # These conservative operation bounds prevent small-leading inputs
        # with enormous dense higher-degree products from bypassing preflight.
        validation_work += (slots[0]*active[1]*active[2]+nonzero[0]*slots[1]*active[2]
                            +nonzero[0]*nonzero[1]*slots[2]
                            +2*(degree_bound+1)*prod(active))
        if validation_work > limits.max_scheme_scalars:
            raise CompilerResourceLimit("formal coefficient operation cap reached; extraction incomplete")
        supports = tuple(frozenset(degree for polynomial in row
                                   for degree, coefficient in enumerate(polynomial.coefficients)
                                   if coefficient and degree <= leading) for row in rows)
        sizes = sorted(len(support) for support in supports)
        inspected_pairs += sizes[0]*sizes[1]
        if inspected_pairs > limits.max_scheme_scalars:
            raise CompilerResourceLimit("degree-pair inspection cap reached; extraction incomplete")
        for _ in _degree_triples(supports, leading):
            output_terms += 1
            limits.scheme(shape, output_terms)
        all_supports.append(supports)
    return tuple(all_supports), output_terms


def extract_leading_coefficient(degeneration, *, limits=DEFAULT_LIMITS):
    """Return a fully verified exact TensorScheme from a TensorDegeneration.

    For r polynomial terms the elementary length bound is
    binomial(leading+2,2)*r. This implementation counts only triples with a
    nonzero coefficient vector on every leg, so it can be much shorter.
    Supports are enumerated directly; no loop runs up to the leading degree.
    Resource caps report incomplete extraction, never an impossibility result.
    """
    if not isinstance(degeneration, TensorDegeneration) or not isinstance(limits, CompilerLimits):
        raise ValueError("supply an actual TensorDegeneration and explicit compiler limits")
    supports, count = _preflight(degeneration, limits)
    # The leading target coefficient alone is insufficient: EVERY coefficient
    # below leading must also vanish, including over finite fields.
    degeneration.require_valid()
    aa, bb, cc = [], [], []
    for rows, support in zip(zip(degeneration.a, degeneration.b, degeneration.c), supports):
        for degrees in _degree_triples(support, degeneration.leading):
            vectors = tuple(tuple(polynomial.coefficient(degree) for polynomial in row)
                            for row, degree in zip(rows, degrees))
            if not all(any(vector) for vector in vectors):
                continue
            aa.append(vectors[0])
            bb.append(vectors[1])
            cc.append(vectors[2])
    if len(aa) != count:
        raise AssertionError("coefficient supports and extracted vector count disagree")
    scheme = TensorScheme(degeneration.target, aa, bb, cc)
    scheme.require_exact()
    return scheme
