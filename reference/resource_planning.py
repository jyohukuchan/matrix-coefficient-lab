"""Integer-only forecasts; a numerical hypothesis is never a rank certificate."""

from dataclasses import asdict, dataclass
from fractions import Fraction

from .gain_planning import _inputs


@dataclass(frozen=True)
class ExpansionBudget:
    """Allocation counts, not peak Python memory or elapsed time estimates."""

    max_tensor_entries: int = 2_000_000
    max_scheme_scalars: int = 2_000_000
    max_terms: int = 100_000
    max_powered_map_entries: int = 2_000_000

    def __post_init__(self):
        if any(type(v) is not int or v < 0 for v in asdict(self).values()):
            raise ValueError("expansion budgets must be nonnegative integers")


@dataclass(frozen=True)
class CandidateCosts:
    """Supplied assumptions only. Use forecast_verified_inputs for witnesses."""

    d: int
    k: int
    m: int
    auxiliary_terms: int
    catalyst_terms: int
    auxiliary_shape: tuple[int, int, int]
    catalyst_shape: tuple[int, int, int]
    characteristic: int = 2
    extension_degree: int = 1

    def __post_init__(self):
        for name in ("d", "k", "m", "auxiliary_terms", "characteristic", "extension_degree"):
            value = getattr(self, name)
            if type(value) is not int or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        if self.d < 2 or self.k <= self.d or self.characteristic < 2:
            raise ValueError("require d>=2, k>d and characteristic>=2")
        if type(self.catalyst_terms) is not int or self.catalyst_terms < 0:
            raise ValueError("catalyst terms must be a nonnegative integer")
        for name in ("auxiliary_shape", "catalyst_shape"):
            shape = tuple(getattr(self, name))
            if len(shape) != 3 or any(type(x) is not int or x < 0 for x in shape):
                raise ValueError("shapes must have three nonnegative integer dimensions")
            object.__setattr__(self, name, shape)


@dataclass(frozen=True)
class BlockForecast:
    block_exponent: int
    block_size: int
    multiplier: int
    powered_gain: int
    catalyst_terms: int
    iterations: int
    matrix_size: int
    extension_terms: int
    prime_terms: int
    strict_gap: bool
    tensor_entries: int
    scheme_scalars: int
    packed_coefficient_bytes: int
    powered_map_entries: int
    budget_failures: tuple[str, ...]

    @property
    def within_expansion_budget(self):
        return not self.budget_failures


@dataclass(frozen=True)
class CostForecast:
    input_status: str
    output_status: str
    target_exponent: Fraction
    blocks: tuple[BlockForecast, ...]

    def to_dict(self):
        return {"input_status": self.input_status, "output_status": self.output_status,
                "target_exponent": str(self.target_exponent),
                "blocks": [asdict(block) for block in self.blocks]}


def forecast_candidate(costs, tau, *, max_b=32, max_iterations=64,
                       budget=ExpansionBudget()):
    """Compare all capped blocks, not just the first sufficient block.

    Each block records its first strict-gap iteration, or the last capped
    iteration. Packed bytes describe the three final coefficient families
    with ideal bit packing; they are NOT Python memory or a runtime forecast.
    Budget checks cover the final dense target/families and powered catalyst
    maps. Intermediate absorption allocations can be larger still.
    """
    if not isinstance(costs, CandidateCosts) or not isinstance(budget, ExpansionBudget):
        raise ValueError("supply candidate costs and an expansion budget")
    if not isinstance(tau, Fraction) or tau <= 0:
        raise ValueError("target exponent must be a positive Fraction")
    if any(type(cap) is not int or cap < 0 for cap in (max_b, max_iterations)):
        raise ValueError("forecast caps must be nonnegative integers")
    d, k, m = costs.d, costs.k, costs.m
    gain, size, copies, C = 0, 1, 1, 0
    Dshape = (0, 0, 0)
    blocks = []
    for b in range(1, max_b+1):
        gain, C = k*gain+m*size, k*C+costs.catalyst_terms*size**3
        Dshape = tuple(k*x+y*size**2 for x, y in zip(Dshape, costs.catalyst_shape))
        size, copies = d*size, k*copies
        Q = copies*costs.auxiliary_terms-gain
        if Q < 1:
            raise ValueError("supplied assumptions give a nonpositive multiplier")
        source = tuple(x+copies*y for x, y in zip(Dshape, costs.auxiliary_shape))
        target = tuple(x+gain+size**2*y for x, y in zip(Dshape, costs.auxiliary_shape))
        map_entries = sum(x*y for x, y in zip(source, target))
        n, r, strict, step = 1, 1, False, 0
        for step in range(1, max_iterations+1):
            n, r = n*size, C+Q*r
            strict = (costs.extension_degree**2*r)**tau.denominator < n**tau.numerator
            if strict:
                break
        prime_terms = costs.extension_degree**2*r
        entries, scalars = n**6, prime_terms*3*n**2
        # Canonical encodings 0..p-1 after one final descent.
        packed = (scalars*(costs.characteristic-1).bit_length()+7)//8
        quantities = {"tensor_entries": entries, "scheme_scalars": scalars,
                      "terms": prime_terms, "powered_map_entries": map_entries}
        failures = tuple(name for name, value in quantities.items()
                         if value > getattr(budget, "max_"+name))
        blocks.append(BlockForecast(b, size, Q, gain, C, step, n, r, prime_terms,
                                    strict, entries, scalars, packed, map_entries, failures))
    return CostForecast("unverified_assumptions", "numerical_only", tau, tuple(blocks))


def forecast_verified_inputs(certificate, auxiliary_scheme, catalyst_scheme, tau, **kwargs):
    """Validate actual inputs; even then the forecast constructs no output."""
    _inputs(certificate, auxiliary_scheme, catalyst_scheme, tau,
            kwargs.get("max_b", 32), kwargs.get("max_iterations", 64))
    costs = CandidateCosts(certificate.d, certificate.k, certificate.m,
                           auxiliary_scheme.terms, catalyst_scheme.terms,
                           certificate.S.shape, certificate.D.shape,
                           certificate.field.p, certificate.field.degree)
    result = forecast_candidate(costs, tau, **kwargs)
    return CostForecast("verified_inputs", result.output_status, result.target_exponent, result.blocks)
