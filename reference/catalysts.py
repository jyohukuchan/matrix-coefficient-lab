"""Finite catalytic restriction certificates and conditional integer planning.

The certificate checks an actual local restriction coefficient by coefficient.
The planner ONLY evaluates a proposed cost recurrence: it constructs neither
matrix decompositions nor a proof of an asymptotic rank bound.
"""

from dataclasses import dataclass
from fractions import Fraction
from functools import cached_property

from .maps import LinearMap, validate_local_maps
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, matrix_multiplication_tensor


def _positive_integer(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _natural_integer(value, name):
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


@dataclass(frozen=True)
class CatalyticCertificate:
    """Certify D ⊕ m*unit ⊕ (M_d tensor S) <= D ⊕ k*S.

    <= here means the supplied three independent LinearMap restrictions.
    Copies are independent full direct-sum blocks. In particular m is a
    positive INTEGER count even in characteristic p dividing m.
    """

    d: int
    k: int
    m: int
    D: FiniteTensor
    S: FiniteTensor
    maps: tuple[LinearMap, LinearMap, LinearMap]

    def __post_init__(self):
        _positive_integer(self.d, "matrix size")
        if self.d < 2:
            raise ValueError("matrix size must be at least two")
        _positive_integer(self.k, "source copy count")
        _positive_integer(self.m, "scalar gain count")
        if not isinstance(self.D, FiniteTensor) or not isinstance(self.S, FiniteTensor):
            raise ValueError("catalyst and auxiliary tensor must be finite tensors")
        if self.D.field != self.S.field:
            raise ValueError("catalyst and auxiliary tensor must use the same field")
        if not any(self.S.coefficients):
            raise ValueError("auxiliary tensor S must be nonzero")
        maps = validate_local_maps(self.field, self.source.shape, self.maps)
        if tuple(mapping.output_size for mapping in maps) != self.target.shape:
            raise ValueError("restriction output axes must match the exact catalytic target")
        object.__setattr__(self, "maps", maps)

    @property
    def field(self):
        return self.S.field

    @cached_property
    def source(self):
        return FiniteTensor.direct_sum((self.D,)+(self.S,)*self.k)

    @cached_property
    def matrix_factor(self):
        return matrix_multiplication_tensor(self.field, self.d, self.d, self.d).tensor_product(self.S)

    @cached_property
    def target(self):
        return FiniteTensor.direct_sum((self.D,)+(FiniteTensor.unit(self.field),)*self.m+(self.matrix_factor,))

    @cached_property
    def restricted(self):
        return self.source.restrict(*self.maps)

    @cached_property
    def _first_mismatch(self):
        for index, (actual, expected) in enumerate(zip(self.restricted.coefficients, self.target.coefficients)):
            if actual != expected:
                ij, z = divmod(index, self.target.shape[2])
                x, y = divmod(ij, self.target.shape[1])
                return (x, y, z), expected, actual
        return None

    def verify(self):
        return self._first_mismatch is None

    def require_valid(self):
        if self._first_mismatch is not None:
            coordinate, expected, actual = self._first_mismatch
            raise ValueError(f"invalid catalytic coefficient {coordinate}: expected {expected}, got {actual}")

    @cached_property
    def core_target(self):
        return FiniteTensor.direct_sum((self.D, self.matrix_factor))

    @cached_property
    def gain_removal_maps(self):
        """Output selections deleting scalar blocks without changing D or M_d*S."""
        return tuple(LinearMap.selection(self.field, output,
                                        tuple(range(catalyst))+tuple(range(catalyst+self.m, output)))
                     for catalyst, output in zip(self.D.shape, self.target.shape))

    def strip_gain(self):
        """Return three SOURCE-to-CORE maps for D ⊕ M_d*S <= D ⊕ k*S.

        These are composed output selections, not subtraction or division by
        the scalar m. Require the original certificate before stripping.
        """
        self.require_valid()
        maps = tuple(selection.compose(mapping) for selection, mapping in zip(self.gain_removal_maps, self.maps))
        if self.source.restrict(*maps) != self.core_target:
            raise ValueError("stripped maps do not preserve the exact core target")
        return maps


@dataclass(frozen=True)
class ConditionalCostPlan:
    """Integer recurrence output, conditional on a separate catalyst compiler.

    status='conditional' means ONLY the requested numerical inequality was
    reached. status='resource_cap' means incomplete search. A failed gap test
    reports a missing sufficient condition, never impossibility of a scheme.
    predicted_rank is not a decomposition certificate or a proven rank bound.
    """

    status: str
    reason: str
    tau: Fraction
    b: int | None
    iterations: int
    predicted_size: int
    predicted_rank: int
    catalyst_cost: int | None
    extension_degree: int

    @property
    def numeric_inequality(self):
        A, B = self.tau.numerator, self.tau.denominator
        return (self.extension_degree**2*self.predicted_rank)**B < self.predicted_size**A


def plan_conditional_costs(d, k, s_rank, catalyst_rank, extension_degree, tau,
                           *, max_b=64, max_iterations=64):
    """Plan a SYNTHETIC/CONDITIONAL cost recurrence using exact integer tests.

    Raw integer costs carry no scheme certificate. For bounds attached to
    verified schemes use plan_catalytic_parameters. Both APIs still lack a
    compiler that realizes the iterated recurrence as a decomposition.
    """
    _positive_integer(d, "matrix size")
    if d < 2:
        raise ValueError("matrix size must be at least two")
    _positive_integer(k, "source copy count")
    _positive_integer(s_rank, "auxiliary rank bound")
    _natural_integer(catalyst_rank, "catalyst rank bound")
    _positive_integer(extension_degree, "extension degree")
    _natural_integer(max_b, "block search cap")
    _natural_integer(max_iterations, "iteration cap")
    if not isinstance(tau, Fraction) or tau <= 0:
        raise ValueError("target exponent must be a positive Fraction")
    A, B = tau.numerator, tau.denominator
    if k**B >= d**A:
        return ConditionalCostPlan("gap_not_certified", "k^B < d^A fails; this sufficient criterion gives no conclusion",
                                   tau, None, 0, 1, 1, None, extension_degree)
    b = next((b for b in range(1, max_b+1) if (s_rank*k**b)**B < d**(b*A)), None)
    if b is None:
        return ConditionalCostPlan("resource_cap", "block search cap reached; search is incomplete",
                                   tau, None, 0, 1, 1, None, extension_degree)
    C_b = catalyst_rank*sum((d**3)**i*k**(b-1-i) for i in range(b))
    u, r = 1, 1
    for step in range(1, max_iterations+1):
        u *= d**b
        r = s_rank*(k**b*r+C_b)
        if (extension_degree**2*r)**B < u**A:
            return ConditionalCostPlan("conditional", "numerical recurrence reaches the strict inequality; decomposition not constructed",
                                       tau, b, step, u, r, C_b, extension_degree)
    return ConditionalCostPlan("resource_cap", "iteration cap reached; search is incomplete",
                               tau, b, max_iterations, u, r, C_b, extension_degree)


def plan_catalytic_parameters(certificate, auxiliary_scheme, catalyst_scheme, tau,
                              *, max_b=64, max_iterations=64):
    """Use supplied exact decompositions as rank upper bounds, then plan costs.

    Verify the actual catalytic restriction and both input schemes first.
    This does not certify the proposed compiler's iterated matrix rank bounds.
    """
    if not isinstance(certificate, CatalyticCertificate):
        raise ValueError("planning requires a catalytic certificate")
    certificate.require_valid()
    for scheme, target in ((auxiliary_scheme, certificate.S), (catalyst_scheme, certificate.D)):
        if not isinstance(scheme, TensorScheme) or scheme.target != target:
            raise ValueError("rank upper bound must come from an exact scheme for the specified tensor")
        scheme.require_exact()
    return plan_conditional_costs(certificate.d, certificate.k, auxiliary_scheme.terms,
                                  catalyst_scheme.terms, certificate.field.degree, tau,
                                  max_b=max_b, max_iterations=max_iterations)
