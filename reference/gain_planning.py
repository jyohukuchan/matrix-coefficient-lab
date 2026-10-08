"""Verified-input gain planning and bounded coefficient compilation.

Planning gives integer costs, not coefficient arrays. Compilation carries one
represented extension throughout all absorption steps and descends only once.
Resource caps explicitly leave the proposed construction incomplete.
"""

from dataclasses import dataclass
from fractions import Fraction

from .catalysts import CatalyticCertificate
from .catalyst_compiler import (DEFAULT_LIMITS, CompilerResourceLimit,
                               compile_gain_absorption, compile_gain_powered_catalyst)
from .schemes import BilinearScheme, descend, naive_scheme
from .tensor_schemes import TensorScheme


@dataclass(frozen=True)
class GainCostPlan:
    status: str
    reason: str
    tau: Fraction
    exponent: int | None
    iterations: int
    predicted_size: int
    predicted_terms: int
    catalyst_terms: int | None
    powered_gain: int | None
    multiplier: int | None
    extension_degree: int

    @property
    def numeric_inequality(self):
        return ((self.extension_degree**2*self.predicted_terms)**self.tau.denominator
                < self.predicted_size**self.tau.numerator)


def _inputs(certificate, auxiliary_scheme, catalyst_scheme, tau, max_b, max_iterations):
    if not isinstance(certificate, CatalyticCertificate):
        raise ValueError("supply an actual positive-gain catalyst")
    if not isinstance(tau, Fraction) or tau <= 0:
        raise ValueError("target exponent must be a positive Fraction")
    if any(type(cap) is not int or cap < 0 for cap in (max_b, max_iterations)):
        raise ValueError("search caps must be nonnegative integers")
    certificate.require_valid()
    for scheme, target in ((auxiliary_scheme, certificate.S), (catalyst_scheme, certificate.D)):
        if not isinstance(scheme, TensorScheme) or scheme.target != target:
            raise ValueError("supply an exact decomposition of the specified tensor")
        scheme.require_exact()


def plan_gain_catalytic_parameters(certificate, auxiliary_scheme, catalyst_scheme, tau,
                                   *, max_b=64, max_iterations=64):
    """Search sufficient fixed-block recurrences using supplied exact schemes.

    Q_b=k^b*R_S-g_b; r_next=C_b+Q_b*r. Unlike the gain-free planner,
    a useful small block is accepted even when k^B<d^A itself fails. No
    exhausted finite cap implies impossibility for other witnesses or bounds.
    """
    _inputs(certificate, auxiliary_scheme, catalyst_scheme, tau, max_b, max_iterations)
    d, k, m = certificate.d, certificate.k, certificate.m
    Rs, CD, e = auxiliary_scheme.terms, catalyst_scheme.terms, certificate.field.degree
    A, B = tau.numerator, tau.denominator
    gain, copies, size, cost = 0, 1, 1, 0
    selected = None
    for b in range(1, max_b+1):
        # D_(b+1)=D*M_(d^b) + k*D_b; matrix gains give m*d^b units.
        gain, cost = k*gain+m*size, k*cost+CD*size**3
        copies, size = k*copies, d*size
        Q = copies*Rs-gain
        if Q < 1:
            raise ValueError("verified catalyst violates the positive supplied-cost invariant")
        if Q**B < size**A:
            selected = (b, size, cost, gain, Q)
            break
    if selected is None:
        return GainCostPlan("resource_cap", "block cap reached; sufficient-block search incomplete",
                            tau, None, 0, 1, 1, None, None, None, e)
    b, block_size, cost, gain, Q = selected
    u, r = 1, 1
    for step in range(1, max_iterations+1):
        u, r = u*block_size, cost+Q*r
        if (e**2*r)**B < u**A:
            return GainCostPlan("planned", "integer costs reach the strict inequality; arrays not yet constructed",
                                tau, b, step, u, r, cost, gain, Q, e)
    return GainCostPlan("resource_cap", "iteration cap reached; construction search incomplete",
                        tau, b, max_iterations, u, r, cost, gain, Q, e)


@dataclass(frozen=True)
class GainCompilationResult:
    status: str
    reason: str
    plan: GainCostPlan
    extension_scheme: BilinearScheme | None = None
    prime_scheme: BilinearScheme | None = None


def compile_gain_target(certificate, auxiliary_scheme, catalyst_scheme, tau,
                        *, max_b=64, max_iterations=64, limits=DEFAULT_LIMITS):
    """Plan, generate all coefficient arrays, descend once, and verify the gap.

    A result with status 'verified' contains both exact schemes. A resource_cap
    result contains a numerical plan only and certifies no final decomposition.
    The current descent targets the represented field's prime subfield.
    """
    plan = plan_gain_catalytic_parameters(certificate, auxiliary_scheme, catalyst_scheme, tau,
                                         max_b=max_b, max_iterations=max_iterations)
    if plan.status != "planned":
        return GainCompilationResult(plan.status, plan.reason, plan)
    try:
        # Preflight the FINAL allocation before constructing any large prefix.
        shape = (plan.predicted_size**2,)*3
        limits.tensor(shape)
        limits.scheme(shape, plan.extension_degree**2*plan.predicted_terms)
        powered = compile_gain_powered_catalyst(certificate, catalyst_scheme, plan.exponent,
                                                limits=limits)
        seed = naive_scheme(certificate.field, 1, 1, 1)
        for _ in range(plan.iterations):
            step = compile_gain_absorption(powered, seed, auxiliary_scheme,
                                           powered.catalyst_scheme, limits=limits)
            seed = step.matrix_scheme
        if (seed.n, seed.terms) != (plan.predicted_size, plan.predicted_terms):
            raise ValueError("generated arrays disagree with the supplied-cost recurrence")
        prime = descend(seed)
        prime.require_exact()
        if prime.terms != plan.extension_degree**2*seed.terms:
            raise ValueError("final descent disagrees with the one-time overhead")
        if prime.terms**tau.denominator >= prime.n**tau.numerator:
            raise ValueError("generated arrays do not meet the strict target inequality")
        return GainCompilationResult("verified", "exact arrays generated and checked after one final descent",
                                     plan, seed, prime)
    except CompilerResourceLimit as error:
        return GainCompilationResult("resource_cap", str(error), plan)
